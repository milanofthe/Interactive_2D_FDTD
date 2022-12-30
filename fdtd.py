import numpy as np

import numba

from misc import (
    timer
    )




## FUNCTIONS FOR SPEEDUP CPU ===============================================


# @numba.njit
def fast_curl_y_minus(A, dirichlet=False):
    
    B = np.copy(A)
    
    B[1:] -= B[:-1]
    
    if not dirichlet:
        B[0] -= A[-1]
        
    return B


# @numba.njit
def fast_curl_y_plus(A, dirichlet=False):
    
    B = np.copy(A)
    
    B[:-1] -= B[1:]
    
    if not dirichlet:
        B[-1] -= A[0]
        
    return B


# @numba.njit
def fast_curl_x_plus(A, dirichlet=False):
    
    B = np.copy(A)
    
    B[:,:-1] -= B[:,1:]
    
    if not dirichlet:
        B[:,-1] -= A[:,0]
        
    return B


# @numba.njit
def fast_curl_x_minus(A, dirichlet=False):
    
    B = np.copy(A)
    
    B[:,1:] -= B[:,:-1]
    
    if not dirichlet:
        B[:,0] -= A[:,-1]
        
    return B


# SIM CPU ===============================================================

class Simulation:
    
    def __init__(self, 
                 bounds=[1,-1,1,-1], 
                 dx=6e-3, 
                 dy=6e-3, 
                 pml_layers=15,
                 boundary_conditions={"left"   :"periodic", 
                                      "right"  :"periodic", 
                                      "bottom" :"pml", 
                                      "top"    :"pml"}
                 ):
        
        
        #natural constants
        self.c_0   = 3e8
        self.eps_0 = 8.854187e-12
        self.mu_0  = 4 * np.pi * 1e-7
        
        
        #simulation properties
        self.total_steps         = 0
        self.bounds              = bounds
        self.dx                  = dx
        self.dy                  = dy
        self.dt                  = 0.5 / (self.c_0 * np.sqrt(1/self.dx**2 + 1/self.dy**2)) #CFL
        self.boundary_conditions = boundary_conditions
        
        
        #sources
        self.sources = []
        
        
        #excitation properties
        self.f_max  = 1 / (100 * self.dt)
        
        
        #grid size
        x_max, x_min, y_max, y_min = self.bounds
        
        self.n_x = int((x_max-x_min)/self.dx)
        self.n_y = int((y_max-y_min)/self.dy)
        
        
        #initialize fields
        self.Ez = np.zeros((self.n_x, self.n_y))
        self.Dz = np.zeros((self.n_x, self.n_y))
        self.Hx = np.zeros((self.n_x, self.n_y))
        self.Hy = np.zeros((self.n_x, self.n_y))
        
        
        #auxiliary terms for integration
        self.I_C_Ex = np.zeros((self.n_x, self.n_y))
        self.I_C_Ey = np.zeros((self.n_x, self.n_y))
        self.I_Dz   = np.zeros((self.n_x, self.n_y))
        
        
        #materials
        self.materials = []
        self.Eps_z = np.ones((self.n_x, self.n_y))
        self.Mu_x  = np.ones((self.n_x, self.n_y))
        self.Mu_y  = np.ones((self.n_x, self.n_y))
        
        
        #pml boundary condition parameters
        self.Nx_lo = pml_layers
        self.Nx_hi = pml_layers
        self.Ny_lo = pml_layers
        self.Ny_hi = pml_layers
        
        
        #setup update coefficients according to boundary conditions
        self._setup_boundary_conditions()
        
        
        #some values for use outside
        self.E_max = 0
        
        
        
    # misc -------------------------------------------------------------------------
    
    def _pos_to_index(self, position):
        
        """
        returns the index in simulation grid
        to a position within bounds
        
        """
        
        x, y = position
        x_max, x_min, y_max, y_min = self.bounds
        
        x_ind = min(int((x-x_min)/self.dx), self.n_x - 1)
        y_ind = min(int((y-y_min)/self.dy), self.n_y - 1)
        
        return x_ind, y_ind
    
    
    
    # sources ----------------------------------------------------------------------
    
    def add_source(self, Source):
        Source.to_grid(self.bounds, self.Dz)
        self.sources.append(Source)
    
    
    # materials -------------------------------------------------------------------
    
    def _apply_material(self, Material):
        
        Material.to_grid(self.bounds, self.Eps_z)
        
        M_eps_r = Material.bitmap * Material.eps_r
        M_mu_r  = Material.bitmap * Material.mu_r
        
        
        self.Eps_z = np.where(M_eps_r >= 1, M_eps_r, self.Eps_z)
        self.Mu_x  = self.Mu_y = np.where(M_mu_r >= 1, M_mu_r, self.Mu_x)
        
        self.reset()
        
    
    def add_material(self, Material):
        
        self.materials.append(Material)
        
        self._apply_material(Material)
        
        
        
        
        
    def _resize_materials(self):
        
        self.Eps_z = np.ones((self.n_x, self.n_y))
        self.Mu_x  = np.ones((self.n_x, self.n_y))
        self.Mu_y  = np.ones((self.n_x, self.n_y))
        
        for Material in self.materials:
            
            #add Eps
            
            self._apply_material(Material)
            
            # M = Material.to_grid(self.bounds, self.Eps_z)
        
            # self.Eps_z  = np.where(M!=0, M, self.Eps_z)
            
            
        
    # pml precomputations ---------------------------------------------------------
    
    def _compute_pml_conductivity(self):
        
        #init pml materials on 2x grid
        self.Sigma_x = np.zeros((2*self.n_x, 2*self.n_y))
        self.Sigma_y = np.zeros((2*self.n_x, 2*self.n_y))
        
        #base pml conductivity
        Sigma_max = self.eps_0 / (2 * self.dt)
        
        
        # sigma x left side
        if self.boundary_conditions["bottom"] == "pml":
            for ny in range(self.Ny_lo*2):
                
                ratio = ny / (self.Ny_lo*2)
                
                self.Sigma_x[:, self.Ny_lo*2 - ny - 1] = Sigma_max * ratio**3
             
        # sigma x right side
        if self.boundary_conditions["top"] == "pml":
            for ny in range(self.Ny_hi*2):
            
                ratio = ny / (self.Ny_hi*2)
                
                self.Sigma_x[:, ny + self.n_y*2 - self.Ny_hi*2] = Sigma_max * ratio**3
            
        # sigma y bottom side
        if self.boundary_conditions["left"] == "pml":
            for nx in range(self.Nx_hi*2):
                
                ratio = nx / (self.Nx_lo*2)
                
                self.Sigma_y[self.Nx_lo*2 - nx - 1, :] = Sigma_max * ratio**3
            
        # sigma y top side
        if self.boundary_conditions["right"] == "pml":
            for nx in range(self.Nx_hi*2):
                
                ratio = nx / (self.Nx_hi*2)
                
                self.Sigma_y[nx + self.n_x*2 - self.Nx_hi*2, :] = Sigma_max * ratio**3
            
            
    def _compute_update_coeffs(self):
        
        #update coefficients for Hx component
        
        Sigma_Hx = self.Sigma_x[1::2, ::2]
        Sigma_Hy = self.Sigma_y[1::2, ::2]
        
        m_Hx_0 = 1 / self.dt + Sigma_Hy / (2 * self.eps_0)
        
        self.m_Hx_1 = (1 / self.dt - Sigma_Hy / (2 * self.eps_0)) / m_Hx_0
        self.m_Hx_2 = -self.c_0 / self.Mu_x / m_Hx_0
        self.m_Hx_3 = -self.c_0 * self.dt / self.eps_0 * Sigma_Hx / self.Mu_x / m_Hx_0
        
        
        #update coefficients for Hy component
        
        Sigma_Hx = self.Sigma_x[::2, 1::2]
        Sigma_Hy = self.Sigma_y[::2, 1::2]
        
        m_Hy_0 = 1 / self.dt + Sigma_Hx / (2 * self.eps_0)
        
        self.m_Hy_1 = (1 / self.dt - Sigma_Hx / (2 * self.eps_0)) / m_Hy_0
        self.m_Hy_2 = -self.c_0 / self.Mu_y / m_Hy_0
        self.m_Hy_3 = -self.c_0 * self.dt / self.eps_0 * Sigma_Hy / self.Mu_y / m_Hy_0
        
        
        #update coefficients for Dz
        
        Sigma_Dx = self.Sigma_x[::2, ::2]
        Sigma_Dy = self.Sigma_y[::2, ::2]
        
        m_Dz_0 = 1 / self.dt + (Sigma_Dx + Sigma_Dy) / (2 * self.eps_0) + Sigma_Dx * Sigma_Dy * self.dt / (4 * self.eps_0**2)
        
        self.m_Dz_1 = (1 / self.dt - (Sigma_Dx + Sigma_Dy) / (2 * self.eps_0) - Sigma_Dx * Sigma_Dy * self.dt / (2 * self.eps_0)**2) / m_Dz_0
        self.m_Dz_2 = self.c_0 / m_Dz_0
        self.m_Dz_4 = -self.dt / self.eps_0**2 * Sigma_Dx * Sigma_Dy / m_Dz_0
        
        
        
    # @timer
    def _setup_boundary_conditions(self):
        
        #precomputations in initialization
        self._compute_pml_conductivity()
        self._compute_update_coeffs()
        
        
        
        
    
    # updating with pml ----------------------------------------------------------
    
    # @timer
    def _update_Hx_from_E(self):
        
        #compute x component of curl of E
        C_Ex = - fast_curl_y_plus(self.Ez, dirichlet=self.boundary_conditions["left"] == "dirichlet") / self.dy
        
        #integration
        self.I_C_Ex += C_Ex
        
        #update Hx
        self.Hx = self.m_Hx_1 * self.Hx + self.m_Hx_2 * C_Ex + self.m_Hx_3 * self.I_C_Ex
        
        
    # @timer
    def _update_Hy_from_E(self):
        
        #compute y component of curl of E
        C_Ey = fast_curl_x_plus(self.Ez, dirichlet=self.boundary_conditions["top"] == "dirichlet") / self.dx
        
        #integration
        self.I_C_Ey += C_Ey
        
        #update Hy
        self.Hy = self.m_Hy_1 * self.Hy + self.m_Hy_2 * C_Ey + self.m_Hy_3 * self.I_C_Ey
        
        
    # @timer
    def _update_Dz_from_H(self):
        
        #compute z component of curl of H
        C_Hz =  fast_curl_x_minus(self.Hy, dirichlet=self.boundary_conditions["bottom"] == "dirichlet") / self.dy
        C_Hz -= fast_curl_y_minus(self.Hx, dirichlet=self.boundary_conditions["right"] == "dirichlet") / self.dx
        
        #integration of D field
        self.I_Dz += self.Dz
        
        #update Dz
        self.Dz = self.m_Dz_1 * self.Dz + self.m_Dz_2 * C_Hz + self.m_Dz_4 * self.I_Dz
        
        
    # @timer
    def _update_E_from_D(self):
        
        #update Ez
        self.Ez = 1/self.Eps_z * self.Dz
        
        
    # @timer
    def _update_sources_D(self):
        
        #evaluate all sources
        for source in self.sources:
            self.Dz += source.evaluate_on_grid(self.dt * self.total_steps)
                
            
    # @timer
    def update(self):
        
        #update H
        self._update_Hx_from_E()
        self._update_Hy_from_E()
        
        #update D
        self._update_Dz_from_H()
        
        #inject source
        self._update_sources_D()
        
        #calculate E
        self._update_E_from_D()
        
        #update steps
        self.total_steps += 1
        
        #update all time maximum
        self.E_max = max(self.E_max, np.max(abs(self.Ez)))
        
        
        
    #reset -----------------------------------------------------------------------
    
    def _reset_sources(self):
        #remove sources
        self.sources = []
        
    def _reset_fields(self):
        #reset fields
        self.Ez = np.zeros((self.n_x, self.n_y))
        self.Dz = np.zeros((self.n_x, self.n_y))
        self.Hx = np.zeros((self.n_x, self.n_y))
        self.Hy = np.zeros((self.n_x, self.n_y))
        
    def _reset_materials(self):
        #reset materials
        self.materials = []
        self.Eps_z = np.ones((self.n_x, self.n_y))
        self.Mu_x  = np.ones((self.n_x, self.n_y))
        self.Mu_y  = np.ones((self.n_x, self.n_y))
        
    def _reset_integration(self):
        #auxiliary terms for integration
        self.I_C_Ex = np.zeros((self.n_x, self.n_y))
        self.I_C_Ey = np.zeros((self.n_x, self.n_y))
        self.I_Dz   = np.zeros((self.n_x, self.n_y))
        
        
    def reset(self):
        
        self._reset_fields()
        self._reset_integration()
        self._reset_sources()
        #self._reset_materials()
        
        #reset steps
        self.total_steps = 0
        
        #reset alltime maximum
        self.E_max = 0
        
        #reset boundary conditions
        self._setup_boundary_conditions()
        
        
    def resize(self, bounds):
        
        #reset bounds
        x_max, x_min, y_max, y_min = self.bounds = bounds
        
        #reset discretization
        self.n_x = int((x_max-x_min)/self.dx)
        self.n_y = int((y_max-y_min)/self.dy)
        
        #resize materials
        self._resize_materials()
        
        #reset simulation
        self.reset()
        
        
        
        #reset boundary conditions
        # self._setup_boundary_conditions()
        
    def get_Ez(self):
        return self.Ez
        
        