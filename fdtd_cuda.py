import numpy as np

from numba import cuda 

from misc import (
    timer
    )



class SimulationCUDA:
    
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
        
        #pml boundary condition parameters
        self.Nx_lo = pml_layers
        self.Nx_hi = pml_layers
        self.Ny_lo = pml_layers
        self.Ny_hi = pml_layers
        
        #some values for use outside
        self.E_max = 0
        
        #materials
        self.materials = []
        
        #setup update coefficients according to boundary conditions
        self._setup()
        
        
        
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
        
        M_eps_r = cuda.to_device( Material.bitmap * Material.eps_r  )
        add_material_cuda[self.blockspergrid, self.threadsperblock](self.Eps_z, M_eps_r)
        
        M_mu_r  = Material.bitmap * Material.mu_r
        self.Mu_x = self.Mu_y = np.where(M_mu_r >= 1, M_mu_r, self.Mu_x)
    
        
    
    def add_material(self, Material):
        
        self.materials.append(Material)
        
        self._apply_material(Material)
        
       
        
        
        
        
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
        
        self.m_Hx_1 = cuda.to_device((1 / self.dt - Sigma_Hy / (2 * self.eps_0)) / m_Hx_0)
        self.m_Hx_2 = cuda.to_device(-self.c_0 / self.Mu_x / m_Hx_0)
        self.m_Hx_3 = cuda.to_device(-self.c_0 * self.dt / self.eps_0 * Sigma_Hx / self.Mu_x / m_Hx_0)
        
        
        #update coefficients for Hy component
        
        Sigma_Hx = self.Sigma_x[::2, 1::2]
        Sigma_Hy = self.Sigma_y[::2, 1::2]
        
        m_Hy_0 = 1 / self.dt + Sigma_Hx / (2 * self.eps_0)
        
        self.m_Hy_1 = cuda.to_device((1 / self.dt - Sigma_Hx / (2 * self.eps_0)) / m_Hy_0)
        self.m_Hy_2 = cuda.to_device(-self.c_0 / self.Mu_y / m_Hy_0)
        self.m_Hy_3 = cuda.to_device(-self.c_0 * self.dt / self.eps_0 * Sigma_Hy / self.Mu_y / m_Hy_0)
        
        
        #update coefficients for Dz
        
        Sigma_Dx = self.Sigma_x[::2, ::2]
        Sigma_Dy = self.Sigma_y[::2, ::2]
        
        m_Dz_0 = 1 / self.dt + (Sigma_Dx + Sigma_Dy) / (2 * self.eps_0) + Sigma_Dx * Sigma_Dy * self.dt / (4 * self.eps_0**2)
        
        self.m_Dz_1 = cuda.to_device((1 / self.dt - (Sigma_Dx + Sigma_Dy) / (2 * self.eps_0) - Sigma_Dx * Sigma_Dy * self.dt / (2 * self.eps_0)**2) / m_Dz_0)
        self.m_Dz_2 = cuda.to_device(self.c_0 / m_Dz_0)
        self.m_Dz_4 = cuda.to_device(-self.dt / self.eps_0**2 * Sigma_Dx * Sigma_Dy / m_Dz_0)
        
        
    # @timer
    def _setup(self):
        
        
        self.total_steps = 0
        self.E_max       = 1e-5
        
        #setup kernel specificatioins
        tpb_x, tpb_y = 32, 32
        self.threadsperblock = (tpb_x, tpb_y)
        
        bpg_x, bpg_y = int(np.ceil(self.n_x / tpb_x)), int(np.ceil(self.n_y / tpb_y))
        self.blockspergrid = (bpg_x, bpg_y)
        
        #dummy variable
        Z = np.zeros((self.n_x, self.n_y))
        
        #initialize fields
        self.Ez = cuda.to_device(Z)
        self.Dz = cuda.to_device(Z)
        self.Hx = cuda.to_device(Z)
        self.Hy = cuda.to_device(Z)
        
        #auxiliary terms for curl calculation
        self.C_Ex  = cuda.to_device(Z)
        self.C_Ey  = cuda.to_device(Z)
        self.C_Hzx = cuda.to_device(Z)
        self.C_Hzy = cuda.to_device(Z)
        
        #auxiliary terms for integration
        self.I_C_Ex = cuda.to_device(Z)
        self.I_C_Ey = cuda.to_device(Z)
        self.I_Dz   = cuda.to_device(Z)
        
        #init materials
        self.Eps_z = cuda.to_device(Z+1)
        self.Mu_x = Z+1
        self.Mu_y = Z+1
        
        #add materials
        for Material in self.materials:
            self._apply_material(Material)
        
        #precomputations in initialization
        self._compute_pml_conductivity()
        self._compute_update_coeffs()
        
        
        
        
        
        
        
    
    # updating with pml ----------------------------------------------------------
    
    def _update_Hx_from_E_cuda(self):
        
        """
        update the x component of the H field 
        from the z component of the E field
        
        functions are implemented as cuda kernels 
        and directly access data on GPU
        """
        
        #compute x component of curl of E
        fast_curl_y_plus_cuda[self.blockspergrid, self.threadsperblock](self.Ez, 
                                                                        self.C_Ex, 
                                                                        self.dy, 
                                                                        self.boundary_conditions["left"] == "dirichlet")
        
        # cuda.synchronize()
        
        #integration
        compute_integration_cuda[self.blockspergrid, self.threadsperblock](self.I_C_Ex, 
                                                                           self.C_Ex)
        
        # cuda.synchronize()
        
        #update Hx
        compute_Hxy_update_cuda[self.blockspergrid, self.threadsperblock](self.Hx, 
                                                                          self.C_Ex, 
                                                                          self.I_C_Ex, 
                                                                          self.m_Hx_1, 
                                                                          self.m_Hx_2, 
                                                                          self.m_Hx_3)
        
        # cuda.synchronize()
        
    
    def _update_Hy_from_E_cuda(self):
        
        """
        update the y component of the H field 
        from the z component of the E field
        
        functions are implemented as cuda kernels 
        and directly access data on GPU
        """
        
        #compute y component of curl of E
        fast_curl_x_plus_cuda[self.blockspergrid, self.threadsperblock](self.Ez, 
                                                                        self.C_Ey, 
                                                                        self.dx, 
                                                                        self.boundary_conditions["top"] == "dirichlet")
        
        # cuda.synchronize()
        
        #integration
        compute_integration_cuda[self.blockspergrid, self.threadsperblock](self.I_C_Ey, 
                                                                           self.C_Ey)
        
        # cuda.synchronize()
        
        #update Hy
        compute_Hxy_update_cuda[self.blockspergrid, self.threadsperblock](self.Hy, 
                                                                          self.C_Ey, 
                                                                          self.I_C_Ey, 
                                                                          self.m_Hy_1, 
                                                                          self.m_Hy_2, 
                                                                          self.m_Hy_3)
        
        # cuda.synchronize()
        
        
        
    
    def _update_Dz_from_H_cuda(self):
        
        """
        update the z component of the D field 
        from the x and y component of the H field
        
        functions are implemented as cuda kernels 
        and directly access data on GPU
        """
        
        #compute z component of curl of Hy
        fast_curl_x_minus_cuda[self.blockspergrid, self.threadsperblock](self.Hy, 
                                                                         self.C_Hzx, 
                                                                         self.dy, 
                                                                         self.boundary_conditions["bottom"]=="dirichlet")
        
        # cuda.synchronize()
        
        #compute z component of curl of Hx
        fast_curl_y_minus_cuda[self.blockspergrid, self.threadsperblock](self.Hx, 
                                                                         self.C_Hzy, 
                                                                         self.dx, 
                                                                         self.boundary_conditions["right"]=="dirichlet")
        
        # cuda.synchronize()
        
        #integration of D field
        compute_integration_cuda[self.blockspergrid, self.threadsperblock](self.I_Dz, 
                                                                           self.Dz)
        
        # cuda.synchronize()
        
        #update Dz
        compute_Dz_update_cuda[self.blockspergrid, self.threadsperblock](self.Dz, 
                                                                         self.C_Hzx, 
                                                                         self.C_Hzy, 
                                                                         self.I_Dz, 
                                                                         self.m_Dz_1, 
                                                                         self.m_Dz_2, 
                                                                         self.m_Dz_4)
        
        # cuda.synchronize()
        
        
    def _update_E_from_D(self):
        
        """
        update D to E with dielectric constant
        
        function is implemented as cuda kernels 
        and directly access data on GPU
        """
        
        #calculate E from D
        update_E_from_D_cuda[self.blockspergrid, self.threadsperblock](self.Dz, 
                                                                       self.Ez, 
                                                                       self.Eps_z
                                                                       )
        
        # cuda.synchronize()
    
    # @timer
    def _update_sources(self):
        
        """
        add sources (soft source) to D field
        
        """
        
        #inject source
        for source in self.sources:
            
            pos_y, pos_x = self._pos_to_index(source.position)
            val = source.evaluate(self.dt * self.total_steps)
            
            add_source_D_cuda[self.blockspergrid, self.threadsperblock](pos_y, pos_x, val, self.Dz)
        
        # cuda.synchronize()
        
    
    @timer
    def update(self):
        
        """
        call all _update functions
        
        """
        
        self._update_sources()
        
        self._update_Hx_from_E_cuda()
        
        self._update_Hy_from_E_cuda()
        
        self._update_Dz_from_H_cuda()
        
        self._update_E_from_D()
        
        #update steps
        self.total_steps += 1
        
        #update all time maximum
        self.E_max = max( np.max( np.abs(self.Ez.copy_to_host())), self.E_max )
        
        # print("grid size = ", self.n_x, self.n_y)
        
        
        
    #reset -----------------------------------------------------------------------
    
    def _reset_sources(self):
        #remove sources
        self.sources = []
        
        
    def _reset_fields(self):
        #reset fields
        set_cuda[self.blockspergrid, self.threadsperblock](self.Ez, 0)
        set_cuda[self.blockspergrid, self.threadsperblock](self.Dz, 0)
        set_cuda[self.blockspergrid, self.threadsperblock](self.Hx, 0)
        set_cuda[self.blockspergrid, self.threadsperblock](self.Hy, 0)
        
        cuda.synchronize()
        
        
    def _reset_materials(self):
        #reset materials
        self.materials = []
        set_cuda[self.blockspergrid, self.threadsperblock](self.Eps_z, 1)
        
        cuda.synchronize()
        
        
    def _reset_integration_curl(self):
        
        #auxiliary terms for curl calculation
        set_cuda[self.blockspergrid, self.threadsperblock](self.C_Ex, 0)
        set_cuda[self.blockspergrid, self.threadsperblock](self.C_Ey, 0)
        set_cuda[self.blockspergrid, self.threadsperblock](self.C_Hzx, 0)
        set_cuda[self.blockspergrid, self.threadsperblock](self.C_Hzy, 0)
        
        cuda.synchronize()
        
        #auxiliary terms for integration
        set_cuda[self.blockspergrid, self.threadsperblock](self.I_C_Ex, 0)
        set_cuda[self.blockspergrid, self.threadsperblock](self.I_C_Ey, 0)
        set_cuda[self.blockspergrid, self.threadsperblock](self.I_Dz, 0)
        
        cuda.synchronize()
        
        
    def reset(self):
        
        self._reset_fields()
        self._reset_integration_curl()
        self._reset_sources()
        
        #reset steps
        self.total_steps = 0
        
        #reset alltime maximum
        self.E_max = 1e-5
        
        #rebuild simulation parameters
        self._setup()
        
        
    def resize(self, bounds):
        
        #reset bounds
        x_max, x_min, y_max, y_min = self.bounds = bounds
        
        #reset grid size
        self.n_x = int((x_max-x_min)/self.dx)
        self.n_y = int((y_max-y_min)/self.dy)
        
        #reset simulation
        self.reset()
        
        #rebuild simulation parameters
        # self._setup()
        
        
    def get_Ez(self):
        return self.Ez.copy_to_host()
        





## CUDA KERNELS ===============================================================


@cuda.jit
def add_material_cuda(A, B):
    
    #determine total grid size
    y , x  = cuda.grid(2)
    ya, xa = A.shape
    
    if x < xa and 0 < y < ya:
        if B[y, x] >= 1:
            A[y, x] = B[y, x]
    
    
@cuda.jit
def set_cuda(A, v):
    
    #determine total grid size
    y , x  = cuda.grid(2)
    ya, xa = A.shape
    
    if x < xa and y < ya:
        A[y, x] = v
    
    
@cuda.jit
def fast_curl_y_minus_cuda(A, B, dy, dirichlet=False):
    
    #determine total grid size
    y , x  = cuda.grid(2)
    ya, xa = A.shape
    
    if x < xa and 0 < y < ya:
        B[y, x] = -(A[y, x] - A[y-1, x]) / dy
        
        if dirichlet:
            B[0, x] = -A[0, x] / dy
        else:
            B[0, x] = -(A[0, x] - A[-1, x]) / dy
    
    
@cuda.jit
def fast_curl_y_plus_cuda(A, B, dy, dirichlet=False):
    
    #determine total grid size
    y , x  = cuda.grid(2)
    ya, xa = A.shape
    
    if x < xa and y < ya-1:
        B[y, x] = -(A[y, x] - A[y+1, x]) / dy
        
        if dirichlet:
            B[-1, x] = -A[-1, x] / dy
        else:
            B[-1, x] = -(A[-1, x] - A[0, x]) / dy
    
    
@cuda.jit
def fast_curl_x_minus_cuda(A, B, dx, dirichlet=False):
    
    #determine total grid size
    y , x  = cuda.grid(2)
    ya, xa = A.shape
    
    if 0 < x < xa and y < ya:
        B[y, x] = (A[y, x] - A[y, x-1]) / dx
        
        if dirichlet:
            B[y, 0] = A[y, 0] / dx
        else:
            B[y, 0] = (A[y, 0] - A[y, -1]) / dx
    
    
@cuda.jit
def fast_curl_x_plus_cuda(A, B, dx, dirichlet=False):
    
    #determine total grid size
    y , x  = cuda.grid(2)
    ya, xa = A.shape
    
    if x < xa-1 and y < ya:
        B[y, x] = (A[y, x] - A[y, x+1]) / dx
        
        if dirichlet:
            B[y, -1] = A[y, -1] / dx
        else:
            B[y, -1] = (A[y, -1] - A[y, 0]) / dx
            
            
@cuda.jit
def compute_integration_cuda(I_C, C):
    
    #determine total grid size
    y , x  = cuda.grid(2)
    ya, xa = I_C.shape
    
    #integration
    if x < xa and y < ya:
        I_C[y, x] += C[y, x]
        
        
@cuda.jit
def compute_Hxy_update_cuda(Hxy, C_Exy, I_C_Exy, m_Hxy_1, m_Hxy_2, m_Hxy_3):
    
    #determine total grid size
    y , x  = cuda.grid(2)
    ya, xa = Hxy.shape
    
    #update Hx
    if x < xa and y < ya:
        Hxy[y, x] = m_Hxy_1[y, x] * Hxy[y, x] + m_Hxy_2[y, x] * C_Exy[y, x] + m_Hxy_3[y, x] * I_C_Exy[y, x]
        
        
@cuda.jit
def compute_Dz_update_cuda(Dz, C_Hzx, C_Hzy, I_Dz, m_Dz_1, m_Dz_2, m_Dz_4):
    
    #determine total grid size
    y , x  = cuda.grid(2)
    ya, xa = Dz.shape
    
    #update Dz
    if x < xa and y < ya:
        Dz[y, x] = m_Dz_1[y, x] * Dz[y, x] + m_Dz_2[y, x] * (C_Hzx[y, x] + C_Hzy[y, x]) + m_Dz_4[y, x] * I_Dz[y, x]
    
    
@cuda.jit
def update_source_to_D_cuda(Source_Grid, Dz):
    
    #determine total grid size
    y , x  = cuda.grid(2)
    ya, xa = Dz.shape
    
    if x < xa and y < ya:
        Dz[y, x] += Source_Grid[y, x]
    
    
@cuda.jit
def update_E_from_D_cuda(Dz, Ez, Eps_z):
    
    #determine total grid size
    y , x  = cuda.grid(2)
    ya, xa = Dz.shape
    
    if x < xa and y < ya:
        Ez[y, x] = Dz[y, x] / Eps_z[y, x]
    
    
@cuda.jit
def add_source_D_cuda(pos_y, pos_x, val, D):
    D[pos_y, pos_x] += val
    
    
    