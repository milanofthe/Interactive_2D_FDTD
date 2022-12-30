
from fdtd import (
    Simulation
    )

from fdtd_cuda import ( 
    SimulationCUDA
    )

from interface import(
    Interface
    )

from materials import (
    ImgMaterial
    )

from configparser import (
    ConfigParser
    )

import sys




# MAIN =======================================================
    
def main():
    
    
    #load settings from config.ini --------------------------- 
    
    P = ConfigParser()
    P.read("config.ini")
    
    fps              = eval(P.get('general', 'fps'))
    resolution       = eval(P.get('general', 'res'))
    color_mode       = eval(P.get('general', 'color_mode'))
    gpu_acceleration = eval(P.get('general', 'gpu_acceleration'))
    
    dx     = eval(P.get('discretization', 'dx'))
    dy     = eval(P.get('discretization', 'dy'))
    bounds = eval(P.get('discretization', 'bounds'))
    
    top_bottom_bc = P.get('boundary_conditions', 'top_bottom')
    left_right_bc = P.get('boundary_conditions', 'left_right')
    pml_layers    = eval(P.get('boundary_conditions', 'pml_layers'))
    
    
    #setup simulation ----------------------------------------
    
    
    
    if gpu_acceleration:
        try:
            S = SimulationCUDA
            print("gpu acceleration enabled")
            
        except:
            S = Simulation
            print("cuda failed -> fallback to cpu")
            
    else:
        
        S = Simulation
    
    
    Sim = S(bounds=bounds, 
            dx=dx, 
            dy=dy, 
            pml_layers=pml_layers,
            boundary_conditions={"left"   : left_right_bc, 
                                 "right"  : left_right_bc, 
                                 "bottom" : top_bottom_bc, 
                                 "top"    : top_bottom_bc
                                 }
            ) 
    
    
    # M = ImgMaterial(eps_r=3000, mu_r=1, path=r"assets\grid02.png")
    
    # Sim.add_material(M)  
    
    
    # setup interface (GUI) ----------------------------------
    
    I = Interface(color_mode=color_mode, 
                  fps=fps, 
                  res=resolution, 
                  steps=2, 
                  Sim=Sim
                  )
    
    I.run()
    
    
if __name__ == '__main__':    
    main()









