#imports
import time 
import numpy as np

# MISC ==============================================================

def timer(func):
    
    """
    This function shows the execution time 
    of the function object passed
    
    """
    
    def wrap_func(*args, **kwargs):
        
        t1 = time.perf_counter()
        result = func(*args, **kwargs)
        t2 = time.perf_counter()
        
        print(f'Function {func.__name__!r} executed in {round((t2-t1)*1e3, 3)}ms')
        
        return result
    
    return wrap_func


def sim_to_screen(x, y, bounds, res):
    
    """
    coordinate trasformation from simulation space to screenspace
    
    """
    
    x_max, x_min, y_max, y_min = bounds
    res_x, res_y = res
    
    x_dsp = (x - x_min) / (x_max - x_min) * res_x
    y_dsp = (y - y_min) / (y_max - y_min) * res_y
    
    return x_dsp, y_dsp


def screen_to_sim(x_dsp, y_dsp, bounds, res):
    
    """
    coordinate trasformation from screenspace to simulation space
    
    """
    
    x_max, x_min, y_max, y_min = bounds
    res_x, res_y = res
    
    x = x_min + x_dsp * (x_max - x_min) / res_x 
    y = y_min + (res_y - y_dsp) * (y_max - y_min) / res_y
    
    return x, y


def pos_in_bounds(position, bounds):
    
    """
    check if position lies within bounding rectangle
    (independent of cordinate system)
    
    """
    
    x_max, x_min, y_max, y_min = bounds
    x, y = position
    
    return x_min < x < x_max and y_min < y < y_max


def pos_to_index(position, bounds, grid):
        
    x, y = position
    
    x_max, x_min, y_max, y_min = bounds
    
    n_x, n_y = np.array(grid).shape
    
    dx = (x_max - x_min) / n_x
    dy = (y_max - y_min) / n_y
        
    #col_ind = max(0, min(int((x-x_min)//dx), n_x - 1))
    #row_ind = max(0, min(int((y-y_min)//dy), n_y - 1))
    
    col_ind = int((x-x_min)//dx)
    row_ind = int((y-y_min)//dy)
        
    return row_ind, col_ind
