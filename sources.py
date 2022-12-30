import numpy as np

from PIL import Image, ImageDraw

import matplotlib.pyplot as plt

from misc import(
    pos_to_index,
    pos_in_bounds
    )

# SOURCE CLASSES ===========================================================


class GaussianSource:
    
    def __init__(self, amplitude=1, f_max=1e9, delay=0, position=(0,0)):
        
        self.delay     = delay
        self.amplitude = amplitude
        self.f_max     = f_max
        self.position  = np.array(position)
        
        
    def to_grid(self, bounds, grid):
        
        pixel_x  , pixel_y  = pos_to_index(self.position, bounds, grid)
        
        self.grid = np.zeros(grid.shape)
        self.grid[pixel_y, pixel_x] = 1 
        
        
    def evaluate(self, t):
        
        tau = 0.5 / self.f_max
        
        val =  np.exp( -( t - 4 * tau - self.delay)**2/tau**2)
                
        return self.amplitude * val
    
    
    def evaluate_on_grid(self, t):
        return self.grid * self.evaluate(t)
        
        
    
    
class SinSource:
    
    def __init__(self, amplitude=1, f_max=1e9, delay=0, position=(0,0)):
        
        self.delay     = delay
        self.amplitude = amplitude
        self.f_max     = f_max
        self.position  = np.array(position)
        
        
    def to_grid(self, bounds, grid):
        
        pixel_x  , pixel_y  = pos_to_index(self.position, bounds, grid)
        
        self.grid = np.zeros(grid.shape)
        self.grid[pixel_y, pixel_x] = 1
        
        
    def evaluate(self, t):
        return self.amplitude * (1 - np.exp(-(t - self.delay) * self.f_max )) * np.sin( 2 * np.pi * self.f_max * (t - self.delay) )
    
    
    def evaluate_on_grid(self, t):
        return self.grid * self.evaluate(t)



class PhasedArray:
    
    def __init__(self, amplitude=1, f_max=1e9, delay=0, position=(0,0)):
        
        self.delay     = delay
        self.amplitude = amplitude
        self.f_max     = f_max
        self.position  = np.array(position)
        
        self.elements = 15
    
    
    def to_grid(self, bounds, grid):
        
        self.grid = np.zeros(grid.shape)
        
        lamb = 3e8 / self.f_max
        
        max_dist = lamb / 4 * self.elements
        
        for d in np.linspace(0, max_dist, self.elements):
            
            d -= max_dist/2
            
            x, y = self.position
            
            if pos_in_bounds((x, y+d), bounds):
                
                pixel_x, pixel_y  = pos_to_index((x, y+d), bounds, grid)
                
                self.grid[pixel_y, pixel_x] = 1
                
                
    def evaluate(self, t):
        return self.amplitude / self.elements * (1 - np.exp(-(t - self.delay) * self.f_max )) * np.sin( 2 * np.pi * self.f_max * (t - self.delay))
                
    
    def evaluate_on_grid(self, t):
        return self.grid * self.evaluate(t)




class PlanarGaussianSource:
    
    def __init__(self, amplitude=1, f_max=1e9, delay=0, position=(0,0)):
        
        self.delay     = delay
        self.amplitude = amplitude
        self.f_max     = f_max
        self.position  = np.array(position)
        
        
    def to_grid(self, bounds, grid):
        
        res_y, res_x = grid.shape
        
        x_max, x_min, y_max, y_min = bounds
        
        x1, y1 = pos_to_index(self.position  , bounds, grid)
        x2, y2 = pos_to_index(((x_max-x_min)/2, (y_max-y_min)/2 ), bounds, grid)
        
        dydx = (y2-y1)/(x2-x1)
        
        y1_new = y1 + dydx * (0 - x1)
        y2_new = y1 + dydx * (res_x - x1)
        
        pixel_p1 = 0    , y1_new
        pixel_p2 = res_x, y2_new
        
        image = Image.new('1', grid.T.shape)
        draw  = ImageDraw.Draw(image)
        
        draw.line([pixel_p1, pixel_p2], fill=1, width=0)
        
        self.grid = np.array(image)
        
        
    def evaluate(self, t):
        
        tau = 0.5 / self.f_max
        
        val =  np.exp( -( t - 4 * tau - self.delay)**2/tau**2)
        
        return self.amplitude * val
    
    
    def evaluate_on_grid(self, t):
        return self.grid * self.evaluate(t)
    
    
    
    
    
    
    
    
    
    
    
    
