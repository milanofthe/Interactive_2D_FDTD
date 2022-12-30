import numpy as np

from PIL import Image, ImageDraw, ImageOps


from skimage import measure

import matplotlib.pyplot as plt

import pygame

from misc import (
    pos_to_index,
    pos_in_bounds,
    sim_to_screen
    )

# MATERIAL CLASSES ====================================================
        
class RectMaterial:
    
    def __init__(self, eps_r=1, mu_r=1, width=0.1, height=0.1, position=(0,0)):
        
        self.mu_r     = mu_r
        self.eps_r    = eps_r
        self.position = np.array(position)
        self.width    = width
        self.height   = height
        
        self.bitmap = None
        
        
    def to_grid(self, bounds, grid):
        
        """
        rasterize material to simulation grid
        
        """
        
        x, y = self.position
        
        pixel_point_1 = pos_to_index(self.position, bounds, grid)
        pixel_point_2 = pos_to_index((x+self.width, y+self.height), bounds, grid)
        
        image = Image.new('1', grid.T.shape)
        draw  = ImageDraw.Draw(image)
        
        draw.rectangle([pixel_point_1,
                        pixel_point_2], 
                       fill=1, 
                       outline=1)
        
        self.bitmap = np.array(image)
    
    
    def draw(self, bounds, font, surf, color=(255,255,255), fill=None):
        
        #get surfece resolution
        res = surf.get_size()
        
        x, y = self.position

        w = self.width
        h = self.height
                    
        x_dsp    , y_dsp      = sim_to_screen(x, y, bounds, res)
        x0_dsp   , y0_dsp     = sim_to_screen(0, 0, bounds, res)
        dsp_width, dsp_height = sim_to_screen(w, h, bounds, res)
        
        dsp_width  -= x0_dsp
        dsp_height -= y0_dsp
                
        
        if fill is not None:
            #draw rectangle
            pygame.draw.rect(surf, fill, (x_dsp, y_dsp, dsp_width, dsp_height), width=0)
        
        #draw rectangle
        pygame.draw.rect(surf, color, (x_dsp, y_dsp, dsp_width, dsp_height), width=2)
                
        #draw text description
        t = f"er={self.eps_r}"
                
        text = font.render(t , True, color)
        surf.blit(pygame.transform.flip(text, False, True), (x_dsp+10, y_dsp+5))
        
        t = f"mur={self.mu_r}"
                
        text = font.render(t , True, color)
        surf.blit(pygame.transform.flip(text, False, True), (x_dsp+10, y_dsp+17))
        
        
    
        
class CircleMaterial:
    
    def __init__(self, eps_r=1, mu_r=1, rad=0.1, position=(0,0)):
        
        self.mu_r     = mu_r
        self.eps_r    = eps_r
        self.position = np.array(position)
        self.rad      = rad
        
        self.bitmap = None
        
        
    def to_grid(self, bounds, grid):
        
        """
        rasterize material to simulation grid
        
        """
        
        pixel_rad, _        = pos_to_index((0, self.rad), bounds, grid)
        pixel_x0 , pixel_y0 = pos_to_index((0,0), bounds, grid)
        pixel_x  , pixel_y  = pos_to_index(self.position, bounds, grid)
        
        pixel_rad -= pixel_x0
        
        image = Image.new('1', grid.T.shape)
        draw  = ImageDraw.Draw(image)
        
        draw.ellipse([(pixel_x - pixel_rad, pixel_y - pixel_rad), 
                      (pixel_x + pixel_rad, pixel_y + pixel_rad)], 
                     fill=1, 
                     outline=1)
        
        self.bitmap = np.array(image)
    
    
    def draw(self, bounds, font, surf, color=(255,255,255), fill=None):
        
        #get surfece resolution
        res = surf.get_size()
        
        x, y = self.position
        
        r = self.rad

        x_dsp, y_dsp = sim_to_screen(x, y, bounds, res)
        r_dsp, _     = sim_to_screen(r, 0, bounds, res)
        x0   , _     = sim_to_screen(0, 0, bounds, res)
        
        r_dsp -= x0
                
        
        if fill is not None:
            #draw circle
            pygame.draw.circle(surf, fill, (x_dsp, y_dsp), r_dsp, width=0)
            
        #draw circle
        pygame.draw.circle(surf, color, (x_dsp, y_dsp), r_dsp, width=2)
                
        #draw text description
        t = f"er={self.eps_r}"
                
        text = font.render(t , True, color)
        surf.blit(pygame.transform.flip(text, False, True), (x_dsp-10, y_dsp-5))
        
        
        t = f"mur={self.mu_r}"
                
        text = font.render(t , True, color)
        surf.blit(pygame.transform.flip(text, False, True), (x_dsp-10, y_dsp+7))
        
        
        
        
        
class PolygonMaterial:
    
    def __init__(self, eps_r=1, mu_r=1, points=[(0, 0), (0.1 ,0), (0, 0.1)]):
        
        self.mu_r     = mu_r
        self.eps_r    = eps_r
        self.points   = points
        
        self.bitmap = None
    
    
    def to_grid(self, bounds, grid):
        
        """
        rasterize material to simulation grid
        
        """
        
        pixel_points = [ pos_to_index(pt, bounds, grid) for pt in self.points ]
        
        image = Image.new('1', grid.T.shape)
        draw  = ImageDraw.Draw(image)
        draw.polygon(xy=pixel_points, fill=1, outline=1) 
        
        self.bitmap = np.array(image)

    
    def draw(self, bounds, font, surf, color=(255,255,255), fill=None):
        
        #get surfece resolution
        res = surf.get_size()
        
        #transform points to screen
        pts = np.array([sim_to_screen(*pt, bounds, res) for pt in self.points])
        
        
        pts_x = pts[:,0]
        pts_y = pts[:,1]
        
        n_pts = len(pts)
        
        #compute center
        center_x, center_y = np.sum(pts_x)/n_pts, np.sum(pts_y)/n_pts
        
        if fill is not None:
            #draw polygon
            pygame.draw.polygon(surf, fill, pts, width=0)
                
        #draw polygon
        pygame.draw.polygon(surf, color, pts, width=2)
        
        #draw text description
        t = f"er={self.eps_r}"
                
        text = font.render(t , True, color)
        surf.blit(pygame.transform.flip(text, False, True), (center_x-10, center_y-5))
        
        t = f"mur={self.mu_r}"
                
        text = font.render(t , True, color)
        surf.blit(pygame.transform.flip(text, False, True), (center_x-10, center_y+7))
    
    
    
class ImgMaterial:
    
    def __init__(self, eps_r=1, mu_r=1, path=""):
        
        self.mu_r     = mu_r
        self.eps_r    = eps_r
        self.path     = path
        self.bitmap   = None
        
        
    def to_grid(self, bounds, grid):
        
        """
        rasterize material to simulation grid
        
        """
        
        image = Image.open(self.path).convert('L')
        
        image = ImageOps.pad(image, grid.shape, color=255)
            
        
        btm = (255 - np.array(image)) / 255
        
        self.bitmap = np.where(btm > 0.1, 1, 0).T
        
        #contours for drawing
        self.contours = measure.find_contours(self.bitmap, 0)
        
    
    
    def draw(self, bounds, font, surf, color=(255,255,255), fill=None):  
        
        #get surfece resolution
        res_x, res_y = surf.get_size()
        rx   , ry    = self.bitmap.shape
        
        scale = np.array([res_x / rx, res_y / ry])
        
        for pts in self.contours:
            
            scale_pts = pts[:]*scale

            if fill is not None:
                #draw polygon
                pygame.draw.polygon(surf, fill, scale_pts, width=0)
                    
            #draw polygon
            pygame.draw.polygon(surf, color, scale_pts, width=2)
        
        
        
    
if __name__ == '__main__':
    
    
    # M1 = CircleMaterial(position=(0.4,0.3), rad=0.25, eps_r=3)
    
    M2 = ImgMaterial(eps_r=3, path=r"C:\Users\Milan\OneDrive\Eigene Projekte\Interactive 2D FDTD\assets\coupler02.png")
    
    # M2 = RectMaterial(position=(0.45,0.4), width = 0.4, height=0.4, eps_r=7)
    
    # M3 = PolygonMaterial(eps_r=4.5, points=[(0.29, 0.74), (0.62, 0.57), (0.99, 0.54), (0.66, 0.71)]) 
    
    
    
    bounds = 1, 0, 1 ,0
    
    n = 200
    
    
    G = np.zeros((200,250))
    
    
    # G1 = M1.to_grid(bounds, G)
    M2.to_grid(bounds, G)
    
    G2 = M2.bitmap
    # G3 = M3.to_grid(bounds, G)
    
    
    
    # G = np.where(G1>1, G1, G)
    G = np.where(G2>0, G2, G)
    # G = np.where(G3>1, G3, G)
    
    contours = measure.find_contours(G, 0)

    # Display the image and plot all contours found
    fig, ax = plt.subplots()
    ax.imshow(G)
    
    for contour in contours:
        print(contour.shape)
        ax.plot(contour[:, 1], contour[:, 0], linewidth=2, color="white")
    
    # plt.imshow(G)
    
    
    
    
    
    
    
    
    
    
    
    
    

