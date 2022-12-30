from misc import ( 
    timer, 
    sim_to_screen, 
    screen_to_sim,
    pos_in_bounds
    )

import pygame

# UI ELEMENTS ===============================================================

class Button:
    
    def __init__(self, text, bounds, is_selected=False, img=None, func=None, args=None, surf=None, font=None):
        
        self.text        = text
        self.bounds      = bounds
        self.is_selected = is_selected
        
        self.has_updated = False
        
        #surface
        self.surf = surf
        self.font = font
        
        #image
        self.img = img
        
        #colors for display
        self.color1 = (0,0,0)
        self.color2 = (55,55,55)
        self.color3 = (255,255,255)
        
        
        self.color      = self.color3 if self.is_selected else self.color1
        self.text_color = self.color1 if self.is_selected else self.color3
        
        #func is function object
        self.func = func
        self.args = args
    
    def interact(self, position=(0,0), pressed=False):
        
        self.has_updated = False
        
        if self.is_selected:
            return
        
        #hovering over button
        if pos_in_bounds(position, self.bounds) :
            
            if pressed:
                
                #call function with args
                if self.func is not None:
                    self.func(*self.args)
                    
                self.select()
                
                self.color = self.color3
                self.text_color = self.color1
                self.has_updated = True
            
            else:
                self.color = self.color2
                self.text_color = self.color3
            
        else:
            self.color = self.color1
            self.text_color = self.color3
    
    
    def draw(self):
        
        x_max, x_min, y_max, y_min = self.bounds
            
        width  = x_max - x_min
        height = y_max - y_min
            
        #draw filled rectangle
        pygame.draw.rect(self.surf, self.color, (x_min,y_min, width, height), width=0)
            
        #draw bounding rectangle
        pygame.draw.rect(self.surf, (255,255,255), (x_min, y_min, width, height), width=2)
            
        if self.img is None:
            #draw text
            text = self.font.render(self.text , True, self.text_color)
            self.surf.blit(text, (x_min + 10 , y_min + height/6 ))
        
                
    def select(self):
        self.is_selected = True
        
        
    def deselect(self):
        self.is_selected = False
    
    
    
    
class ButtonArray:
    
    def __init__(self, buttons=[]):
        
        self.buttons = buttons
        
        self.currently_selected = None
        
        self.has_updated = False
        
        for btn in self.buttons:
            if btn.is_selected:
                self.currently_selected = btn
        
    
    def interact(self, position=(0,0), pressed=False):
        
        self.has_updated = False
        
        #handle button interactions
        for btn in self.buttons:
            btn.interact(position, pressed)
            if btn.has_updated:
                self.currently_selected = btn
                self.has_updated = True
                
        #resolve selections
        for btn in self.buttons:
            if btn != self.currently_selected:
                btn.deselect()
                
                
    def draw(self):
        
        #draw every button
        for btn in self.buttons:
            btn.draw()
    

## KEYBOARD INPUT ==================================================================

class KeyboardInput:
        
    def __init__(self, text, value, bounds, surf=None, font=None):
        
        self.text        = text
        self.bounds      = bounds
        
        #surface
        self.surf = surf
        self.font = font
        
        self.is_selected = False
        
        #colors for display
        self.color1 = (0,0,0)
        self.color2 = (55,55,55)
        self.color3 = (255,255,255)
        
        self.color      = self.color1
        self.text_color = self.color3
        
        #value
        self.value_str = str(value)
        self.value = value
        
        
        self.key_event_dict = {pygame.K_1      : "1",
                               pygame.K_2      : "2",
                               pygame.K_3      : "3",
                               pygame.K_4      : "4",
                               pygame.K_5      : "5",
                               pygame.K_6      : "6",
                               pygame.K_7      : "7",
                               pygame.K_8      : "8",
                               pygame.K_9      : "9",
                               pygame.K_0      : "0",
                               pygame.K_e      : "e",
                               pygame.K_PERIOD : ".",
                               pygame.K_MINUS  : "-",
                               }
        
        
        
    def select(self):
        self.is_selected = True
        
    def deselect(self):
        self.is_selected = False
        
        
    def interact(self, position=(0,0), pressed=False, events=[]):
        
        
        if self.is_selected:
            self.color      = self.color3
            self.text_color = self.color1
            
            if pressed:
                self.deselect()
                self.color      = self.color1
                self.text_color = self.color3
                
            for event in events:
            
                if event in self.key_event_dict:
                    self.value_str += self.key_event_dict[event]
                
                if event == pygame.K_BACKSPACE and len(self.value_str) > 0:
                    self.value_str = self.value_str[:-1]
                    
        else:
            self.color      = self.color1
            self.text_color = self.color3
                    
        if pos_in_bounds(position, self.bounds):
            
            if pressed:
                self.select()
                self.color      = self.color3
                self.text_color = self.color1
                
            if not self.is_selected:
                self.color      = self.color2
                self.text_color = self.color3
                
                    
                    
    def get_value(self):
        
        if len(self.value_str) > 0:
            try:
                self.value = eval(self.value_str)
                return self.value
            except:
                return self.value
                
        else:
            return self.value
        
        
    def draw(self):
        
        x_max, x_min, y_max, y_min = self.bounds
            
        width  = x_max - x_min
        height = y_max - y_min
            
        #draw filled rectangle
        pygame.draw.rect(self.surf, self.color, (x_min,y_min, width, height), width=0)
            
        #draw bounding rectangle
        pygame.draw.rect(self.surf, (255,255,255), (x_min, y_min, width, height), width=2)
        
        #draw text
        text = self.font.render(self.text + self.value_str , True, self.text_color)
        self.surf.blit(text, (x_min + 10 , y_min + height/6 ))
        
        

class Menu:
    
    def __init__(self, interactables, descriptions, title, background):
        
        self.interactables = interactables
        self.descriptions  = descriptions
        self.title         = title
        self.background    = background
        
    





class Grid:
    
    
    """
    draw grid for cartesian coordinate system
    
    """
    
    def __init__(self, sim_bounds, nx, ny, surf, font):
        
        self.sim_bounds = sim_bounds
        
        self.nx = nx
        self.ny = ny
        
        self.surf = surf
        self.font = font
        
        self.color = (55, 55, 55)
        
        
    def draw(self):
        
        x_max, x_min, y_max, y_min = self.sim_bounds
        
        res_x, res_y = self.surf.get_size()
        
        for n in range(self.nx):
            
            x        = (n + 1) / (self.nx + 1) * res_x
            x_sim ,_ = screen_to_sim(x, 0, self.sim_bounds, (res_x, res_y))

            #draw line
            pygame.draw.line(self.surf, self.color, (x, 0), (x, res_y), 2)
            
            #draw text
            text = self.font.render(f"{round(x_sim, 3)}" , True, self.color)
            self.surf.blit(text, (x+5, 0))
            
            
            
        for n in range(self.ny):
            
            y        = (n + 1) / (self.ny + 1) * res_y
            _, y_sim = screen_to_sim(0, y, self.sim_bounds, (res_x, res_y))
            
            
            #draw line
            pygame.draw.line(self.surf, self.color, (0, y), (res_x, y), 2)
            
            #draw text
            text = self.font.render(f"{round(y_sim, 3)}" , True, self.color)
            self.surf.blit(text, (5, y))
            
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        