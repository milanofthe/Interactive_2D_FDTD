#imports
import numpy as np
import pygame

# import skimage as sk

# import matplotlib.pyplot as plt

# import sys

import os

#modules
from sources import (
    GaussianSource, 
    SinSource,
    PhasedArray
    )

from materials import (
    RectMaterial,
    CircleMaterial,
    PolygonMaterial,
    ImgMaterial
    )

from misc import ( 
    timer, 
    sim_to_screen, 
    screen_to_sim
    )

from ui import (
    Button,
    ButtonArray,
    KeyboardInput,
    Grid
    )



class Interface:
    
    
    """
    class for handling interaction with and display of the FDTD simulation
    
    """
    
    def __init__(self, 
                 fps=30, 
                 steps=2, 
                 res=(750,750), 
                 color_mode=True, 
                 Sim = None
                 ):
        
        #some display properties
        self.fps = fps
        self.frametime = 1/self.fps
        
        #resolution of window
        self.res = res
        
        #substeps for simulation
        self.steps = steps
        
        #simulation time
        self.sim_time = 0
        
        #color mode
        self.color_mode = color_mode
        
        #setup simulation
        self.Sim = Sim
        
        #min and max resolution of main window
        self.res_min = (460 , 460)
        self.res_max = (1920, 1080)
        
        
        
        
        
        
        
    ## EDITOR =================================================================
    
    def _set_material_selection(self, mat):
        self.material_selecion = mat
        self.new_material_coords = []
        
        
        
    def _setup_editor(self):
        
        #size of surface
        res_x, res_y = self.Editor_surf.get_size()
        
        #coords for adding new material
        self.new_material_coords = []
        
        
        #material dielectric constant for editor
        self.new_material_er = 3
        self.new_material_mur = 1
        
        
        #create grid
        self.editor_grid = Grid(sim_bounds=self.Sim.bounds, nx=res_x//120, ny=res_y//120, surf=self.Editor_surf, font=self.smallfont)
        
        
        #button location
        btn_h = 40
        btn_w = 60
        
        btn_sep = btn_w + 20
        
        #button sizes
        btn_x = 20
        btn_y = res_y - 20 - btn_h
        
        
        self.editor_buttons_material_selection = ButtonArray([Button(text="rect", 
                                                                     bounds=[btn_x + btn_w, btn_x, btn_y + btn_h, btn_y], 
                                                                     is_selected=(self.material_selecion==RectMaterial),
                                                                     func=self._set_material_selection,
                                                                     args=[RectMaterial],
                                                                     surf=self.Editor_surf,
                                                                     font=self.font
                                                                     ),
                                                              Button(text="circ", 
                                                                     bounds=[btn_x + btn_w + btn_sep, btn_x + btn_sep, btn_y + btn_h, btn_y], 
                                                                     is_selected=(self.material_selecion==CircleMaterial),
                                                                     func=self._set_material_selection,
                                                                     args=[CircleMaterial],
                                                                     surf=self.Editor_surf,
                                                                     font=self.font
                                                                     ),
                                                              Button(text="poly", 
                                                                     bounds=[btn_x + btn_w + 2*btn_sep, btn_x + 2*btn_sep, btn_y + btn_h, btn_y], 
                                                                     is_selected=(self.material_selecion==PolygonMaterial),
                                                                     func=self._set_material_selection,
                                                                     args=[PolygonMaterial],
                                                                     surf=self.Editor_surf,
                                                                     font=self.font
                                                                     )
                                                              ])
        
        #keyboard input for material parameters
        
        btn_x = res_x - 20 - 2*btn_w
            
        self.editor_keyboard_input_er = KeyboardInput(text="er = ", 
                                                      value=self.new_material_er,
                                                      bounds=[btn_x + 2*btn_w, btn_x, btn_y + btn_h, btn_y], 
                                                      surf=self.Editor_surf, 
                                                      font=self.font
                                                      )
        
        btn_y = res_y - 80 - btn_h
        
        self.editor_keyboard_input_mur = KeyboardInput(text="mur = ", 
                                                      value=self.new_material_mur,
                                                      bounds=[btn_x + 2*btn_w, btn_x, btn_y + btn_h, btn_y ], 
                                                      surf=self.Editor_surf, 
                                                      font=self.font
                                                      )
        
        
        
        
    
    
    def _place_material(self):
        
        
        if self.leftclick:
            point = screen_to_sim(*self.mouse_position, self.Sim.bounds, self.res)
            self.new_material_coords.append( point )
            
            
        if len(self.new_material_coords) == 2:
            
            (x1, y1), (x2, y2) = self.new_material_coords
            
            if self.material_selecion == RectMaterial:
                
                x = min(x1, x2)
                y = min(y1, y2)
                
                w = abs(x1 - x2)
                h = abs(y1 - y2)
                
                self.Sim.add_material(RectMaterial(eps_r=self.new_material_er, mu_r=self.new_material_mur, width=w, height=h, position=(x,y)))
                
                self.new_material_coords = []
            
            elif self.material_selecion == CircleMaterial:
                
                r = np.sqrt((x1-x2)**2 + (y1-y2)**2)
                
                self.Sim.add_material(CircleMaterial(eps_r=self.new_material_er, mu_r=self.new_material_mur, rad=r, position=(x1,y1)))
                
                self.new_material_coords = []
                
            else:
                pass
            
        if len(self.new_material_coords) >= 5:
            
            if self.material_selecion == PolygonMaterial:
                
                self.Sim.add_material(PolygonMaterial(eps_r=self.new_material_er, mu_r=self.new_material_mur, points=self.new_material_coords))
            
                self.new_material_coords = []
        
    def _remove_last_material(self):
        
        materials = self.Sim.materials[:-1]
                        
        self.Sim._reset_materials()
        self.new_material_coords = []
                            
        for m in materials:
            self.Sim.add_material(m)
        
    
    def _update_editor(self):
        
        
        #update keyboard input
        self.editor_keyboard_input_er.interact(position=self.mouse_position, pressed=self.leftclick, events=self.keyboard_input_buffer )
        self.new_material_er = self.editor_keyboard_input_er.get_value()
        
        self.editor_keyboard_input_mur.interact(position=self.mouse_position, pressed=self.leftclick, events=self.keyboard_input_buffer )
        self.new_material_mur = self.editor_keyboard_input_mur.get_value()
        
        #update buttons for shape selection
        self.editor_buttons_material_selection.interact(position=self.mouse_position, pressed=self.leftclick)
        
        
        if (self.leftclick 
            and not self.editor_buttons_material_selection.has_updated 
            and not self.editor_keyboard_input_er.is_selected
            and not self.editor_keyboard_input_mur.is_selected):
            
            #place material
            self._place_material()
        
        
        if (pygame.K_BACKSPACE in self.keyboard_input_buffer 
            and not self.editor_keyboard_input_er.is_selected 
            and not self.editor_keyboard_input_mur.is_selected):
            
            #remove last material
            self._remove_last_material()
        
        
        
    def _draw_editor_new_material_coords(self):
        
        res = self.Editor_surf.get_size()
        
        
        cds = pygame.Surface(res)
        cds.set_colorkey((0,0,0))
        
        for (x, y) in self.new_material_coords:
            
            x_dsp, y_dsp = sim_to_screen(x, y, self.Sim.bounds, res)
            
            rad = 6
            d = rad / np.sqrt(2)
            
            pygame.draw.line(cds, (155,155,155), (x_dsp+d, y_dsp+d), (x_dsp-d, y_dsp-d), width=3)
            pygame.draw.line(cds, (155,155,155), (x_dsp+d, y_dsp-d), (x_dsp-d, y_dsp+d), width=3)
            
            
        self.Editor_surf.blit(pygame.transform.flip(cds, False, True), (0,0))
        
        
    def _draw_editor_grid(self):
        self.editor_grid.draw()
        
        
        
    def _draw_editor_materials(self):
        
        res = self.Editor_surf.get_size()
        
        #temporary surface
        mats = pygame.Surface(res)
        mats.set_colorkey((0,0,0))
        
        for M in self.Sim.materials:
            M.draw(bounds=self.Sim.bounds, surf=mats, font=self.smallfont, color=(155,155,155), fill=(55,55,55))
            
        self.Editor_surf.blit(pygame.transform.flip(mats, False, True), (0,0))
    
    
    def _draw_editor_ui(self):
        
        
        #text
        text = self.largefont.render("E D I T - M O D E" , True, (255,255,255))
        self.Editor_surf.blit(text, (20, 5))
        
        #buttons
        self.editor_buttons_material_selection.draw()
        
        #keyboard input
        self.editor_keyboard_input_er.draw()
        self.editor_keyboard_input_mur.draw()
        
        
    
    def _render_editor(self):
        
        
        self.Editor_surf.fill((0,0,0))
        
        self._draw_editor_grid()
        
        self._draw_editor_materials()
        
        self._draw_editor_new_material_coords()
        
        self._draw_editor_ui()
        
        #blit to screen
        self.SCREEN.blit(self.Editor_surf, (0,0))
        
        #update screen
        pygame.display.update()
    
    
    def _run_editor(self):
        
        self.sub_active = True
        
        while self.sub_active:
            
            #update general state
            self._update_state()
            
            #check static key events
            for event in pygame.event.get():
                
                if event.type == pygame.QUIT:
                    self.sub_active = False
                    self.active = False
                    
                if event.type == pygame.KEYDOWN:
                    
                    self.keyboard_input_buffer.append(event.key)
                    
                    if event.key == pygame.K_ESCAPE:
                        self.sub_active = False
                        self.what_to_run_next = "mainmenu"
                        
                    if event.key == pygame.K_r:
                        self.Sim._reset_materials()
                        self.new_material_coords = []
                        
                        
                #resize event
                if event.type == pygame.VIDEORESIZE:
                    self._resize(event.w, event.h)
                    
                        
                #leftclick
                if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
                    self.leftclick = True

            #update menu state
            self._update_editor()
            
            #draw menu
            self._render_editor()
            
            #limit fps
            self.clock.tick(self.fps)
        
        
        
        
        
        
    ## MAIN MENU ==============================================================
        
        
    # menu logic --------------------------------------------------------------
    
    
    def _set_what_to_run_next(self, what):
        self.what_to_run_next = what
        #to break current subloop
        self.sub_active = False
    
    def _setup_main_menu(self):
        
        
        #size of surface
        res_x, res_y = self.MainMenu_surf.get_size()
        
        
        
        #load image assets
        if self.color_mode:
            img_fdtd = pygame.image.load(os.path.join("assets", "main_menu.png"))
        
        else:
            img_fdtd = pygame.image.load(os.path.join("assets", "main_menu_gr.png"))
        
        img_wh = res_x - 40
        
        self.img_fdtd = pygame.transform.smoothscale(img_fdtd, (img_wh, img_wh))
        self.img_fdtd = self.img_fdtd.subsurface(0, 0, img_wh,  min(img_wh, res_y - 60))
        
        
        
        n_btns = 2
        
        
        
        #button sizes
        btn_x = 40
        
        btn_h = 40
        btn_w = (res_x - 40 - (n_btns + 1) * 20)//n_btns
        
        btn_sep = btn_w + 20
        
        
        
        btn_y = res_y - btn_h - 40 - 60
        
        self.menu_button_sim = Button(text="Simulator", 
                                      bounds=[btn_x + btn_w, btn_x, btn_y + btn_h, btn_y], 
                                      is_selected=False,
                                      func=self._set_what_to_run_next,
                                      args=["sim"],
                                      surf=self.MainMenu_surf,
                                      font=self.font
                                      )
        
        
        self.menu_button_edt = Button(text="Material Editor", 
                                      bounds=[btn_x + btn_w + btn_sep, btn_x + btn_sep, btn_y + btn_h, btn_y], 
                                      is_selected=False,
                                      func=self._set_what_to_run_next,
                                      args=["edit"],
                                      surf=self.MainMenu_surf,
                                      font=self.font
                                      )
        
        
        btn_y = res_y - btn_h - 40
        
        self.menu_button_sim_settings = Button(text="Simulation Settings", 
                                               bounds=[btn_x + btn_w, btn_x, btn_y + btn_h, btn_y], 
                                               is_selected=False,
                                               func=self._set_what_to_run_next,
                                               args=["sim_settings"],
                                               surf=self.MainMenu_surf,
                                               font=self.font
                                               )

        
        self.menu_button_general_settings = Button(text="General Settings", 
                                                   bounds=[btn_x + btn_w + btn_sep, btn_x + btn_sep, btn_y + btn_h, btn_y], 
                                                   is_selected=False,
                                                   func=self._set_what_to_run_next,
                                                   args=["general_settings"],
                                                   surf=self.MainMenu_surf,
                                                   font=self.font
                                                   )
        
        
    def _update_main_menu(self):
        
        self.menu_button_sim.interact(position=self.mouse_position, pressed=self.leftclick)
        self.menu_button_sim.deselect()
        
        self.menu_button_edt.interact(position=self.mouse_position, pressed=self.leftclick)
        self.menu_button_edt.deselect()
        
        self.menu_button_sim_settings.interact(position=self.mouse_position, pressed=self.leftclick)
        self.menu_button_sim_settings.deselect()
        
        self.menu_button_general_settings.interact(position=self.mouse_position, pressed=self.leftclick)
        self.menu_button_general_settings.deselect()
        
        
    def _draw_main_menu_buttons(self):
        
        self.menu_button_sim.draw()
        
        self.menu_button_edt.draw()
        
        self.menu_button_sim_settings.draw()
        
        self.menu_button_general_settings.draw()
        
        
    def _draw_main_menu_exterior(self):
        
        #size of surface
        res_x, res_y = self.MainMenu_surf.get_size()
        
        #image
        self.MainMenu_surf.blit(self.img_fdtd, (20, 40))
        
        
        #draw bounding rectangle
        pygame.draw.rect(self.MainMenu_surf, (255,255,255), (20,40, res_x-40, res_y-60), width=2)
        
        #text
        text = self.largefont.render("2 D - F D T D" , True, (255,255,255))
        self.MainMenu_surf.blit(text, (20, 5))
        
        
        
        
    def _render_main_menu(self):
        
        
        #menu surface
        self.MainMenu_surf.fill((0,0,0))
        
        #menu background
        self._draw_main_menu_exterior()
        
        #menu buttons
        self._draw_main_menu_buttons()
        
        #blit to screen
        self.SCREEN.blit(self.MainMenu_surf, (0,0))
        
        #update screen
        pygame.display.update()
        
        
        
        
    def _run_main_menu(self):
        
        self.sub_active = True
        
        while self.sub_active:
            
            #update general state
            self._update_state()
            
            #check static key events
            for event in pygame.event.get():
                
                if event.type == pygame.QUIT:
                    self.sub_active = False
                    self.active = False
                    
                if event.type == pygame.KEYDOWN:
                    
                    if event.key == pygame.K_ESCAPE:
                        self.sub_active = False
                        self.active = False
                        
                        
                #resize event
                if event.type == pygame.VIDEORESIZE:
                    self._resize(event.w, event.h)
                    
                        
                #leftclick
                if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
                    self.leftclick = True

            #update menu state
            self._update_main_menu()
            
            #draw menu
            self._render_main_menu()
            
            #limit fps
            self.clock.tick(self.fps)
        
    
    
    ## SIM SETTINGS MENU ======================================================
    
    
    def _set_boundary_conditions_TB(self, bc):
        self.Sim.boundary_conditions["top"]    = bc
        self.Sim.boundary_conditions["bottom"] = bc
        
        # self.Sim._setup_boundary_conditions()
        self.Sim.reset()
    
    
    def _set_boundary_conditions_LR(self, bc):
        self.Sim.boundary_conditions["left"]  = bc
        self.Sim.boundary_conditions["right"] = bc
        
        # self.Sim._setup_boundary_conditions()
        self.Sim.reset()
        
    def _set_source_selection(self, Source):
        self.source_selection = Source
    
    
    #setup --------------------------------------------------------------------
    
    def _setup_settings_sim(self):
        
        #size of surface
        res_x, res_y = self.SettingsSim_surf.get_size()
        
        
        
        
    
        #button sizes
        btn_x = 40
        btn_y = 80
        btn_h = 40
        
        #top and bottom boundary conditions
        
        n_btns = 3
        
        btn_w = (res_x - 40 - (n_btns + 1) * 20)//n_btns
        btn_sep = btn_w + 20
        
        
        self.options_buttons_BC_TB = ButtonArray([Button(text="periodic", 
                                                      bounds=[btn_x + btn_w, btn_x, btn_y + btn_h, btn_y], 
                                                      is_selected=(self.Sim.boundary_conditions["top"] == "periodic"),
                                                      func=self._set_boundary_conditions_TB,
                                                      args=["periodic"],
                                                      surf=self.SettingsSim_surf,
                                                      font=self.font
                                                      ),
                                              Button(text="dirichlet", 
                                                      bounds=[btn_x + btn_w + btn_sep, btn_x + btn_sep, btn_y + btn_h, btn_y], 
                                                      is_selected=(self.Sim.boundary_conditions["top"] == "dirichlet"),
                                                      func=self._set_boundary_conditions_TB,
                                                      args=["dirichlet"],
                                                      surf=self.SettingsSim_surf,
                                                      font=self.font
                                                      ),
                                              Button(text="pml", 
                                                      bounds=[btn_x + btn_w + 2*btn_sep, btn_x + 2*btn_sep, btn_y + btn_h, btn_y], 
                                                      is_selected=(self.Sim.boundary_conditions["top"] == "pml"),
                                                      func=self._set_boundary_conditions_TB,
                                                      args=["pml"],
                                                      surf=self.SettingsSim_surf,
                                                      font=self.font
                                                     ),
                                              ])
        
        
        
        #left and right boundary conditions
        btn_y = 160
        
        n_btns = 3
        
        btn_w = (res_x - 40 - (n_btns + 1) * 20)//n_btns
        btn_sep = btn_w + 20
        
        self.options_buttons_BC_LR = ButtonArray([Button(text="periodic", 
                                                      bounds=[btn_x + btn_w, btn_x, btn_y + btn_h, btn_y], 
                                                      is_selected=(self.Sim.boundary_conditions["left"] == "periodic"),
                                                      func=self._set_boundary_conditions_LR,
                                                      args=["periodic"],
                                                      surf=self.SettingsSim_surf,
                                                      font=self.font
                                                      ),
                                               Button(text="dirichlet", 
                                                      bounds=[btn_x + btn_w + btn_sep, btn_x + btn_sep, btn_y + btn_h, btn_y], 
                                                      is_selected=(self.Sim.boundary_conditions["left"] == "dirichlet"),
                                                      func=self._set_boundary_conditions_LR,
                                                      args=["dirichlet"],
                                                      surf=self.SettingsSim_surf,
                                                      font=self.font
                                                      ),
                                               Button(text="pml", 
                                                      bounds=[btn_x + btn_w + 2*btn_sep, btn_x + 2*btn_sep, btn_y + btn_h, btn_y], 
                                                      is_selected=(self.Sim.boundary_conditions["left"] == "pml"),
                                                      func=self._set_boundary_conditions_LR,
                                                      args=["pml"],
                                                      surf=self.SettingsSim_surf,
                                                      font=self.font
                                                      ),
                                               ])
        
        #self.img_dirichlet  = pygame.transform.scale(img_dirichlet, (btn_w, btn_w))
        #self.img_periodic   = pygame.transform.scale(img_periodic, (btn_w, btn_w))
        #self.img_pml        = pygame.transform.scale(img_pml, (btn_w, btn_w))
    
        #source selection
        btn_y = 240
        
        n_btns = 2
        
        btn_w = (res_x - 40 - (n_btns + 1) * 20)//n_btns
        btn_sep = btn_w + 20
        
        self.options_buttons_src_selection = ButtonArray([Button(text="Gaussian", 
                                                              bounds=[btn_x + btn_w, btn_x, btn_y + btn_h, btn_y], 
                                                              is_selected=self.source_selection==GaussianSource,
                                                              func=self._set_source_selection,
                                                              args=[GaussianSource],
                                                              surf=self.SettingsSim_surf,
                                                              font=self.font
                                                              ),
                                                       Button(text="Sinusoidal", 
                                                              bounds=[btn_x + btn_w + btn_sep, btn_x + btn_sep, btn_y + btn_h, btn_y], 
                                                              is_selected=self.source_selection==SinSource,
                                                              func=self._set_source_selection,
                                                              args=[SinSource],
                                                              surf=self.SettingsSim_surf,
                                                              font=self.font
                                                              ),
                                                        # Button(text="Array", 
                                                        #        bounds=[btn_x + btn_w + 2*btn_sep, btn_x + 2*btn_sep, btn_y + btn_h, btn_y], 
                                                        #        is_selected=self.source_selection==PhasedArray,
                                                        #        func=self._set_source_selection,
                                                        #        args=[PhasedArray],
                                                        #        surf=self.SettingsSim_surf,
                                                        #        font=self.font
                                                        #        )
                                                       ])
        
        #frequency selection
        
        btn_y = 320
        
        n_btns = 1
        
        btn_w = (res_x - 40 - (n_btns + 1) * 20)//n_btns
        btn_sep = btn_w + 20
        
        self.options_keyboard_input_freq = KeyboardInput(text="f = ", 
                                                         value=self.max_sim_frequency,
                                                         bounds=[btn_x + btn_w, btn_x, btn_y + btn_h, btn_y], 
                                                         surf=self.SettingsSim_surf, 
                                                         font=self.font
                                                         )
    
    
    
    
    
    
    
    # updating ------------------------------------------------------------------
    
    def _update_settings_sim(self):
        
        self.options_buttons_BC_TB.interact(position=self.mouse_position, pressed=self.leftclick)
        
        self.options_buttons_BC_LR.interact(position=self.mouse_position, pressed=self.leftclick)
        
        self.options_buttons_src_selection.interact(position=self.mouse_position, pressed=self.leftclick)
        
        self.options_keyboard_input_freq.interact(position=self.mouse_position, pressed=self.leftclick, events=self.keyboard_input_buffer )
        self.max_sim_frequency = self.options_keyboard_input_freq.get_value()
        
        
    # drawing -------------------------------------------------------------------
    
    def _draw_settings_sim_buttons(self):
        
        self.options_buttons_BC_TB.draw()
        
        self.options_buttons_BC_LR.draw()
        
        self.options_buttons_src_selection.draw()
        
        self.options_keyboard_input_freq.draw()
        
        
    
    def _draw_settings_sim_exterior(self):
        
        #size of surface
        res_x, res_y = self.SettingsSim_surf.get_size()
        
        #draw bounding rectangle
        pygame.draw.rect(self.SettingsSim_surf, (255,255,255), (20,40, res_x-40, res_y-60), width=2)
        
        #text
        text = self.largefont.render("S I M U L A T I O N - S E T T I N G S" , True, (255,255,255))
        self.SettingsSim_surf.blit(text, (20, 5))
        
        #text
        text = self.font.render("TOP and BOTTOM boundary conditions" , True, (255,255,255))
        self.SettingsSim_surf.blit(text, (40, 50))
        
        #text
        text = self.font.render("LEFT and RIGHT boundary conditions" , True, (255,255,255))
        self.SettingsSim_surf.blit(text, (40, 130))
        
        #text
        text = self.font.render("select SOURCE type" , True, (255,255,255))
        self.SettingsSim_surf.blit(text, (40, 210))
        
        #text
        text = self.font.render("set source FREQUENCY in Hz" , True, (255,255,255))
        self.SettingsSim_surf.blit(text, (40, 290))
        
    
    def _render_settings_sim(self):
        
         #menu surface
        self.SettingsSim_surf.fill((0,0,0))
        
        #menu background
        self._draw_settings_sim_exterior()
        
        #menu buttons
        self._draw_settings_sim_buttons()
        
        
        #blit to screen
        self.SCREEN.blit(self.SettingsSim_surf, (0,0))
        
        #update screen
        pygame.display.update()
    
    # menu loop ---------------------------------------------------------------------------
    
    def _run_settings_sim(self):
        
        self.sub_active = True
        
        while self.sub_active:
            
            #update general state
            self._update_state()
            
            #check static key events
            for event in pygame.event.get():
                
                if event.type == pygame.QUIT:
                    self.sub_active = False
                    self.active = False
                    
                if event.type == pygame.KEYDOWN:
                    
                    self.keyboard_input_buffer.append(event.key)
                    
                    if event.key == pygame.K_ESCAPE:
                        self.sub_active = False
                        self.what_to_run_next = "mainmenu"
                        
                        
                #resize event
                if event.type == pygame.VIDEORESIZE:
                    self._resize(event.w, event.h)
                        
                        
                #leftclick
                if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
                    self.leftclick = True
                
            #update menu state
            self._update_settings_sim()
            
            #draw menu
            self._render_settings_sim()
            
            #limit fps
            self.clock.tick(self.fps)
        
        
    ## GENERAL SETTINGS MENU ==================================================
        
    # menu logic --------------------------------------------------------------

    
    
    def _set_target_fps(self, fps):
        self.fps = fps        
            
    def _setup_settings_general(self):
        
        #size of surface
        res_x, res_y = self.SettingsGeneral_surf.get_size()
        
        #button sizes
        btn_x = 40
        btn_y = 80
        btn_h = 40
        
        #top and bottom boundary conditions
        
        n_btns = 3
        
        btn_w = (res_x - 40 - (n_btns + 1) * 20)//n_btns
        btn_sep = btn_w + 20
        
        
        

        #target framerate
        btn_y = 80
        
        n_btns = 3
        
        btn_w = (res_x - 40 - (n_btns + 1) * 20)//n_btns
        btn_sep = btn_w + 20
        
        self.options_buttons_FPS = ButtonArray([Button(text="30", 
                                                    bounds=[btn_x + btn_w, btn_x, btn_y + btn_h, btn_y], 
                                                    is_selected=(self.fps == 30),
                                                    func=self._set_target_fps,
                                                    args=[30],
                                                    surf=self.SettingsGeneral_surf,
                                                    font=self.font
                                             ),
                                             Button(text="60", 
                                                    bounds=[btn_x + btn_w + btn_sep, btn_x + btn_sep, btn_y + btn_h, btn_y], 
                                                    is_selected=(self.fps == 60),
                                                    func=self._set_target_fps,
                                                    args=[60],
                                                    surf=self.SettingsGeneral_surf,
                                                    font=self.font
                                             ),
                                             Button(text="90", 
                                                    bounds=[btn_x + btn_w + 2*btn_sep, btn_x + 2*btn_sep, btn_y + btn_h, btn_y], 
                                                    is_selected=(self.fps == 90),
                                                    func=self._set_target_fps,
                                                    args=[90],
                                                    surf=self.SettingsGeneral_surf,
                                                    font=self.font
                                                    )
                                             ])
        
        

        
        #self.img_gaussian   = pygame.transform.scale(img_gaussian, (btn_w, btn_w))
        #self.img_sinusoidal = pygame.transform.scale(img_sinusoidal, (btn_w, btn_w))
        
        
 
    def _update_settings_general(self):
        
        self.options_buttons_FPS.interact(position=self.mouse_position, pressed=self.leftclick)
    
    # menu drawing ------------------------------------------------------------
    
    def _draw_settings_general_buttons(self):
        
        self.options_buttons_FPS.draw()
        
    
    def _draw_settings_general_exterior(self):
        
        #size of surface
        res_x, res_y = self.SettingsGeneral_surf.get_size()
        
        #draw bounding rectangle
        pygame.draw.rect(self.SettingsGeneral_surf, (255,255,255), (20,40, res_x-40, res_y-60), width=2)
        
        #text
        text = self.largefont.render("G E N E R A L - S E T T I N G S" , True, (255,255,255))
        self.SettingsGeneral_surf.blit(text, (20, 5))
        
        #text
        text = self.font.render("target FPS" , True, (255,255,255))
        self.SettingsGeneral_surf.blit(text, (40, 50))
        
      
        
        
    def _render_settings_general(self):
        
        #menu surface
        self.SettingsGeneral_surf.fill((0,0,0))
        
        #menu background
        self._draw_settings_general_exterior()
        
        #menu buttons
        self._draw_settings_general_buttons()
        
        #blit to screen
        self.SCREEN.blit(self.SettingsGeneral_surf, (0,0))
        
        #update screen
        pygame.display.update()
        
        
        
    # menu loop ---------------------------------------------------------------------------
    
    def _run_settings_general(self):
        
        self.sub_active = True
        
        while self.sub_active:
            
            #update general state
            self._update_state()
            
            #check static key events
            for event in pygame.event.get():
                
                if event.type == pygame.QUIT:
                    self.sub_active = False
                    self.active = False
                    
                if event.type == pygame.KEYDOWN:                    
                    
                    if event.key == pygame.K_ESCAPE:
                        self.sub_active = False
                        self.what_to_run_next = "mainmenu"
                        
                        
                #resize event
                if event.type == pygame.VIDEORESIZE:
                    self._resize(event.w, event.h)
                        
                        
                #leftclick
                if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
                    self.leftclick = True
                
            #update menu state
            self._update_settings_general()
            
            #draw menu
            self._render_settings_general()
            
            #limit fps
            self.clock.tick(self.fps)
            
            
            
            
            
            
    
    ## SIMULATION =================================================================
        
    # sim drawing -----------------------------------------------------------------
    
    @timer
    def _draw_Ez(self):
        
        res = self.Sim_surf.get_size()
        
        Ez = np.flip(self.Sim.get_Ez(), axis=1)
        E0 = 0.05 * self.Sim.E_max
        
        Ez_img = np.zeros((*Ez.shape,3), dtype=np.uint8)
        
        
        
        if self.color_mode:
            Ez_p = abs(np.where(Ez > 0, Ez, 0))
            Ez_m = abs(np.where(Ez < 0, Ez, 0))
            
            Ez_img[:,:,0] = 255 * (1 - np.exp(- Ez_p / E0))
            Ez_img[:,:,2] = 255 * (1 - np.exp(- Ez_m / E0))
        
        else:
            
            Ez_abs = abs(Ez)
            value = 255 * (1 - np.exp(- Ez_abs / E0))
            
            Ez_img[:,:,0] = value
            Ez_img[:,:,1] = value
            Ez_img[:,:,2] = value
            
            
        
        pygame.surfarray.blit_array(self.Ez_surf, Ez_img)
        
        self.Sim_surf.blit( pygame.transform.smoothscale(self.Ez_surf, res), (0, 0))
        
        
    # @timer
    def _draw_materials(self):
        
        #fill surface
        self.Mat_surf.fill((0,0,0))
        
        for M in self.Sim.materials:
            M.draw(bounds=self.Sim.bounds, surf=self.Mat_surf, font=self.smallfont )
            
        self.Sim_surf.blit(pygame.transform.flip(self.Mat_surf, False, True), (0,0))
        
    # @timer
    def _draw_real_simulation_time(self):
                
        t  = self.Sim.total_steps * self.Sim.dt
        
        if t < 1e-9:
            t *= 1e12
            unit = "fs"
            
        elif t < 1e-6:
            t *= 1e9
            unit = "ns"
            
        elif t < 1e-3:
            t *=1e6
            unit = "us"
            
        else:
            t *=1e3
            unit = "ms"
        
        
        #display text
        text = self.font.render(f"simulation time [{unit}] : {round(t,2)}" , True, (255,255,255))
        self.Sim_surf.blit(text, (5, 0))
        
        # text = self.font.render(unit , True, (255,255,255))
        # self.Sim_surf.blit(text, (70,20))
        
        
    # @timer
    def _draw_real_framerate(self):
        
        res_x, res_y = self.res
        
        #display text
        text = self.font.render(f"FPS : {int(1/self.real_frametime)}" , True, (255,255,255))
        self.Sim_surf.blit(text, (res_x - 75, 0))
        
        # text = self.font.render("FPS" , True, (255,255,255))
        # self.Sim_surf.blit(text, (70,0))
        
        
    # @timer
    def _render_sim(self):
        
        #fill black
        self.Sim_surf.fill((255,255,255))
        
        #draw fields
        self._draw_Ez()
        
        #draw materials
        self._draw_materials()
        
        #draw texts
        self._draw_real_framerate()
        self._draw_real_simulation_time()
        
        
        #blit to screen
        self.SCREEN.blit(self.Sim_surf, (0,0))
        
        #update screen
        pygame.display.update()
        
        
    
        
        
        
        
        
    # simulation loop ---------------------------------------------------------------------
        
    def _run_sim(self):
        
        self.sub_active = True
        
        while self.sub_active:
            
            #update general state
            self._update_state()
            
            #check static key events
            for event in pygame.event.get():
                
                if event.type == pygame.QUIT:
                    self.sub_active = False
                    self.active = False
                    
                if event.type == pygame.KEYDOWN:
                    
                    if event.key == pygame.K_ESCAPE:
                        self.sub_active = False
                        self.what_to_run_next = "mainmenu"
                        
                    if event.key == pygame.K_r:
                        self.Sim.reset()
                        
                    if event.key == pygame.K_BACKSPACE:
                        self.Sim.sources = self.Sim.sources[:-1]
                        
                    if event.key == pygame.K_SPACE:
                        self.pause = False if self.pause else True
                        
                        
                #resize event
                if event.type == pygame.VIDEORESIZE:
                    self._resize(event.w, event.h)
                
                
                #leftclick -> add source
                if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
                    self.leftclick = True
                    
                    #mouse position in simulation
                    x_sim, y_sim = screen_to_sim( *self.mouse_position, self.Sim.bounds, self.res)
                    
                    self.Sim.add_source(self.source_selection(amplitude = 1, 
                                                              f_max     = self.max_sim_frequency, 
                                                              delay     = self.Sim.dt * self.Sim.total_steps, 
                                                              position  = (x_sim, y_sim)))
                    
            
            #update Sim
            if not self.pause:
                for _ in range(self.steps):
                    self.Sim.update()
                
                
            #draw simulation
            self._render_sim()
            
            #limit fps
            self.clock.tick(self.fps)
    
    
    
    ## GENERAL ===========================================================================
    
    # setup and ending ------------------------------------------------------------------
        
    def _setup(self):
        
        #get resolution
        res_x, res_y = self.res
        
        
        #initialize simulation frequency
        self.max_sim_frequency = self.Sim.f_max * 2
        
        
        #init pygame
        pygame.init()
        
        #init fonts
        pygame.font.init()
        self.hugefont  = pygame.font.Font(os.path.join("fonts", "OpenSans-Semibold.ttf"), 36)
        self.largefont = pygame.font.Font(os.path.join("fonts", "OpenSans-Semibold.ttf"), 22)
        self.font      = pygame.font.Font(os.path.join("fonts", "OpenSans-Semibold.ttf"), 18)
        self.smallfont = pygame.font.Font(os.path.join("fonts", "OpenSans-Semibold.ttf"), 14)
        
        #icon
        programIcon = pygame.image.load(os.path.join("assets", "icon.png"))
        programIcon  = pygame.transform.scale(programIcon, (32, 32))
        pygame.display.set_icon(programIcon)
        
        #caption
        pygame.display.set_caption("2D FDTD Simulation")
        
        #init display
        self.SCREEN = pygame.display.set_mode(self.res, pygame.RESIZABLE)
        
        
        #current choice for source
        self.source_selection = GaussianSource
        
        
        #current choice of material
        self.material_selecion = RectMaterial
        
        
        #surface for simulation
        self.Sim_surf = pygame.Surface(self.res)
        
        
        #material surface
        self.Mat_surf = pygame.Surface(self.res)
        self.Mat_surf.set_colorkey((0,0,0))
        self.Mat_surf.set_alpha(150)
        
        
        #general settings menu surface
        self.SettingsGeneral_surf =  pygame.Surface(self.res)
        self._setup_settings_general()
        
        
        #general settings menu surface
        self.SettingsSim_surf = pygame.Surface(self.res)
        self._setup_settings_sim()
        
        
        #main menu surface
        self.MainMenu_surf =  pygame.Surface(self.res)
        self._setup_main_menu()
        
        
        #editor surface
        self.Editor_surf =  pygame.Surface(self.res)
        self._setup_editor()
        
        
        #init display data
        self.total_frames   = 0
        self.real_frametime = 1
        self.total_time     = 0
        
        
        #init clock
        self.clock = pygame.time.Clock()
        
        
        #loop conditions
        self.active = True
        self.sub_active = True
        
        
        #pause simulation
        self.pause = False
        
        
        #initial simulation surface
        self.Ez_surf = pygame.Surface((self.Sim.n_x, self.Sim.n_y))
        
        
    def _end(self):
        
        pygame.font.quit()
        pygame.quit()
        #sys.exit()

        
    def _update_state(self):
        
        #count frames
        self.total_frames += 1
            
        #calculate the real frametime in seconds
        self.real_frametime = (pygame.time.get_ticks() - self.total_time) * 1e-3
        
        #update total runtime
        self.total_time = pygame.time.get_ticks()
        
        #get current mouse position
        self.mouse_position = pygame.mouse.get_pos()
        
        #get continuous mouse events
        #self.leftclick, _, self.rightclick = pygame.mouse.get_pressed()
        self.leftclick = False
        
        #buffer keyboard input
        self.keyboard_input_buffer = []
        


    def _resize(self, w, h):
        
        res_x_min, res_y_min = self.res_min
        res_x_max, res_y_max = self.res_max
                    
        new_res = max( min(res_x_max, w), res_x_min ), max( min(res_y_max, h), res_y_min )

        self.SCREEN = pygame.display.set_mode(new_res, pygame.RESIZABLE)
        
        new_res_x, new_res_y = new_res
        res_x    , res_y     = self.res
        
        self.res = new_res
        
        #calculate new bounds for simulation
        x_max, x_min, y_max, y_min = self.Sim.bounds
        
        x_max = (x_max - x_min) * new_res_x / res_x + x_min
        y_max = (y_max - y_min) * new_res_y / res_y + y_min
        
        bounds = x_max , x_min, y_max, y_min
        
        #update simulation
        self.Sim.resize(bounds)
        
        
        #surface for simulation
        self.Sim_surf = pygame.Surface(self.res)
        
        
        #material surface
        self.Mat_surf = pygame.Surface(self.res)
        self.Mat_surf.set_colorkey((0,0,0))
        self.Mat_surf.set_alpha(150)
        
        #general settings menu surface
        self.SettingsGeneral_surf = pygame.Surface(self.res)
        self._setup_settings_general()
        
        
        #general settings menu surface
        self.SettingsSim_surf = pygame.Surface(self.res)
        self._setup_settings_sim()
        
        
        #main menu surface
        self.MainMenu_surf = pygame.Surface(self.res)
        self._setup_main_menu()
        
        #editor surface
        self.Editor_surf = pygame.Surface(self.res)
        self._setup_editor()
        
        #surface for simulation values
        self.Ez_surf = pygame.Surface((self.Sim.n_x, self.Sim.n_y))
        
        
        
    # main loop ----------------------------------------------------------------
    
    def run(self):
        
        """
        main method to get called from outside
        
        """
        
        self._setup()
        
        self.what_to_run_next = "mainmenu"
        
        run_dict = {"sim"              : self._run_sim, 
                    "general_settings" : self._run_settings_general,
                    "sim_settings"     : self._run_settings_sim,
                    "mainmenu"         : self._run_main_menu,
                    "edit"             : self._run_editor}
        
        # main loop that switches between different loops
        while self.active:
            run_dict[self.what_to_run_next]()
            
        
        self._end()
        
        
        
        