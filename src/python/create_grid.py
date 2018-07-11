# -*- coding: utf-8 -*-
'''
Create soil column grid

1D vertical grid types
- FVGrid1D: 1D finite-volume grid with exponentially increasing node depth, face-centered
- QuadGrid1D: 1D collocated grid, with quadratically increasing node depth

2D axi-symmetrical cylindrical grid types
- EquilGridCyl2D: 2D cylindrical grid with equal grid sizes 
- FVGridCyl2D: 2D cylindrical finite-volume grid, with reduced numbers of grid nodes
- FVGridCyl2D_Exp: 2D cylindrical finite-volume grid, with exponentially increasing node depth, not used
'''
import numpy as np

class FVGrid1D:
    def __init__(self):
        n_level = 25
        self.n_level = n_level
        self.min_depth = 0.
        self.max_depth = 1.
        self.level = np.array(range(self.n_level)) + 1
        self.grid_node = np.exp(self.level * 0.2 - 5)
        self.cv_size = np.zeros(self.n_level)
        self.cv_face = np.zeros(self.n_level)
        
        self.cv_size[0] = (self.grid_node[0] + self.grid_node[1]) / 2.
        self.cv_size[1:n_level-1] = (self.grid_node[2:n_level] - 
            self.grid_node[0:n_level-2]) / 2.
        self.cv_size[n_level-1] = self.grid_node[n_level-1] - \
            self.grid_node[n_level-2]
        
        self.cv_face[0:n_level-1] = (self.grid_node[0:n_level-1] + 
            self.grid_node[1:n_level]) / 2.
        self.cv_face[n_level-1] = self.grid_node[n_level-1] + \
            self.cv_size[n_level-1] / 2.

class QuadGrid1D:
    def __init__(self):
        n_level = 26
        self.n_level = n_level
        self.min_depth = 0.
        self.max_depth = 1.
        self.level = np.array(range(self.n_level))
        self.grid_node = np.linspace(self.min_depth, self.max_depth, 26)**2

class EquilGridCyl2D:
    def __init__(self):
        N_r = 50
        N_z = 50
        self.N_r = N_r
        self.N_z = N_z
        self.dr = 1./N_r
        self.dz = 1./N_z
        self.r = np.arange(N_r+1) * self.dr
        self.z = np.arange(N_z+1) * self.dz
        self.r_grid, self.z_grid = np.meshgrid(self.r, self.z)

class FVGridCyl2D: 
    def __init__(self):
        n_vert_level = 25
        n_rad_level = 20
        self.n_vert_level = n_vert_level
        self.n_rad_level = n_rad_level
        self.vert_level = np.arange(n_vert_level)
        self.rad_level = np.arange(n_rad_level)
        self.vert_node = np.concatenate((np.arange(10), np.round(np.arange(0,15)**1.4+10))) * 0.02
        self.rad_node = np.concatenate((np.arange(10), np.round(np.arange(4,14)**1.5+3))) * 0.02
        
        # FV grid nodes can be extracted from the equil-spacing grid nodes with these indices
        self.vert_index = np.concatenate((np.arange(10), np.round(np.arange(0,15)**1.4+10))).astype('int')
        self.rad_index = np.concatenate((np.arange(10), np.round(np.arange(4,14)**1.5+3))).astype('int')

        self.vert_cv_size = np.zeros(self.n_vert_level)
        self.vert_cv_face = np.zeros(self.n_vert_level)
        self.rad_cv_size = np.zeros(self.n_rad_level)
        self.rad_cv_face = np.zeros(self.n_rad_level)

        self.vert_cv_size[0] = np.nan  # the topmost layer is the atmosphere with infinite thickness 
        self.vert_cv_size[1] = (self.vert_node[1] + self.vert_node[2]) / 2.
        self.vert_cv_size[2:n_vert_level-1] = (self.vert_node[3:n_vert_level] - 
            self.vert_node[1:n_vert_level-2]) / 2.
        self.vert_cv_size[n_vert_level-1] = self.vert_node[n_vert_level-1] - \
            self.vert_node[n_vert_level-2]
        
        self.vert_cv_face[0] = 0. 
        self.vert_cv_face[1:n_vert_level-1] = (self.vert_node[1:n_vert_level-1] + 
            self.vert_node[2:n_vert_level]) / 2.
        self.vert_cv_face[n_vert_level-1] = self.vert_node[n_vert_level-1] + \
            self.vert_cv_size[n_vert_level-1] / 2.

        self.rad_cv_size[0] = self.rad_node[1]
        self.rad_cv_size[1:n_rad_level-1] = (self.rad_node[2:n_rad_level] - 
            self.rad_node[0:n_rad_level-2]) / 2.
        self.rad_cv_size[n_rad_level-1] = self.rad_node[n_rad_level-1] - \
            self.rad_node[n_rad_level-2]
        
        self.rad_cv_face[0] = self.rad_node[1]/2.
        self.rad_cv_face[1:n_rad_level-1] = (self.rad_node[1:n_rad_level-1] + 
            self.rad_node[2:n_rad_level]) / 2.
        self.rad_cv_face[n_rad_level-1] = self.rad_node[n_rad_level-1] + \
            self.rad_cv_size[n_rad_level-1] / 2.

        self.r_grid, self.z_grid = np.meshgrid(self.rad_node, self.vert_node)

class FVGridCyl2D_Exp:
    def __init__(self):
        n_vert_level = 26
        n_rad_level = 21
        self.n_vert_level = n_vert_level
        self.n_rad_level = n_rad_level
        self.vert_level = np.arange(n_vert_level)
        self.rad_level = np.arange(n_rad_level)
        self.vert_node = np.exp(self.vert_level * 0.2 - 5)
        self.rad_node = np.exp(self.rad_level * 0.2 - 4)
        self.vert_node[0] = 0.
        self.rad_node[0] = 0.
        self.vert_cv_size = np.zeros(self.n_vert_level)
        self.vert_cv_face = np.zeros(self.n_vert_level)
        self.rad_cv_size = np.zeros(self.n_rad_level)
        self.rad_cv_face = np.zeros(self.n_rad_level)

        self.vert_cv_size[0] = np.nan  # the topmost layer is the atmosphere with infinite thickness 
        self.vert_cv_size[1] = (self.vert_node[1] + self.vert_node[2]) / 2.
        self.vert_cv_size[2:n_vert_level-1] = (self.vert_node[3:n_vert_level] - 
            self.vert_node[1:n_vert_level-2]) / 2.
        self.vert_cv_size[n_vert_level-1] = self.vert_node[n_vert_level-1] - \
            self.vert_node[n_vert_level-2]
        
        self.vert_cv_face[0] = 0. 
        self.vert_cv_face[1:n_vert_level-1] = (self.vert_node[1:n_vert_level-1] + 
            self.vert_node[2:n_vert_level]) / 2.
        self.vert_cv_face[n_vert_level-1] = self.vert_node[n_vert_level-1] + \
            self.vert_cv_size[n_vert_level-1] / 2.

        self.rad_cv_size[0] = self.rad_node[1]
        self.rad_cv_size[1:n_rad_level-1] = (self.rad_node[2:n_rad_level] - 
            self.rad_node[0:n_rad_level-2]) / 2.
        self.rad_cv_size[n_rad_level-1] = self.rad_node[n_rad_level-1] - \
            self.rad_node[n_rad_level-2]
        
        self.rad_cv_face[0] = self.rad_node[1]/2.
        self.rad_cv_face[1:n_rad_level-1] = (self.rad_node[1:n_rad_level-1] + 
            self.rad_node[2:n_rad_level]) / 2.
        self.rad_cv_face[n_rad_level-1] = self.rad_node[n_rad_level-1] + \
            self.rad_cv_size[n_rad_level-1] / 2.
        
        self.r_grid, self.z_grid = np.meshgrid(self.rad_node, self.vert_node)
