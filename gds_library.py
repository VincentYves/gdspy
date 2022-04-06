# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:46:41 2021

@author: pelgrin-adm
"""

import os
import numpy as np
import gdspy
import math

def pair(N):
    m = int(N/2)
    N = m * 2
    return N

def ispair(N):
    m = pair(N)
    if m-N ==0:
        p = 1
    else:
        p = 0
    return p
# print('Using gdspy module version ' + gdspy.__version__)

#--------------------------------------------------------------------------------------------------#

#FILE DEFINITIONS

print('Using gdspy module version ' + gdspy.__version__)
precision = 2.0e-9
lib = gdspy.GdsLibrary(name = 'First_gds', precision=precision)

# =============================================================================
# Fontions utiles
# =============================================================================
    
def check_radius_Flexpath(radius,points):
    xs = []
    ys = []
    for i,tup in enumerate (points):
        xs.append(tup[0])
        ys.append(tup[1])
    
    r = abs(min(min(np.diff(np.array(xs))),min(np.diff(np.array(ys))))/2)    
    if r < radius:
        print("radius is to big, correcting...")
        radius = r
    return radius

def add_zero_rectangle(Cell,xmin,xmax,ymin,ymax,Background = 100):
      
    poly = gdspy.Rectangle((xmin,ymin),(xmax,ymax),layer = Background)
    Background += 1 
    Cell.add(poly)
    return Background

# =============================================================================
# Straight waveguide
# =============================================================================

class Straight_wg():

    def __init__(self,width, initial_point,length, direction,etch_width = 3, layer = 1,**kwargs):
        """
        Parameters
        ----------
        width : float
            width of the waveguide (um)
        initial_point : tuple
            Coordinates of the origin (um, um)
        length : float
            length of the waveguide (um)
        direction : str
            '+x', '-x', '+y', '-y'
        layer : int
            order of the layer

        """
        self.width = width
        self.initial_point = initial_point
        self.length = length 
        self.direction = direction 
        self.etch_width = etch_width
        self.layer = layer
        
        self.Ox = self.initial_point[0]
        self.Oy = self.initial_point[1]
        
        if self.direction == '+x':
            self.Ex = self.Ox + self.length
            self.Ey = self.Oy

        elif self.direction == '-x':
            self.Ex = self.Ox - self.length
            self.Ey = self.Oy

        elif self.direction == '+y':
            self.Ex = self.Ox 
            self.Ey = self.Oy + self.length
            
        else:
            self.Ex = self.Ox 
            self.Ey = self.Oy - self.length
            
        self.name = "straight_waveguide"+str(width)+str(length)+direction+str(self.Ox)+str(self.Oy)

    def add2cell(self,Cell):
        
        # path = gdspy.Path(3,self.initial_point, number_of_paths=2, distance=self.width+3)
        path = gdspy.Path(self.etch_width,self.initial_point, number_of_paths=2, distance=self.width+self.etch_width)
        path.segment(self.length,self.direction,layer = self.layer)
        Cell.add(path)


# =============================================================================
# Round spiral
# =============================================================================

class Round_spiral():
    
    def __init__(self,width, initial_point,roundt,rot = np.pi/2,
                 turn_radius = 40,etch_width = 3,layer = 1,**kwargs):
        """
        Parameters
        ----------
        width : float
            width of the waveguide (um)
        initial_point : tuple
            Coordinates of the origin (um, um)
        roundt : int
            Number of full round trip.
        rot : float
            rotation angle in radians
        turn_radius : float, optional
            radius of the half circles in the center of the spiral. The default is 8.
        layer : int, optional
            order of the layer. The default is 0.

        """
        
        self.width = width
        self.initial_point = initial_point
        self.roundt = roundt
        self.etch_width = etch_width
        self.layer = layer
        self.turn_radius = turn_radius
        self.rot = rot
        self.Ox = self.initial_point[0]
        self.Oy = self.initial_point[1]
        
        self.waveguide_spacing = 0.5* (self.etch_width+ self.width)
        self.enr = self.roundt*2       
        self.turn_radius = turn_radius
        self.Iw = 4 * self.turn_radius
        self.dw = 2 * self.waveguide_spacing
        self.space = self.enr * self.dw
        self.size = (self.Iw + 2 * self.space)/2
        self.name = "Spiral"+str(self.width)+str(self.Ox)+str(self.Oy)+str(self.roundt)
        self.Ex = self.Ox
        self.Ey = self.Oy-self.size*2

        
    def add2cell(self,Cell):    
        
        def spiral(u):
            r = self.size - self.space * u
            theta = self.enr * u * np.pi
            x = r * np.cos(theta) - self.size
            y = -r * np.sin(theta) 
            return (x, y)
        
        def spiral2(u):
            u = 1-u
            r = self.size - self.space * u
            theta = self.enr * u * np.pi
            x = - r * np.cos(theta) - self.size 
            y =  r * np.sin(theta) 
            return (x-2*self.size+self.space, y)
        
        def dspiral_dt(u):
            theta = self.enr * u * np.pi
            dx_dt = -np.sin(theta)
            dy_dt = -np.cos(theta)
            return (dx_dt, dy_dt)
        
        def dspiral_dt2(u):
            u = 1-u
            theta = self.enr * u * np.pi
            dx_dt = np.sin(theta)
            dy_dt = np.cos(theta)
            return (dx_dt, dy_dt)
            
        
        path1 = gdspy.Path(self.etch_width,(0,0), number_of_paths=2, distance=self.width+self.etch_width)
        path1.parametric(spiral, dspiral_dt,layer = self.layer)        
        path1.turn(self.turn_radius,'rr',layer = self.layer,tolerance= 0.005)
        path1.turn(self.turn_radius,'ll',layer = self.layer,tolerance= 0.005)        
        path2 = gdspy.Path(self.etch_width,(-path1.x,-path1.y), number_of_paths=2, distance=self.width+self.etch_width)
        path2.parametric(spiral2, dspiral_dt2,layer = self.layer)
        
        spiral = gdspy.boolean(path1, path2, 'or',layer = self.layer)
        spiral.rotate(self.rot)
        spiral.translate(self.Ox, self.Oy)

        Cell.add(spiral)
        
    def adjust_output(self,Cell,L, DeltaY):
        
        S = abs(self.Oy-self.Ey) - DeltaY
        path = gdspy.Path(self.etch_width,(self.Ex,self.Ey), number_of_paths=2, distance=self.width+self.etch_width)
        path.bezier([(L/2,0),(L/2,S),(L,S)],layer = self.layer)
        Cell.add(path)
        self.Ex = self.Ex + L
        self.Ey = self.Ey + S
        return S

# =============================================================================
# Ring
# =============================================================================

class Simple_Ring():
    
    def __init__(self,width,center,radius,gap,updown = -1,etch_width = 3,layer = 1,**kwargs):
        """
        Parameters
        ----------
        width : TYPE
            DESCRIPTION.
        center : TYPE
            DESCRIPTION.
        radius : TYPE
            DESCRIPTION.
        gap : TYPE
            DESCRIPTION.
        layer : TYPE, optional
            DESCRIPTION. The default is 1.

        """
        
        self.width = width
        self.center = center 
        self.radius = radius 
        self.gap = gap 
        self.etch_width = etch_width
        self.layer = layer
        self.updown = updown
        
        self.Ox = self.center[0]-(self.radius+10)
        self.Oy = self.center[1]+self.updown*(self.radius+self.gap+self.width)
        self.Ex = self.Ox + 2*self.radius+20
        self.Ey = self.Oy
        
    def add2cell(self,Cell):
        
        circle_etch = gdspy.Round(self.center, self.radius+(self.etch_width+self.width/2),inner_radius=self.radius-(self.etch_width+self.width/2),tolerance= 0.005,layer = self.layer)
        circle_ring = gdspy.Round(self.center, self.radius+(self.width/2),inner_radius=self.radius-(self.width/2),tolerance= 0.005,layer = self.layer)
       
        bus_etch = gdspy.Path(2*self.etch_width+self.width,initial_point=(self.center[0]-(self.radius+10),self.center[1]+self.updown*(self.radius+self.gap+self.width)),number_of_paths=1)
        bus_etch.segment(2*self.radius+20,direction = '+x',layer = self.layer)
        bus = gdspy.Path(self.width,initial_point=(self.center[0]-(self.radius+10),self.center[1]+self.updown*(self.radius+self.gap+self.width)),number_of_paths=1)
        bus.segment(2*self.radius+20,direction = '+x',layer = self.layer)
        
        ring = gdspy.boolean(gdspy.boolean(bus, circle_ring, 'or',max_points=5,layer = self.layer),
                     gdspy.boolean(bus_etch,circle_etch, 'or',layer = self.layer),'xor',layer = self.layer)

        Cell.add(ring)
        

# =============================================================================
# Bend_guides
# =============================================================================
class Bend_guides():
    
    def __init__(self,width, initial_point,radius,Nturn, rot = 0, 
                 etch_width = 3, layer = 1,Lin = 50,Lout = 50,**kwargs):
        """
        Parameters
        ----------
        width : float
            width of the waveguide (um)
        initial_point : tuple
            Coordinates of the origin (um, um)
        radius : float
            radius of the turns.
        Nturn : int
            number of turns.
        rot : float, optional
            rotation angle in radians. The default is 0.
        layer : int, optional
            order of the layer. The default is 0.
        Lin : float, optional
            length of the input waveguide. The default is 50.
        Lout : float, optional
            length of the output waveguide. The default is 50.

        """
        self.width = width 
        self.initial_point = initial_point
        self.Ox = self.initial_point[0]
        self.Oy = self.initial_point[1]
        self.radius = radius
        self.Nturn = 2*Nturn
        self.rot = rot
        self.etch_width = etch_width
        self.layer = layer
        self.Lin = Lin
        self.Lout = Lout
        
        self.span = 2 * self.radius
        
        self.Ex = self.Ox+self.Nturn*self.span +self.Lin+self.Lout 
        self.Ey = self.Oy     
        self.length = abs(self.Ex-self.Ox)
        
    def add2cell(self,Cell):
        
        path = gdspy.Path(self.etch_width,self.initial_point, number_of_paths=2, distance=self.width+self.etch_width)
        path.segment(self.Lin,'+x',layer = self.layer)
        path.turn(self.radius, 'l',layer = self.layer)
        for n in range (int(self.Nturn/2)-1):    
            path.turn(self.radius, 'rr',layer = self.layer)
            path.segment(self.span,'-y',layer = self.layer)
            path.turn(self.radius, 'll',layer = self.layer)
            path.segment(self.span,'+y',layer = self.layer)
        
        path.turn(self.radius, 'rr',layer = self.layer)  
        path.turn(self.radius, 'l',layer = self.layer)
        path.segment(self.Lout,'+x',layer = self.layer)
        Cell.add(path)
    
    
# =============================================================================
# Draw Elipse
# =============================================================================

class DrawElipseLine(): 
    
    def __init__(self,initial_point, Np_i, GrArc_ini, GrArc_end, Gr_Lsi, Gr_Lgap,
                   lambda_o, FocalLength, nSlab, nc, theta_in, NPoints, layer,rot = 0):
        # self.Cell = Cell
        self.initial_point = initial_point
        self.Ox = self.initial_point[0]
        self.Oy = self.initial_point[1]
        self.Np_i = Np_i
        self.GrArc_ini = GrArc_ini
        self.GrArc_end = GrArc_end
        self.Gr_Lsi = Gr_Lsi
        self.Gr_Lgap = Gr_Lgap
        self.lambda_o = lambda_o
        self.FocalLength = FocalLength
        self.nSlab = nSlab
        self.nc = nc
        self.theta_in = theta_in 
        self.NPoints = NPoints 
        self.layer = layer
        self.rot = rot 
        
        self.Gr_Pitch_L = self.Gr_Lsi + self.Gr_Lgap
        self.Gr_DC = self.Gr_Lsi / self.Gr_Pitch_L 
        self.nBF = self.nc*np.sin(self.theta_in)+self.lambda_o/self.Gr_Pitch_L
        self.k = (self.nSlab-self.nc*np.sin(self.theta_in))/(self.nBF-self.nc*np.sin(self.theta_in)) 
        self.Np_ini = int(np.floor(self.FocalLength/self.Gr_Pitch_L))

    def add2cell(self,Cell):
        def elipse(t):
            phi = (self.GrArc_end-self.GrArc_ini) * t + self.GrArc_ini
            r=(self.Np_i*self.lambda_o)/(self.nBF-self.nc*np.cos(phi)*np.sin(self.theta_in))+(self.k*self.Np_ini*self.lambda_o)/(self.nSlab-self.nc*np.cos(phi)*np.sin(self.theta_in))
            x = r * np.cos(phi)
            y = r * np.sin(phi)
            return (x, y)
    	
        def delipse_dt(t):
            phi = (self.GrArc_end-self.GrArc_ini) * t + self.GrArc_ini
            r=(self.Np_i*self.lambda_o)/(self.nBF-self.nc*np.cos(phi)*np.sin(self.theta_in))+(self.k*self.Np_ini*self.lambda_o)/(self.nSlab-self.nc*np.cos(phi)*np.sin(self.theta_in))
            x = r * -np.sin(phi)
            y = r * np.cos(phi)
            return (x, y)	
    
        path1 = gdspy.Path(0.5, (0,0))	
        path1.parametric(elipse, delipse_dt,final_width=lambda t: (self.lambda_o)/(self.nBF-self.nc*np.cos((self.GrArc_end-self.GrArc_ini) * t + self.GrArc_ini)*np.sin(self.theta_in))*(1-self.Gr_DC),number_of_evaluations=self.NPoints, layer=self.layer)
        path1.rotate(self.rot)
        path1.translate(self.Ox, self.Oy)
        Cell.add(path1)

# =============================================================================
# Butt Coupling Taper 
# =============================================================================

class ButtCouplingTaper():
    
    def __init__(self,width,L_taper_to_wg = 100,Final_width = 3,L_multimode = 500,
                 etch_width = 3,layer = 1, ID = 'Butt_taper',**kwargs):
        self.width = width
        self.L_taper_to_wg = L_taper_to_wg
        self.Final_width = Final_width
        self.L_multimode = L_multimode
        self.etch_width = etch_width
        self.layer = layer
        self.ID = ID
        
        self.name = ID +str(self.width)+str(self.L_taper_to_wg)+str(self.Final_width)+str(self.L_multimode)
        
        
    def add2cell(self,Cell,Ox,Oy,rot):
        
        self.Ox = Ox 
        self.Oy = Oy
        self.rot = rot 
        # self.name = self.name+str(self.Ox)+str(self.Oy)+str(rot)
        
        if self.name not in lib.cells:              
            coupler = lib.new_cell(self.name)
                
            path1 = gdspy.Path(self.etch_width, (-self.L_taper_to_wg, 0),number_of_paths = 2,distance = self.width+self.etch_width) 
            path1.segment(self.L_taper_to_wg, '+x',final_distance =  self.Final_width+self.etch_width,layer = self.layer)
            path1.segment(self.L_multimode, '+x',layer=self.layer)
            coupler.add(path1)
            path1.rotate(self.rot,(-self.L_taper_to_wg/2,0))
            Cell.add(gdspy.copy(path1,self.Ox,self.Oy))
            

        else: 
            ref = lib.cells[self.name]
            coupler = ref.copy(self.name+str(Ox)+str(Oy)+str(rot),translation=(Ox,Oy),rotation=rot+np.pi)
            Cell.add(gdspy.copy(coupler,self.Ox,self.Oy))
            
            
# =============================================================================
# Grating
# =============================================================================

    
class FocusingGratingCoupler_SWG(DrawElipseLine): 
    
    def __init__(self,Lsi,Lgap,Wsi,Wgap,
                               Nperiod = 50, FocalLength = 45,
								GrArc = 0.18, TaperArc = 0.1651,
								theta_in = 30*np.pi/180, nc = 1.45,
                                nSlab = 2.83, lambda_o = 1.55,
								W_EndTaper = 0.55,L_taper_to_wg = 10, W_wgInterconnect = 0.45,
                                etch_width = 3, NPoints = 3, layer = 1,ID = 'FocusSWGC'):

        self.Lsi = Lsi
        self.Lgap = Lgap
        self.Wsi = Wsi 
        self.Wgap = Wgap 
        self.Nperiod = Nperiod 
        self.FocalLength = FocalLength
        self.GrArc = GrArc 
        self.TaperArc = TaperArc						
        self.theta_in = theta_in          
        self.nc = nc								 
        self.nSlab = nSlab
        self.lambda_o =lambda_o
        self.W_EndTaper = W_EndTaper 
        self.W_wgInterconnect = W_wgInterconnect 
        self.etch_width = etch_width 
        self.NPoints = NPoints
        self.layer = layer                  

        self.L_taper_to_wg = L_taper_to_wg 
        self.Ox,self.Oy = (self.L_taper_to_wg,0)
        self.Np_ini = np.floor(self.FocalLength/(self.Lsi+self.Lgap))
        self.nBF = self.nc*np.sin(self.theta_in)+self.lambda_o/(self.Lsi+self.Lgap)
        self.k=(self.nSlab-self.nc*np.sin(self.theta_in))/(self.nBF-self.nc*np.sin(self.theta_in)) 
        self.ID = ID
        self.name = ID+str(self.Lsi)+str(self.Lgap)+str(self.Wsi)+str(self.Wgap)+str(self.Nperiod)+str(self.FocalLength)+str(self.GrArc)+str(self.TaperArc)+str(self.theta_in)+str(self.nc)+str(self.nSlab)+str(self.lambda_o)+str(self.W_EndTaper)+str(self.W_wgInterconnect)
    
    def add2cell(self,Cell,Ox,Oy,rot):
                        
        self.Ox = Ox + self.L_taper_to_wg
        self.Oy = Oy
        self.rot = rot
        self.name = self.name+str(self.Ox)+str(self.Oy)+str(rot)
        
        if self.name not in lib.cells:                   
            coupler = lib.new_cell(self.name)
            
            for iL in range(0,self.Nperiod):
                Np_i     = iL+1
                Radius   = (Np_i*self.lambda_o)/(self.nBF-self.nc*np.cos(0)*np.sin(self.theta_in))+(self.k*self.Np_ini*self.lambda_o)/(self.nSlab-self.nc*np.cos(0)*np.sin(self.theta_in))
                PitchArc = (self.Wsi+self.Wgap)/Radius
                GapArc   = PitchArc*self.Wgap/(self.Wsi+self.Wgap)
                SiArc    = PitchArc*self.Wsi/(self.Wsi+self.Wgap)
                Gr_NumPeriodos_T = np.floor(self.GrArc/PitchArc)
                Gr_NumPeriodos_T = int(Gr_NumPeriodos_T)+1
                Swg_GrArc_ini=0
                
                for i_T in range(0, Gr_NumPeriodos_T):
                    E1 = DrawElipseLine((0,0),Np_i,Swg_GrArc_ini,
                                   Swg_GrArc_ini+GapArc,self.Lsi,self.Lgap,self.lambda_o,
                                   self.FocalLength,self.nSlab,self.nc,self.theta_in,self.NPoints,self.layer, rot = 0)#self.rot)
                    E1.add2cell(coupler)
                    E2 = DrawElipseLine((0,0),Np_i,-Swg_GrArc_ini-SiArc,
                                   -Swg_GrArc_ini-GapArc-SiArc,self.Lsi,self.Lgap,self.lambda_o,
                                   self.FocalLength,self.nSlab,self.nc,self.theta_in,self.NPoints,self.layer, rot = 0)#self.rot)
                    E2.add2cell(coupler)
                    Swg_GrArc_ini=Swg_GrArc_ini+PitchArc
        
            self.L_Taper = self.FocalLength+(self.Lsi+self.Lgap)*self.Nperiod         
            Grat_FinalWidth = self.L_Taper*np.tan(2*self.TaperArc)
            
            path1 = gdspy.Path(self.etch_width, (-self.L_taper_to_wg, 0),number_of_paths = 2,distance = self.W_wgInterconnect+self.etch_width) 
            path1.segment(self.L_taper_to_wg, '+x',final_distance = self.W_EndTaper+self.etch_width,layer = self.layer)
            path1.segment(self.L_Taper, '+x',final_distance = Grat_FinalWidth+self.etch_width,layer=self.layer)

            grat = gdspy.boolean(path1, coupler, 'or',layer = self.layer)
            grat.rotate(rot,(-self.L_taper_to_wg,0))
            Cell.add(gdspy.copy(grat,self.Ox,self.Oy))

        else: 
            ref = lib.cells[self.name]
            coupler = ref.copy(self.name+str(Ox)+str(Oy)+str(rot),translation=(Ox,Oy),rotation=rot)
            path1 = gdspy.Path(self.etch_width, (-self.L_taper_to_wg, 0),number_of_paths = 2,distance = self.W_wgInterconnect+self.etch_width) 
            path1.segment(self.L_taper_to_wg, '+x',final_distance = self.W_EndTaper+self.etch_width,layer = self.layer)
            path1.segment(self.L_Taper, '+x',final_distance = Grat_FinalWidth+self.etch_width,layer=self.layer)

            grat = gdspy.boolean(path1, coupler, 'or',layer = self.layer)
            grat.rotate(rot,(-self.L_taper_to_wg,0))
            Cell.add(gdspy.copy(grat,self.Ox,self.Oy))
         
# =============================================================================
# PhC    
# =============================================================================


def maille_tri(Ox, Oz, d, Nx = 2,Nz = 9, Inv_xs = False):
    """
    Parameters
    ----------
    Ox : TYPE
        DESCRIPTION.
    Oz : TYPE
        DESCRIPTION.
    d : TYPE
        DESCRIPTION.
    Nx : TYPE, optional
        DESCRIPTION. The default is 2.
    Nz : TYPE, optional
        DESCRIPTION. The default is 9.
    Inv_xs : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    xx1 : TYPE
        DESCRIPTION.
    xx2 : TYPE
        DESCRIPTION.
    zz1 : TYPE
        DESCRIPTION.
    zz2 : TYPE
        DESCRIPTION.
    Nx1 : TYPE
        DESCRIPTION.
    Nz1 : TYPE
        DESCRIPTION.
    width : TYPE
        DESCRIPTION.
    length : TYPE
        DESCRIPTION.

    """
   
    # Fixes the proper size for Lattices 1 & 2
    if ispair(Nx) == 1:
        Nx1 = int(Nx/2)
    else:
        Nx1 = int(Nx/2)+1
    if ispair(Nz) == 1:
        Nz1 = int(Nz/2)
    else:
        Nz1 = int(Nz/2)+1
        
    # Lattice 1
    Ox1 = Ox
    Oz1 = Oz
    dz = d*np.sqrt(3)
    
    xx1 = np.arange(-(Nx1/2)*d+d/2,(1+Nx1/2)*d-d/2,d)+Ox1
    zz1 = np.arange(-(Nz1/2)*dz+dz/2,(0.1+Nz1/2)*dz-dz/2,dz)+Oz1
    
    # Lattice 2
    
    xx2 = xx1[:Nx-Nx1] + d/2
    zz2 = zz1[:Nz-Nz1] + d*np.sqrt(3)/2
    
    if Inv_xs == True:     
        xx1,xx2 = -xx1,-xx2

    # Size of the whole PhC
    width_PhC = abs(min(np.concatenate((zz1,zz2)))-max(np.concatenate((zz1,zz2))))
    length = abs(min(np.concatenate((xx1,xx2)))-max(np.concatenate((xx1,xx2))))
    
    return xx1,xx2,zz1,zz2,Nx1,Nz1,width_PhC,length

 

class PhC_W1():
    
    def __init__(self,width,initial_point,length,d,radius, r1 = 1, x1 = 0, s1 = 0, Nx = 2,Ny = 9,
                 etch_width = 3, layer = 1,name = 'PhC_W1',**kwargs):

        self.etch_width = etch_width
        self.layer = layer
        self.width = width
        self.initial_point = initial_point
        self.length = length
        self.d = d 
        self.Ox = self.initial_point[0]
        self.Oy = self.initial_point[1]        
        
        self.radius = radius
        self.Nx = Nx
        self.Ny = Ny 
        self.width_PhC = 2 * np.sqrt(3)*d
        self.xx1,self.xx2,self.yy1,self.yy2,self.Nx1,self.Ny1,self.span_y,_ = maille_tri(self.Ox, self.Oy, self.d, self.Nx,self.Ny)                
        self.span_x = self.d/2               
        self.Ncycl = int(self.length/(2*self.span_x))                          # Number of time the base cell is repeated
        self.taper = 0
        
        
        # Slow waveguide parameters
        self.r1 = r1
        self.x1 = x1 * self.d
        self.s1 = s1 * self.d
        
        # Coordinates of the end of the W1 waveguide 
        self.Ex = self.Ox+2*self.Ncycl*self.span_x
        self.Ey = self.Oy 
        
        self.name = name + str(round(np.random.randn()*1e6))
        # Library of all components
        self.lib_PhC = gdspy.GdsLibrary(name = "lib_PhC"+name+str(self.width_PhC)+str(self.length)+str(self.Ox)+str(self.Oy))
        # self.name = str(round(10000*self.d))+str(self.width_PhC)+str(self.length)+str(self.Ox)+str(self.Oy)+str(self.radius)
        

    def add2cell(self,Cell):

        # self.name = str(round(10000*self.d))+str(self.width_PhC)+str(self.length)+str(self.Ox)+str(self.Oy)+str(self.radius)
        #Construct the base Cell   
        # if 'base_W1'+self.name not in lib.cells:
        base_W1 = lib.new_cell('base_W1'+self.name)#+str(np.random.randn()))       
    
        for x in self.xx2:
            for y in self.yy1:
                circle = gdspy.Round((x, y), self.radius, tolerance=0.0005,layer = self.layer)
                base_W1.add(circle)
        for x in self.xx1:
            for y in self.yy2:
                circle = gdspy.Round((x, y), self.radius, tolerance=0.0005,layer = self.layer)
                base_W1.add(circle) 
                
        circle = gdspy.Round((self.initial_point[0] +self.x1+0*self.d/2, self.initial_point[1]+self.s1-self.span_y/2-self.d*np.sqrt(3)/2), self.r1*self.radius, tolerance=0.0005,layer = self.layer)
        base_W1.add(circle)    
            
        # Duplicates the base Cell in an Array
        # PhC_array = gdspy.CellArray(base_W1,self.Ncycl,1,(self.d,self.d),(0,0))
        PhC_array = gdspy.CellArray(base_W1,self.Ncycl,1,(self.d,self.d),(0,0))
        PhC = gdspy.Cell('PhC'+self.name)
        PhC.add(PhC_array)
        
        
        # Add the last hole for symetrization
        circle = gdspy.Round((self.initial_point[0] +self.x1+(2*self.Ncycl)*self.span_x, self.initial_point[1]+self.s1-self.span_y/2-self.d*np.sqrt(3)/2), self.r1*self.radius, tolerance=0.0005,layer = self.layer)        
        PhC.add(circle)

        
        # Makes copy of the PhCs in order to form a W1 waveguide
        PhC_up = PhC.copy("PhC_up"+self.name,translation=(0,(self.span_y+self.width_PhC)/2))
        PhC_down = PhC.copy("PhC_down"+self.name,translation=(2*self.initial_point[0]+2*self.span_x*self.Ncycl,2*self.initial_point[1]-(self.span_y+self.width_PhC)/2),rotation=(np.pi))
        Cell.add(PhC_up)
        Cell.add(PhC_down) 
        
        for cell in self.lib_PhC.cells : 
            element = self.lib_PhC.cells[cell]
            Cell.add(element)
            
        # lib.remove(base_W1,remove_references=False)
               
    def set_taper(self,d,N):       
        
        if self.taper == 0:
            Nxt = 2*N
        else:
            Nxt = 2*N#+1
        name = str(round(1000*d))+str(self.width_PhC)+str(self.length)+str(self.Ox)+str(self.Oy)+str(self.taper)+self.name
        taper = self.lib_PhC.new_cell('taper'+name) 
        
        # Input taper
        xx1_in,xx2_in,yy1_in,yy2_in,Nx1_in,Ny1_in,span_y_in,span_x_in = maille_tri(0,0,d,Nxt,self.Ny,Inv_xs = False)                 
        for x in xx1_in:
            for y in yy1_in:
                circle = gdspy.Round((self.Ox+x-(N+1)*d/2, y+self.Oy), self.radius, tolerance=0.0005,layer = self.layer)
                taper.add(circle)
        for x in xx2_in: 
            for y in yy2_in:
                circle = gdspy.Round((self.Ox+x-(N+1)*d/2, y+self.Oy), self.radius, tolerance=0.0005,layer = self.layer)
                taper.add(circle) 
                
        
        # Output taper  
        xx1_out,xx2_out,yy1_out,yy2_out,Nx1_out,Ny1_out,span_y_out,span_x_out = maille_tri(0,0,d,Nxt,self.Ny,Inv_xs = True)     
        for x in xx1_out:
            for y in yy1_out:
                circle = gdspy.Round((self.Ex+x+(N+1)*d/2, y+self.Oy), self.radius, tolerance=0.0005,layer = self.layer)
                taper.add(circle)
        for x in xx2_out: 
            for y in yy2_out:
                circle = gdspy.Round((self.Ex+x+(N+1)*d/2, y+self.Oy), self.radius, tolerance=0.0005,layer = self.layer)
                taper.add(circle)    
         
        # Sets the tapers 
        tapers =  self.lib_PhC.new_cell("tapers"+name)
        taper_up = taper.copy("taper_up"+name,translation=(0,(span_y_out+self.width_PhC/2)/2))
        taper_down = taper.copy("taper_down"+name,translation=(0,-(span_y_out+self.width_PhC/2)/2))
        tapers.add(taper_up)
        tapers.add(taper_down) 
        self.lib_PhC.remove(taper)              
        
        self.taper += 1

        # Actualize the injection points
        self.Ox = self.Ox - span_x_in - d/2
        self.Ex = self.Ex + span_x_out + d/2       
        self.length = abs(self.Ex-self.Ox)
        
    
    def set_outputs(self,width,L_ino=50):

        # Create the waveguides
        distance = self.width_PhC/2+(self.etch_width-(self.width_PhC/2-width))
        wg_in = gdspy.Path(self.etch_width,(self.Ox,self.Oy),2,distance = self.width_PhC/2+self.etch_width)
        wg_out = gdspy.Path(self.etch_width,(self.Ex,self.Ey),2,distance = self.width_PhC/2+self.etch_width)       
        wg_in.segment(L_ino,direction='-x',final_distance=self.width+self.etch_width,layer = self.layer)
        wg_out.segment(L_ino,direction='+x',final_distance=self.width+self.etch_width,layer = self.layer)

        # Add to the Cell 

        strip_ino = self.lib_PhC.new_cell("strip_ino"+self.name)
        strip_ino.add(wg_in)
        strip_ino.add(wg_out)

        
        # Actualize the injection points
        self.Ox = self.Ox - L_ino
        self.Ex = self.Ex + L_ino
        self.length = abs(self.Ex-self.Ox)
    # def __init__(self,width,initial_point,length,d,radius, Nx = 2,Ny = 9,
    #              etch_width = 3, layer = 1,**kwargs):
    #     """
    #     Parameters
    #     ----------
    #     width_PhC : TYPE
    #         DESCRIPTION.
    #     length : TYPE
    #         DESCRIPTION.
    #     Ox : TYPE
    #         DESCRIPTION.
    #     Oy : TYPE
    #         DESCRIPTION.
    #     d : TYPE
    #         DESCRIPTION.
    #     radius : TYPE
    #         DESCRIPTION.
    #     Nx : TYPE, optional
    #         DESCRIPTION. The default is 2.
    #     Ny : TYPE, optional
    #         DESCRIPTION. The default is 9.
    #     """
        
    #     self.etch_width = etch_width
    #     self.layer = layer
    #     self.width = width
    #     self.initial_point = initial_point
    #     self.length = length
    #     self.Ox = self.initial_point[0] 
    #     self.Oy = self.initial_point[1]        
    #     self.d = d 
    #     self.radius = radius
    #     self.Nx = Nx
    #     self.Ny = Ny 
    #     self.width_PhC = np.sqrt(3)*d
    #     self.xx1,self.xx2,self.yy1,self.yy2,self.Nx1,self.Ny1,self.span_y,_ = maille_tri(self.Ox, self.Oy, self.d, self.Nx,self.Ny,)                
    #     self.span_x = self.d/2               
    #     self.Ncycl = int(self.length/(2*self.span_x))                          # Number of time the base cell is repeated
        
    #     # Coordinates of the end of the W1 waveguide 
    #     self.Ex = self.Ox+2*self.Ncycl*self.span_x 
    #     self.Ey = self.Oy 
        
    #     # Library of all components
    #     self.lib_PhC = gdspy.GdsLibrary(name = "lib_PhC"+str(self.width_PhC)+str(self.length)+str(self.Ox)+str(self.Oy))
    #     self.name = str(round(10000*self.d))+str(self.width_PhC)+str(self.length)+str(self.Ox)+str(self.Oy)+str(self.radius)

    # def add2cell(self,Cell):
    
    #     #Construct the base Cell   
    #     base_W1 = lib.new_cell('base_W1'+self.name)       
        
    #     for x in self.xx1:
    #         for y in self.yy1:
    #             circle = gdspy.Round((x, y), self.radius, tolerance=0.0005,layer = self.layer)
    #             base_W1.add(circle)
    #     for x in self.xx2:
    #         for y in self.yy2:
    #             circle = gdspy.Round((x, y), self.radius, tolerance=0.0005,layer = self.layer)
    #             base_W1.add(circle) 

    #     # Duplicates the base Cell in an Array
    #     PhC_array = gdspy.CellArray(base_W1,self.Ncycl,1,(self.d,self.d),(0,0))
    #     PhC = gdspy.Cell('PhC'+self.name)
    #     PhC.add(PhC_array)

    #     # Add the last row for symetrization
    #     zz = np.arange(-(self.Ny1/2)*self.d*np.sqrt(3) +self.d*np.sqrt(3) /2,(0.1+self.Ny1/2)*self.d*np.sqrt(3) -self.d*np.sqrt(3) /2,self.d*np.sqrt(3) ) 
    #     for z in zz : 
    #         circle = gdspy.Round((self.initial_point[0]+2*self.Ncycl*self.span_x, self.initial_point[1]+z), self.radius, tolerance=0.0005,layer = self.layer)
    #         PhC.add(circle)
        
    #     # Makes copy of the PhCs in order to form a W1 waveguide
    #     PhC_up = PhC.copy("PhC_up"+self.name,translation=(0,(self.span_y+self.width_PhC)/2))
    #     PhC_down = PhC.copy("PhC_down"+self.name,translation=(0,-(self.span_y+self.width_PhC)/2))
    #     Cell.add(PhC_up)
    #     Cell.add(PhC_down)      
        
    #     for cell in self.lib_PhC.cells : 
    #         element = self.lib_PhC.cells[cell]
    #         Cell.add(element)
               
    # def set_taper(self,d,N):
    #     """
    #     Parameters
    #     ----------
    #     d : TYPE
    #         DESCRIPTION.
    #     N : TYPE
    #         DESCRIPTION.

    #     """
            
    #     Nxt = 2 * N
    #     name = str(round(1000*d))+str(self.width_PhC)+str(self.length)+str(self.Ox)+str(self.Oy)+self.name
    #     taper = self.lib_PhC.new_cell('taper'+name) 
        
    #     # Input taper
    #     xx1_in,xx2_in,yy1_in,yy2_in,Nx1_in,Ny1_in,span_y_in,span_x_in = maille_tri(0,0,d,Nxt,self.Ny,Inv_xs = False)                 
    #     for x in xx1_in:
    #         for y in yy1_in:
    #             circle = gdspy.Round((self.Ox+x-(N+1)*d/2, y+self.Oy), self.radius, tolerance=0.0005,layer = self.layer)
    #             taper.add(circle)
    #     for x in xx2_in:
    #         for y in yy2_in:
    #             circle = gdspy.Round((self.Ox+x-(N+1)*d/2, y+self.Oy), self.radius, tolerance=0.0005,layer = self.layer)
    #             taper.add(circle) 
        
    #     # Output taper  
    #     xx1_out,xx2_out,yy1_out,yy2_out,Nx1_out,Ny1_out,span_y_out,span_x_out = maille_tri(0,0,d,Nxt,self.Ny,Inv_xs = True)     
    #     for x in xx1_out:
    #         for y in yy1_out:
    #             circle = gdspy.Round((self.Ex+x+(N+1)*d/2, y+self.Oy), self.radius, tolerance=0.0005,layer = self.layer)
    #             taper.add(circle)
    #     for x in xx2_out:
    #         for y in yy2_out:
    #             circle = gdspy.Round((self.Ex+x+(N+1)*d/2, y+self.Oy), self.radius, tolerance=0.0005,layer = self.layer)
    #             taper.add(circle)    
         
    #     # Sets the tapers 
    #     tapers =  self.lib_PhC.new_cell("tapers"+name)
    #     taper_up = taper.copy("taper_up"+name,translation=(0,(span_y_out+self.width_PhC)/2))
    #     taper_down = taper.copy("taper_down"+name,translation=(0,-(span_y_out+self.width_PhC)/2))
    #     tapers.add(taper_up)
    #     tapers.add(taper_down) 
    #     self.lib_PhC.remove(taper)              

    #     # Actualize the injection points
    #     self.Ox = self.Ox - span_x_in - d/2
    #     self.Ex = self.Ex + span_x_out + d/2       
    #     self.length = abs(self.Ex-self.Ox)
        
    
    # def set_outputs(self,width,L_ino):
    #     """
    #     Parameters
    #     ----------
    #     width : TYPE
    #         DESCRIPTION.
    #     L_ino : TYPE
    #         DESCRIPTION.

    #     """
    #     # Create the waveguides
    #     distance = self.width_PhC+(self.etch_width-(self.width_PhC-width))
    #     wg_in = gdspy.Path(self.etch_width,(self.Ox,self.Oy),2,distance = self.width_PhC+self.etch_width)
    #     wg_out = gdspy.Path(self.etch_width,(self.Ex,self.Ey),2,distance = self.width_PhC+self.etch_width)       
    #     wg_in.segment(L_ino,direction='-x',final_distance=self.width+self.etch_width,layer = self.layer)
    #     wg_out.segment(L_ino,direction='+x',final_distance=self.width+self.etch_width,layer = self.layer)

    #     # Add to the Cell 
    #     strip_ino = self.lib_PhC.new_cell("strip_ino"+self.name)
    #     strip_ino.add(wg_in)
    #     strip_ino.add(wg_out)
        
    #     # Actualize the injection points
    #     self.Ox = self.Ox - L_ino
    #     self.Ex = self.Ex + L_ino
    #     self.length = abs(self.Ex-self.Ox)
 
# =============================================================================
# MZI    
# =============================================================================

class MZI(Straight_wg,Bend_guides,PhC_W1):
    
    def __init__(self,initial_point,potato1,potato2,etch_width = 3,layer = 1,**kwargs):
        """
        Parameters
        ----------
        initial_point : TYPE
            DESCRIPTION.
        potato1 : TYPE
            DESCRIPTION.
        potato2 : TYPE
            DESCRIPTION.

        Raises
        ------
        Exception
            DESCRIPTION.

        """
        self.etch_width = etch_width
        self.layer = layer  
        self.initial_point = initial_point
        self.potato1 = potato1
        self.potato2 = potato2
        self.Ox = self.initial_point[0]
        self.Oy = self.initial_point[1]
        self.width = self.potato1.width
        
        if self.potato1.Oy == self.potato2.Oy :
            raise Exception("Potatoes are on top of each other, should have diffrent 'y'")
        
        if self.potato1.Ox != self.potato2.Ox :
            raise Exception("Potatoes are not in line, should have the same 'x'")
            
        if abs(self.potato1.width-self.potato2.width)>0:
            raise Exception("Sorry, unmacthing waveguide widths") 
        
        self.space = abs(self.potato1.Oy-self.potato2.Oy)
        self.ShiftOx = min(self.potato1.Ox,self.potato2.Ox)#+1.4*self.space
        self.ShiftOy = min(self.potato1.Oy,self.potato2.Oy)#+self.space/2            
        self.Ex = self.Ox+2*4*self.space + max(self.potato1.length,self.potato2.length)
        self.Ey = self.Oy      
        self.length = self.Ex-self.Ox
        self.name = self.potato1.name + self.potato2.name + str(self.Ox) + str(self.Oy) + str(self.length)
        
    def add2cell(self,Cell):
        
        patch = gdspy.Cell("patch"+self.name)

        self.potato1.add2cell(patch)
        self.potato2.add2cell(patch)

        paste = patch.copy("paste"+self.name,translation=(-self.ShiftOx + self.Ox+4*self.space,
                                                -self.ShiftOy + self.Oy-self.space/2))
       
        if self.potato1.length > self.potato2.length:
            addon = gdspy.Path(self.etch_width,initial_point=(-self.ShiftOx + self.Ox+4*self.space + self.potato2.Ex,
                                                -self.ShiftOy + self.Oy-self.space/2 + self.potato2.Ey),
                               number_of_paths = 2,distance=self.etch_width+self.width )
            addon.segment(abs(self.potato1.length - self.potato2.length),layer = self.layer)
            Cell.add(addon)

        elif self.potato2.length > self.potato1.length:
            addon = gdspy.Path(self.etch_width,initial_point=(-self.ShiftOx + self.Ox+4*self.space + self.potato1.Ex,
                                                -self.ShiftOy + self.Oy-self.space/2 + self.potato1.Ey),
                               number_of_paths = 2,distance=self.etch_width+self.width )
            addon.segment(abs(self.potato1.length - self.potato2.length),layer = self.layer)
            Cell.add(addon)
            
        path_up_etch = gdspy.Path(2*self.etch_width+self.width,(self.Ox,self.Oy),number_of_paths=1)
        path_up_etch.bezier([(4*self.space/2,0),(4*self.space/2,self.space/2),(4*self.space,self.space/2)],layer = self.layer)
        
        path_down_etch = gdspy.Path(2*self.etch_width+self.width,(self.Ox,self.Oy),number_of_paths=1)
        path_down_etch.bezier([(4*self.space/2,0),(4*self.space/2,-self.space/2),(4*self.space,-self.space/2)],layer = self.layer)
        
        path_up = gdspy.Path(self.width,(self.Ox,self.Oy),number_of_paths=1)
        path_up.bezier([(4*self.space/2,0),(4*self.space/2,self.space/2),(4*self.space,self.space/2)],layer = self.layer)
        
        path_down = gdspy.Path(self.width,(self.Ox,self.Oy),number_of_paths=1)
        path_down.bezier([(4*self.space/2,0),(4*self.space/2,-self.space/2),(4*self.space,-self.space/2)],layer = self.layer)
        
        why_l = gdspy.boolean(gdspy.boolean(path_up, path_down, 'or',layer = self.layer), 
                            gdspy.boolean(path_up_etch, path_down_etch, 'or',layer = self.layer),'xor',layer = self.layer)
        
        why_r = gdspy.copy(why_l,self.Ox,self.Oy)#
        why_r.rotate(np.pi,(self.Ox,self.Oy))
        why_r.translate(self.Ex, self.Oy)
        
        Cell.add(why_l)
        Cell.add(why_r)
        Cell.add(paste) 
    
# =============================================================================
# in_n_out_waveguide
# =============================================================================
        

class In_n_out_waveguides_bend(Straight_wg,Round_spiral,Bend_guides,
                          FocusingGratingCoupler_SWG,PhC_W1,Simple_Ring):
    
    def __init__(self,inject,potato,Linject,Rinject,
                          Deltas = (100,20),**kwargs):
        """

        Parameters
        ----------
        inject : TYPE
            DESCRIPTION.
        potato : TYPE
            DESCRIPTION.
        Linject : TYPE
            DESCRIPTION.


        """
        
        self.inject = inject
        self.potato = potato     
        self.Linject = Linject
        self.Rinject = Rinject
        self.width = self.potato.width
        self.etch_width = self.potato.etch_width
        self.layer = self.potato.layer
        self.Deltas = Deltas
        self.Dx = self.Deltas[0]
        self.Dy = self.Deltas[1]
        
    def add2cell(self,Cell):
        global template

        path_in = gdspy.Path(self.etch_width ,(self.potato.Ox,self.potato.Oy),
                              number_of_paths=2,distance = self.etch_width  + self.potato.width)
        path_in.segment((self.Linject-self.Dx)/2,'-x',layer = self.layer)
        path_in.bezier([(-self.Dx/2,0),(-self.Dx/2,self.Dy),(-self.Dx,self.Dy)],layer = self.layer)
        path_in.segment((self.Linject-self.Dx)/2,'-x',layer = self.layer)
        Xtip = self.potato.Ox-self.Linject
        Ytip = self.potato.Oy+self.Dy
        Cell.add(path_in)
        rot = np.pi
        
        if self.inject.name not in lib.cells:   
            template = lib.new_cell('template_'+self.inject.name)  
            self.inject.add2cell(Cell,Xtip, Ytip, rot)
            self.inject.add2cell(template,0, 0, np.pi)
        else:
            Cell.add(gdspy.CellReference(template, (Xtip,Ytip)))
        
        path_out = gdspy.Path(self.etch_width ,(self.potato.Ex,self.potato.Ey),
                              number_of_paths=2,distance = self.etch_width  + self.potato.width)
        path_out.segment((self.Rinject-self.Dx)/2,'+x',layer = self.layer)
        path_out.bezier([(self.Dx/2,0),(self.Dx/2,-self.Dy),(self.Dx,-self.Dy)],layer = self.layer)
        path_out.segment((self.Rinject-self.Dx)/2,'+x',layer = self.layer)
        Xtip = self.potato.Ex+self.Rinject
        Ytip = self.potato.Ey-self.Dy
        Cell.add(path_out)
         
        if self.inject.name not in lib.cells: 
            template = lib.new_cell('template_'+self.inject.name)  
            self.inject.add2cell(Cell,Xtip, Ytip, rot = rot + np.pi)
            self.inject.add2cell(template,0, 0, rot = np.pi)
        else:
            Cell.add(gdspy.CellReference(template, (Xtip,Ytip), rotation = (rot)*180/np.pi ))



class In_n_out_waveguides(Straight_wg,Round_spiral,Bend_guides,
                          FocusingGratingCoupler_SWG,PhC_W1,Simple_Ring):
    
    def __init__(self,inject,potato,Linject,Rinject,**kwargs):
        """

        Parameters
        ----------
        inject : TYPE
            DESCRIPTION.
        potato : TYPE
            DESCRIPTION.
        Linject : TYPE
            DESCRIPTION.


        """
        
        self.inject = inject
        self.potato = potato     
        self.Linject = Linject
        self.Rinject = Rinject
        self.etch_width = self.potato.etch_width
        self.width = self.potato.width
        self.layer = self.potato.layer
        
    def add2cell(self,Cell):
        global template
        path_in = gdspy.Path(self.etch_width,(self.potato.Ox,self.potato.Oy),
                              number_of_paths=2,distance = self.etch_width + self.potato.width)
        path_in.segment(self.Linject,'-x',layer = self.layer)
        Xtip = self.potato.Ox-self.Linject
        Ytip = self.potato.Oy
        Cell.add(path_in)
        rot = np.pi
        
        if self.inject.name not in lib.cells:   
            template = lib.new_cell('template_'+self.inject.name)  
            self.inject.add2cell(Cell,Xtip, Ytip, rot)
            self.inject.add2cell(template,0, 0, np.pi)
        else:         
            Cell.add(gdspy.CellReference(template, (Xtip,Ytip)))
        
        path_out = gdspy.Path(self.etch_width,(self.potato.Ex,self.potato.Ey),number_of_paths=2,distance = self.etch_width + self.potato.width)
        path_out.segment(self.Rinject,'+x',layer = self.layer)
        Xtip = self.potato.Ex+self.Rinject
        Ytip = self.potato.Ey
        Cell.add(path_out)
        
        if self.inject.name not in lib.cells: 
            template = lib.new_cell('template_'+self.inject.name)  
            self.inject.add2cell(Cell,Xtip, Ytip, rot = rot + np.pi)
            self.inject.add2cell(template,0, 0, rot = np.pi)
        else:
            Cell.add(gdspy.CellReference(template, (Xtip,Ytip), rotation = (rot)*180/np.pi ))


# =============================================================================
# Oh hi Marc ! (allignment)
# =============================================================================

class MarkAlign():
    
    def __init__(self,X,Y,layer = 1,**kwargs):
        
        self.X = X
        self.Y = Y
        self.layer = layer
        
    def add2cell(self,Cell):
        
        R1 = gdspy.Rectangle((self.X-150,self.Y-150), (self.X-50,self.Y+150))
        R2 = gdspy.Rectangle((self.X+50,self.Y-150), (self.X+150,self.Y+150))
        
        Cell.add(R1)
        Cell.add(R2)

        # R3 = gdspy.Rectangle((self.X-150,self.Ys[1]-150), (self.X-50,self.Ys[1]+150))
        # R4 = gdspy.Rectangle((self.X+50,self.Ys[1]-150), (self.X+150,self.Ys[1]+150))

        # Cell.add(R3)
        # Cell.add(R4)

        

class MarkAlign_double():
    
    def __init__(self,X,Ys,layer = 1,**kwargs):
        
        self.X = X
        self.Ys = Ys
        self.layer = layer
        
    def add2cell(self,Cell):
        
        R1 = gdspy.Rectangle((self.X-150,self.Ys[0]-150), (self.X-50,self.Ys[0]+150))
        R2 = gdspy.Rectangle((self.X+50,self.Ys[0]-150), (self.X+150,self.Ys[0]+150))
        
        Cell.add(R1)
        Cell.add(R2)

        R3 = gdspy.Rectangle((self.X-150,self.Ys[1]-150), (self.X-50,self.Ys[1]+150))
        R4 = gdspy.Rectangle((self.X+50,self.Ys[1]-150), (self.X+150,self.Ys[1]+150))

        Cell.add(R3)
        Cell.add(R4)

# =============================================================================
# Control Etch
# =============================================================================
    
class TestEtch():
    def __init__(self,initial_point,d,radius,Nx,Ny,coeff = 0.25,layer = 1):
        
        self.Ox = initial_point[0]
        self.Oy = initial_point[1]
        self.d = d
        self.radius = radius
        self.Nx = Nx
        self.Ny = Ny
        self.coeff =  coeff        
        self.layer = layer
        
        self.xx = self.d*np.arange(1,Nx+1,1) + self.Ox - self.Nx*self.d/2
        self.yy = self.d*np.arange(1,Ny+1,1) + self.Oy - self.Ny*self.d/2
        
    def add2cell(self,Cell):
                        
        shift = 0
        
        for y in self.yy:                
            shift += self.coeff * self.d
            for x in self.xx:
                circle = gdspy.Round((x+shift,y), self.radius, tolerance=0.0005,layer = self.layer)
                Cell.add(circle)
        
    
    
