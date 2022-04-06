# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 10:57:32 2021

@author: pelgrin-adm
"""

import os
import numpy as np 
import gdspy 
import matplotlib.pyplot as plt

from gds_library import *




wref= 0.40
width = wref
a = 0.350 

start = 1800
end = 2900
taper_to_wg = 300
multimode = 2200#1800

# Gratings
Lsi         = 0.355
Lgap        = 0.345
Wsi         = 0.3
Wgap        = 0.15
# gc = FocusingGratingCoupler_SWG(Lsi, Lgap, Wsi, Wgap, Nperiod = 60,
#                                 FocalLength = 25.6, GrArc = (0.27-0.019),
#                                 TaperArc = (0.234-0.018), theta_in = 20 * np.pi / 180,
#                                 nc = 1.0, nSlab = 3.05,lambda_o = 1.55, W_EndTaper = 0.5,
#                                  W_wgInterconnect = width, NPoints = 3,ID = ID)    

gc = ButtCouplingTaper(width,L_taper_to_wg = taper_to_wg,Final_width = 3,L_multimode = multimode,
             etch_width = 3,layer = 1, ID = 'Butt_taper')

# =============================================================================
# Rings
# =============================================================================

def rings(Background,xstart,ystart,width,Gaps,R,Cell,**kwargs):
    """
    Parameters
    ----------
    Background : TYPE
        DESCRIPTION.
    xstart : TYPE
        DESCRIPTION.
    ystart : TYPE
        DESCRIPTION.
    width : TYPE
        DESCRIPTION.
    Gaps : TYPE
        DESCRIPTION.
    R : TYPE
        DESCRIPTION.
    Cell : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    Background : TYPE
        DESCRIPTION.
    ystart : TYPE
        DESCRIPTION.

    """
    
    
    for ii,gap in enumerate(Gaps):
        
        ud = (-1)**(ii+1)
        center = (xstart+4*R+ud*3*R,ystart)
        # ring = Simple_Ring(width,center,R,gap,updown = ud ,layer = 1,**kwargs)
        ring = Simple_Ring(width,center,R,gap,updown = ud ,**kwargs)
        ring.add2cell(Cell)
                
        ralong = abs(xstart+end-ring.Ex)
        ralong2 = start+ring.Ox-xstart
        
        gc = ButtCouplingTaper(ring.width,L_taper_to_wg = taper_to_wg,
                               L_multimode = multimode,
                               ID = 'Butt_taper',**kwargs)
        
        innout = In_n_out_waveguides(gc, ring, start+ring.Ox-xstart,ralong,**kwargs)#1450-shift_x)
        innout.add2cell(Cell) 
        poly1 = gdspy.Rectangle((center[0]-R-20,center[1]-R-20),
                              (center[0]+R+20,center[1]+R+20),layer = Background)
        Background += 1 
        
        poly2 = gdspy.Rectangle((center[0]-ralong2-R-20-taper_to_wg-multimode,ring.Oy-11),
                              (center[0]-R-20,ring.Oy+11),layer = Background)
        Background += 1 
        
        poly3 = gdspy.Rectangle((ring.Ex+10,ring.Ey-10),
                      (ring.Ex+taper_to_wg+multimode+10+ralong,ring.Ey+11),layer = Background)
        Background += 1 
       
        Cell.add(poly1)
        Cell.add(poly2)
        Cell.add(poly3) 
        
        title = "R-"+str(R)+"_Gap-"+str(round(1000*gap)/1000)+"_W-"+str(ring.width)
        text = gdspy.Text(title, 30, (xstart-1500,ring.Oy - 15 - ud * (30 + width + ring.etch_width)))
        Cell.add(text)

        text2 = gdspy.Text(title, 30, (xstart+2000,ring.Oy - 15 - ud * (30 + width + ring.etch_width)))
        Cell.add(text2)
         
        ystart += (150+ud*100)
        
    return Background,ystart

# =============================================================================
# W1
# =============================================================================

def W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L,Cell):
    
    """
    Parameters
    ----------
    Background : TYPE
        DESCRIPTION.
    xstart : TYPE
        DESCRIPTION.
    ystart : TYPE
        DESCRIPTION.
    radius : TYPE
        DESCRIPTION.
    a : TYPE
        DESCRIPTION.
    ap : TYPE
        DESCRIPTION.
    pp : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    Cell : TYPE
        DESCRIPTION.

    Returns
    -------
    Background : TYPE
        DESCRIPTION.
    y_final : TYPE
        DESCRIPTION.

    """
    

    d  = a
    
    width_PhC = np.sqrt(3)*d
    mono_width = width
    L_ino = 55
    PhC = PhC_W1(width,(xstart,ystart),L,a,radius,)
    
    title = "Taper_study_"
    
    for ii in range (len(ap)):
        d = ap[ii]
        p = pp[ii] 
        PhC.set_taper(d,p)
        title = title + "T-" + str(ii) + "_d-" + str(d)+ "-p-" + str(p) + "__"

    PhC.set_outputs(mono_width,L_ino)
    PhC.add2cell(Cell)
    device = PhC
    ralong = abs(xstart+end-PhC.Ex)
    innout = In_n_out_waveguides(gc, device, start+PhC.Ox-xstart,ralong)#1500-PhC.length)
    innout.add2cell(Cell)
    y_final = ystart + 50
        
    xmin = xstart-20 - 2*gc.L_taper_to_wg - gc.L_multimode
    xmax = xstart + end + 20 + gc.L_taper_to_wg + gc.L_multimode
    ymin = ystart -10
    ymax = ystart +10
    Background = add_zero_rectangle(Cell,xmin,xmax,ymin,ymax,Background = Background)
           
    # title = ''.join(['taper-','_PhC_','L-',str(round(PhC.length))])
    text = gdspy.Text(title, 10, (xstart-1000,ystart +6))
    Cell.add(text)
        
    start_taper = PhC.Ox
    
    return Background,y_final


# =============================================================================
# Spiral
# =============================================================================
def spirals(Background,xstart,ystart,RR,Cell,**kwargs):  
    """
    

    Parameters
    ----------
    Background : TYPE
        DESCRIPTION.
    xstart : TYPE
        DESCRIPTION.
    ystart : TYPE
        DESCRIPTION.
    RR : TYPE
        DESCRIPTION.
    Cell : TYPE
        DESCRIPTION.

    Returns
    -------
    Background : TYPE
        DESCRIPTION.
    y_final : TYPE
        DESCRIPTION.

    """
    
    
    shift_x = 0
    shift_y = 0
    
    y_init = ystart
    y_final = ystart
    
    for rr in RR:
        
        spi = Round_spiral(width, (xstart+shift_x, 50+y_init+rr*250), roundt = rr,rot = np.pi/2,**kwargs)
        S = spi.adjust_output(Cell, 500, 4*spi.turn_radius+4*spi.etch_width)
        spi.add2cell(Cell)
        device = spi
        ralong = abs(xstart+end-spi.Ex)

        gc = ButtCouplingTaper(width,L_taper_to_wg = taper_to_wg,
                       L_multimode = multimode,
                       ID = 'Butt_taper',**kwargs)

        innout = In_n_out_waveguides(gc, device, start+spi.Ox-xstart,ralong,**kwargs)#1450-shift_x)
        innout.add2cell(Cell) 
        
 
        # title = ''.join(['spiral-','W-',str(round(1000*width)/1000),'R-',str(rr),'L-'])
        # text = gdspy.Text(title, 30, (xstart+2500,  y_init+rr*65 + 6))
        # Cell.add(text)
        
        poly1 = gdspy.Rectangle((spi.initial_point[0]-(spi.size)-10,spi.initial_point[1]-2*spi.size-11),
                              (spi.initial_point[0]+(spi.size)+10,spi.initial_point[1]+11),layer = Background)
        Background += 1 
        
        poly2 = gdspy.Rectangle((-start-20 - gc.L_taper_to_wg - gc.L_multimode-shift_x+spi.Ox,spi.Oy-9),
                              (spi.initial_point[0]-spi.size-10,spi.initial_point[1]+9),layer = Background)
        Background += 1 
        
        poly3 = gdspy.Rectangle((spi.initial_point[0]+spi.size+10,spi.initial_point[1]-2*spi.size-9),
                      (spi.initial_point[0]+ end + 20 + gc.L_taper_to_wg + gc.L_multimode - shift_x,spi.initial_point[1]-2*spi.size+9+S),layer = Background)
        Background += 1 
        

        Cell.add(poly1)
        Cell.add(poly2)
        Cell.add(poly3)
      
        shift_x += (spi.space+ 2 * spi.turn_radius) + 200
        shift_y += -(spi.space+ 2 * spi.turn_radius)/2  
        
    y_final = y_init +50+ rr*(250)
    
    return Background,y_final




def str_waveguide(Background,xstart,ystart,L,Cell,**kwargs):
    strwg = Straight_wg(width,(xstart,ystart),L,'+x',**kwargs)
    strwg.add2cell(Cell)
    
    device = strwg
    ralong = abs(xstart+end-strwg.Ex)
    
    gc = ButtCouplingTaper(width,L_taper_to_wg = taper_to_wg,
               L_multimode = multimode,
               ID = 'Butt_taper',**kwargs)

    innout = In_n_out_waveguides(gc, device, start,ralong,**kwargs)#350)
    innout.add2cell(Cell) 
    # print(strwg.name)
    y_final = ystart + 50
    
    xmin = xstart-20 - gc.L_taper_to_wg - gc.L_multimode-start
    xmax = xstart + end + 20 + gc.L_taper_to_wg + gc.L_multimode
    ymin = ystart -10
    ymax = ystart +10
    Background = add_zero_rectangle(Cell,xmin,xmax,ymin,ymax,Background = Background)
    
    return Background,y_final    

# Geometry must be placed in cells.
cell1 = lib.new_cell('ref')

name = "SiN_sample_1"
lib.write_gds(name+'.gds')

EW = 5#3#
FW = 4# 3#
# A #############################################################################
xstart = 0 
ystart = 0
Background = 100
width = 1.0#0.4#
L = 1000
Background,ystart = str_waveguide(Background,xstart,ystart,L,cell1,etch_width = EW,layer = 1,Final_width = FW)
L = 2000
Background,ystart = str_waveguide(Background,xstart,ystart+100,L,cell1,etch_width = EW,layer = 1,Final_width = FW)
L = 2000
Background,ystart = str_waveguide(Background,xstart,ystart+100,L,cell1,etch_width = EW,layer = 1,Final_width = FW)

ystart += 100
R = 100#20#
Gaps =  np.arange(0.3,1.1,0.02)#np.arange(0.1,0.350,0.005)#
Background,ystart = rings(Background,xstart,ystart+100,width,Gaps,R,cell1,etch_width = EW,layer = 1,Final_width = FW)

RR = [1,2,3,4,5,6]
Background,ystart = spirals(Background,xstart,ystart+100,RR,cell1,turn_radius = 40,etch_width = EW,layer = 1,Final_width = FW)

ymark = ystart+400
Mark1 = MarkAlign_double(-3500,(-400,ymark),layer = 2)
Mark1.add2cell(cell1)
Background = add_zero_rectangle(cell1,Mark1.X - 175,Mark1.X + 175,Mark1.Ys[0] - 175,Mark1.Ys[0] + 175,Background = Background)
Background = add_zero_rectangle(cell1,Mark1.X - 175,Mark1.X + 175,Mark1.Ys[1] - 175,Mark1.Ys[1] + 175,Background = Background)


# B #############################################################################
xstart = 8000
ystart = 75
width = 1.5#0.5#
L = 1000
Background,ystart = str_waveguide(Background,xstart,ystart,L,cell1,etch_width = EW,layer = 1,Final_width = FW)
L = 2000
Background,ystart = str_waveguide(Background,xstart,ystart+100,L,cell1,etch_width = EW,layer = 1,Final_width = FW)
L = 2000
Background,ystart = str_waveguide(Background,xstart,ystart+100,L,cell1,etch_width = EW,layer = 1,Final_width = FW)

ystart += 100

R = 100#20#
Background,ystart = rings(Background,xstart,ystart+175,width,Gaps,R,cell1,etch_width = EW,layer = 1,Final_width = FW)

RR = [1,2,3,4,5,6]
Background,ystart = spirals(Background,xstart,ystart-175,RR,cell1,turn_radius = 40,etch_width = EW,layer = 1,Final_width = FW)

Mark2 = MarkAlign_double(xstart-3500,(-400,ymark),layer = 2)
Mark2.add2cell(cell1)
Background = add_zero_rectangle(cell1,Mark2.X - 175,Mark2.X + 175,Mark2.Ys[0] - 175,Mark2.Ys[0] + 175,Background = Background)
Background = add_zero_rectangle(cell1,Mark2.X - 175,Mark2.X + 175,Mark2.Ys[1] - 175,Mark2.Ys[1] + 175,Background = Background)


# C #############################################################################
xstart = 16000
ystart = 0
width = 2.0#0.6#
L = 1000
Background,ystart = str_waveguide(Background,xstart,ystart,L,cell1,etch_width = EW,layer = 1,Final_width = FW)
L = 2000
Background,ystart = str_waveguide(Background,xstart,ystart+100,L,cell1,etch_width = EW,layer = 1,Final_width = FW)
L = 2000
Background,ystart = str_waveguide(Background,xstart,ystart+100,L,cell1,etch_width = EW,layer = 1,Final_width = FW)

ystart += 100
R = 100#20#
Background,ystart = rings(Background,xstart,ystart+100,width,Gaps,R,cell1,etch_width = EW,layer = 1,Final_width = FW)

RR = [1,2,3,4,5,6]
Background,ystart = spirals(Background,xstart,ystart-150,RR,cell1,turn_radius = 40,etch_width = EW,layer = 1,Final_width = FW)

Mark3 = MarkAlign_double(xstart-3500,(-400,ymark),layer = 2)
Mark3.add2cell(cell1)
Background = add_zero_rectangle(cell1,Mark3.X - 175,Mark3.X + 175,Mark3.Ys[0] - 175,Mark3.Ys[0] + 175,Background = Background)
Background = add_zero_rectangle(cell1,Mark3.X - 175,Mark3.X + 175,Mark3.Ys[1] - 175,Mark3.Ys[1] + 175,Background = Background)


Mark4 = MarkAlign_double(xstart+4500,(-400,ymark),layer = 2)
Mark4.add2cell(cell1)
Background = add_zero_rectangle(cell1,Mark4.X - 175,Mark4.X + 175,Mark4.Ys[0] - 175,Mark4.Ys[0] + 175,Background = Background)
Background = add_zero_rectangle(cell1,Mark4.X - 175,Mark4.X + 175,Mark4.Ys[1] - 175,Mark4.Ys[1] + 175,Background = Background)



lib.write_gds(name+'.gds')
