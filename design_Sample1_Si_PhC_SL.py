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
# Lsi         = 0.355
# Lgap        = 0.345
# Wsi         = 0.3
# Wgap        = 0.15
# gc = FocusingGratingCoupler_SWG(Lsi, Lgap, Wsi, Wgap, Nperiod = 60,
#                                 FocalLength = 25.6, GrArc = (0.27-0.019),
#                                 TaperArc = (0.234-0.018), theta_in = 20 * np.pi / 180,
#                                 nc = 1.0, nSlab = 3.05,lambda_o = 1.55, W_EndTaper = 0.5,
#                                  W_wgInterconnect = width, NPoints = 3,ID = ID)    

# =============================================================================
# 
# =============================================================================

def full_cell(width,a,radius,Cell,xstart,ystart,Background = 100,**kwargs):
   
    LL = [1400,1400]
    
    y_init = ystart
    
    for ii,L in enumerate(LL): 
        Background,ystart = str_waveguide(Background,xstart,ystart,L,cell1)
    
    ### 1 W1
    ## 150
    PhC150 = PhC_W1(width,(0,0),100,a,radius,**kwargs)
    # PhC150.set_taper(a, 1)
    PhC150.set_outputs(width,50)
    strwg1 = Straight_wg(width,(PhC150.Ox,20),20,direction = '+x',layer = 1)
    Background,y_final = set_MZI(Background,PhC150,strwg1,xstart,ystart,cell1,etch_width = 3) 
    
    # 250
    PhC250 = PhC_W1(width,(0,0),250,a,radius,**kwargs)
    PhC250.set_outputs(width,50)
    strwg1 = Straight_wg(width,(PhC250.Ox,20),20,direction = '+x',layer = 1)
    ystart = y_final
    Background,y_final = set_MZI(Background,PhC250,strwg1,xstart,ystart,cell1,etch_width = 3) 
    
    ## 350
    PhC350 = PhC_W1(width,(0,0),350,a,radius,**kwargs)
    PhC350.set_outputs(width,50)
    strwg1 = Straight_wg(width,(PhC350.Ox,20),20,direction = '+x',layer = 1)
    ystart = y_final
    Background,y_final = set_MZI(Background,PhC350,strwg1,xstart,ystart,cell1,etch_width = 3) 
    
    ### 2 W1
    ## 150
    ystart = y_final
    PhC50 = PhC_W1(width,(0,20),50,a,radius,**kwargs)
    PhC50.set_outputs(width,50)
    PhC150 = PhC_W1(width,(0,0),150,a,radius,**kwargs)
    PhC150.set_outputs(width,50)
    Background,y_final = set_MZI(Background,PhC150,PhC50,xstart,ystart,cell1,etch_width = 3) 
    ## 250
    
    ystart = y_final
    PhC50 = PhC_W1(width,(0,20),50,a,radius,**kwargs)
    PhC50.set_outputs(width,50)
    PhC250 = PhC_W1(width,(0,0),250,a,radius,**kwargs)
    PhC250.set_outputs(width,50)
    Background,y_final = set_MZI(Background,PhC250,PhC50,xstart,ystart,cell1,etch_width = 3) 
    ## 350
    
    ystart = y_final
    PhC50 = PhC_W1(width,(0,20),50,a,radius,**kwargs)
    PhC50.set_outputs(width,50)
    PhC350 = PhC_W1(width,(0,0),350,a,radius,**kwargs)
    PhC350.set_outputs(width,50)
    Background,ystart = set_MZI(Background,PhC350,PhC50,xstart,ystart,cell1,etch_width = 3) 
    L = [100,150,200,250,800,850,900]
    Background,y_final = W1_PhC(Background,xstart,ystart,L,cell1,**kwargs)
  
    return Background,y_final
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

def W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L,Cell,**kwargs):
    
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
       
    gc = ButtCouplingTaper(width,L_taper_to_wg = taper_to_wg,
                   L_multimode = multimode,
                   ID = 'Butt_taper',**kwargs)
    
    innout = In_n_out_waveguides(gc, device, start+PhC.Ox-xstart,ralong)#1500-PhC.length)
    innout.add2cell(Cell)
    y_final = ystart + 50
        
    xmin = xstart-20 - gc.L_taper_to_wg - gc.L_multimode-start
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


def W1_PhC(Background,xstart,ystart,L,Cell,**kwargs):
    for ii, l in enumerate(L):
        
        # Construction of the W1 waveguide and printing on the cell
        PhC = PhC_W1(width,(xstart,ystart + ii*50),l,a,radius,**kwargs)
        try:
            PhC.set_taper(d2,20)
        except NameError:
            pass
        try:
            PhC.set_taper(d3,10) 
        except NameError:
            pass
            
        PhC.set_outputs(width)
        PhC.add2cell(Cell)
    
        # Addition of the injection device
        gc = ButtCouplingTaper(width,L_taper_to_wg = taper_to_wg,
               L_multimode = multimode,
               ID = 'Butt_taper',**kwargs)
        
        innout = In_n_out_waveguides(gc, PhC, start+PhC.Ox-xstart,abs(xstart+end-PhC.Ex))#1500-PhC.length)
        innout.add2cell(Cell)
        xmin = xstart-20 - gc.L_taper_to_wg - gc.L_multimode-start
        xmax = xstart + end + 20 + gc.L_taper_to_wg + gc.L_multimode
        ymin = ystart +  ii*50-10
        ymax = ystart +  ii*50+10
        Background = add_zero_rectangle(Cell,xmin,xmax,ymin,ymax,Background = Background)
           
        title = ''.join(['taper-','W-',str(round(1000*PhC.width_PhC)/1000),' -W-',str(round(1000*width)/1000),'L-',str(round(PhC.length))])
        text = gdspy.Text(title, 10, (xstart-1000,ystart + ii*50+6))
        Cell.add(text)
        
        start_taper = PhC.Ox
        y_final = ystart + ii*50
        
    return Background,y_final  

def set_MZI(Background,potato1,potato2,xstart,ystart,Cell,**kwargs): 
    
    # # MZI 1 ___________________________________________________________________

    mzi1 = MZI((xstart,ystart),potato1,potato2,**kwargs)
    mzi1.add2cell(Cell)
    
    gc = ButtCouplingTaper(width,L_taper_to_wg = taper_to_wg,
               L_multimode = multimode,
               ID = 'Butt_taper',**kwargs)
    ralong = abs(xstart+end-mzi1.Ex)
    innout = In_n_out_waveguides(gc, mzi1, start+mzi1.Ox-xstart, abs(xstart+end-mzi1.Ex))#1450-mzi1.length)
    innout.add2cell(Cell)
      
    # Background cell
    poly1 = gdspy.Rectangle((mzi1.Ox-10,mzi1.Oy-15),
                          (mzi1.Ex+10,mzi1.Ey+15),layer = Background)
    Background += 1   
    poly2 = gdspy.Rectangle((xstart-10 - start -gc.L_taper_to_wg - gc.L_multimode,mzi1.Oy-10),
                          (mzi1.Ox-10,mzi1.Oy+10),layer = Background)
    Background += 1 
    poly3 = gdspy.Rectangle((mzi1.Ex+10,mzi1.Ey-10),
                          (mzi1.Ex+10+taper_to_wg+multimode+ralong,mzi1.Ey+10),layer = Background)
    Background += 1 

    Cell.add(poly1)
    Cell.add(poly2)
    Cell.add(poly3)
    
    y_final = ystart + 50
    return Background,y_final  
    

        
## W1 w/o taper ###############################################################


# Geometry must be placed in cells.
cell1 = lib.new_cell('ref')

name = "SL_sample1"
lib.write_gds(name+'.gds')

Background = 100



xstart = 0 
ystart = 0


radius = 0.090
a = 0.360
width = 0.45

Background, ystart = full_cell(width, a, radius, cell1, xstart, ystart,r1 =  0.84, s1=0.12)
ystart += 150
Background, ystart = full_cell(width, a, radius, cell1, xstart, ystart,r1 =  0.86, s1=0.12)
ystart += 150
Background, ystart = full_cell(width, a, radius, cell1, xstart, ystart,r1 =  0.88, s1=0.12)
ystart += 150
Background, ystart = full_cell(width, a, radius, cell1, xstart, ystart,r1 =  0.90, s1=0.12)
ystart += 150
Background, ystart = full_cell(width, a, radius, cell1, xstart, ystart)


radius = 0.080
a = 0.360
width = 0.45

ystart += 150
Background, ystart = full_cell(width, a, radius, cell1, xstart, ystart,r1 =  0.92, s1=0.16)
ystart += 150
Background, ystart = full_cell(width, a, radius, cell1, xstart, ystart,r1 =  0.94, s1=0.16)
ystart += 150
Background, ystart = full_cell(width, a, radius, cell1, xstart, ystart,r1 =  0.96, s1=0.16)
ystart += 150
Background, ystart = full_cell(width, a, radius, cell1, xstart, ystart)

RR = [1,2,3,4,5,6]
Background,ystart = spirals(Background,xstart,ystart+175,RR,cell1,turn_radius = 40,etch_width = 3,layer = 1)

ymark = ystart+400
Mark1 = MarkAlign_double(-3500,(-400,ymark),layer = 2)
Mark1.add2cell(cell1)
Background = add_zero_rectangle(cell1,Mark1.X - 175,Mark1.X + 175,Mark1.Ys[0] - 175,Mark1.Ys[0] + 175,Background = Background)
Background = add_zero_rectangle(cell1,Mark1.X - 175,Mark1.X + 175,Mark1.Ys[1] - 175,Mark1.Ys[1] + 175,Background = Background)

lib.write_gds(name+'.gds')

# # # =============================================================================
# # # Taper Study
# # # =============================================================================

# Mark2 = MarkAlign(xstich- start -gc.L_taper_to_wg - gc.L_multimode -950/2,(-500,y_final + 450),layer = 2)
# Mark2.add2cell(cell1)
# Background = add_zero_rectangle(cell1,Mark2.X - 175,Mark2.X + 175,Mark2.Ys[0] - 175,Mark2.Ys[0] + 175,Background = Background)
# Background = add_zero_rectangle(cell1,Mark2.X - 175,Mark2.X + 175,Mark2.Ys[1] - 175,Mark2.Ys[1] + 175,Background = Background)


# test = TestEtch((xstich- start -gc.L_taper_to_wg - gc.L_multimode -950/2, y_final + 200),a,0.09,Nx = 800,Ny = 50,coeff = 0.1,layer = 3)
# test.add2cell(cell1)

# xmin = test.Ox- 20 - test.Nx*test.d/2
# xmax = test.Ox+ 20 + test.Nx*test.d/2
# ymin = test.Oy - test.Ny*test.d/2-20
# ymax = test.Oy + test.Ny*test.d/2+20
# Background = add_zero_rectangle(cell1,xmin,xmax,ymin,ymax,Background = Background)


a = 0.360
radius = 0.090
L1 = 100
L2 = 200
L3 = 700
L4 = 800
xstart = 8000
ystart = -25

Mark2 = MarkAlign_double(xstart-3500,(-400,ymark),layer = 2)
Mark2.add2cell(cell1)
Background = add_zero_rectangle(cell1,Mark2.X - 175,Mark2.X + 175,Mark2.Ys[0] - 175,Mark2.Ys[0] + 175,Background = Background)
Background = add_zero_rectangle(cell1,Mark2.X - 175,Mark2.X + 175,Mark2.Ys[1] - 175,Mark2.Ys[1] + 175,Background = Background)

Mark4 = MarkAlign_double(xstart+4500,(-400,ymark),layer = 2)
Mark4.add2cell(cell1)
Background = add_zero_rectangle(cell1,Mark4.X - 175,Mark4.X + 175,Mark4.Ys[0] - 175,Mark4.Ys[0] + 175,Background = Background)
Background = add_zero_rectangle(cell1,Mark4.X - 175,Mark4.X + 175,Mark4.Ys[1] - 175,Mark4.Ys[1] + 175,Background = Background)


####

ap = [a+0.010,a+0.020]
pp = [5,5]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.010,a+0.020]
pp = [10,10]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.010,a+0.020]
pp = [15,15]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.010,a+0.020]
pp = [18,18]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 200

####

ap = [a+0.010,a+0.020,a+0.030]
pp = [5,5,5]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.010,a+0.020,a+0.030]
pp = [10,10,10]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.005,a+0.010,a+0.015,a+0.020]
pp = [10,10,10,10]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.010,a+0.020,a+0.030]
pp = [18,18,18]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 200

####

ap = [a+0.005,a+0.010]
pp = [6,6]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.005,a+0.010]
pp = [10,10]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.005,a+0.010]
pp = [15,15]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.005,a+0.010]
pp = [18,18]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 200

####

ap = [a+0.005,a+0.010,a+0.015]
pp = [8,8,8]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.005,a+0.010,a+0.015]
pp = [10,10,10]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.005,a+0.010,a+0.015]
pp = [15,15,15]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.005,a+0.010,a+0.015]
pp = [18,18,18]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 200

####

ap = [a+0.005,a+0.010,a+0.015,a+0.020]
pp = [8,8,8,8]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.005,a+0.010,a+0.015,a+0.020]
pp = [10,10,10,10]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.005,a+0.010,a+0.015,a+0.020]
pp = [15,15,15,15]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 100

ap = [a+0.005,a+0.010,a+0.015,a+0.020]
pp = [18,18,18,18]
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L1,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L2,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L3,cell1)
Background,ystart = W1_taper_test(Background,xstart,ystart,radius,a,ap,pp,L4,cell1)

ystart += 50

L = [L1,L2,L3,L4]

## a sweep
a = 0.350
r = 0.090
Background,ystart= W1_PhC(Background,xstart,ystart,L,cell1)
ystart += 100

a = 0.355
r = 0.090
Background,ystart= W1_PhC(Background,xstart,ystart,L,cell1)
ystart += 100

a = 0.360
r = 0.090
Background,ystart= W1_PhC(Background,xstart,ystart,L,cell1)
ystart += 100

a = 0.365
r = 0.090
Background,ystart= W1_PhC(Background,xstart,ystart,L,cell1)
ystart += 100

a = 0.370
r = 0.090
Background,ystart= W1_PhC(Background,xstart,ystart,L,cell1)
ystart += 100

## R sweep
a = 0.360
r = 0.086
Background,ystart= W1_PhC(Background,xstart,ystart,L,cell1)
ystart += 100

a = 0.360
r = 0.088
Background,ystart= W1_PhC(Background,xstart,ystart,L,cell1)
ystart += 100

a = 0.360
r = 0.090
Background,ystart= W1_PhC(Background,xstart,ystart,L,cell1)
ystart += 100

a = 0.360
r = 0.092
Background,ystart= W1_PhC(Background,xstart,ystart,L,cell1)
ystart += 100

a = 0.360
r = 0.094
Background,ystart= W1_PhC(Background,xstart,ystart,L,cell1)
ystart += 100
####



lib.write_gds(name+'.gds')

# # R = 100
# # Gaps = np.arange(0.1,0.29,0.01)

# # Background,ystart = rings(Background,xstart,ystart,Gaps,R,cell1)

# xmax = xstart + end  + gc.L_taper_to_wg + gc.L_multimode
# Mark3 = MarkAlign(xmax - 950/2 ,(-500,y_final + 450),layer = 2)
# Mark3.add2cell(cell1)
# Background = add_zero_rectangle(cell1,Mark3.X - 175,Mark3.X + 175,Mark3.Ys[0] - 175,Mark3.Ys[0] + 175,Background = Background)
# Background = add_zero_rectangle(cell1,Mark3.X - 175,Mark3.X + 175,Mark3.Ys[1] - 175,Mark3.Ys[1] + 175,Background = Background)

# test = TestEtch((xmax-950/2, y_final + 200),a,radius,Nx = 800,Ny = 50,coeff = 0.1,layer = 3)
# test.add2cell(cell1)

# xmin = test.Ox- 20 - test.Nx*test.d/2
# xmax = test.Ox+ 20 + test.Nx*test.d/2
# ymin = test.Oy - test.Ny*test.d/2 -20
# ymax = test.Oy + test.Ny*test.d/2 +20
# Background = add_zero_rectangle(cell1,xmin,xmax,ymin,ymax,Background = Background)

# lib.write_gds(name+'.gds')
