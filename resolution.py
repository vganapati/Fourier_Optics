#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 11:30:24 2021

@author: vganapa1
"""

from functions import make_meshgrid, F, Ft, make_circle, plot
from propagators import NAfilter
import numpy as np
import matplotlib.pyplot as plt
####
### 4F imaging system and impulse response ###
####

# units are microns
wavelength = .6
focal_length = 1000
dx = 0.3
dy = 0.3
Lx_max = 100
Ly_max = 100
              
x,y,xm,ym,extent,F_extent = make_meshgrid(dx = dx,
                                          dy = dy,
                                          Lx_max = Lx_max,
                                          Ly_max = Ly_max)    

dots = np.zeros_like(xm)

dots[len(x)//2-5,len(y)//2]=1
dots[len(x)//2+5,len(y)//2]=1

# first lens
output_field = F(dots)
extent2 = np.array(F_extent)*wavelength*focal_length
dx2 = extent2[1]*2/output_field.shape[0]
dy2 = extent2[3]*2/output_field.shape[1]

# circular aperture in Fourier plane multiplies the field
aperture, extent, F_mat2D, F_extent = make_circle(dx = dx2,
                                                  dy = dy2,
                                                  Lx_max = extent2[1],
                                                  Ly_max = extent2[3],
                                                  a = .005, # stretch in x
                                                  b=.005, # stretch in y
                                                  )
plot(aperture,extent)
output_field = output_field*aperture

# second lens
output_field = F(output_field)
extent2 = np.array(F_extent)*wavelength*focal_length

plot(output_field,extent2)


####
# Imaging system without aberrations, circular OTF with radius defined by 
# numerical aperture (NA), magnification M
####

# units are microns
wavelength = .6
dx = 0.1
dy = 0.1
Lx_max = 100
Ly_max = 100
            
NA = 0.1
M = 1.0
  
x,y,xm,ym,extent,F_extent = make_meshgrid(dx = dx,
                                          dy = dy,
                                          Lx_max = Lx_max,
                                          Ly_max = Ly_max)    

dots = np.zeros_like(xm)

dots[len(x)//2-50,len(y)//2]=1
dots[len(x)//2+50,len(y)//2]=1


 
OTF = NAfilter(len(x), len(y), Lx_max*2, Ly_max*2, wavelength, NA)
plt.imshow(OTF)
output_field = Ft(F(dots)*OTF)

plot(output_field,np.array(extent)*M)


# coherent light and resolution
plot(np.abs(output_field)**2,np.array(extent)*M)


# incoherent light and resolution
OTF_incoherent = F(np.abs(Ft(OTF))**2)
output_field = Ft(F(dots)*OTF_incoherent)
plot(np.abs(output_field)**2,np.array(extent)*M)

# space invariant aberrations can be represented in the OTF