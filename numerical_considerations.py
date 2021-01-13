#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 18:24:53 2021

@author: vganapati

Numerical considerations in applying the scalar or Fresnel propagation kernel
are zero padding and sampling in real space (i.e. sampling and zero padding 
respectively in frequency space). 

"""

from functions import make_circle, plot, make_meshgrid
from propagators import scalar_prop
import matplotlib.pyplot as plt
import numpy as np

# units are microns
dx = 0.3 # Decrease to increase sampling
dy = 0.3 # Decrease to increase sampling
Lx_max = 100 # Increase to increase zero padding
Ly_max = 100 # Increase to increase zero padding
wavelength = .6
z = 200


 
aperture, extent, F_aperture, F_extent = make_circle(dx = dx,
                                                  dy = dy,
                                                  Lx_max = Lx_max,
                                                  Ly_max = Ly_max,
                                                  a = .05, # stretch in x
                                                  b=.05, # stretch in y
                                                  )

output_field = scalar_prop(aperture,dx,dy,z,wavelength)

plot(aperture, extent, figsize_x = Lx_max/10, figsize_y = Ly_max/10)
plot(F_aperture, F_extent, figsize_x = Lx_max/10, figsize_y = Ly_max/10)
plot(output_field, extent, figsize_x = Lx_max/10, figsize_y = Ly_max/10)


# Loop different sampling
dx_vec = [0.1, 1.0]
Lx_max = 100 # Increase to increase zero padding
Ly_max = 100 # Increase to increase zero padding
wavelength = .6
z = 200

plt.figure()
for dx in dx_vec:

    
    x,y,xm,ym,extent,F_extent = make_meshgrid(dx = dx,
                                              dy = dx,
                                              Lx_max = Lx_max,
                                              Ly_max = Ly_max)
    
    aperture, extent, F_aperture, F_extent = make_circle(dx = dx,
                                                         dy = dx,
                                                         Lx_max = Lx_max,
                                                         Ly_max = Ly_max,
                                                         a = .05, # stretch in x
                                                         b=.05, # stretch in y
                                                         )

    output_field = scalar_prop(aperture,dx,dx,z,wavelength)
    plt.plot(y, np.abs(output_field[len(x)//2,:]),label="%.2f" % dx)
plt.legend()

# Loop different zero padding
Lx_max_vec = [50, 100]
dx = 0.6 # Decrease to increase sampling
dy = 0.6 # Decrease to increase sampling
wavelength = .6
z = 200

plt.figure()
for Lx_max in Lx_max_vec:

    
    x,y,xm,ym,extent,F_extent = make_meshgrid(dx = dx,
                                              dy = dy,
                                              Lx_max = Lx_max,
                                              Ly_max = Lx_max)
    
    aperture, extent, F_aperture, F_extent = make_circle(dx = dx,
                                                         dy = dy,
                                                         Lx_max = Lx_max,
                                                         Ly_max = Lx_max,
                                                         a = .05, # stretch in x
                                                         b=.05, # stretch in y
                                                         )

    output_field = scalar_prop(aperture,dx,dy,z,wavelength)
    plt.plot(y, np.abs(output_field[len(x)//2,:]),label="%.0f" % Lx_max)
plt.legend()