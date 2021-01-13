#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 19:04:03 2021

@author: vganapa1
"""
import numpy as np
from functions import make_meshgrid, plot, make_circle, F
from propagators import scalar_prop, fresnel_prop
import matplotlib.pyplot as plt

# units are microns
wavelength = .6
focal_length = 10000
dx = 1.0
dy = 1.0
Lx_max = 1000
Ly_max = 1000
                  
x,y,xm,ym,extent,F_extent = make_meshgrid(dx = dx,
                                          dy = dy,
                                          Lx_max = Lx_max,
                                          Ly_max = Ly_max)

# Create the transmission function of a lens
lens_phase = -np.pi/(wavelength*focal_length)*(xm**2+ym**2) 
lens_transmission = np.exp(1j*lens_phase) # Equation 6-10

# Simulate a plane wave passing through a lens and propagating the focal length

output_field = scalar_prop(lens_transmission,dx,dy,
                            focal_length, # z
                            wavelength)

plt.imshow(lens_phase)

plot(lens_transmission, extent)
plot(output_field, extent)



# Simulate an input field at the front focal plane propagating to the lens,
# passing the lens, and propagating the focal length propagating with scalar transform

input_field, extent, F_mat2D, F_extent = make_circle(dx = dx,
                                                     dy = dy,
                                                     Lx_max = Lx_max,
                                                     Ly_max = Ly_max,
                                                     a = .05, # stretch in x
                                                     b=.05, # stretch in y
                                                     )
output_field = scalar_prop(input_field,dx,dy,
                           focal_length, # z
                           wavelength)

output_field = output_field*lens_transmission

output_field = scalar_prop(output_field,dx,dy,
                           focal_length, # z
                           wavelength)

plot(output_field, extent)

# Lens performs a Fourier transform

F_input_field = F(input_field)
extent2 = np.array(F_extent)*wavelength*focal_length
plot(F_input_field, extent2)

# Lens law
# 1/z_1 + 1/z_2 = 1/f # 6-34
# Magnification M = -z_2/z_1

z_1 = focal_length*2
z_2 = (1/focal_length - 1/z_1)**(-1)
print('Magnification is: ' + str(-z_2/z_1))
# Propagate from z_1 to focal plane 1

input_field, extent, F_mat2D, F_extent = make_circle(dx = dx,
                                                     dy = dy,
                                                     Lx_max = Lx_max,
                                                     Ly_max = Ly_max,
                                                     a = .005, # stretch in x
                                                     b=.005, # stretch in y
                                                     )

plot(input_field, extent)
output_field = scalar_prop(input_field,dx,dy,
                           z_1 - focal_length, # z
                           wavelength)

# Lens performs a Fourier transform
output_field = F(output_field)
extent2 = np.array(F_extent)*wavelength*focal_length

# Propagate from focal plane 2 to z_2
dx = (extent2[1] - extent2[0])/output_field.shape[0]
dy = (extent2[3] - extent2[2])/output_field.shape[1]
output_field = scalar_prop(output_field,dx,dy,
                           z_2 - focal_length, # z
                           wavelength)
plot(output_field, extent2)

# Exercise: Simulate an optical system that creates an image of a pyramidal 
# aperture that is double in size