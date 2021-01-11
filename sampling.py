#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 15:27:46 2021

@author: vganapati
"""

from propagators import create_plane_wave
from functions import plot

# units are microns
dx = 0.3
dy = 0.3
Lx_max = 100
Ly_max = 100
wavelength = 0.6


extent = [-Lx_max, Lx_max-dx, -Ly_max, Ly_max-dy]
alpha_max = wavelength/(2*dx)
beta_max = wavelength/(2*dy)


# No aliasing
alpha =  alpha_max*1
beta = 0
plane_wave = create_plane_wave(dx, dy,
                               Lx_max, Ly_max, 
                               alpha, beta, wavelength)

    
# Aliasing example 1
plot(plane_wave, extent)
alpha =  alpha_max*2
beta = 0
plane_wave = create_plane_wave(dx, dy,
                               Lx_max, Ly_max, 
                               alpha, beta, wavelength)
plot(plane_wave, extent)


# Aliasing example 2
alpha =  alpha_max*7
beta = 0
plane_wave = create_plane_wave(dx, dy,
                               Lx_max, Ly_max, 
                               alpha, beta, wavelength)
plot(plane_wave, extent)