#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 18:02:56 2021

@author: vganapati
"""

from functions import plot, F, Ft, make_meshgrid, make_circle, make_2D, rect
import numpy as np
from scipy.signal import convolve2d

dx = 0.01
dy = 0.01
Lx_max = 1
Ly_max = 1
                  
x,y,xm,ym,extent,F_extent = make_meshgrid(dx = dx,
                                          dy = dy,
                                          Lx_max = Lx_max,
                                          Ly_max = Ly_max)

dots = np.zeros_like(xm)

dots[len(x)//3,len(y)//3]=1
dots[len(x)*2//3, len(y)*2//3]=1

plot(dots, extent)

# shape_0, extent, F_circle, F_extent = make_circle(dx = dx,
#                                                   dy = dy,
#                                                   Lx_max = Lx_max,
#                                                   Ly_max = Ly_max,
#                                                   a = 10, # stretch in x
#                                                   b=10, # stretch in y
#                                                   )

shape_0, extent, F_circle, F_extent = make_2D(dx = dx,
                                              dy = dy,
                                              Lx_max = Lx_max,
                                              Ly_max = Ly_max,
                                              a = 10,
                                              b=10,
                                              func=rect # options are rect, triangle, sgn, comb, or exp_func_n
                                              )
plot(shape_0, extent)

# Convolve dots and shape_0
convolved_0 = convolve2d(dots,shape_0,'same')
plot(convolved_0, extent)

# Convolve with FFT
convolved_0_fft = Ft(F(dots)*F_circle)
plot(convolved_0_fft, extent)
