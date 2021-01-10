#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 10:39:10 2021

@author: vganapa1
"""

from functions import F, Ft, make_meshgrid, plot
import numpy as np

def scalar_prop(u0,dx,dy,z,wavelength): # x is the row coordinate, y is the column coordinate
    Nx = u0.shape[0]
    Ny = u0.shape[1]
    Lx = Nx*dx
    Ly = Ny*dy
    
    U0 = F(u0)
    
    fx=np.linspace(-1/(2*dx),1/(2*dx)-1/Lx,Nx) #freq coords
    fy=np.linspace(-1/(2*dy),1/(2*dy)-1/Ly,Ny) #freq coords
    
    FX,FY=np.meshgrid(fx,fy, indexing = 'ij')

    H = np.exp(1j*2*np.pi*(1./wavelength)*z*np.sqrt(1-(wavelength*FX)**2-(wavelength*FY)**2))
    H[np.nonzero( np.sqrt(FX**2+FY**2) > (1./wavelength) )] = 0
    
    U1 = U0*H
    u1 = Ft(U1)
    
    return(u1)

def fresnel_prop(u0,dx,dy,z,wavelength): # x is the row coordinate, y is the column coordinate
    Nx = u0.shape[0]
    Ny = u0.shape[1]
    Lx = Nx*dx
    Ly = Ny*dy
    
    U0 = F(u0)
    
    fx=np.linspace(-1/(2*dx),1/(2*dx)-1/Lx,Nx) #freq coords
    fy=np.linspace(-1/(2*dy),1/(2*dy)-1/Ly,Ny) #freq coords
    
    FX,FY=np.meshgrid(fx,fy, indexing = 'ij')

    H_fresnel=np.exp(-1j*np.pi*wavelength*z*(FX**2+FY**2))
    H_fresnel[np.nonzero( np.sqrt(FX**2+FY**2) > (1./wavelength) )] = 0
    
    U1 = U0*H_fresnel
    u1 = Ft(U1)
    
    return(u1)

def create_plane_wave(dx, dy,
                      Lx_max, Ly_max, 
                      alpha, beta, wavelength):
    
    fx = alpha/wavelength
    fy = beta/wavelength
    
    x,y,xm,ym,extent,F_extent = make_meshgrid(dx,
                                              dy,
                                              Lx_max,
                                              Ly_max)
    
    plane_wave = np.exp(1j*2*np.pi*(fx*xm + fy*ym))
    return(plane_wave)

if __name__ == "__main__":
    from functions import make_circle
    
    # units are microns
    dx = 0.6
    dy = 0.6
    Lx_max = 1000
    Ly_max = 1000
    wavelength = .6
    alpha =  .1
    beta = .1
    z = 1000
    
    plane_wave = create_plane_wave(dx, dy,
                                   Lx_max, Ly_max, 
                                   alpha, beta, wavelength)
   
    
    mat2D, extent, F_mat2D, F_extent = make_circle(dx = dx,
                                                   dy = dy,
                                                   Lx_max = Lx_max,
                                                   Ly_max = Ly_max,
                                                   a = .005, # stretch in x
                                                   b=.005, # stretch in y
                                                   )
    plot(plane_wave, extent)
    plot(mat2D, extent)
    output_field = scalar_prop(mat2D*plane_wave,dx,dy,z,wavelength)
    
    
    plot(output_field, extent)
    