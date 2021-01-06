#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 10:12:44 2021

@author: vganapati
"""

import numpy as np
import matplotlib.pyplot as plt

def F(mat2D): # Fourier transform assuming center pixel as center
    return np.fft.fftshift(np.fft.fft2(mat2D))

def Ft(mat2D): # Inverse Fourier transform assuming center pixel as center
    return np.fft.ifft2(np.fft.ifftshift(mat2D))

def make_meshgrid(dx = 0.01,
                  dy = 0.01,
                  Lx_max = 1,
                  Ly_max = 1):
    
    Lx_min = -Lx_max
    Ly_min = -Ly_max
    x = np.arange(Lx_min,Lx_max,dx)
    y = np.arange(Ly_min,Ly_max,dy)
    xm, ym = np.meshgrid(x,y, indexing='ij')
    extent = [Lx_min, Lx_max-dx, Ly_min, Ly_max-dy]
    
    # Fourier transform dimensions
    F_extent = [len(x)/(4*Lx_min), len(x)/(4*Lx_max) - 1/(Lx_max*2), 
                len(y)/(4*Ly_min), len(y)/(4*Ly_max) - 1/(Ly_max*2)]
    
    return(x,y,xm,ym,extent,F_extent)

def make_circle(dx = 0.01,
                dy = 0.01,
                Lx_max = 1,
                Ly_max = 1,
                a = 1, # stretch in x
                b=1, # stretch in y
                ):
    
    x,y,xm,ym,extent,F_extent = make_meshgrid(dx,
                                              dy,
                                              Lx_max,
                                              Ly_max)

    circle = np.zeros_like(xm)
    circle[((a*xm)**2+(b*ym)**2) < 1] = 1
    circle[((a*xm)**2+(b*ym)**2) == 1] = 0.5

    F_circle = F(circle)
    
    return(circle, extent, F_circle, F_extent)

def plot(img, extent):

    plt.figure()
    plt.title('Magnitude')
    plt.imshow(np.abs(img), interpolation='none', extent = extent, cmap='gray')
    plt.colorbar()

    plt.figure()
    plt.title('Phase')
    plt.imshow(np.angle(img), interpolation='none', extent = extent, cmap='hsv',
               vmin = -np.pi, vmax = np.pi)
    plt.colorbar()
        
def rect(mesh):
    output = np.zeros_like(mesh)
    output[np.abs(mesh)<0.5] = 1
    output[np.abs(mesh)==0.5] = 0.5
    return(output)
    
def triangle(mesh):
    output = np.zeros_like(mesh)
    output[np.abs(mesh)<=1] = 1 - np.abs(mesh[np.abs(mesh)<=1])
    return(output)

def sgn(mesh):
    output = np.zeros_like(mesh)
    output[mesh>0] = 1
    output[mesh<0] = -1
    return(output)

def comb(mesh):
    output = np.zeros_like(mesh)
    output[np.abs(mesh % 1) <= 1e-10] = 1
    return(output)

exp_func_0 = lambda mesh: np.exp(-np.pi*mesh**2)
exp_func_1 = lambda mesh: np.exp(1j*np.pi*mesh)
exp_func_2 = lambda mesh: np.exp(1j*np.pi*mesh**2)
exp_func_3 = lambda mesh: np.exp(-(np.abs(mesh)))

def make_2D(dx = 0.01,
            dy = 0.01,
            Lx_max = 1,
            Ly_max = 1,
            a = 1,
            b=1,
            func=rect # options are rect, triangle, sgn, comb, or exp_func_n
            ):
    
    x,y,xm,ym,extent,F_extent = make_meshgrid(dx,
                                              dy,
                                              Lx_max,
                                              Ly_max)

    mat2D = func(a*xm)*func(b*ym)
        
    return(mat2D, extent, F(mat2D), F_extent)


def make_2D_delta(dx = 0.01,
                  dy = 0.01,
                  Lx_max = 1,
                  Ly_max = 1,
                  a = 1,
                  b=1,
                  ):
    
    x,y,xm,ym,extent,F_extent = make_meshgrid(dx,
                                              dy,
                                              Lx_max,
                                              Ly_max)
    delta_x = np.zeros_like(xm)
    delta_y = np.zeros_like(xm)
    delta_x[np.abs(xm)<1e-10]=1/np.abs(a)
    delta_y[np.abs(ym)<1e-10]=1/np.abs(b)
    
    mat2D = delta_x*delta_y
        
    return(mat2D, extent, F(mat2D), F_extent)

if __name__ == "__main__":
    func = exp_func_0
    mat2D, extent, F_mat2D, F_extent = make_2D(func=func)
    plot(mat2D, extent)

    mat2D, extent, F_mat2D, F_extent = make_2D_delta()
    plot(mat2D, extent)