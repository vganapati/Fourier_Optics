#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 15:26:09 2021

@author: vganapa1
"""

from functions import make_meshgrid, plot, Ft
import numpy as np

def make_annulus(dfx = 0.01,
                 dfy = 0.01,
                 fx_max = 1,
                 fy_max = 1,
                 inner_radius = 0.5, 
                 d_radius = .04, # stretch in y
                 ):
    
    fx,fy,fxm,fym,F_extent,_ = make_meshgrid(dfx,
                                           dfy,
                                           fx_max,
                                           fy_max)

    annulus = np.zeros_like(fxm)
    annulus[(fxm**2+fym**2) <= (inner_radius+d_radius)] = 1
    annulus[(fxm**2+fym**2) < inner_radius] = 0


    return(annulus, F_extent)



def scalar_prop_freq(U0,
                     dfx,
                     dfy,
                     fx_max,
                     fy_max,
                     z,wavelength): # x is the row coordinate, y is the column coordinate

    fx,fy,fxm,fym,F_extent,_ = make_meshgrid(dfx,
                                             dfy,
                                             fx_max,
                                             fy_max)
    
    FX,FY=np.meshgrid(fx,fy, indexing = 'ij')

    H = np.exp(1j*2*np.pi*(1./wavelength)*z*np.sqrt(1-(wavelength*FX)**2-(wavelength*FY)**2))
    H[np.nonzero( np.sqrt(FX**2+FY**2) > (1./wavelength) )] = 0
    
    U1 = U0*H
    u1 = Ft(U1)

    
    return(u1)



if __name__ == "__main__":
    
    
    dfx = 0.0005
    dfy = 0.0005
    fx_max = 1
    fy_max = 1
    wavelength = 0.5
    z = 0
              
    F_annulus, F_extent = make_annulus(dfx = dfx,
                                       dfy = dfy,
                                       fx_max = fx_max,
                                       fy_max = fy_max,
                                       d_radius = dfx*3)

    plot(F_annulus, F_extent)
    
    u1 = scalar_prop_freq(F_annulus,
                          dfx,
                          dfy,
                          fx_max,
                          fy_max,
                          z,wavelength)
    plot(u1, F_extent, log=True)

    z = 500
    u1 = scalar_prop_freq(F_annulus,
                         dfx,
                         dfy,
                         fx_max,
                         fy_max,
                         z,wavelength)
    plot(u1, F_extent, log=True)