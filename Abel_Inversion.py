# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 16:48:50 2020

@author: samia
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import abel
from scipy.optimize import curve_fit
import random


def smooth(x,window_len=11,window='hanning'): #smoothing function

    
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y
    


def parabolic(x,A,c): #infinite parabolic waveguide function to retrieve the matched spot size
    return A * x**2 + c 


"""
Abel Inversion function, parameters are:
    data: .txt file with vertical strip phase profile
    sign: predetermined sign of the phase profile
    offset: Offset from axis of symetry (predetermined using an optimising loop) 
    smoohing: size of the smoothing Hamming window
"""
def abel_inversion(data, sign, offset, smoothing): 
    
    df0 = pd.read_csv(data) #read datat from .txt file
    
    micperpix = 0.34 #camera resolution caibration
    r = 2.85e-15
    y = 1e-6
    K = 1/(r*y)
    w_tot = []
    
 
    tally = np.ndarray(0) #find axis of symmetry 
    
    if sign == -1:
        for i in range(df0.shape[0]):
            if df0.values[i,0] <=0 :
                tally = np.append(tally, [i])
    else :
        for i in range(df0.shape[0]):
            if df0.values[i,0] >= sign :
                tally = np.append(tally, [i])
             
    R = 1
    for j in range(0, R): #change range R of for loop if you want to perform a Monte Carlo analysis

        sym = np.min(tally)+((np.max(tally) - np.min(tally))/2) + offset 
        rand_sym = random.randint(0, 5)  
        #sym = sym + rand_sym   #random variations for Monte Carlo analysis
        print ("axis of symmetry", sym)
             
        right = df0.values[int(sym):]  #split data in half
        left = df0.values[:int(sym)]
       
        right1D =right[:,2]         
        left1D = left[:,2]
        
        z = np.ndarray(0)
        
        if left1D.size < right1D.size: #adjust length of right and left side
            for i in range(left1D.size, right1D.size):
                z = np.append(z, right1D[i])
            left1D = np.append(z, left1D)
                
        else:
            for i in range(right1D.size, left1D.size):
                z = np.append(z, left1D[i])
            right1D = np.append(z, right1D)
        
                
        avr0 = (right1D + left1D[::-1])/2 #average both halves and set background to zero
        avr0 = avr0 - avr0[-1]
        
        if sign == -1:
            for i in range(avr0.size):
                if avr0[i] > 0:
                     avr0[i:] = 0
        else:
            for i in range(avr0.size):
                if avr0[i] < 0:
                    avr0[i:] = 0
        
        rand_smooth = random.randint(0, 10)  
        #smoothing = smoothing + rand_smoothing 
        avr = smooth(avr0, smoothing, "hamming") #smooth pahse
                    
        inverse_avr = abel.hansenlaw.hansenlaw_transform(avr,
             dr=0.34e-6, direction = 'inverse',
             hold_order = 0) #perform abel inversion on the half phase using the Hansenlaw method
        
        bavr = ((inverse_avr)*K) #include physical constants
        
        if sign == -1: 
            for i in range(bavr.size):
                if bavr[i] > 0:
                    bavr[i:] = 0
        else:
            for i in range(bavr.size):
                if bavr[i] < 0:
                    bavr[i:] = 0
    
        bavrflip = bavr[::-1]
        totavr = np.append(bavrflip, bavr) #make full density by symmetry, can plot it if needed
        
        if sign == -1: #adjust the sign (make density positive)
            bavr = -bavr
            
        sizeavr = avr.size*micperpix
        lengthavr = np.linspace(0, sizeavr, num = avr.size) #convert pixel length to microns
    
        ma = bavr.argmax() 
        d = bavr[:ma]/1e24
        l = lengthavr[:ma]/1e6
        
        coef, cov = curve_fit(parabolic, l, d) #fit parabolic curve
        
        fit = parabolic(l, coef[0], coef[1])
        
        """ plot individual density lineout with parabolic fit if needed
        d = d* 1e18
        fit = fit*1e18
        plt.plot(l, d)
        plt.plot(l, fit)
        """
        
        w = np.sqrt(np.sqrt(1/(np.pi * coef[0] * 2.85e-15 * 1e24))) #get matched spot size
        w_tot.append(w)
        
        if j == R-1:
            w_avr = np.mean(w_tot) #get matched spot size average
            w_std = np.std(w_tot) #standard deviation
            print data
            print"matched spot size: average:", w_avr, " std: ", w_std
            print "----------------------"
        

        
    return bavr #return half density lineout


    
dens1a3 = abel_inversion("1-air-3.txt", 0, -8, 65)/1e6
dens1a9 = abel_inversion("1-air-9.txt", 5, 9, 65)/1e6
dens2a2 = abel_inversion("2-air-2.txt", -1, -8, 45)/1e6 
dens2a3 = abel_inversion("2-air-3.txt", -1, 0, 25) /1e6
dens3a8 = abel_inversion("3-air-8.txt", 2, -21, 50)/1e6
densvlow4 = abel_inversion("vlow4.txt", -1, 2, 85)/1e6


dens4h19 = abel_inversion("4-he-19.txt", 5, 12, 55)/1e6
dens3h3 = abel_inversion("3-he-3.txt", -1, 18, 110)/1e6
dens2h7 = abel_inversion("2-he-7.txt", -1, 24, 115)/1e6
dens2h10 = abel_inversion("2-he-10.txt", 5, 20, 100)/1e6
dens1h3 = abel_inversion("1-he-3.txt", -1, -14, 150)/1e6

"""
Function to get size of density lineout in microns. 
Parameters are the density lineout and the camera reolsuion in microns per pixels
"""
def size(dens, micperpix): 
    size = dens.size*micperpix
    length = np.linspace(0, size, num = dens.size)
    return length

lenght1a3 = size(dens1a3)
lenght1a9 = size(dens1a9)
lenght2a2 = size(dens2a2)
lenght2a3 = size(dens2a3)
lenght3a8 = size(dens3a8)
lenghtvlow4 = size(densvlow4)

lenght4h19 = size(dens4h19)
lenght3h3 = size(dens3h3)
lenght2h10 = size(dens2h10)
lenght1h3 = size(dens1h3)
    
    
fig, axs = plt.subplots(1, 2, figsize=(14, 7))
    
axs[0].plot(lenght3a8, dens3a8, label = 'Pressure 3', linewidth = 3)
axs[0].plot(lenght2a3, dens2a3, label = 'Pressure 2', linewidth = 3)
axs[0].plot(lenght1a3, dens1a3, label = 'Pressure 1', linewidth = 3)

axs[0].legend(prop = {'size':18})

axs[1].plot(lenght4h19, dens4h19, label = 'Pressure 4', linewidth = 3)
axs[1].plot(lenght3h3, dens3h3, label = 'Pressure 3', linewidth = 3)
axs[1].plot(lenght2h10, dens2h10, label = 'Pressure 2', linewidth = 3)
axs[1].plot(lenght1h3, dens1h3, label = 'Pressure 1', linewidth = 3)


axs[1].legend(prop = {'size':20})
  
axs[0].set_title('Air', fontsize = 32)
axs[0].tick_params(labelsize = 20)
axs[0].set_xlabel('x [microns]', fontsize = 22)
axs[0].set_ylabel('Electron Number Density [cm$^{-3}$]', fontsize = 22)

axs[1].set_title('Helium', fontsize = 32)
axs[1].tick_params(labelsize = 20)
axs[1].set_ylabel('Electron Number Density [cm$^{-3}$]', fontsize = 22)
axs[1].set_xlabel('x [microns]', fontsize = 22)


    
plt.show()



