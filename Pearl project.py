# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 15:02:24 2016

@author: Pearl
"""
#Wastewater treatment plant design

import scipy as sp
import numpy as np
from scipy import integrate
import pylab as pl

print 'Press the "Return" button to choose the Default values. \n '

mum = raw_input ('mum: ' '\n')
if mum == '': mum = 0.35
else: mum = float(mum)

Ks = raw_input('Ks:    ' '\n')
if Ks == '': Ks = 7.1
else: Ks = float(Ks)

Kd = raw_input('Kd:    ' '\n')
if Kd == '': Kd = 0.02
else: Kd = float(Kd)

Ki = raw_input('Ki:    ' '\n')
if Ki == '': Ki = 0.5
else: Ki = float(Ki)

Yx = raw_input('Yx:    ' '\n')
if Yx == '': Yx = 0.3
else: Yx = float(Yx)

Ys = raw_input('Ys:    ' '\n')
if Ys == '': Ys = 0.04
else: Ys = float(Ys)

Yp = raw_input('Yp:    ' '\n')
if Yp == '': Yp = 4.35
else: Yp = float(Yp)

Kmx = raw_input('Kmx:  ' '\n')
if Kmx == '': Kmx = 0.4
else: Kmx = float(Kmx)

Ksx = raw_input('Ksx:  ' '\n')
if Ksx == '': Ksx = 0.983
else: Ksx = float(Ksx)

D = raw_input('Dilution Rate:  ' '\n')
if D == '': D = 0.025
else: D = float(D)

Xo = raw_input('Initial MO conc:  ' '\n')
if Xo == '': Xo = 3
else: Xo = float(Xo)

So = raw_input('Initial substarte conc:   ' '\n')
if So == '': So = 6
else: So = float(So)


'''to design the anaerobic deigestor volume'''


Zo=0
mu=mum
uo = np.array([Xo,So,Zo])

Si=So 
Xi=Xo
Xe=2 #g/L
Se=1.5 #g/L
iters = np.arange(0.0,10.0,0.001)
'''
eta=(So-Se)/So
Q=300 
thetac=((Ks+Se)/(mu*Se))
theta=((thetac*Yx(So-Se))/(Xo(1+thetac*Kd)))
Vol=Q*theta

print(Vol)
'''
def sludge(u,t):
    x = (D*(Xi-Xe)+mu*u[0]-Kd*u[0])
    z = (Yp*mu*u[0])
    y = (D*(Si-Se)-((mu*u[0])/Yx)-Ksx*u[0]*mu-Kmx*u[0]*(u[1]/(Ks+u[1])))
	
	
   
    return np.array([x,y,z], dtype = float)


# The equations are integrated in vector form using odeint

results = integrate.odeint(sludge,uo,iters)


pl.xlabel('Time (days)') # set x-axis label
pl.ylabel('Conc (g/L)') # set y-axis label
#pl.title('')# set plot title.plot(iters,np.log10(results)) 
pl.title('\n' + '\n' + 'Anaerobic Treatment Of Sludge Conc Profiles \n'+ '\n'+'Blue = MO, Red = Substrate, Green = Methane')
pl.plot(iters,results[:,0],iters,results[:,1],iters,results[:,2])

pl.show()

