# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 22:57:27 2015

@author: manu
"""
from __future__ import division
from pylab import *
from numpy import *
from pylab import rcParams

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})

f=load("Exstirrer1source.npz")["f"]
Nx1l=load("Exlinear.npz")["Nx"]
Nmx1l=load("Exlinear.npz")["Nmodesx"]
Nx1b=load("Exbrownian.npz")["Nx"]

figure(1,figsize=(8,4))
plot(f[0:960:5]/1e6,Nx1l[0:960:5],label="$N_{x_{eff}}$")
plot(f[0:960:5]/1e6,Nmx1l[0:960:5],label="$N_{x_{modes}}$")
grid()
xlabel("$f$/MHz")
legend(loc=2)
savefig("fig7.pdf",bbox_inches="tight")




Nx1=load("Exstirrer1source.npz")["Nx"]
Nx100=load("Exstirrer100sources.npz")["Nx"]
Nx6=load("Exstirrer.npz")["Nx"]


def runningMeanFast(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]

def runningmean(x,N):
    a=zeros(len(x)-N)
    for i in range(len(x)-N):
        a[i]=sum(x[i:i+N])/N
    return a

def runningmin(x,N):
    a=zeros(len(x)-N)
    for i in range(len(x)-N):
        a[i]=min(x[i:i+N])
    return a

def runningmax(x,N):
    a=zeros(len(x)-N)
    for i in range(len(x)-N):
        a[i]=max(x[i:i+N])
    return a    
    
figure(2)   
plot(f/1e6,Nx1l,"gx",ms=5,alpha=0.3,label="$N_{x_{eff}}$, 1 source linear path")    
plot(f[:-50]/1e6,runningmean(Nx1l, 50),"g",lw=2,label="Moving av. 1 source linear path")

plot(f/1e6,Nx1b,"b+",ms=5,alpha=0.3,label="$N_{x_{eff}}$, 1 source random path")    
plot(f[:-50]/1e6,runningmean(Nx1b, 50),"b",lw=2,label="Moving av. 1 source random path")

legend(loc=2)
xlim(40,950)
ylim(1,15)
grid()
xlabel("$f$/MHz")
savefig("fig_N1s.pdf",bbox_inches="tight")




#figure(2,figsize=(14,6))
##plot(f/1e6,Nx1l,"k*")
#plot(f/1e6,runningMeanFast(Nx1l, 50),"g",label="1 source linear path")
#
#plot(f/1e6,runningMeanFast(Nx1b, 50),"b",label="1 source random path")
#
##plot(f/1e6,Nx1,"b+")
#plot(f/1e6,runningMeanFast(Nx1, 50),color='y',label="$N_s=1$ source rotating")
##plot(f/1e6,Nx6,"g.")
#plot(f/1e6,runningMeanFast(Nx6, 50),'#EF5915',label="$N_s=6$ sources rotating")
##plot(f/1e6,Nx100,"rs")
#plot(f/1e6,runningMeanFast(Nx100, 50),'k',label="$N_s=100$ sources rotating")
#legend(loc=2)
#xlim(40,950)
#ylim(1,7)
#grid()
#xlabel("$f$/MHz")
#savefig("fig_N.pdf",bbox_inches="tight")


figure(3)
plot(f/1e6,runningMeanFast(Nx1, 50),color='y',label="1 source rotating")
plot(f/1e6,runningMeanFast(Nx6, 50),'#EF5915',label="6 sources rotating")
plot(f/1e6,runningMeanFast(Nx100, 50),'k',label="100 sources rotating")
plot(f/1e6,runningMeanFast(Nx1l, 50),"g",label="1 source linear path")
plot(f/1e6,runningMeanFast(Nx1b, 50),"b",label="1 source random path")
legend(loc=2)
xlim(40,950)
ylim(1,7)
grid()
xlabel("$f$/MHz")
savefig("fig_Nm.pdf",bbox_inches="tight")


