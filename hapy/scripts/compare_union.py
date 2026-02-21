#!/usr/bin/env python
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table
# read in nsa
nsa = fits.getdata('NRGs27_nsa.fits')
hdi_pixelscale = 0.43
# read in my data
sdat =  fits.getdata('halpha-data-rfinn-2019-Sep-17.fits')
nsadict = dict((a,b) for a,b in zip(sdat.NSAID,np.arange(len(sdat.NSAID))))
fitflag = (sdat.GAL_XC != 0) 
# read in matched union catalog
udat = fits.getdata('union_measurements/N027_Union_nonan.fits',1)
# matched array
mdat=np.zeros(len(sdat),dtype=udat.dtype)
#mdat = np.zeros((sdat.shape[0],len(udat[0])),dtype=udat.dtype)
                
# match union to siena catalog based on nsaid
for i,n in enumerate(udat.NSAID):
    sindex = nsadict[n]
    #print(i,n,sindex)
    #print(mdat[sindex])
    #print(udat[i])
    for j in range(len(udat[i])):
        mdat[nsadict[n]][j] = udat[i][j]

mdat = Table(mdat)
temp = mdat.colnames
for i,c in enumerate(temp):
    temp[i] = c.strip(',')
mdat.names = temp
# cull tables
nsa = nsa[fitflag]
sdat = sdat[fitflag]
mdat = mdat[fitflag]
def plotxy(x,y):
    plt.plot(x,y,'bo')
    xl = np.linspace(min(x),max(x),20)
    plt.plot(xl,xl,'k--')

####  Compare NSA and galfit Sersic fits
# n, Re, PA, BA, mag

def compare_sersic():
    xvars = ['SERSICFLUX','SERSIC_N','SERSIC_BA','SERSIC_PHI','SERSIC_TH50']
    yvars = ['GAL_MAG','GAL_N','GAL_BA','GAL_PA','GAL_RE']
    x = nsa
    y = sdat
    plt.figure(figsize=(10,5))
    plt.subplots_adjust(hspace=.5,wspace=.35)
    for i in range(len(xvars)):
        plt.subplot(2,3,i+1)
        if i == 0:
            plt.plot(22.5-2.5*np.log10(x[xvars[i]][:,4]),y[yvars[i]],'bo')
        elif i == 4:
            plt.plot(x[xvars[i]],y[yvars[i]]*hdi_pixelscale,'bo')
        else:
            plt.plot(x[xvars[i]],y[yvars[i]],'bo')
        plt.xlabel(xvars[i])
        plt.ylabel(yvars[i])
        x1,x2 = plt.xlim()
        xl = np.linspace(x1,x2,100)
        plt.plot(xl,xl,'k--')

def nsasfr():
    xvars = ['HAFLUX','D4000','ABSMAG','B300']
    yvars = ['LOG_SFR_HA','LOG_SFR_HA','LOG_SFR_HA','LOG_SFR_HA']
    yvars2 = ['GAL_LOG_SFR_HA','GAL_LOG_SFR_HA','GAL_LOG_SFR_HA','GAL_LOG_SFR_HA']
    x = nsa
    y = sdat
    plt.figure(figsize=(8,6))
    plt.subplots_adjust(hspace=.5,wspace=.35)
    for i in range(len(xvars)):
        plt.subplot(2,2,i+1)
        
        if (i == 2): 
            plt.plot(x[xvars[i]][:,1],y[yvars[i]],'bo',label='PHOTUTIL')
            plt.plot(x[xvars[i]][:,1],y[yvars2[i]],'cs', label='GALFIT')
        else:
            plt.plot(x[xvars[i]],y[yvars[i]],'bo', label='PHOTUTIL')
            plt.plot(x[xvars[i]],y[yvars2[i]],'cs', label='GALFIT')
        if (i == 0) | (i == 3):
            #pass
            plt.gca().set_xscale('log')
        plt.xlabel(xvars[i])
        plt.ylabel(yvars[i])
        x1,x2 = plt.xlim()
        xl = np.linspace(x1,x2,100)
        plt.legend()
        #plt.plot(xl,xl,'k--')
def sfr():
    xvars = ['gnhaf','inhaf','onhaf']
    yvars = ['LOG_SFR_HA','SSFR_IN','SSFR_OUT',]
    yvars2 = ['GAL_LOG_SFR_HA','GAL_SSFR_IN','GAL_SSFR_OUT']
    x = mdat
    y = sdat
    plt.figure(figsize=(10,3.5))
    plt.subplots_adjust(hspace=.5,wspace=.35,bottom=.15)
    for i in range(len(xvars)):
        plt.subplot(1,3,i+1)
        if i == 0:
            yn = np.log10(sdat['HF_R24']/sdat['F_R24'])
            plt.plot(x[xvars[i]],yn,'bo', label='PHOTUTIL')
            yn = np.log10(sdat['GAL_HF_R24']/sdat['GAL_F_R24'])
            plt.plot(x[xvars[i]],yn,'cs', label='GALFIT')
        else:
            plt.plot(x[xvars[i]],np.log10(y[yvars[i]]),'bo', label='PHOTUTIL')
            plt.plot(x[xvars[i]],np.log10(y[yvars2[i]]),'cs', label='GALFIT')
        plt.xlabel(xvars[i])
        plt.ylabel(yvars[i])
        plt.legend()
        x1,x2 = plt.xlim()
        xl = np.linspace(x1,x2,100)
        plt.legend()
        plt.plot(xl,xl,'k--')
    
def morphology():
    xvars = ['CLUMPY','ASYMMETRY','SERSIC_N','SERSIC_N','CLUMPY','ASYMMETRY','SERSIC_N','SERSIC_N',]
    yvars = ['ELLIP_GINI','ELLIP_ASYM','GAL_C30','C30','ELLIP_GINI2','ELLIP_ASYM2','GAL_HC30','HC30',]
    x = nsa
    y = sdat

    plt.figure(figsize=(10,5))
    plt.subplots_adjust(hspace=.5,wspace=.35)
    for i in range(len(xvars)):
        plt.subplot(2,4,i+1)
        if i < 2:
            plt.plot(x[xvars[i]][:,4],y[yvars[i]],'bo')
        elif (i == 4) | (i == 5):
            plt.plot(x[xvars[i]][:,1],y[yvars[i]],'bo')
        else:
            plt.plot(x[xvars[i]],y[yvars[i]],'bo')
        plt.xlabel(xvars[i])
        plt.ylabel(yvars[i])
        #plt.legend()
        x1,x2 = plt.xlim()
        xl = np.linspace(x1,x2,100)
        #plt.plot(xl,xl,'k--')

def size():
    ####  SIZE
    yvars = ['R24','R25','R26','R_F50','HR17','HR_F50']
    yvars2 = ['GAL_R24','GAL_R25','GAL_R26','GAL_R_F50','GAL_HR17','GAL_HR_F50']
    xvars = ['r24','r25','r26','halfr_1','rad17','halfr_2']
    x = mdat
    y = sdat
    plt.figure(figsize=(10,8))
    plt.subplots_adjust(hspace=.5,wspace=.35)
    for i in range(len(xvars)):
        plt.subplot(2,3,i+1)
        if (i == 3) or (i == 5):
            plt.plot(x[xvars[i]],y[yvars[i]]/y[yvars[0]],'bo',label='PHOTUTILS')
            plt.plot(x[xvars[i]],y[yvars2[i]]/y[yvars2[0]],'cs',label='GALFIT')
        else:
            plt.plot(x[xvars[i]],y[yvars[i]],'bo',label='PHOTUTILS')
            plt.plot(x[xvars[i]],y[yvars2[i]],'cs',label='GALFIT')
            
        plt.xlabel(xvars[i])
        plt.ylabel(yvars[i])
        x1,x2 = plt.xlim()
        plt.legend()
        xl = np.linspace(x1,x2,100)

        plt.plot(xl,xl,'k--')
        if i < 3:        
            plt.ylim(x1,x2)
        else:
            plt.ylim(x1,1.3*x2)
def oldsize():
    # compare radii
    # nsa Re on X, galfit Re, ellip R50
    x = nsa.PETROTH50
    y = sdat.GAL_RE*hdi_pixelscale
    plt.figure()
    plotxy(x,y)
    plt.xlabel('NSA PETRO THETA 50')
    plt.ylabel('GALFIT GALFIT Re')
def mag():
    ####  SIZE
    yvars = ['M24','M25','HM17','HF_30R24']
    yvars2 = ['GAL_M24','GAL_M25','GAL_HM17','GAL_H
        F_30R24']    

    xvars = ['mag24','mag25','tot17flux','tot30r24flux,']
    x = mdat
    y = sdat
    plt.figure(figsize=(8,6))
    plt.subplots_adjust(hspace=.5,wspace=.35)
    for i in range(len(xvars)):
        plt.subplot(2,2,i+1)
        if (i == 5):
            plt.plot(x[xvars[i]],y[yvars[i]]/y[yvars[0]],'bo',label='PHOTUTILS')
            plt.plot(x[xvars[i]],y[yvars2[i]]/y[yvars2[0]],'cs',label='GALFIT')
        elif i == 3:
            plt.plot(x[xvars[i]]*1.e-18,y[yvars[i]],'bo',label='PHOTUTILS')
            plt.plot(x[xvars[i]]*1.e-18,y[yvars2[i]],'cs',label='GALFIT')
            plt.gca().set_xscale('log')
            plt.gca().set_yscale('log')

        else:
            plt.plot(x[xvars[i]],y[yvars[i]],'bo',label='PHOTUTILS')
            plt.plot(x[xvars[i]],y[yvars2[i]],'cs',label='GALFIT')
        plt.xlabel(xvars[i])
        plt.ylabel(yvars[i])
        x1,x2 = plt.xlim()
        plt.legend()
        xl = np.linspace(x1,x2,100)
        if (i < 2) | (i==3):
            plt.plot(xl,xl,'k--')
        
            plt.ylim(x1,x2)
        elif i == 2:
            plt.gca().set_xscale('log')
        #else:
        #    plt.ylim(x1,1.3*x2)
def concentration():
    # this is a mess
    yvars = ['C30','HC30']
    yvars2 = ['GAL_C30','GAL_HC30']

    xvars = ['c30','conc30,']
    x = mdat
    y = sdat
    plt.figure(figsize=(8,6))
    plt.subplots_adjust(hspace=.5,wspace=.35)
    for i in range(len(xvars)):
        plt.subplot(2,2,i+1)
        if (i == 3) or (i == 5):
            plt.plot(x[xvars[i]],y[yvars[i]]/y[yvars[0]],'bo',label='PHOTUTILS')
            plt.plot(x[xvars[i]],y[yvars2[i]]/y[yvars2[0]],'cs',label='GALFIT')
        else:
            plt.plot(x[xvars[i]],y[yvars[i]],'bo',label='PHOTUTILS')
            plt.plot(x[xvars[i]],y[yvars2[i]],'cs',label='GALFIT')
            
        plt.xlabel(xvars[i])
        plt.ylabel(yvars[i])
        x1,x2 = plt.xlim()
        plt.legend()
        xl = np.linspace(x1,x2,100)
        if i < 2:
            plt.plot(xl,xl,'k--')
        
            #plt.ylim(x1,x2)
        elif i == 2:
            plt.gca().set_xscale('log')
        #else:
    
def oldmagnitude():
    ###  MAGNITUDE
    # nsa
    x = sdat.NSA_RMAG
    y = sdat.GAL_MAG
    plt.figure()
    plotxy(x,y)
    plt.show()
    plt.xlabel('NSA r mag')
    plt.ylabel('GALFIT r mag')
    
    # isophotal magnitudes
    x = mdat['mag24']
    y = sdat.M24
    plt.figure()
    
    
    x = mdat['mag24']
    y = sdat.M24
    plt.figure()
    plt.plot(x,y,'bo',label='GALFIT')

    y = sdat.GAL_M24
    plt.plot(x,y,'cs',label='PHOTUTILS')
    plt.legend()
    plt.xlabel('Union m24')
    plt.ylabel('Siena m24')
    xl = np.linspace(min(x),max(x),20)
    plt.plot(xl,xl,'k--')
    plt.show()
###  SHAPE
def BA():
    x = nsa.SERSIC_BA
    y = sdat.GAL_BA
    plt.figure()
    plt.plot(x,y,'bo',label='GALFIT')
    x = nsa.SERSIC_BA
    y = 1-sdat.ELLIP_EPS
    plt.plot(x,y,'cs',label='PHOTUTILS')
    plt.legend()
    plt.xlabel('NSA B/A')
    plt.ylabel('Siena B/A')
    xl = np.linspace(min(x),max(x),20)
    plt.plot(xl,xl,'k--')
    plt.show()
