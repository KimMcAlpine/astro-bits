#! /usr/bin/env

import math 
import numpy
import scipy.integrate

H0=70.0 # H0 in units of  km s-1 Mpc-1
c= 299792.458 # km / s
Mpctom=3.085677581e+22


def lum_dist(limz=numpy.array([4.0])): #luminosity distance in Mpc
    """calculate luminosity distance at redshift limz in Mpc"""
    if type(limz)==float:
    	lin=scipy.integrate.quad(lambda z: numpy.power(numpy.power((0.3*numpy.power(z+1,3)+0.7),0.5),-1),0,limz)
    	DL=(1+limz)*(lin)[0]*c/H0
    if type(limz)==numpy.ndarray:
	zeros=numpy.zeros([len(limz)])
	for i in range(0,len(limz)):
            lin=scipy.integrate.quad(lambda z: numpy.power(numpy.power((0.3*numpy.power(z+1,3)+0.7),0.5),-1),0,limz[i])
    	    zeros[i]=(1+limz[i])*(lin)[0]*c/H0
            DL=zeros
    return DL 


def trans_dist(limz=4.0): # transverse comoving distance in Mpc
	"""calculate radial comoving distance at redshift limz in Mpc"""
	if type(limz)==float:
    		lin=scipy.integrate.quad(lambda z: numpy.power(numpy.power((0.3*numpy.power(z+1,3)+0.7),0.5),-1),0,limz)
    		Dm=(lin)[0]*c/H0
    	if type(limz)==numpy.ndarray:
	  zeros=numpy.zeros([len(limz)])
	  for i in range(0,len(limz)):
            lin=scipy.integrate.quad(lambda z: numpy.power(numpy.power((0.3*numpy.power(z+1,3)+0.7),0.5),-1),0,limz[i])
    	    zeros[i]=(lin)[0]*c/H0
            Dm=zeros
    	return Dm 

def com_volume(limz=4.0): # comoving volume within z
        """ calculate comoving volume at redshift limz in cubic Gpc"""
	if type(limz)==float:
    		lin=scipy.integrate.quad(lambda z: numpy.power(numpy.power((0.3*numpy.power(z+1,3)+0.7),0.5),-1),0,limz)
    		Dm=(lin)[0]*c/H0
    	if type(limz)==numpy.ndarray:
	  zeros=numpy.zeros([len(limz)])
	  for i in range(0,len(limz)):
            lin=scipy.integrate.quad(lambda z: numpy.power(numpy.power((0.3*numpy.power(z+1,3)+0.7),0.5),-1),0,limz[i])
    	    zeros[i]=(lin)[0]*c/H0
            Dm=zeros
	V=4.*math.pi/3.0*(Dm**3)*1e-9 # gives output in cubic Gpc
	return V 


def radio_lum(fluxin=1.0,obsnu=1.4,emnu=1.4,alpha=-0.7,z=2):
        """ luminosity of a radio source, fluxin must be in mJy, alpha=-0.7 ie specI must include the negative sign """
	fluxobs=fluxin*(numpy.power(numpy.divide(emnu,obsnu),alpha))*(numpy.power((1+z),-(alpha+1))) #convert from obs to emitted flux at appropriate frequency
	lumobs=fluxobs*4*math.pi*(lum_dist(limz=z)*Mpctom)**2 # convert from flux to luminosity
	lumobsW=lumobs*1e-29 # convert mJy to W 
	return lumobsW


