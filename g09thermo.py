#!/usr/bin/env python3
import numpy as np

# CONSTANTS
from scipy.constants import N_A as CONST_NA # 1/mol
from scipy.constants import k as CONST_kB  # J/K
from scipy.constants import Planck as CONST_h # Js
from scipy.constants import speed_of_light as CONST_c # m/s
CONST_R = CONST_NA*CONST_kB
kcalmol2cm= 349.75 # kcal/mol to cm^-1
cm2kcalmol= 1/349.75

class Thermo(object):
	def __init__(self,freqs,T):
		self.T     = T
		self.freqs = np.asarray(freqs)
		self.Nfreq = self.freqs.shape # numbers of frequencies
		self.Ei     = 10**2*CONST_c*CONST_h*self.freqs # J
		self.E      = np.sum(self.Ei) # J
		self.Θi_vib = self.Ei/CONST_kB # K^-1
		# self.Ei_zpe = self.Θi_vib/2
		self.Ei_zpe = self.Ei/2
		self.E_zpe  = np.sum(self.Ei_zpe)
		# self.Ui_kcalmol[:]  = self.Ei_kcalmol[:]+self.Ei_zpe[:]*
		self.partition()
		self.pratition_molar()
		# self.thermodynamics()
		self.atomicunits()
		return None

	def atomicunits(self):
		kj2kcalmol = 0.23901*CONST_NA
		self.Ei_kcalmol     = self.Ei*10**-3*kj2kcalmol # kcal/mol
		self.E_kcalmol      = np.sum(self.Ei_kcalmol)
		self.Ei_zpe_kcalmol = kj2kcalmol*10**-3*self.Ei_zpe
		self.E_zpe_kcalmol  = np.sum(self.Ei_zpe_kcalmol)
		self.Ui_kcalmol     = self.Ei_zpe_kcalmol+self.Ei_kcalmol

	def partition(self):
		self.Qi    = np.zeros((self.Nfreq))
		self.Qi[:] = (1/(1-np.exp(-self.Θi_vib[:]/self.T)))
		self.Q     = np.prod(self.Qi)
		self.lnQ   = np.log(self.Q)
		self.lnQi  = np.log(self.Qi)
		return self.Qi,self.Q

	def pratition_molar(self):
		self.Qi_m  = CONST_NA * self.Qi
		self.Q_m   = CONST_NA * self.Q
		return self.Qi_m,self.Q_m

	def thermodynamics(self):
		self.Ui    = self.Ei_zpe + self.Ei # J
		self.U     = np.sum(self.Ui) # J
		# self.S     = CONST_kB + CONST_kB*
		self.H     = None
		self.Cv    = None
		self.A     = None
		self.G     = None
	def thermodynamics_molar(self):
		self.S_m   = None
		self.H_m   = None
		self.Cv_m  = None
		self.A_m   = None
		self.G_m   = None

ipfb2dabco = Thermo([ 5.2901, 7.0478, 7.2467, 14.8743, 16.5557],298.15)
print("Vibrational frequencies: ",ipfb2dabco.freqs)
print("Vibrational Temperatures: ",ipfb2dabco.Θi_vib)
# ipfb2dabco.partition()
print("Partial Partition Functions: ",ipfb2dabco.Qi)
print("Partial Molar Partition Functions: ",ipfb2dabco.Qi_m)
print("Molar Internal Thermal Energy: ",ipfb2dabco.Ei_kcalmol)
print("E(ZPE): ",ipfb2dabco.E_zpe)
print("Internal thermal energy Ezpe+E(T): ",ipfb2dabco.Ui_kcalmol)