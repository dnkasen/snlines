###############################################
# Simple script to calculate a single 
# supernova line in the Sobolev approximation
# Given some 3D distribution
# of optical depth and a source function
#
# Note that everything is done in velocity 
# coordinates with units km/s
# The conversion to radius is r = v*t
# where t is the time since explosion
#
# Dan Kasen (kasen@berkeley.edu)
###########################################

import numpy as np
import matplotlib.pyplot as plt 

##########################################

##### parameters of the calculation

## wavelength grid (in angstroms)
l_start  = 4300;    # starting wavelength of spectrum
l_stop   = 5500;    # ending wavelength of spectrum
l_del    = 10;      # spectral bin size

## photospheric properties
v_phot  = 10000.0;  # velocity of the photosphere
T_ph    = 10000.0   # temperature of the photosphere

## line properties
lambda_0 = 5000.0;  # rest wavelength of the line (in angstroms)
tau_0    = 5.0      # Sobolev optical depth at photosphere

## spatial grid
n_x     = 80;       # number of zones per dimension
v_max   = 20000.0;  # velocity edge of the grid, in km/s

  
##########################################

# physical constants
c_light_kms   = 2.99e5
c_light_cms   = 2.99e10
h_planck      = 6.6260755e-27 
k_boltz       = 1.380658e-16

lambda_grid   = np.arange(l_start,l_stop,l_del)
spatial_grid  = np.arange(-1*v_max,v_max,2*v_max/(1.0*n_x))


###########################################
# CALCULATE SPECTRUM
###########################################

# direction to the observer is in the negative z direction
# At each wavelength, wwe integrate the intensity over
# the perendicular directions x and y 

# array to hold output spectrum
# zero it out to start
L_observed = np.zeros(len(lambda_grid))

# loop over all wavelengths
i = 0
for lam in lambda_grid:

	# z coordinate of the common direction plane
	# for this wavelength 
	vz = (lam-lambda_0)/lambda_0*c_light_kms;

	# set photospheric intensity to a blackbody
	lam_cm = lam*1e-8
	I_ph = 2.0*h_planck*c_light_cms**2.0/(lam_cm)**5.0
	I_ph = I_ph/(np.exp(h_planck*c_light_cms/lam_cm/k_boltz/T_ph))

	# loop over x and y dimensions to integrate
	for vx in spatial_grid:
		for vy in spatial_grid:

			# radius at this point
			vr = (vx*vx + vy*vy + vz*vz)**0.5
			# impact parameter at this point
			vp = (vx*vx + vy*vy)**0.5

			# get the line optical depth at this point
			# uses a simple power-law for now
			if (vr < v_phot): etau = 1
			else:
				tau = tau_0*(vr/v_phot)**(-8.0)
				etau = np.exp(-1.0*tau)

			# get the line source function at this point
			# uses simple pure scattering for now
			if (vr < v_phot): S = 0
			else:
				S  = I_ph*0.5*(1 - (1 - (v_phot/vr)**2.0)**0.5)	

			# calculate intensity from this ray
			this_I = 0

			# see if this position is inside the photosphere
			# If so, just we just see the photospheric intensity
			if (vr < v_phot):
				this_I = I_ph

			# see if this position is behind the photosphere
			# If so, just we just see the photospheric intensity
			elif (vp < v_phot and vz > 0): 
				this_I = I_ph

			# see if this position is in front of the photosphere
			# if so, add in attenuated I_ph plus line emission
			elif (vp < v_phot):
				this_I = I_ph*etau + S*(1-etau)
	
			# see if this region is to the side of the photosphere
			# if so, just add in line emission
			else:
				this_I += S*(1 - etau);

			# add in this intensity ray
			L_observed[i] += this_I
	
	i = i + 1

plt.plot(lambda_grid,L_observed)
plt.ion()
plt.show()
j = raw_input()
