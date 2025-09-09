#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#        Physical Dimensions and Material Properties            #
#################################################################

import numpy as np

#--------------------Physical dimensions-----------------------# 
# Pipe dimensions #
D = 0.192 # diameter of the pipe, in m
L = 16.2 # length of the pipe, in m

# Particles dimensions #
# The input is an array
dp = [64e-6] # in m

#--------------------Material properties----------------------#
# Pipe #
# Pipe's material in Rautiainen experiment is Plexiglass
# The used plexiglass roughness is based on average rougness (Sa) 
# from Lepore et al paper (2008). doi:10.1155/2008/194524

rough = 0.033e-6 # roughness, in m

# Particles #
# Density of particles based on Rautiainen's paper (glass beads).
rho_p = 2450 # density of the particle, in kg/m3

# Fluid
rho_f = 1.17 # density of the fluid (air), in kg/m3
mu = 1.789e-5 # dynamic viscosity of fluid (air), in kg/m/s

#---------------------Flow parameters-----------------------#
# Fluid
# The input is an array
Vf = np.linspace(3.5,13,num=100) # fluid velocity, in m/s.
# The range of this fluid velocity is based on Rautiainen's paper, may be changed

# Solid
# The input is an array
mflux_part = np.array([5.1, 10.5, 19.6, 31.8, 49.1, 59, 80.7, 141]) # solid mass flux, in kg/m2/s

#----------------------Intial Guesses-----------------------#
inguess_frict = 0.01 # initial guess of friction factor
inguess_phi = 0.3 # initial guess of phi

#-------Switch to include particle friction at wall--------#
partFrict = 0 #1: on, 0: off, default is 0

# Compulsory additional input if partFrict = 1
partFrictCoeff = 0.01 # constant to be fitted on the experiments
mu_p = 0.5 # Coulomb dynamic coefficient of friction

#-------------Fitting coefficient, optional----------------#
bodyForceCoeff = 1 # If not used, just keep it as 1, default is 1!!
fluidFrictCoeff = 1 # If not used, just keep it as 1, default is 1!! 

#-------------------Saving file settings--------------------#
output_filename = "pneumatic_transport_results.csv"
# output_filename = "pneumatic_transport_results_withPartFrict0.01.csv"