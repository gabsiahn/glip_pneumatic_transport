#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#        Physical Dimensions and Material Properties            #
#################################################################

g = 9.81 # graitational acceleration, in m/s2 

#--------------------Physical dimensions-----------------------#

# Pipe dimensions #
D = 0.192 # diameter of the pipe, in m
L = 16.2 # length of the pipe, in m

# Particles dimensions #
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