#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#                  Written for GLiP project                     #
#                           by                                  # 
#  Gabriele SIAHAAN, Maral EBADI, Ehab bin Nauman & Wajih UDDIN #
#################################################################

#-----------------------Initialization-----------------------#
# Python packages #
import numpy as np
import matplotlib.pyplot as plot
import scipy.optimize as sci

#-----------------------Input parameters----------------------#
import physical_props as phys

# Initial guesses #
inguess_frict = 0.001 # initial guess of friction factor based on Colebrook formula
inguess_phi = 0.3 # initial guess of phi

# Flow parameters #
# Pipe
A = np.pi/4*phys.D**2 # cross-section area of the pipe, in m2

# Fluid
Vf = np.linspace(3.5,13,num=100) # fluid velocity, in m/s
# The range of this fluid velocity is based on Rautiainen's paper

# Solid
mf_part = [5.1, 10.5, 19.6, 49.1, 80.7] # solid mass flux, in kg/m2/s
m_part = mf_part*A # solid mass flowrate, in kg/s, vector
Vp = m_part/phys.rho_p/A # solid flow rate, in m/s, vector

vt = np.zeros(len(phys.dp)) # terminal velocity by Stokes' drag, in m/s, vector
for i in range(0, len(phys.dp)):
    vt[i] = 2/9*(phys.rho_p-phys.rho_f)*phys.g*phys.dp**2/4/phys.mu

#---------------------Dimensionless number------------------#
Re_pipe = np.zeros(len(Vf)) # Reynolds number of the fluid in the pipe, vector
Re_termvel = np.zeros(len(phys.dp)) # Reynolds number of the fluid to the particle at its terminal vel

for i in range(0, len(Vf)):
    Re_pipe[i] = phys.rho_f*Vf[i]*phys.D/phys.mu 

for i in range(0, len(phys.dp)):
    Re_termvel[i] = phys.rho_f*vt[i]*phys.dp/phys.mu 

#----------Calculation of 'n' for Buoyancy Calculation------#
n = np.zeros(len(phys.dp))

for i in range (0, len(phys.dp)):
    if Re_termvel[i] < 0.2:
        n[i] = 4.65
    elif 0.2 < Re_termvel[i] < 1:
        n[i] = 4.35*Re_termvel[i]**(-0.03)
    elif 1 < Re_termvel[i] < 500:
        n[i] = 4.45*Re_termvel[i]**(-0.1)
    else:
        n[i] = 2.39

#---------------Find friction factor------------------------#
# Calculation using Colebrook-White formula

def Colebrook(f):
    return 1/np.sqrt(f) + 2*np.log10(phys.rough/3.7/phys.D + 2.51/Re_pipe/np.sqrt(f))


frict = sci.fsolve(Colebrook, inguess_frict)
print("Pipe friction factor:", frict[0])

#--------------------Finding phi---------------------------#
wall_frict = -frict[0]*phys.rho_f*Vf**2

def eps(phi_dummy):
    return 1-phi_dummy

def beta(phi_dummy):
    return (phys.rho_p-phys.rho_f)*phi_dummy*phys.g / vt / (1-phi_dummy)**(n-2)

def dpdz(phi_dummy):
    return (
        -phi_dummy*phys.rho_p*phys.g - frict[0]*eps(phi_dummy)*4/phys.D*phys.rho_f*Vf**2
    ) 

def NS_fluid(phi_dummy):
    return (
        -eps(phi_dummy)*dpdz(phi_dummy) + 
        eps(phi_dummy)*4/phys.D*(-frict[0]*phys.rho_f*Vf**2) - 
        beta(phi_dummy)*(Vf/(1-phi_dummy)-Vp/phi_dummy)
    )

def NS_particle(phi_dummy): # assume no particle friction to the wall
    return (
        -phi_dummy*dpdz(phi_dummy) - phi_dummy*phys.rho_p*phys.g + 
        beta(phi_dummy)*(Vf/(1-phi_dummy)-Vp/phi_dummy)
    )

def NS(phi_dummy):
    return NS_fluid(phi_dummy) - NS_particle(phi_dummy)

phi = sci.fsolve(NS, inguess_phi)
print("Fluid volume fraction:", phi[0])

print("Pressure gradient:", dpdz(phi))
<<<<<<< HEAD

#branchh
=======
#new
>>>>>>> a4f8f7031af340d10c8724c6a59af737da06c5d6









