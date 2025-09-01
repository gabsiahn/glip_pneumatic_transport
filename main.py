#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#                  Written for GLiP project                     #
#                           by                                  # 
#  Gabriele SIAHAAN, Maral EBADI, Ehab bin Nauman & Wajih UDDIN #
#################################################################

#-----------------------Initialization-----------------------#
# Python packages #
import numpy as np
import pandas as pd
import os

# Module functions
from reynolds import reynolds
from computeBuoyancyFactor import computeBuoyancyFactor
from termvel import termvel
from colebrook import colebrook
from computePhi import computePhi


# Postprocessing by Wajih!

import input

g = 9.81 # gravitational acceleration, in m/s2 

#---------------------Flow parameters-----------------------#
# Pipe
A = np.pi/4*input.D**2 # cross-section area of the pipe, in m2

# Solid
m_part = input.mflux_part*A # solid mass flowrate, in kg/s, vector
Vp = m_part/input.rho_p/A # solid flow rate, in m/s, vector

lenVf = len(input.Vf)
lendp = len(input.dp)
lenVp = len(Vp)

# Terminal velocity of the particle
vt = np.zeros(lendp) # terminal velocity by Stokes' drag, in m/s, vector

for i in range(0, lendp):
    vt[i] = termvel(input.rho_p, input.rho_f, input.dp[i], input.mu)

#---------------------Dimensionless number------------------#
Re_pipe = np.zeros(lenVf) # Reynolds number of the fluid in the pipe, vector
Re_termvel = np.zeros(lendp) # Reynolds number of the fluid to the particle at its terminal vel

for i in range(0, lenVf):
    Re_pipe[i] = reynolds(input.rho_f, input.Vf[i], input.D, input.mu)

for i in range(0, lendp):
    Re_termvel[i] = reynolds(input.rho_f, vt[i], input.dp[i], input.mu)

#----------Calculation of 'n' for Buoyancy Calculation------#
n = np.zeros(lendp)

for i in range (0, lendp):
    n[i] = computeBuoyancyFactor(Re_termvel[i])

#---------------Find friction factor------------------------#
# Calculation using Colebrook-White formula

f = np.zeros(lenVf)

for i in range (0, lenVf):
    f[i] = colebrook(input.rough, input.D, Re_pipe[i], input.inguess_frict)

#--------------------Finding phi---------------------------#
data_results = []

for i, Vf_i in enumerate(input.Vf):
    f_i = f[i]

    for j, mflux_part_j in enumerate(input.mflux_part):
        Vp_j = Vp[j]
        
        for k, dp_k in enumerate(input.dp):
            n_k = n[k]
            vt_k = vt[k]
            
            result = computePhi(Vf_i,Vp_j,vt_k,n_k,f_i,input.rho_p,input.rho_f,input.D, input.inguess_phi)
            
            phi, dpdz, epsilon, beta = result["phi"], result["dpdz"], result["epsilon"], result["beta"]
            
            data_results.append([Vf_i, mflux_part_j, Vp_j, dp_k, vt_k, n_k, f_i, phi, dpdz, epsilon, beta])

results = pd.DataFrame(
    data_results,
    columns = ["Vf","mflux_part","Vp","dp","vt","n","f","phi","dpdz","epsilon","beta"]
    )

#-------------------Store the result-----------------------#
def save_results(results, base_filename):
    filename = base_filename
    name, ext = os.path.splitext(base_filename)
    counter = 1

    while os.path.exists(filename):
        filename = f"{name}_{counter}{ext}"
        counter += 1
        
    results.to_csv(filename, index=False)

save_results(results, input.output_filename)










