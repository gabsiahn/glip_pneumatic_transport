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
from functions.reynolds import reynolds
from functions.computeBuoyancyFactor import computeBuoyancyFactor
from functions.termvel import termvel
from functions.colebrook import colebrook
from functions.computePhi import computePhi

# Insert inputs
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
            
            result = computePhi(Vf_i,Vp_j,vt_k,n_k,f_i,input.rho_p,input.rho_f,input.D, 
                                input.inguess_phi, input.partFrict, input.mu_p, input.partFrictConst)
            
       
            phi, dpdz, epsilon, beta, fluidWallStress = result["phi"], result["dpdz"], result["epsilon"], result["beta"], result["fluidWallStress"]

            dpdz_inv = -dpdz
            fluidWallStress_inv = -fluidWallStress
            
            particleWallStress = result.get("particleWallStress", None)
            particleWallStress_inv = -particleWallStress if particleWallStress is not None else None
            
            data_results.append([Vf_i, mflux_part_j, Vp_j, dp_k, vt_k, n_k, f_i, phi, dpdz, epsilon, beta, fluidWallStress, particleWallStress,
                                 dpdz_inv, fluidWallStress_inv, particleWallStress_inv])

results = pd.DataFrame(
    data_results,
    columns = ["Vf","mflux_part","Vp","dp","vt","n","f","phi","dpdz","epsilon","beta", "fluidWallStress", "particleWallStress", "dpdz_inv",
               "fluidWallStress_inv", "particleWallStress_inv"]
    )

#-------------------Store the result-----------------------#
results_folder = os.path.join("results", "model")
os.makedirs(results_folder, exist_ok=True)


def save_results(results, base_filename):
    # filename = base_filename
    name, ext = os.path.splitext(base_filename)
    filename = os.path.join(results_folder, base_filename)
    counter = 1

    while os.path.exists(filename):
        filename = os.path.join(results_folder, f"{name}_{counter}{ext}")
        counter += 1
        
    results.to_csv(filename, index=False)
    print(f"Saved results to: {filename}")

save_results(results, input.output_filename)


#branchh
#new










