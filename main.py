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
from functions.termvel import termvel, termVel
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
# lenVp = len(Vp)

# Calculate Reynolds number of particle
Re_particle = reynolds(input.rho_f, Vp, input.dp, input.mu)

# Terminal velocity of the particle, flow regimes, Re at terminal velocity and buoyancy factor ('n')
vt = np.zeros((lenVf, lendp))
Re_termvel = np.zeros((lenVf, lendp))
regime = np.empty((lenVf, lendp), dtype=object)

n = np.zeros(lendp)

for i in range(lenVf):
    for j in range(lendp):
        solnTermVel = termVel(input.Vf[i], input.rho_p, input.rho_f, input.dp[j], input.mu)

        vt[i, j] = solnTermVel[0]
        Re_termvel[i, j] = solnTermVel[1]
        regime[i, j] = solnTermVel[2]

        n[j] = computeBuoyancyFactor(Re_termvel[i,j])

#---------------------Dimensionless number------------------#
Re_pipe = np.zeros(lenVf) # Reynolds number of the fluid in the pipe, vector

for i in range(0, lenVf):
    Re_pipe[i] = reynolds(input.rho_f, input.Vf[i], input.D, input.mu)

#----------Calculation of 'n' for Buoyancy Calculation------#
# n = np.zeros(lendp)

# for i in range (0, lendp):
#     n[i] = computeBuoyancyFactor(Re_termvel[i])

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

            vt_ik = vt[i,k]
            Re_vt_ik = Re_termvel[i,k]
            regime_ik = regime[i,k]
            
            result = computePhi(Vf_i,Vp_j,vt_ik,n_k,f_i,input.rho_p,input.rho_f,input.D, 
                                input.inguess_phi, input.bodyForceCoeff, input.fluidFrictCoeff,
                                input.partFrict, input.mu_p, input.partFrictCoeff)
            
       
            phi, epsilon, beta = result["phi"], result["epsilon"], result["beta"]

            dpdz, lossPartBodyForce, lossFluidWallFrict = -result["dpdz"], -result["lossPartBodyForce"], -result["lossFluidWallFrict"]
            
            lossPartWallFrict = result.get("lossPartWallFrict", None)
            lossPartWallFrict = -lossPartWallFrict if lossPartWallFrict is not None else None

            lossTotalFrict = lossFluidWallFrict + lossPartWallFrict if lossPartWallFrict is not None else None

            Vp_phi = Vp_j/phi

            Vsl = Vf_i-Vp_phi
            
            data_results.append([Vf_i, Re_pipe[i], mflux_part_j, Vp_phi, Vsl,dp_k, vt_ik, Re_vt_ik, regime_ik, n_k, f_i, phi, epsilon, beta, 
                                 dpdz, lossPartBodyForce, lossFluidWallFrict, lossPartWallFrict, lossTotalFrict])

results = pd.DataFrame(
    data_results,
    columns = ["Vf","Re_pipe","mflux_part","Vp","Vsl","dp","vt","Re_vt","flowRegime","n","f","phi","epsilon",
               "beta","dpdz","lossPartBodyForce","lossFluidWallFrict","lossPartWallFrict", "lossTotalFrict"]
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










