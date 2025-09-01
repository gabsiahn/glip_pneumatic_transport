#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#                  Written for GLiP project                     #
#                           by                                  # 
#  Gabriele SIAHAAN, Maral EBADI, Ehab bin Nauman & Wajih UDDIN #
#################################################################

# ****** Solver to find phi using root-finding solver ********* #

import scipy.optimize as sci
import scipy.constants as const

# Description of the function
def computePhi(Vf,Vp,vt,n,frictCoeff,rho_p, rho_f, D, inguess_phi):
    
    wallFrictionFluid = -frictCoeff*rho_f*Vf**2
    
    def epsilon(phi_dummy):
        return 1-phi_dummy
    
    def beta(phi_dummy):
        return (rho_p-rho_f)*phi_dummy*const.g / vt / (1-phi_dummy)**(n-2)
    
    def dpdz(phi_dummy):
        return (
        -phi_dummy*rho_p*const.g - frictCoeff*epsilon(phi_dummy)*4/D*rho_f*Vf**2
        ) 
    
    def NS_fluid(phi_dummy):
        return (
        -epsilon(phi_dummy)*dpdz(phi_dummy) + 
        epsilon(phi_dummy)*4/D*wallFrictionFluid - 
        beta(phi_dummy)*(Vf/(1-phi_dummy)-Vp/phi_dummy)
        )
    
    # def NS_particle(phi_dummy): # assume no particle friction to the wall
    #     return (
    #     -phi_dummy*dpdz(phi_dummy) - phi_dummy*rho_p*9.81 + 
    #     beta(phi_dummy)*(Vf/(1-phi_dummy)-Vp/phi_dummy)
    #     )
    
    # def NS(phi_dummy):
    #     return NS_fluid(phi_dummy) - NS_particle(phi_dummy)
    
    phi_soln = sci.fsolve(NS_fluid, inguess_phi)
    phi = phi_soln[0]
    return {
        "phi" : phi,
        "epsilon" : epsilon(phi),
        "dpdz" : dpdz(phi),
        "beta" : beta(phi)
    }
