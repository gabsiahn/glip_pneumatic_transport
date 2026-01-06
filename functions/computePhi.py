#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#                  Written for GLiP project                     #
#                           by                                  # 
#  Gabriele SIAHAAN, Maral EBADI, Ehab bin Nauman & Wajih UDDIN #
#################################################################

# ****** Solver to find phi using root-finding solver ********* #

# PHI IS SOLID VOLUME FRACTION AND EPSILON IS FLUID VOLUME FRACTION
# BASED ON JACKSON'S BOOK

import scipy.optimize as sci
import scipy.constants as const

# Description of the function
def computePhi(Vf,Vp,vt,n,frictCoeff,rho_p, rho_f, D, inguess_phi, bodyForceCoeff=1, fluidFrictCoeff=1, 
               partFrict=0, mu_p=None, partFrictCoeff=None):
    
    Vp = max(Vp, 1e-6)
    
    if partFrict == 1:
        if mu_p is None or partFrictCoeff is None:
            raise ValueError("Less input than expected when partFrict=1")
        
        def lossPartWallFrict(phi_dummy):

            particleWallStress = -partFrictCoeff*mu_p*rho_p*phi_dummy*(Vp/phi_dummy)**2

            return 4/D*particleWallStress

    def lossFluidWallFrict(phi_dummy):

        fluidWallStress = -frictCoeff*rho_f*Vf**2*fluidFrictCoeff/8 #test

        return epsilon(phi_dummy)*4/D*fluidWallStress
    
    def lossPartBodyForce(phi_dummy):
        return -phi_dummy*rho_p*const.g*bodyForceCoeff
    
    def epsilon(phi_dummy):
        return 1-phi_dummy
    
    def beta(phi_dummy):
        return (rho_p-rho_f)*phi_dummy*const.g / vt / (1-phi_dummy)**(n-2)
    
    def dpdz(phi_dummy):
        dpdz_total = lossPartBodyForce(phi_dummy) + lossFluidWallFrict(phi_dummy)

        if partFrict == 1: 
            dpdz_total += lossPartWallFrict(phi_dummy)

        return dpdz_total
    
    def NS_fluid(phi_dummy):
        return (
        -epsilon(phi_dummy)*dpdz(phi_dummy) + 
        lossFluidWallFrict(phi_dummy) - 
        beta(phi_dummy)*(Vf/(1-phi_dummy)-Vp/phi_dummy)
        )
    
    phi_soln = sci.fsolve(NS_fluid, inguess_phi)
    phi = phi_soln[0]
    resultsComputePhi =  {
        "phi" : phi,
        "epsilon" : epsilon(phi),
        "beta" : beta(phi),
        "dpdz" : dpdz(phi),
        "lossPartBodyForce" : lossPartBodyForce(phi),
        "lossFluidWallFrict" : lossFluidWallFrict(phi),
    }

    if partFrict == 1:
        resultsComputePhi["lossPartWallFrict"] = lossPartWallFrict(phi)

    return resultsComputePhi
