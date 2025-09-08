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
def computePhi(Vf,Vp,vt,n,frictCoeff,rho_p, rho_f, D, inguess_phi, partFrict=0, mu_p=None, partFrictConst=None):
    
    if partFrict == 1:
        if mu_p is None or partFrictConst is None:
            raise ValueError("Less input than expected when partFrict=1")
        
        def particleWallStress(phi_dummy):
            return -partFrictConst*mu_p*rho_p*phi_dummy*(Vp/phi_dummy)**2
    
    def fluidWallStress(phi_dummy):
        return -epsilon(phi_dummy)*frictCoeff*rho_f*Vf**2
    
    def epsilon(phi_dummy):
        return 1-phi_dummy
    
    def beta(phi_dummy):
        return (rho_p-rho_f)*phi_dummy*const.g / vt / (1-phi_dummy)**(n-2)
    
    def dpdz(phi_dummy):
        dpdz_total = -phi_dummy*rho_p*const.g + 4/D*fluidWallStress(phi_dummy)

        if partFrict == 1: 
            dpdz_total += 4/D*particleWallStress(phi_dummy)

        return dpdz_total
    
    def NS_fluid(phi_dummy):
        return (
        -epsilon(phi_dummy)*dpdz(phi_dummy) + 
        4/D*fluidWallStress(phi_dummy) - 
        beta(phi_dummy)*(Vf/(1-phi_dummy)-Vp/phi_dummy)
        )
    
    # def NS_particle(phi_dummy): # assume no particle friction to the wall
    #     return (
    #     -phi_dummy*dpdz(phi_dummy) - phi_dummy*rho_p*9.81 + 
    #     beta(phi_dummy)*(Vf/(1-phi_dummy)-Vp/phi_dummy)
    #     )
    
    phi_soln = sci.fsolve(NS_fluid, inguess_phi)
    phi = phi_soln[0]
    resultsComputePhi =  {
        "phi" : phi,
        "epsilon" : epsilon(phi),
        "dpdz" : dpdz(phi),
        "beta" : beta(phi),
        "fluidWallStress" : fluidWallStress(phi),
    }

    if partFrict == 1:
        resultsComputePhi["particleWallStress"] = particleWallStress(phi)

    return resultsComputePhi
