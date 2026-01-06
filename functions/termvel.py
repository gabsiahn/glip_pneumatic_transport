#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#                  Written for GLiP project                     #
#                           by                                  # 
#  Gabriele SIAHAAN, Maral EBADI, Ehab bin Nauman & Wajih UDDIN #
#################################################################

# **************** Compute terminal velocity *******************#

# Computer terminal velocity under Stoke's regime

import scipy.constants as const
import numpy as np
from functions.reynolds import reynolds
from scipy.optimize import fsolve

def termvelStokes(Vf, rho_p, rho_f, dp, mu):
    stokes = 2/9*(rho_p-rho_f)*const.g*dp**2/4/mu # + Vf

    return stokes

def termvelAllen(Vf, rho_p, rho_f, dp, mu_f):
    
    # We consider 1D velocity in y-direction. (+) represents an upward direction and (-) represents the opposite.
    
    vp = np.pi*dp**3/6
    Ap = np.pi*dp**2/4

    allen = ( (2*(rho_p-rho_f)*vp*const.g*dp**0.6) / (18.5*rho_f**0.4*Ap*mu_f**0.6) )**(5/7) # + Vf

    return allen

def termvelNewton(Vf, rho_p, rho_f, dp):

    vp = np.pi*dp**3/6
    Ap = np.pi*dp**2/4

    newton = (2*(rho_p-rho_f)*vp*const.g / 0.44 / rho_f / Ap)**0.5 # + Vf

    return newton

def termVel(Vf, rho_p, rho_f, dp, mu_f):

    Vt = np.zeros(3)
    Re = np.zeros(3)
    
    Vt[0] = termvelStokes(Vf, rho_p, rho_f, dp, mu_f)
    Vt[1] = termvelAllen(Vf, rho_p, rho_f, dp, mu_f)
    Vt[2] = termvelNewton(Vf, rho_p, rho_f, dp)

    flowRegimes = ["Stokes", "Allen", "Newton"]
    Re_ranges = [(1e-4, 1), (1, 800), (1e3, 5e5)]
    
    for i, regime in enumerate(flowRegimes):
        # V = abs(Vt[i]-Vf)
        Re[i] = reynolds(rho_f, Vt[i], dp, mu_f)
        low, up = Re_ranges[i]

        if (low <= Re[i] < up):
            soln = (Vt[i], Re[i], regime)
            break

    if soln is None:
        raise ValueError(f"No valid regime found. Re = {Re}, Vt = {Vt}")

    return soln       
    


def termvel(Vf, rho_p, rho_f, dp, mu_f, inguessVp):
    
    # We consider 1D velocity in y-direction. (+) represents an upward direction and (-) represents the opposite.
    
    def Vrel(Vp):
        return Vp - Vf

    def Vmag(Vp):
        return max(((Vrel(Vp))**2)**0.5, 1e-10)
    
    def Re_p(Vp):

        Re_min = 1e-10 # to avoid Re = 0
        Re_calc = reynolds(rho_f, Vp, dp, mu_f)

        Re = max (Re_min, Re_calc)

        return Re
    
    def Cd_Stokes(Vp): # Stokes regime
        return 24/Re_p(Vp) 
    
    def Cd_Schiller(Vp): #Schiller-Nauman
        return 24/Re_p(Vp) * (1+0.15*Re_p(Vp)**0.687) 
    
    def Cd_Newton(Vp): # Newton regime
        return 0.44
    
    flowRegimes = {"Stokes": (Cd_Stokes,(1e-4,1)), 
                   "Schiller-Nauman": (Cd_Schiller,(1,1e3)),
                   "Newton": (Cd_Newton,(1e3,5e5))
                   }

    vp = np.pi*dp**3/6
    Ap = np.pi*dp**2/4

    soln = None

    for regime, (Cd_func, Re_range) in flowRegimes.items():
        
        def eq(Vp):
            return (2*(rho_p-rho_f)*vp*const.g/rho_f/Ap) - (Vmag(Vp)**2*Cd_func(Vp))
        
        try:
            Vp_sol = fsolve(eq, inguessVp)[0]
            Re_sol = Re_p(Vp_sol)

            # Check the proper regime
            if (Re_range[0] <= Re_sol < Re_range[1]):
                soln = (Vp_sol, Re_sol, regime)
                break

        except Exception:
            continue

    if soln is None:
        raise RuntimeError("No valid terminal velocity found in any of defined regimes.")

    return soln




    

        
    
    











