#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#                  Written for GLiP project                     #
#                           by                                  # 
#  Gabriele SIAHAAN, Maral EBADI, Ehab bin Nauman & Wajih UDDIN #
#################################################################

# *************** Calculation for Reynolds number ************* #

import numpy as np

def reynolds(rho, velocity, D, mu):
    Re = rho*velocity*D/mu

    return Re

