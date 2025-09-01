#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#                  Written for GLiP project                     #
#                           by                                  # 
#  Gabriele SIAHAAN, Maral EBADI, Ehab bin Nauman & Wajih UDDIN #
#################################################################

# ************ Function to find friction factor ****************# 
# ************** using Colebrook-White formula **************** #

import numpy as np
import scipy.optimize as sci

def colebrook(roughness,D,Re_pipe,inguess_frict):

    def computeColebrook(f):
        return 1/np.sqrt(f) + 2*np.log10(roughness/3.7/D + 2.51/Re_pipe/np.sqrt(f))
    
    frict = sci.fsolve(computeColebrook, inguess_frict)
    return frict[0]

