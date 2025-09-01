#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#                  Written for GLiP project                     #
#                           by                                  # 
#  Gabriele SIAHAAN, Maral EBADI, Ehab bin Nauman & Wajih UDDIN #
#################################################################

# **************** Compute terminal velocity *******************#

import scipy.constants as const

def termvel(rho_p, rho_f, dp, mu):
    u_terminal = 2/9*(rho_p-rho_f)*const.g*dp**2/4/mu

    return u_terminal




