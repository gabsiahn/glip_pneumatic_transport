#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#                  Written for GLiP project                     #
#                           by                                  # 
#  Gabriele SIAHAAN, Maral EBADI, Ehab bin Nauman & Wajih UDDIN #
#################################################################

# ******** Calculation of 'n' for Buoyancy Calculation ******** #

def computeBuoyancyFactor(Re_termvel):
    if Re_termvel < 0.2:
        n = 4.65 
    elif 0.2 < Re_termvel < 1:
        n = 4.35*Re_termvel**(-0.03)
    elif 1 < Re_termvel < 500:
        n = 4.45*Re_termvel**(-0.1)
    else:
        n = 2.39
              
    return n


