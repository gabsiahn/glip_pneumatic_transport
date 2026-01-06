# test

import numpy as np
import input  # renamed from input
from functions.termvel import termvel, termvelAllen, termVel
import pandas as pd

# Pipe and particle properties
A = np.pi / 4 * input.D**2  # cross-section area
m_part = input.mflux_part * A  # solid mass flowrate
Vp = m_part / input.rho_p / A  # solid flow velocity

lenVf = len(input.Vf)
lendp = len(input.dp)

# Initialize array to store terminal velocities
vt = np.zeros((lenVf, lendp))

rows = []

# Loop over fluid velocities and particle diameters
for i in range(lenVf):
    for j in range(lendp):
        soln = termVel(
            input.Vf[i],
            input.rho_p,
            input.rho_f,
            input.dp[j],
            input.mu
        )
        
        rows.append({
            "Vf": input.Vf[i],
            "dp": input.dp[j],
            "Vt": soln[0],
            "Re": soln[1],
            "Regime": soln[2]
        })

# Convert to DataFrame
df = pd.DataFrame(rows)

# Save to CSV
df.to_csv("terminal_velocities_test_new_2.csv", index=False)






