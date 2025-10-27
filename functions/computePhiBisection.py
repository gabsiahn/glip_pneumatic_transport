import scipy.constants as const

def computePhi(Vf, Vp, vt, n, frictCoeff, rho_p, rho_f, D, inguess_phi,
                        tol=1e-6, max_iter=100):
    
    wallFrictionFluid = -frictCoeff*rho_f*Vf**2

    def epsilon(phi_dummy):
        return 1 - phi_dummy

    def beta(phi_dummy):
        return (rho_p - rho_f) * phi_dummy * const.g / vt / (1 - phi_dummy) ** (n - 2)

    def dpdz(phi_dummy):
        return -phi_dummy * rho_p * const.g - frictCoeff * epsilon(phi_dummy) * 4 / D * rho_f * Vf**2

    def NS_fluid(phi_dummy):
        return -epsilon(phi_dummy) * dpdz(phi_dummy) + epsilon(phi_dummy) * 4 / D * wallFrictionFluid - beta(phi_dummy) * (Vf / (1 - phi_dummy) - Vp / phi_dummy)

    # --- Shooting solver using bisection ---
    phi_low = 0.0 + 1e-8  # avoid division by zero
    phi_high = 0.99       # avoid phi = 1
    phi_mid = inguess_phi

    for iteration in range(max_iter):
        NS_mid = NS_fluid(phi_mid)

        # Check if solution is close enough
        if abs(NS_mid) < tol:
            break

        NS_low = NS_fluid(phi_low)

        # Decide which half to keep
        if NS_low * NS_mid < 0:
            phi_high = phi_mid
        else:
            phi_low = phi_mid

        phi_mid = 0.5 * (phi_low + phi_high)

    phi = phi_mid

    return {
        "phi": phi,
        "epsilon": epsilon(phi),
        "dpdz": dpdz(phi),
        "beta": beta(phi)
    }
