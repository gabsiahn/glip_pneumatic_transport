import numpy as np

def beta_drag(phi, Vf, Vp, mu, rho, d_p, C_D=None, phi_switch=0.8, tiny=1e-12):
    """
    Piecewise drag/momentum-exchange coefficient beta.

    Parameters
    ----------
    phi : float or array
        Void fraction (gas volume fraction) ϕ.
    Vf : (3,) or (...,3) array
        Fluid velocity vector(s).
    Vp : (3,) or (...,3) array
        Particle/solid velocity vector(s).
    mu : float or array broadcastable
        Fluid dynamic viscosity μ [Pa·s].
    rho : float or array broadcastable
        Fluid density ρ [kg/m^3].
    d_p : float or array broadcastable
        Particle diameter d_p [m].
    C_D : float or array broadcastable (required when phi >= phi_switch)
        Drag coefficient.
    phi_switch : float
        Switching value (0.8 in your model).
    tiny : float
        Small number for numerical safety.

    Returns
    -------
    beta : array
        Momentum-exchange coefficient β.
    """
    phi = np.asarray(phi, dtype=float)
    Vf  = np.asarray(Vf, dtype=float)
    Vp  = np.asarray(Vp, dtype=float)

    slip = Vp - Vf
    slip_mag = np.linalg.norm(slip, axis=-1)

    phi_safe = np.maximum(phi, tiny)
    dp_safe  = np.maximum(d_p, tiny)

    mask = phi >= phi_switch

    # φ < 0.8 branch (Ergun-type)
    beta_low = (150.0 * mu * (1.0 - phi_safe)**2) / (phi_safe * dp_safe**2) \
               + (1.75 * (1.0 - phi_safe) * rho / dp_safe) * slip_mag

    # φ ≥ 0.8 branch (dilute suspension)
    if np.any(mask):
        if C_D is None:
            raise ValueError("C_D must be provided for phi >= 0.8.")
        beta_high = (3.0 / (4.0 * dp_safe)) * C_D * (1.0 - phi_safe) \
                    * (phi_safe**(-1.65)) * rho * slip_mag
        beta = np.where(mask, beta_high, beta_low)
    else:
        beta = beta_low

    return beta
import numpy as np

def total_force_drag_pressure(
    Vf, Vp, rho, mu, d_p, C_D, A_perp, epsilon=1.0
):
    """
    F = 0.5 * rho * |Vf - Vp| * C_D * A_perp * epsilon^(1-chi) * (Vf - Vp)

    with:
      Re_p = rho * d_p * |Vf - Vp| / mu
      chi  = 3.7 - 0.65 * exp( - (1.5 - log(Re_p))^2 / 2 )

    Parameters
    ----------
    Vf : (3,) or (...,3) array
        Fluid velocity vector(s) [m/s]
    Vp : (3,) or (...,3) array
        Particle velocity vector(s) [m/s]
    rho : float or array broadcastable
        Fluid density [kg/m^3]
    mu : float or array broadcastable
        Fluid dynamic viscosity [Pa·s]
    d_p : float or array broadcastable
        Particle diameter [m]
    C_D : float or array broadcastable
        Drag coefficient
    A_perp : float or array broadcastable
        Projected area normal to slip direction [m^2]
    epsilon : float or array broadcastable
        ε (must be > 0)

    Returns
    -------
    F : (...,3) array
        Total force vector(s) [N]
    chi : (...) array
        Chi value(s) used (handy for debugging/plotting)
    Re_p : (...) array
        Particle Reynolds number(s)
    """
    Vf = np.asarray(Vf, dtype=float)
    Vp = np.asarray(Vp, dtype=float)

    slip = Vf - Vp
    slip_mag = np.linalg.norm(slip, axis=-1)

    rho = np.asarray(rho, dtype=float)
    mu = np.asarray(mu, dtype=float)
    d_p = np.asarray(d_p, dtype=float)
    epsilon = np.asarray(epsilon, dtype=float)

    if np.any(mu <= 0.0):
        raise ValueError("mu must be > 0.")
    if np.any(d_p <= 0.0):
        raise ValueError("d_p must be > 0.")
    if np.any(epsilon <= 0.0):
        raise ValueError("epsilon must be > 0 (required for fractional powers).")

    # Re_p (particle Reynolds number)
    Re_p = rho * d_p * slip_mag / mu

    # chi(Re_p) exactly as given
    if np.any(Re_p <= 0.0):
        raise ValueError("Re_p must be > 0 (check rho, d_p, mu, or velocities).")
    chi = 3.7 - 0.65 * np.exp(-((1.5 - np.log(Re_p))**2) / 2.0)

    pref = 0.5 * rho * slip_mag * C_D * A_perp * (epsilon ** (1.0 - chi))
    F = pref[..., None] * slip
    return F, chi, Re_p
