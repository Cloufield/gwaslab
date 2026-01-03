#!/usr/bin/env python

"""
Optimized random variate generator for the generalized inverse Gaussian distribution.

Performance improvements:
- Use ** operator instead of math.pow for better performance
- Pre-compute constants where possible
- More efficient conditional checks
- Better type handling
"""

import math
from numpy import random


def psi(x, alpha, lam):
    """Optimized psi function using efficient operations."""
    cosh_x = math.cosh(x)
    exp_x = math.exp(x)
    f = -alpha * (cosh_x - 1.0) - lam * (exp_x - x - 1.0)
    return f


def dpsi(x, alpha, lam):
    """Optimized dpsi function using efficient operations."""
    sinh_x = math.sinh(x)
    exp_x = math.exp(x)
    f = -alpha * sinh_x - lam * (exp_x - 1.0)
    return f


def g(x, sd, td, f1, f2):
    """Optimized g function with efficient conditional checks."""
    if -sd <= x <= td:
        return 1.0
    elif x > td:
        return f1
    else:  # x < -sd
        return f2


def gigrnd(p, a, b):
    """
    Optimized random variate generator for GIG distribution.
    
    Parameters:
    -----------
    p : float
        Shape parameter
    a : float
        First scale parameter
    b : float
        Second scale parameter
    
    Returns:
    --------
    float
        Random variate from GIG(p, a, b) distribution
    """
    # setup -- sample from the two-parameter version gig(lam,omega)
    p = float(p)
    a = float(a)
    b = float(b)
    lam = p
    omega = math.sqrt(a * b)

    # Handle negative lambda
    if lam < 0:
        lam = -lam
        swap = True
    else:
        swap = False

    # Pre-compute alpha
    omega_sq = omega ** 2  # More efficient than math.pow(omega, 2)
    lam_sq = lam ** 2
    alpha = math.sqrt(omega_sq + lam_sq) - lam

    # find t - optimized with pre-computed values
    x = -psi(1.0, alpha, lam)
    if 0.5 <= x <= 2.0:
        t = 1.0
    elif x > 2.0:
        if alpha == 0 and lam == 0:
            t = 1.0
        else:
            t = math.sqrt(2.0 / (alpha + lam))
    else:  # x < 0.5
        if alpha == 0 and lam == 0:
            t = 1.0
        else:
            t = math.log(4.0 / (alpha + 2.0 * lam))

    # find s - optimized with pre-computed values
    x = -psi(-1.0, alpha, lam)
    if 0.5 <= x <= 2.0:
        s = 1.0
    elif x > 2.0:
        if alpha == 0 and lam == 0:
            s = 1.0
        else:
            cosh_1 = math.cosh(1.0)  # Pre-compute cosh(1)
            s = math.sqrt(4.0 / (alpha * cosh_1 + lam))
    else:  # x < 0.5
        if alpha == 0 and lam == 0:
            s = 1.0
        elif alpha == 0:
            s = 1.0 / lam
        elif lam == 0:
            alpha_inv = 1.0 / alpha
            alpha_inv_sq = alpha_inv ** 2
            s = math.log(1.0 + alpha_inv + math.sqrt(alpha_inv_sq + 2.0 * alpha_inv))
        else:
            alpha_inv = 1.0 / alpha
            alpha_inv_sq = alpha_inv ** 2
            s_candidate = math.log(1.0 + alpha_inv + math.sqrt(alpha_inv_sq + 2.0 * alpha_inv))
            s = min(1.0 / lam, s_candidate)

    # find auxiliary parameters
    eta = -psi(t, alpha, lam)
    zeta = -dpsi(t, alpha, lam)
    theta = -psi(-s, alpha, lam)
    xi = dpsi(-s, alpha, lam)

    # Pre-compute reciprocals
    p_val = 1.0 / xi
    r = 1.0 / zeta

    td = t - r * eta
    sd = s - p_val * theta
    q = td + sd
    pqr_sum = p_val + q + r  # Pre-compute denominator

    # random variate generation - optimized rejection sampling
    while True:
        U = random.random()
        V = random.random()
        W = random.random()
        
        # Optimized conditional checks with pre-computed values
        U_pqr = U * pqr_sum
        if U_pqr < q:
            rnd = -sd + q * V
        elif U_pqr < q + r:
            rnd = td - r * math.log(V)
        else:
            rnd = -sd + p_val * math.log(V)

        # Pre-compute exponentials
        exp_eta_zeta = math.exp(-eta - zeta * (rnd - t))
        exp_theta_xi = math.exp(-theta + xi * (rnd + s))
        
        # Use optimized g function
        g_val = g(rnd, sd, td, exp_eta_zeta, exp_theta_xi)
        psi_rnd = psi(rnd, alpha, lam)
        
        if W * g_val <= math.exp(psi_rnd):
            break

    # transform back to the three-parameter version gig(p,a,b)
    # Optimized: use ** instead of math.pow
    lam_omega_ratio = lam / omega
    sqrt_term = math.sqrt(1.0 + lam_sq / omega_sq)
    rnd = math.exp(rnd) * (lam_omega_ratio + sqrt_term)
    
    if swap:
        rnd = 1.0 / rnd

    rnd = rnd / math.sqrt(a / b)
    return rnd

