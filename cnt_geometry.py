"""
cnt_geometry.py

Basic CNT geometry and analytic bandgap model based on
the graphene zone-folding picture (Lambin / Saito style).

Eg ≈ 2 * gamma0 * d_CC / d
"""

import numpy as np

# Graphene / CNT constants (literature values)
D_CC = 0.142  # nm, C–C bond length
GAMMA0 = 2.9  # eV, nearest-neighbor π-orbital hopping


def cnt_diameter(n: int, m: int) -> float:
    """
    Return CNT diameter in nm for chirality (n, m).

    d = (a_cc / pi) * sqrt(n^2 + m^2 + n*m)
    """
    return (D_CC / np.pi) * np.sqrt(n**2 + m**2 + n * m)


def cnt_bandgap_nm(
    n: int,
    m: int,
    gamma0: float = GAMMA0,
    d_cc: float = D_CC,
) -> float:
    """
    Approximate semiconducting CNT bandgap Eg (eV) using:

        Eg ≈ 2 * gamma0 * d_CC / d

    This is the standard 1/d scaling from the zone-folding model.
    """
    d = cnt_diameter(n, m)
    return 2.0 * gamma0 * d_cc / d
