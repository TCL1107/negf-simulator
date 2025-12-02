# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 04:23:26 2025

@author: USER
"""

# surface_gf.py
import numpy as np

def surface_gf_1d(E, eps0, t, eta=1e-6):
    """
    半無限一維最近鄰 tight-binding 鏈的 surface Green's function g(E).
    E: float，能量 (eV)
    eps0: float，on-site energy
    t: float，最近鄰 hopping (例如 -3.0 eV)
    eta: float，小的 +i*eta，避免數值發散
    
    回傳 complex scalar g(E)
    """
    # 加一個 i*eta 穩定
    z = E + 1j*eta - eps0  # (E + iη - ε0)
    # 導帶寬度相關項
    rad = z**2 - (2*t)**2  # = (E-ε0 + iη)^2 - 4 t^2

    # 我們需要選正確的根號支，讓 Im(g) <= 0 對 retarded GF
    root = np.sqrt(rad)

    # retarded convention: choose branch with negative imaginary part in g(E)
    # 我們可以透過這個 trick：如果 imag(root) > 0，就取 -root
    if np.imag(root) > 0:
        root = -root

    g = (z - root) / (2 * t**2)
    return g
