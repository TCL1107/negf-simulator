# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 17:08:33 2025

@author: USER
"""

"""
examples/cnt_zonefold_bandstructure.py

Zone-folding bandstructure for CNT (n,m) using the Dirac approximation.

This script:
    - computes the low-energy subbands E_m(k_parallel)
    - estimates the minimal bandgap Eg
    - plots E(k) in a style similar to Lambin/Saito figures
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# 讓 examples 可以匯入 src/negf
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(ROOT, "src"))

from negf.cnt_zonefold import (
    zonefold_subbands,
    cnt_min_gap,
    cnt_diameter,
)


def plot_cnt_zonefold(n, m, Nk=201, m_range=6, k_max=0.30):
    k_par, E_plus, E_minus, Eg_min = zonefold_subbands(
        n, m, Nk=Nk, m_range=m_range, k_max=k_max
    )

    d = cnt_diameter(n, m)
    Eg_simple = cnt_min_gap(n, m)

    print(f"CNT ({n},{m})")
    print(f"  diameter d ≈ {d:.3f} nm")
    print(f"  Eg (zone-folding low-energy) ≈ {Eg_min:.3f} eV")
    print(f"  Eg (simple 2 γ0 a_cc / d ) ≈ {Eg_simple:.3f} eV")

    # 只畫能量比較低的部分，避免太亂
    Emax_plot = 2.5  # eV
    Emin_plot = -2.5

    plt.figure(figsize=(4, 6))

    # 畫所有子能帶
    for j in range(E_plus.shape[1]):
        plt.plot(k_par, E_plus[:, j], lw=0.8)
        plt.plot(k_par, E_minus[:, j], lw=0.8)

    plt.ylim(Emin_plot, Emax_plot)
    plt.xlabel(r"$k_{\parallel}$ (1/nm)")
    plt.ylabel("Energy (eV)")
    plt.title(f"CNT ({n},{m}) zone-folding\nEg ≈ {Eg_min:.3f} eV")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # 你可以在這裡換不同的 (n,m) 試：
    #   (17,0), (13,6), (10,10) ...
    for (n, m) in [(17, 0), (13, 6), (10, 10)]:
        plot_cnt_zonefold(n, m, Nk=301, m_range=4, k_max=0.25)
