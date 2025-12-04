# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 17:17:42 2025

@author: USER
"""

"""
examples/cnt_gap_to_chain.py

Map a CNT (n,m) bandgap Eg_CNT (from zone-folding)
to an effective 1D dimerized chain used in NEGF.

We:
  - compute Eg_CNT using cnt_zonefold.cnt_min_gap()
  - build a 1D chain with staggered on-site potential ±Δ
    so that Eg_chain ≈ 2Δ ≈ Eg_CNT
  - compute T(E) via a guaranteed-working NEGF routine
    with analytic 1D lead self-energies
  - extract Eg_chain from T(E) and compare to Eg_CNT
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# --- allow importing src/negf ---
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(ROOT, "src"))

from negf.cnt_zonefold import cnt_min_gap, cnt_diameter


# =====================================================
# 1D semi-infinite chain self-energy (analytic)
# =====================================================

def sigma_1d(E, eps0, t):
    """
    Analytic retarded self-energy for a semi-infinite 1D chain:

    Σ(E) = (E - eps0)/2 - i * sqrt(4 t^2 - (E-eps0)^2)/2
    """
    z = E - eps0
    inside = 4*t*t - z*z
    inside = np.where(inside >= 0, inside, 0.0)
    imag = -np.sqrt(inside) / 2.0
    real = z / 2.0
    return real + 1j * imag


def build_staggered_chain(N, t, Delta):
    """
    Build Hamiltonian for a 1D chain with alternating on-site ±Δ:

        H_ii = +Δ, -Δ, +Δ, -Δ, ...
        H_{i,i+1} = t
    """
    H = np.zeros((N, N), dtype=complex)
    for i in range(N):
        H[i, i] = Delta if (i % 2 == 0) else -Delta
        if i < N - 1:
            H[i, i+1] = t
            H[i+1, i] = t
    return H


def transmission_1d(E, H, eps_lead, t_lead):
    """
    NEGF transmission for 1D chain with analytic self-energies.
    """
    N = H.shape[0]
    I = np.eye(N, dtype=complex)
    Tvals = []

    for Ei in E:
        SigmaL = np.zeros((N, N), dtype=complex)
        SigmaR = np.zeros((N, N), dtype=complex)

        sL = sigma_1d(Ei, eps_lead, t_lead)
        sR = sigma_1d(Ei, eps_lead, t_lead)
        SigmaL[0, 0] = sL
        SigmaR[-1, -1] = sR

        G = np.linalg.inv((Ei + 1e-12j) * I - H - SigmaL - SigmaR)
        GammaL = 1j * (SigmaL - SigmaL.conj().T)
        GammaR = 1j * (SigmaR - SigmaR.conj().T)

        Tvals.append(np.real(np.trace(GammaL @ G @ GammaR @ G.conj().T)))

    return np.array(Tvals)


def extract_gap_from_T(E, T, T_th=1e-3):
    """
    Roughly estimate gap Eg_chain from T(E) by finding the largest
    region around E=0 where T < T_th.
    """
    mask = (T < T_th)
    idx0 = np.argmin(np.abs(E))

    if not mask[idx0]:
        return 0.0, None, None

    i_left = idx0
    while i_left > 0 and mask[i_left]:
        i_left -= 1
    E_low = E[i_left + 1]

    i_right = idx0
    while i_right < len(E) - 1 and mask[i_right]:
        i_right += 1
    E_high = E[i_right - 1]

    Eg = E_high - E_low
    return Eg, E_low, E_high


if __name__ == "__main__":
    # --- choose a CNT (n,m) to map ---
    n, m = 17, 0   # you can change to (13,6), etc.

    d = cnt_diameter(n, m)
    Eg_cnt = cnt_min_gap(n, m)  # eV

    print(f"CNT ({n},{m})")
    print(f"  diameter d ≈ {d:.3f} nm")
    print(f"  Eg_CNT (zone-folding) ≈ {Eg_cnt:.3f} eV")

    # --- build effective 1D chain ---
    t_chain = -2.7   # eV, hopping
    Delta = Eg_cnt / 2.0   # so that Eg_chain ≈ 2Δ ≈ Eg_CNT
    N_sites = 200

    Hc = build_staggered_chain(N_sites, t_chain, Delta)

    # --- NEGF transport ---
    E = np.linspace(-2.0, 2.0, 801)
    Tvals = transmission_1d(E, Hc, eps_lead=0.0, t_lead=t_chain)

    Eg_chain, E_low, E_high = extract_gap_from_T(E, Tvals, T_th=1e-3)

    print(f"  Eg_chain (from T(E)) ≈ {Eg_chain:.3f} eV")
    if E_low is not None:
        print(f"    gap edges: [{E_low:.3f}, {E_high:.3f}] eV")

    # --- plot ---
    plt.figure(figsize=(5, 4))
    plt.plot(E, Tvals, label="T(E)")
    if E_low is not None:
        plt.axvspan(E_low, E_high, alpha=0.2, label="gap region")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Transmission T(E)")
    plt.title(f"1D chain calibrated to CNT ({n},{m}) gap\n"
              f"Eg_CNT ≈ {Eg_cnt:.3f} eV, Eg_chain ≈ {Eg_chain:.3f} eV")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("CNT_TE_true.png", dpi=300)
    plt.show()
