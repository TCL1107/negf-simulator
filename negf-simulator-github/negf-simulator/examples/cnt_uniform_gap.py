"""
cnt_uniform_gap.py

Toy example:
    - Build a uniform 1D chain with a staggered on-site potential
      to create a band gap (like a simple dimerized semiconductor).
    - Use the existing NEGF core:
          from negf.negf import transmission
          from negf.surface_gf import surface_gf_1d
      to compute T(E).
    - Automatically extract the central energy gap Eg_num from T(E).

This is NOT a full CNT model, but a minimal working example that
demonstrates how to go from T(E) to an estimated band gap, using
only the public NEGF API (transmission, surface_gf_1d).
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# ---- Make sure "negf" package is importable when run from examples/ ----
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(ROOT, "src"))

from negf.negf import transmission
from negf.surface_gf import surface_gf_1d


# ======================
# 1. Physical parameters
# ======================

# Base hopping (CNT-like scale)
t = -2.9  # eV

# Staggered on-site potential: (+Delta, -Delta, +Delta, -Delta, ...)
# This opens a band gap ~ 2*Delta in a dimerized 1D chain.
Delta = 0.3  # eV

Nx = 200  # number of sites
a = 0.246  # nm, lattice spacing (for length L = Nx * a if needed)

# Energy grid
E_min, E_max = -2.0, 2.0  # eV
nE = 801
E = np.linspace(E_min, E_max, nE)

# "Almost zero" transmission threshold
T_th = 1e-4


# ============================
# 2. Build the Hamiltonian Hc
# ============================

Hc = np.zeros((Nx, Nx), dtype=complex)

for i in range(Nx):
    # Staggered on-site: +Delta, -Delta, +Delta, ...
    Hc[i, i] = Delta if (i % 2 == 0) else -Delta
    if i < Nx - 1:
        Hc[i, i + 1] = t
        Hc[i + 1, i] = t


# ==============================
# 3. Build lead self-energies
# ==============================

# For now we use simple 1D semi-infinite leads with the same eps0, t,
# and do NOT include the staggering in the leads (just a simple metal-like lead).
eps_lead = 0.0

SigmaL = surface_gf_1d(E, eps_lead, t)
SigmaR = surface_gf_1d(E, eps_lead, t)

# ============================
# 4. Compute T(E)
# ============================

T = transmission(E, Hc, SigmaL, SigmaR)  # (nE,)


# ============================
# 5. Extract central band gap
# ============================

mask_gap = (T < T_th)
idx0 = np.argmin(np.abs(E))  # index closest to E = 0

Eg_num = None
E_low = None
E_high = None

if not mask_gap[idx0]:
    print("[Warning] T(E=0) is not in a low-transmission region.")
    print("         This likely means the toy model is not gapped at E=0.")
else:
    # Walk left to find the start of the gap
    i_left = idx0
    while i_left > 0 and mask_gap[i_left]:
        i_left -= 1
    E_low = E[i_left + 1]

    # Walk right to find the end of the gap
    i_right = idx0
    while i_right < nE - 1 and mask_gap[i_right]:
        i_right += 1
    E_high = E[i_right - 1]

    Eg_num = E_high - E_low

    print("=======================================")
    print("Toy dimerized 1D chain: band gap from T(E)")
    print("---------------------------------------")
    print(f"E_low  ≈ {E_low:.4f} eV")
    print(f"E_high ≈ {E_high:.4f} eV")
    print(f"Eg_num ≈ {Eg_num:.4f} eV")
    print("=======================================")


# ============================
# 6. Plot T(E) and gap region
# ============================

plt.figure()
plt.plot(E, T, label="T(E)")
plt.axhline(T_th, ls="--", alpha=0.5, label="T_th")

if E_low is not None and E_high is not None:
    plt.axvspan(E_low, E_high, alpha=0.2, label="gap region")

plt.xlabel("Energy E (eV)")
plt.ylabel("Transmission T(E)")
plt.title("Toy dimerized chain: T(E) and band gap region")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
