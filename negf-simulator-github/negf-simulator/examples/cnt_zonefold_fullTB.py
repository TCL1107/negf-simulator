"""
examples/cnt_zonefold_fullTB.py

Full-band style CNT (n,m) bandstructure using graphene π-band TB
and a simple zone-folding construction.

This is qualitative: it reproduces the 1/d gap scaling and
metallic vs semiconducting behavior, and gives Lambin-style
multi-band panels over ~[-3, 3] eV.

It relies only on numpy + matplotlib and the diameter/gap formulas
already in cnt_zonefold.py.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# --- allow importing src/negf ---
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(ROOT, "src"))

from negf.cnt_zonefold import cnt_diameter, chiral_index_nu, gamma0, a_cc


# ===== graphene π-band dispersion (full TB) =====

def graphene_E(kx, ky, a):
    """
    Full nearest-neighbor π-band dispersion of graphene:

    E(k) = ± γ0 * sqrt(1 + 4 cos(√3 kx a / 2) cos(ky a / 2)
                        + 4 cos^2(ky a / 2))

    Here kx, ky in 1/nm, a is lattice constant ≈ 0.246 nm.
    """
    term1 = 1.0
    term2 = 4.0 * np.cos(np.sqrt(3) * kx * a / 2.0) * np.cos(ky * a / 2.0)
    term3 = 4.0 * np.cos(ky * a / 2.0) ** 2
    inside = term1 + term2 + term3
    inside = np.clip(inside, 0.0, None)
    Ek = gamma0 * np.sqrt(inside)
    return Ek, -Ek  # conduction, valence


# ===== simple zone-fold mapping =====

def cnt_zonefold_fullTB(n, m, Nk=301, m_range=15, k_par_max=2.5):
    """
    Qualitative full-TB zone-folding:

    - We treat k_parallel along y-axis, k_perp along x-axis in graphene k-space.
    - Quantization: k_perp(q) = 2π (q + ν/3) / |C_h|
      where |C_h| = π d, d = diameter, ν = chiral index.

    Parameters
    ----------
    n, m : chirality
    Nk : number of k_parallel points
    m_range : draw subbands for q = -m_range..+m_range
    k_par_max : max |k_parallel| (1/nm) to plot

    Returns
    -------
    k_par, E_plus_list, E_minus_list
    """
    d = cnt_diameter(n, m)
    C_len = np.pi * d
    nu = chiral_index_nu(n, m)
    if nu == 2:
        nu_eff = -1
    else:
        nu_eff = nu

    a = np.sqrt(3) * a_cc  # graphene lattice constant

    k_par = np.linspace(-k_par_max, k_par_max, Nk)  # along tube axis (treated as ky)
    q_vals = np.arange(-m_range, m_range + 1)
    E_plus_list = []
    E_minus_list = []

    for q in q_vals:
        k_perp = (2.0 * np.pi / C_len) * (q + nu_eff / 3.0)  # 1/nm
        kx = np.full_like(k_par, k_perp)
        ky = k_par

        Ep, Em = graphene_E(kx, ky, a)
        E_plus_list.append(Ep)
        E_minus_list.append(Em)

    E_plus = np.stack(E_plus_list, axis=1)
    E_minus = np.stack(E_minus_list, axis=1)

    return k_par, E_plus, E_minus


def plot_fullTB_panel(n, m, Nk=301, m_range=15, k_par_max=2.5,
                      Emax_plot=3.0):
    d = cnt_diameter(n, m)
    k_par, E_plus, E_minus = cnt_zonefold_fullTB(
        n, m, Nk=Nk, m_range=m_range, k_par_max=k_par_max
    )

    plt.figure(figsize=(4, 6))
    for j in range(E_plus.shape[1]):
        plt.plot(k_par, E_plus[:, j], lw=0.5)
        plt.plot(k_par, E_minus[:, j], lw=0.5)

    plt.ylim(-Emax_plot, Emax_plot)
    plt.xlabel(r"$k_{\parallel}$ (1/nm)")
    plt.ylabel("Energy (eV)")
    plt.title(f"CNT ({n},{m}) full-TB zone-folding\n d ≈ {d:.3f} nm")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # 重現 Lambin 圖裡幾顆代表性 CNT：
    panels = [(17, 0), (13, 6), (10, 10)]
    for (n, m) in panels:
        plot_fullTB_panel(n, m,
                          Nk=401,
                          m_range=20,
                          k_par_max=3.0,
                          Emax_plot=3.0)
