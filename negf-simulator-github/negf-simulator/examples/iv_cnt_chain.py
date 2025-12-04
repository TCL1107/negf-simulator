"""
iv_cnt_chain_example.py

Simple Landauer I–V example for a 1D dimerized chain calibrated
to a CNT bandgap (Eg ≈ 0.62 eV → Δ ≈ Eg/2).

This script is self-contained and does NOT depend on src/negf.
You can drop it into examples/ and run it with your Python env.

It uses:
 - wide-band limit self-energies (Σ_L/R = -i Γ/2 on first/last site)
 - finite chain Hamiltonian with ±Δ on-site energies
 - Landauer formula with Fermi functions at 300 K
"""

import numpy as np
import matplotlib.pyplot as plt

e_charge = 1.0      # work in eV units (e = 1)
h_bar = 1.0         # set ħ = 1
h_planck = 2 * np.pi * h_bar
kBT = 0.02585       # ~ 300 K in eV

def fermi(E, mu):
    return 1.0 / (1.0 + np.exp((E - mu) / kBT))

def make_chain_hamiltonian(N, Delta, t, V_bias=0.0):
    """
    Build a 1D dimerized chain:
      on-site: +Δ, -Δ alternating,
      hopping: t between neighbors,
      optional linear potential drop across the chain for bias.
    """
    H = np.zeros((N, N), dtype=complex)

    # simple linear potential drop from +V/2 (left) to -V/2 (right)
    if N > 1:
        x = np.linspace(+V_bias/2, -V_bias/2, N)
    else:
        x = np.array([0.0])

    for i in range(N):
        onsite_sign = +1 if (i % 2 == 0) else -1
        H[i, i] = onsite_sign * Delta + x[i]
        if i < N - 1:
            H[i, i+1] = t
            H[i+1, i] = t

    return H

def self_energies_wbl(N, Gamma):
    """
    Wide-band limit self-energies:
      Σ_L on site 0, Σ_R on site N-1 = -i Γ/2
    """
    Sigma_L = np.zeros((N, N), dtype=complex)
    Sigma_R = np.zeros((N, N), dtype=complex)

    Sigma_L[0, 0] = -1j * Gamma / 2.0
    Sigma_R[-1, -1] = -1j * Gamma / 2.0

    return Sigma_L, Sigma_R

def transmission(E, H, Sigma_L, Sigma_R, eta=1e-6):
    """
    Compute transmission T(E) = Tr[Γ_L G Γ_R G^†]
    """
    I = np.eye(H.shape[0], dtype=complex)
    z = (E + 1j * eta) * I
    G = np.linalg.inv(z - H - Sigma_L - Sigma_R)

    Gamma_L = 1j * (Sigma_L - Sigma_L.conjugate().T)
    Gamma_R = 1j * (Sigma_R - Sigma_R.conjugate().T)

    GR = G
    GA = G.conjugate().T

    return np.real(np.trace(Gamma_L @ GR @ Gamma_R @ GA))

def current_vs_voltage(N=40, Delta=0.31, t=-2.7, Gamma=0.2,
                       E_min=-2.0, E_max=2.0, nE=2000,
                       V_points=np.linspace(0.0, 0.8, 17)):
    """
    Compute I(V) using the Landauer formula:
      I(V) = (2e/h) ∫ dE T(E,V) [f(E, μ_L) - f(E, μ_R)]

    Here we work in eV units with e = 1, so 2e/h becomes a constant factor.
    """
    E_grid = np.linspace(E_min, E_max, nE)
    dE = E_grid[1] - E_grid[0]

    currents = []

    # Prefactor in natural units (e=1). For plotting shape, constant scale is fine.
    prefactor = 2.0 / h_planck

    for V in V_points:
        # Build Hamiltonian at this bias
        H = make_chain_hamiltonian(N, Delta, t, V_bias=V)
        Sigma_L, Sigma_R = self_energies_wbl(N, Gamma)

        T_vals = np.zeros_like(E_grid)
        for i, E in enumerate(E_grid):
            T_vals[i] = transmission(E, H, Sigma_L, Sigma_R)

        mu_L = +V / 2.0
        mu_R = -V / 2.0
        fL = fermi(E_grid, mu_L)
        fR = fermi(E_grid, mu_R)

        integrand = T_vals * (fL - fR)
        I = prefactor * np.trapz(integrand, E_grid)
        currents.append(I)

    return V_points, np.array(currents)

if __name__ == "__main__":
    # Parameters roughly calibrated to a CNT with Eg ~ 0.62 eV
    N = 40
    Delta = 0.31   # Eg ≈ 2Δ → Eg ~ 0.62 eV
    t = -2.7
    Gamma = 0.3

    V_array, I_array = current_vs_voltage(
        N=N,
        Delta=Delta,
        t=t,
        Gamma=Gamma,
        V_points=np.linspace(0.0, 0.8, 17)
    )

    plt.figure(figsize=(4, 3))
    plt.plot(V_array, I_array, marker="o")
    plt.xlabel("Bias V (V)")
    plt.ylabel("Current I (arb. units)")
    plt.title("I–V of CNT-calibrated 1D chain (Landauer, WBL)")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("IV_cnt_chain_example.png", dpi=300)
    plt.show()
