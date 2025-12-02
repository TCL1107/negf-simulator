"""
Created on Tue Oct 28 06:52:45 2025

@author: USER
"""

import numpy as np
import matplotlib.pyplot as plt
from surface_gf import surface_gf_1d
from negf import transmission

# -----------------------
# graphene-like 2-atom chain model
# -----------------------

t = -3.0
eps0 = 0.0
Ncell = 5  # 散射區的原胞數

# 單一原胞 (2 atoms)
H0 = np.array([[eps0, t],
               [t, eps0]], dtype=complex)

# 原胞間耦合
V = np.array([[0, 0],
              [t, 0]], dtype=complex)

# 組成散射區 Hc (2*Ncell x 2*Ncell)
dim = 2 * Ncell
Hc = np.zeros((dim, dim), dtype=complex)

for i in range(Ncell):
    # cell 內部
    Hc[2*i:2*i+2, 2*i:2*i+2] = H0
    # cell 間
    if i < Ncell - 1:
        Hc[2*i:2*i+2, 2*(i+1):2*(i+1)+2] = V
        Hc[2*(i+1):2*(i+1)+2, 2*i:2*i+2] = V.conj().T

# 耦合矩陣 VL, VR （只連最左、最右一個cell）
VL = np.zeros((dim, 2), dtype=complex)
VR = np.zeros((dim, 2), dtype=complex)
VL[0:2, :] = V
VR[-2:, :] = V.conj().T

# 掃能量
energies = np.linspace(-10, 10, 400)
Tvals = []

for E in energies:
    gL = surface_gf_1d(E, eps0, t)
    gR = surface_gf_1d(E, eps0, t)
    T_E = transmission(E, Hc, VL, VR, gL, gR)
    Tvals.append(T_E)

Tvals = np.array(Tvals)

plt.figure(figsize=(7,4))
plt.plot(energies, Tvals, lw=1.5)
plt.xlabel("Energy (eV)")
plt.ylabel("Transmission T(E)")
plt.title("Graphene-like 2-atom chain (A-B alternating)")
plt.grid(True, ls="--", alpha=0.6)
plt.tight_layout()

# 存成圖片
plt.savefig("graphene_2atom_TE.png", dpi=300)
plt.savefig("graphene_2atom_TE.pdf")
plt.close()

# 存成 CSV
np.savetxt(
    "graphene_2atom_TE.csv",
    np.column_stack([energies, Tvals]),
    delimiter=",",
    header="E(eV),T(E)",
    comments=""
)

