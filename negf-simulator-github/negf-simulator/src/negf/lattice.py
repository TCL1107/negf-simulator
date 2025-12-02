# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 07:16:24 2025

@author: USER
"""

import numpy as np
import matplotlib.pyplot as plt

# 物理參數
a_cc = 1.42    # C-C distance (Å)
t1 = -3.0      # 最近鄰 hopping (eV)
Ny = 4         # ribbon 寬度 (幾層原子)
Nx = 6         # ribbon 長度 (原胞數)

# 生成原子坐標
atoms = []
sublattice = []  # A 或 B
for ix in range(Nx):
    for iy in range(Ny):
        # A 原子
        xA = 3/2 * a_cc * ix
        yA = np.sqrt(3) * a_cc * iy
        atoms.append([xA, yA])
        sublattice.append('A')
        # B 原子
        xB = xA + a_cc / 2
        yB = yA + np.sqrt(3) * a_cc / 2
        atoms.append([xB, yB])
        sublattice.append('B')

atoms = np.array(atoms)
Natom = len(atoms)

# 建 H 矩陣
H = np.zeros((Natom, Natom), dtype=complex)
for i in range(Natom):
    for j in range(i+1, Natom):
        dist = np.linalg.norm(atoms[i] - atoms[j])
        if np.isclose(dist, a_cc, atol=0.1):
            H[i,j] = t1
            H[j,i] = t1

# 畫 lattice
plt.figure(figsize=(6,4))
for i, (x,y) in enumerate(atoms):
    color = 'royalblue' if sublattice[i]=='A' else 'tomato'
    plt.scatter(x, y, color=color, s=60)
    plt.text(x+0.05, y+0.05, str(i), fontsize=8)
plt.axis('equal')
plt.title(f"Zigzag Graphene Nanoribbon (Ny={Ny}, Nx={Nx})")
plt.show()

# 看看 H 矩陣的大小
print(f"H matrix size: {H.shape}")
np.save("ZGNR_H.npy", H)
np.save("ZGNR_atoms.npy", atoms)
print("Graphene Hamiltonian saved as ZGNR_H.npy")

