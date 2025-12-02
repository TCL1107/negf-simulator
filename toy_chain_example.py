# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 04:53:56 2025

@author: USER
"""
# toy_chain_example.py
import numpy as np
import matplotlib.pyplot as plt

from surface_gf import surface_gf_1d
from negf import transmission

# 物理參數
t = -3.0      # eV, 最近鄰 hopping
eps0 = 0.0    # eV, on-site
N = 5         # 散射區 sites 數量

# 建立散射區哈密頓量 Hc (N x N)
Hc = np.zeros((N, N), dtype=complex)
for i in range(N):
    Hc[i, i] = eps0
    if i < N-1:
        Hc[i, i+1] = t
        Hc[i+1, i] = t
        Hc[2, 2] = 1.0  # eV，隨便放一個位能抬高


# 左右耦合矩陣 VL, VR
# 在這個最簡版本裡，我們假裝每個 lead 只有 "1-site surface block"
VL = np.zeros((N, 1), dtype=complex)
VR = np.zeros((N, 1), dtype=complex)
VL[0,0] = t   # 左 lead 連到散射區的第0個site
VR[-1,0] = t  # 右 lead 連到散射區的最後一個site

# 掃能量
energies = np.linspace(-10, 10, 400)
Tvals = []

for E in energies:
    # 算左右 lead 的 surface GF
    gL = surface_gf_1d(E, eps0, t)
    gR = surface_gf_1d(E, eps0, t)
    # 算傳輸
    T_E = transmission(E, Hc, VL, VR, gL, gR)
    Tvals.append(T_E)

# 畫圖
plt.figure()
plt.plot(energies, Tvals)
plt.xlabel("Energy (eV)")
plt.ylabel("Transmission T(E)")
plt.title("1D chain transmission")
plt.grid(True)
plt.show()

