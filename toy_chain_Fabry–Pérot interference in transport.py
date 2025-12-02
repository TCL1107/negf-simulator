# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 05:57:37 2025

@author: USER
"""

import numpy as np
import matplotlib.pyplot as plt

from surface_gf import surface_gf_1d
from negf import transmission

# -----------------------
# 物理參數設定
# -----------------------

# 中央散射區 (device/channel) 的 hopping
t_device = -3.0   # eV  (類似我們前面的值)

# leads 的 hopping 我們故意改掉，做出 mismatch
t_lead   = -2.0   # eV  (不同於 t_device，造成介面反射)

# on-site energy
eps0 = 0.0        # eV

# 散射區長度：N 個 site
N = 5

# -----------------------
# 建立散射區的哈密頓量 Hc
# -----------------------
Hc = np.zeros((N, N), dtype=complex)

for i in range(N):
    Hc[i, i] = eps0
    if i < N-1:
        Hc[i, i+1] = t_device
        Hc[i+1, i] = t_device

# 你也可以在這邊加障礙來玩，例如打開下面這行：
# Hc[2,2] = 1.0  # eV, 在中間 site 加一個位能抬高，模擬雜質/缺陷/閘極勢壘

# -----------------------
# 裝置與 leads 的耦合矩陣 VL, VR
# -----------------------
# 還是最簡case：lead 的表面只有一個 site，耦合到裝置最邊邊的 site
VL = np.zeros((N, 1), dtype=complex)
VR = np.zeros((N, 1), dtype=complex)

# 我們這裡用的是 device 端的 hopping 當界面耦合
VL[0,0]  = t_device
VR[-1,0] = t_device

# -----------------------
# 掃能量，計算 T(E)
# -----------------------
energies = np.linspace(-10, 10, 400)
Tvals = []

for E in energies:
    # 左右 lead 的 surface green function：注意！用 t_lead
    gL = surface_gf_1d(E, eps0, t_lead)
    gR = surface_gf_1d(E, eps0, t_lead)

    T_E = transmission(E, Hc, VL, VR, gL, gR)
    Tvals.append(T_E)

# -----------------------
# 畫圖
# -----------------------
plt.figure()
plt.plot(energies, Tvals)
plt.xlabel("Energy (eV)")
plt.ylabel("Transmission T(E)")
plt.title("1D chain with contact mismatch (Fabry–Pérot-like)")
plt.grid(True)
plt.ylim(-0.05, 1.05)  # 把y軸拉在[0,1]上下，方便看
plt.show()
