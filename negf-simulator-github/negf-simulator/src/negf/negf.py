# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 04:40:20 2025

@author: USER
"""

# negf.py
import numpy as np

def self_energy(g_surf, V_couple):
    """
    Σ = V^\dagger * g_surf * V
    g_surf: complex scalar (surface GF of the lead)
    V_couple: (Nc x 1) coupling between device and lead surface site

    return Σ: (Nc x Nc)
    """
    # 保證 g 是 scalar
    if not np.isscalar(g_surf):
        # 如果之後你把 lead 做成多軌域，g_surf 會是矩陣
        # 我們那時候再升級這段
        raise ValueError("For now, g_surf should be a scalar for 1D single-orbital lead.")

    # V^\dagger g V
    # 先做 g * V 會變成 (Nc x 1)
    gV = g_surf * V_couple        # (Nc x 1)
    # 再做 V^\dagger * (gV) → (1 x Nc) @ (Nc x 1) = scalar
    # 但我們不想要 scalar，我們想要整個 Σ 矩陣 (Nc x Nc)
    # 正確形式其實是：Σ = V * g * V^\dagger  (注意轉置的位置!)
    # 這裡要小心，我們來寫成外積：
    Sigma = V_couple @ (g_surf * V_couple.conj().T)
    return Sigma


def retarded_green(E, Hc, SigmaL, SigmaR, eta=1e-6):
    """
    G^r(E) = [ (E+iη)I - Hc - ΣL - ΣR ]^{-1}
    """
    dim = Hc.shape[0]
    I = np.eye(dim, dtype=complex)
    M = (E + 1j*eta)*I - Hc - SigmaL - SigmaR
    G = np.linalg.inv(M)
    return G

def gamma_matrix(Sigma):
    """
    Γ = i (Σ - Σ^\dagger)
    """
    return 1j * (Sigma - Sigma.conj().T)

def transmission(E, Hc, VL, VR, gL, gR, eta=1e-6):
    """
    計算單一能量點的傳輸係數 T(E)
    Hc: 散射區哈密頓量 (Nc x Nc)
    VL: 左 lead 和散射區的耦合 (Nc x 1) 在1D玩具情況
    VR: 右 lead 和散射區的耦合 (Nc x 1)
    gL, gR: 左右 lead 的 surface GF (scalar 或 1x1)
    """
    SigmaL = self_energy(gL, VL)
    SigmaR = self_energy(gR, VR)

    G = retarded_green(E, Hc, SigmaL, SigmaR, eta=eta)

    GammaL = gamma_matrix(SigmaL)
    GammaR = gamma_matrix(SigmaR)

    # T = Tr[GammaL * G * GammaR * G^\dagger]
    T = np.trace(GammaL @ G @ GammaR @ G.conj().T).real
    return T
