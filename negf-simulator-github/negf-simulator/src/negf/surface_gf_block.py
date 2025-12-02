# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 17:26:41 2025

@author: USER
"""

import numpy as np

def surface_gf_block(E, H0, V, eta=1e-6, max_iter=200, tol=1e-10):
    """
    求半無限 periodic lead 的 surface Green's function g(E)
    對象是「多軌域」原胞：
    - H0: 單元胞內哈密頓量 (n x n)
    - V:  該單元胞 與 下一個單元胞 的耦合 (n x n)

    我們使用一個常見的 iterative decimation：
    不斷把更深層的單元胞向表面折疊，最後取得有效 self-energy，
    最後 g = [ (E+iη)I - H0_eff ]^{-1}
    """

    n = H0.shape[0]
    I = np.eye(n, dtype=complex)

    # 初始 "effective" 參數
    alpha = V.copy()        # forward coupling
    beta  = V.conj().T.copy()  # backward coupling
    eps   = H0.copy()       # onsite block that we'll renormalize

    z = (E + 1j*eta)

    for _ in range(max_iter):
        # 解 (zI - eps)
        M = z*I - eps
        M_inv = np.linalg.inv(M)

        # renormalize
        eps_new   = eps + alpha @ M_inv @ beta
        alpha_new = alpha @ M_inv @ alpha
        beta_new  = beta  @ M_inv @ beta

        # 收斂檢查
        if (np.linalg.norm(eps_new - eps) < tol and
            np.linalg.norm(alpha_new - alpha) < tol and
            np.linalg.norm(beta_new  - beta)  < tol):
            eps = eps_new
            break

        eps   = eps_new
        alpha = alpha_new
        beta  = beta_new

    # surface GF
    gsurf = np.linalg.inv(z*I - eps)
    return gsurf
