# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 18:11:46 2025

@author: USER
"""

"""
make_cnt_pipeline_figure.py
產生最終作品集 / README 會用到的三圖合一 Pipeline 圖

需求：
 - cnt_zonefold_bandstructure.py 可呼叫 zonefold_subbands
 - cnt_gap_to_chain.py 已能產生 T(E) 數據（或你可以 inline 計算）
 - matplotlib

輸出：
 - CNT_bandstructure.png
 - CNT_gap_mapping.png
 - CNT_T_E_gap.png
 - CNT_pipeline_composite.png
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(ROOT, "src"))

from negf.cnt_zonefold import zonefold_subbands, cnt_min_gap, cnt_diameter

#############################################
# (A) 1. CNT zone-folding bandstructure 圖
#############################################

def plot_cnt_bandstructure(n, m, save_path="CNT_bandstructure.png"):
    k, E_p, E_m, Eg = zonefold_subbands(n, m, Nk=300, m_range=5, k_max=0.30)

    plt.figure(figsize=(4, 6))
    for j in range(E_p.shape[1]):
        plt.plot(k, E_p[:, j], 'b-', lw=0.6)
        plt.plot(k, E_m[:, j], 'b-', lw=0.6)

    plt.axhline(Eg/2, color='r', ls='--')
    plt.axhline(-Eg/2, color='r', ls='--')
    plt.title(f"CNT ({n},{m}) zone-folding\nEg = {Eg:.3f} eV")
    plt.xlabel(r"$k_{\parallel}$ (1/nm)")
    plt.ylabel("Energy (eV)")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    return Eg


#############################################
# (B) 2. CNT → Δ Mapping (簡易示意圖)
#############################################

def plot_gap_mapping(Eg, save_path="CNT_gap_mapping.png"):
    Δ = Eg / 2

    plt.figure(figsize=(4, 3))
    plt.plot([0, 1], [Δ, Δ], 'r-', lw=3, label=f"+Δ = +{Δ:.3f} eV")
    plt.plot([0, 1], [-Δ, -Δ], 'b-', lw=3, label=f"-Δ = -{Δ:.3f} eV")
    plt.title(f"Gap Mapping: Eg = 2Δ = {Eg:.3f} eV")
    plt.legend()
    plt.xticks([])
    plt.ylabel("Energy (eV)")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()


#############################################
# (C) 3. NEGF T(E) with Gap Highlight
# 你可以換成你自己跑的 T(E)
#############################################

def fake_TE(E, Eg):
    # 一個簡易模型（之後可換成你的真實 T(E)）
    Δ = Eg/2
    T = np.where(np.abs(E) < Δ, 0.0, 1.0 / (1 + (E/2)**2))
    return T

def plot_T_E_gap(E, T, Eg, save_path="CNT_T_E_gap.png"):
    Δ = Eg / 2

    plt.figure(figsize=(5, 3))
    plt.plot(E, T, lw=2, label="T(E)")
    plt.axvspan(-Δ, Δ, color='red', alpha=0.15, label="Gap region")
    plt.title(f"NEGF Transmission (Eg ≈ {Eg:.3f} eV)")
    plt.xlabel("Energy (eV)")
    plt.ylabel("T(E)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()


#############################################
# (D) 4. 三圖合一 Composite Figure
#############################################

def make_composite(output="CNT_pipeline_composite.png"):
    # Load images
    import matplotlib.image as mpimg

    fig, axs = plt.subplots(1, 3, figsize=(12, 4))

    imgs = [
        ("CNT_bandstructure.png", "CNT bandstructure"),
        ("CNT_gap_mapping.png", "Gap Mapping"),
        ("CNT_T_E_gap.png", "NEGF T(E)"),
    ]

    for ax, (path, title) in zip(axs, imgs):
        img = mpimg.imread(path)
        ax.imshow(img)
        ax.set_title(title)
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


#############################################
# Main Execution
#############################################

if __name__ == "__main__":
    n, m = 17, 0

    Eg = plot_cnt_bandstructure(n, m)
    plot_gap_mapping(Eg)

    E = np.linspace(-2, 2, 600)
    T = fake_TE(E, Eg)
    plot_T_E_gap(E, T, Eg)

    make_composite()
    print("Pipeline figures generated:")
    print("- CNT_bandstructure.png")
    print("- CNT_gap_mapping.png")
    print("- CNT_T_E_gap.png")
    print("- CNT_pipeline_composite.png")
