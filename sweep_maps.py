# -*- coding: utf-8 -*-
"""
Created on Sun Nov 30 20:03:44 2025

@author: USER
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ---------------- Core: Fermi, SGF, NEGF kernels ----------------
def fermi(E, mu, kT=0.025):
    x = np.clip((E - mu)/kT, -60, 60)
    return 1.0/(1.0 + np.exp(x))

def surface_gf_1d(E, eps0, t, eta=1e-3):
    z = E + 1j*eta
    delta = z - eps0
    band_edge = 2.0*abs(t)
    if np.abs(delta.real) <= band_edge:
        root = np.sqrt((2*t)**2 - delta**2)
        g = (delta - 1j*root) / (2*(t**2))
    else:
        root = np.sqrt(delta**2 - (2*t)**2)
        sgn = np.sign(delta.real) if delta.real!=0 else 1.0
        g = (delta - sgn*root) / (2*(t**2))
    return g

def self_energy_scalar(gsurf, t_couple): 
    return t_couple*gsurf*np.conjugate(t_couple)

def gamma_from_sigma(Sigma): 
    return 1j*(Sigma - np.conjugate(Sigma))

def build_hetero_device(epsL, tL, epsR, tR, NL, NR,
                        barrier_shift=0.0,
                        eta_deph=0.0,
                        deph_region="none"):
    N = NL + NR
    Hc = np.zeros((N, N), dtype=complex)

    for i in range(NL): 
        Hc[i,i] = epsL
    for j in range(NR): 
        Hc[NL+j, NL+j] = epsR + barrier_shift

    for i in range(NL-1):
        Hc[i,i+1]=tL; Hc[i+1,i]=tL
    for j in range(NR-1):
        a=NL+j; b=NL+j+1
        Hc[a,b]=tR; Hc[b,a]=tR

    Hc[NL-1, NL] = 0.5*(tL+tR)
    Hc[NL, NL-1] = 0.5*(tL+tR)

    if eta_deph>0.0:
        if deph_region=="right":
            for j in range(NR): Hc[NL+j, NL+j] += -1j*eta_deph
        elif deph_region=="all":
            for i in range(N): Hc[i,i] += -1j*eta_deph

    return Hc, tL, tR, {"N":NL+NR,"NL":NL,"NR":NR}

def retarded_green(E, Hc, SigmaL, SigmaR, eta=1e-9):
    N = Hc.shape[0]
    I = np.eye(N, dtype=complex)
    M = (E + 1j*eta)*I - Hc - SigmaL - SigmaR
    return np.linalg.inv(M)

def transmission_at_E(E, Hc, epsL, tL, epsR, tR, tLc, tRc):
    N = Hc.shape[0]
    gL = surface_gf_1d(E, epsL, tL)
    gR = surface_gf_1d(E, epsR, tR)

    SigmaL = np.zeros((N,N), dtype=complex)
    SigmaR = np.zeros((N,N), dtype=complex)
    SigmaL[0,0]   = self_energy_scalar(gL, tLc)
    SigmaR[-1,-1] = self_energy_scalar(gR, tRc)

    G = retarded_green(E, Hc, SigmaL, SigmaR)

    GammaL = np.zeros((N,N), dtype=complex)
    GammaR = np.zeros((N,N), dtype=complex)
    GammaL[0,0]   = gamma_from_sigma(SigmaL[0,0])
    GammaR[-1,-1] = gamma_from_sigma(SigmaR[-1,-1])

    T = np.trace(GammaL @ G @ GammaR @ G.conjugate().T).real
    return max(T, 0.0)

def compute_TE(Hc, epsL,tL, epsR,tR, tLc,tRc, E_grid):
    return np.array([transmission_at_E(E,Hc,epsL,tL,epsR,tR,tLc,tRc) for E in E_grid])

def current_from_TE(E, T_E, muL, muR, kT=0.025):
    fL, fR = fermi(E, muL, kT), fermi(E, muR, kT)
    return float(np.trapz(T_E*(fL - fR), E))

def TE_window_stats(E, T_E, muL, muR):
    a, b = sorted([muL, muR]); m = (E>=a) & (E<=b)
    if not np.any(m): return 0.0, 0.0
    Tmax = float(np.max(T_E[m]))
    Tavg = float(np.trapz(T_E[m], E[m])/(b-a)) if b>a else 0.0
    return Tmax, Tavg

def rectification_ratio_area(V, I, Vmax=1.5, floor=1e-6):
    m_f = (V>=0) & (V<=Vmax)
    m_r = (V<=0) & (V>=-Vmax)
    If = np.trapz(np.maximum(I[m_f],0), V[m_f])
    Ir = np.trapz(np.abs(np.minimum(I[m_r],0)), V[m_r])
    denom = max(Ir, floor)          # 量測下限 floor：若反向小於此值，RR 用下限計
    RR = If/denom
    saturated = (Ir < floor)        # True → 代表「>RR」的下限報告
    return RR, saturated

def turn_on_voltage_abs(V, I, Ith=0.05):
    Vp, Ip = V[V>=0], np.abs(I[V>=0])
    if np.all(Ip < Ith): return np.nan
    k = np.argmax(Ip >= Ith)
    if k == 0: return float(Vp[0])
    v0, v1 = Vp[k-1], Vp[k]; i0, i1 = Ip[k-1], Ip[k]
    return float(v0 + (Ith - i0) * (v1 - v0) / max(i1 - i0, 1e-12))

def sweep_metrics(NL, NR, epsL, tL, epsR, tR, alpha, muL0, muR0,
                  E_grid, V_list, kT=0.025, eta_deph=0.0, deph_region="right"):
    I_list, Tavg_list = [], []
    for V in V_list:
        muL = muL0 + 0.5*V
        muR = muR0 - 0.5*V
        barrier_shift = alpha * V

        Hc, tLc, tRc, _ = build_hetero_device(
            epsL, tL, epsR, tR, NL, NR,
            barrier_shift=barrier_shift,
            eta_deph=eta_deph,
            deph_region=deph_region
        )
        T_E = compute_TE(Hc, epsL,tL, epsR,tR, tLc,tRc, E_grid)
        I = current_from_TE(E_grid, T_E, muL, muR, kT=kT)
        _, Tavg = TE_window_stats(E_grid, T_E, muL, muR)
        I_list.append(I); Tavg_list.append(Tavg)

    I = np.array(I_list); Tavg = np.array(Tavg_list)
    RR_area, sat = rectification_ratio_area(V_list, I, Vmax=1.5, floor=1e-6)
    Von = turn_on_voltage_abs(V_list, I, Ith=0.05)
    return RR_area, sat, Von, float(np.mean(Tavg))

# ---------------- Global settings ----------------
NL_default, NR_default = 12, 24
epsL_default, tL_default = 0.0, -3.0
alpha_default = -0.4
muL0_default, muR0_default = 0.9, -1.1
kT_default = 0.025
E_grid = np.linspace(-3, 3, 300)      # 調細可更精準（計算較慢）
V_list = np.linspace(-1.5, 1.5, 51)   # 同上

# ---------------- Grid 1: (epsR, tR) ----------------
epsR_vals = np.linspace(2.0, 2.8, 6)
tR_vals   = np.linspace(-1.4, -0.6, 6)

RR_mat_1  = np.zeros((len(epsR_vals), len(tR_vals)))
sat_mat_1 = np.zeros_like(RR_mat_1, dtype=bool)
Von_mat_1 = np.zeros_like(RR_mat_1)
Tav_mat_1 = np.zeros_like(RR_mat_1)

for i, epsR in enumerate(epsR_vals):
    for j, tR in enumerate(tR_vals):
        RR, sat, Von, Tav = sweep_metrics(
            NL_default, NR_default,
            epsL_default, tL_default,
            epsR, tR,
            alpha_default, muL0_default, muR0_default,
            E_grid, V_list, kT=kT_default,
            eta_deph=0.0, deph_region="right"
        )
        RR_mat_1[i,j]  = RR
        sat_mat_1[i,j] = sat
        Von_mat_1[i,j] = Von
        Tav_mat_1[i,j] = Tav

df1 = pd.DataFrame({
    "epsR": np.repeat(epsR_vals, len(tR_vals)),
    "tR":   np.tile(tR_vals, len(epsR_vals)),
    "RR_area": RR_mat_1.flatten(),
    "sat(IR_floor)": sat_mat_1.flatten(),
    "V_on": Von_mat_1.flatten(),
    "<T>":  Tav_mat_1.flatten()
})
df1.to_csv("sweep_epsR_tR.csv", index=False)

plt.figure(figsize=(6,5))
im = plt.imshow(RR_mat_1, origin="lower",
                extent=[tR_vals[0], tR_vals[-1], epsR_vals[0], epsR_vals[-1]],
                aspect="auto")
plt.colorbar(im, label="RR_area (±1.5V)")
plt.xlabel("tR (eV)"); plt.ylabel("epsR (eV)")
plt.title("Rectification ratio (area) vs (epsR, tR)")
plt.tight_layout(); plt.savefig("heatmap_rr_epsR_tR.png", dpi=150); plt.show()

plt.figure(figsize=(6,5))
im = plt.imshow(Von_mat_1, origin="lower",
                extent=[tR_vals[0], tR_vals[-1], epsR_vals[0], epsR_vals[-1]],
                aspect="auto")
plt.colorbar(im, label="V_on (arb. V)")
plt.xlabel("tR (eV)"); plt.ylabel("epsR (eV)")
plt.title("Turn-on voltage vs (epsR, tR)")
plt.tight_layout(); plt.savefig("heatmap_von_epsR_tR.png", dpi=150); plt.show()

plt.figure(figsize=(6,5))
im = plt.imshow(Tav_mat_1, origin="lower",
                extent=[tR_vals[0], tR_vals[-1], epsR_vals[0], epsR_vals[-1]],
                aspect="auto")
plt.colorbar(im, label="<T> in bias window")
plt.xlabel("tR (eV)"); plt.ylabel("epsR (eV)")
plt.title("Average transmission <T> vs (epsR, tR)")
plt.tight_layout(); plt.savefig("heatmap_tav_epsR_tR.png", dpi=150); plt.show()

# ---------------- Grid 2: (NR, alpha) ----------------
NR_vals    = np.array([10, 14, 18, 22, 26, 30])
alpha_vals = np.linspace(-0.7, -0.2, 6)

RR_mat_2  = np.zeros((len(NR_vals), len(alpha_vals)))
sat_mat_2 = np.zeros_like(RR_mat_2, dtype=bool)
Von_mat_2 = np.zeros_like(RR_mat_2)
Tav_mat_2 = np.zeros_like(RR_mat_2)

epsR_ref, tR_ref = 2.6, -0.8

for i, NR in enumerate(NR_vals):
    for j, alpha in enumerate(alpha_vals):
        RR, sat, Von, Tav = sweep_metrics(
            NL_default, int(NR),
            epsL_default, tL_default,
            epsR_ref, tR_ref,
            alpha, muL0_default, muR0_default,
            E_grid, V_list, kT=kT_default,
            eta_deph=0.0, deph_region="right"
        )
        RR_mat_2[i,j]  = RR
        sat_mat_2[i,j] = sat
        Von_mat_2[i,j] = Von
        Tav_mat_2[i,j] = Tav

df2 = pd.DataFrame({
    "NR": np.repeat(NR_vals, len(alpha_vals)),
    "alpha": np.tile(alpha_vals, len(NR_vals)),
    "RR_area": RR_mat_2.flatten(),
    "sat(IR_floor)": sat_mat_2.flatten(),
    "V_on": Von_mat_2.flatten(),
    "<T>":  Tav_mat_2.flatten()
})
df2.to_csv("sweep_NR_alpha.csv", index=False)

plt.figure(figsize=(6,5))
im = plt.imshow(RR_mat_2, origin="lower",
                extent=[alpha_vals[0], alpha_vals[-1], NR_vals[0], NR_vals[-1]],
                aspect="auto")
plt.colorbar(im, label="RR_area (±1.5V)")
plt.xlabel("alpha"); plt.ylabel("NR (sites)")
plt.title("Rectification ratio (area) vs (NR, alpha)")
plt.tight_layout(); plt.savefig("heatmap_rr_NR_alpha.png", dpi=150); plt.show()

plt.figure(figsize=(6,5))
im = plt.imshow(Von_mat_2, origin="lower",
                extent=[alpha_vals[0], alpha_vals[-1], NR_vals[0], NR_vals[-1]],
                aspect="auto")
plt.colorbar(im, label="V_on (arb. V)")
plt.xlabel("alpha"); plt.ylabel("NR (sites)")
plt.title("Turn-on voltage vs (NR, alpha)")
plt.tight_layout(); plt.savefig("heatmap_von_NR_alpha.png", dpi=150); plt.show()

plt.figure(figsize=(6,5))
im = plt.imshow(Tav_mat_2, origin="lower",
                extent=[alpha_vals[0], alpha_vals[-1], NR_vals[0], NR_vals[-1]],
                aspect="auto")
plt.colorbar(im, label="<T> in bias window")
plt.xlabel("alpha"); plt.ylabel("NR (sites)")
plt.title("Average transmission <T> vs (NR, alpha)")
plt.tight_layout(); plt.savefig("heatmap_tav_NR_alpha.png", dpi=150); plt.show()

print("Done. Files written:",
      "heatmap_rr_epsR_tR.png, heatmap_von_epsR_tR.png, heatmap_tav_epsR_tR.png, sweep_epsR_tR.csv,",
      "heatmap_rr_NR_alpha.png, heatmap_von_NR_alpha.png, heatmap_tav_NR_alpha.png, sweep_NR_alpha.csv")
