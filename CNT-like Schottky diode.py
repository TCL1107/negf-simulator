import numpy as np
import matplotlib.pyplot as plt

# =========================
# 0) 全域參數「旋鈕」
# =========================
# 幾何
NL, NR = 12, 24

# 左段（金屬側）
epsL, tL = 0.0, -3.0

# 右段（障壁側）— 高障壁/窄帶寬
epsR, tR = +2.6, -0.8

# 偏壓耦合：ΔU = alpha * V（alpha<0 → 正偏為順向）
alpha = -0.4

# 能量/偏壓網格與溫度
E_grid = np.linspace(-3, 3, 900)
V_list = np.linspace(-2.0, 2.0, 81)
kT = 0.025  # ~300K

# 接觸功函數差（零偏）
muL0, muR0 = 0.9, -1.1

# 量測口徑
V_probe = 1.5       # RR 評估用的 ±V
Vmax_RR = 1.5       # 面積 RR 積分範圍
I_floor  = 1e-6     # 面積/點狀 RR 的「量測下限」
Ith_abs  = 0.05     # turn-on 的絕對門檻（相對單位）

# 物理/工程開關
deph_region = "right"   # "none" | "right" | "all"
eta_deph    = 0.008     # eV，0 代表無去相干
contact_scale = 1.0     # <1 降接觸透明度（兩端同時乘上）
T_bg = 0.0              # 背景漏電（並聯路徑），建議 1e-6~1e-5 起試

# =========================
# 1) 工具函式（Fermi、KPI）
# =========================
def fermi(E, mu, kT=0.025):
    x = np.clip((E - mu)/kT, -60, 60)
    return 1.0/(1.0 + np.exp(x))

def rectification_ratio_point_safe(V, I, V_probe=1.5, floor=1e-6):
    vpos = V[np.argmin(np.abs(V - (+V_probe)))]
    vneg = V[np.argmin(np.abs(V - (-V_probe)))]
    Ipos = float(I[V == vpos][0]); Ineg = float(I[V == vneg][0])
    denom = max(abs(Ineg), floor)
    RR = abs(Ipos)/denom
    saturated = (abs(Ineg) < floor)
    return RR, saturated

def rectification_ratio_area_safe(V, I, Vmax=1.5, floor=1e-6):
    m_f = (V>=0) & (V<=Vmax)
    m_r = (V<=0) & (V>=-Vmax)
    If = np.trapz(np.maximum(I[m_f],0), V[m_f])
    Ir = np.trapz(np.abs(np.minimum(I[m_r],0)), V[m_r])
    denom = max(Ir, floor)
    RR = If/denom
    saturated = (Ir < floor)
    return RR, saturated

def turn_on_voltage_abs(V, I, Ith=0.05):
    Vp, Ip = V[V>=0], np.abs(I[V>=0])
    if np.all(Ip < Ith): return np.nan
    k = np.argmax(Ip >= Ith)
    if k == 0: return float(Vp[0])
    v0, v1 = Vp[k-1], Vp[k]; i0, i1 = Ip[k-1], Ip[k]
    return float(v0 + (Ith - i0)*(v1 - v0)/max(i1 - i0, 1e-12))

def TE_window_stats(E, T_E, muL, muR):
    a, b = sorted([muL, muR]); m = (E>=a) & (E<=b)
    if not np.any(m): return 0.0, 0.0
    Tmax = float(np.max(T_E[m]))
    Tavg = float(np.trapz(T_E[m], E[m])/(b-a)) if b>a else 0.0
    return Tmax, Tavg

# =========================
# 2) 半無限 1D 導線 surface GF
# =========================
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
    return t_couple * gsurf * np.conjugate(t_couple)

def gamma_from_sigma(Sigma):
    return 1j*(Sigma - np.conjugate(Sigma))

# =========================
# 3) 裝置建構（含局部/全域去相干）
# =========================
def build_hetero_device(epsL, tL, epsR, tR, NL, NR,
                        barrier_shift=0.0,
                        eta_deph=0.0,
                        deph_region="none"):
    N = NL + NR
    Hc = np.zeros((N, N), dtype=complex)

    # 左/右 on-site（右段含偏壓障壁）
    for i in range(NL):
        Hc[i,i] = epsL
    for j in range(NR):
        Hc[NL+j, NL+j] = epsR + barrier_shift

    # 左/右 hopping
    for i in range(NL-1):
        Hc[i,i+1] = tL; Hc[i+1,i] = tL
    for j in range(NR-1):
        a = NL + j; b = NL + j + 1
        Hc[a,b] = tR; Hc[b,a] = tR

    # 界面耦合（平均）
    Hc[NL-1, NL] = 0.5*(tL + tR)
    Hc[NL,   NL-1] = 0.5*(tL + tR)

    # 去相干（把 -i*eta_deph 加到指定區域的對角）
    if eta_deph > 0.0:
        if deph_region == "right":
            for j in range(NR):
                Hc[NL+j, NL+j] += -1j*eta_deph
        elif deph_region == "all":
            for i in range(N):
                Hc[i,i] += -1j*eta_deph

    # 接觸透明度調整
    t_left_contact  = contact_scale * tL
    t_right_contact = contact_scale * tR

    return Hc, t_left_contact, t_right_contact, {"N": N, "NL": NL, "NR": NR}

# =========================
# 4) NEGF 核心
# =========================
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

def compute_TE(Hc, epsL, tL, epsR, tR, tLc, tRc, E_grid):
    T = np.array([transmission_at_E(E, Hc, epsL, tL, epsR, tR, tLc, tRc) for E in E_grid])
    if T_bg > 0.0:
        T = T + T_bg  # 背景漏電
    return T

def current_from_TE(E, T_E, muL, muR, kT=0.025):
    fL, fR = fermi(E, muL, kT), fermi(E, muR, kT)
    return float(np.trapz(T_E*(fL - fR), E))

# =========================
# 5) 掃偏壓 & 報告
# =========================
def sweep_IV(eta_deph=0.0, deph_region="right"):
    I_list, Tavg_list = [], []
    for V in V_list:
        # 化學勢 & 障壁位移
        muL = muL0 + 0.5*V
        muR = muR0 - 0.5*V
        barrier_shift = alpha * V

        # 重建 Hc（含去相干）+ 計 T(E,V)
        Hc, tLc, tRc, _ = build_hetero_device(
            epsL, tL, epsR, tR, NL, NR,
            barrier_shift=barrier_shift,
            eta_deph=eta_deph,
            deph_region=deph_region
        )
        T_E = compute_TE(Hc, epsL, tL, epsR, tR, tLc, tRc, E_grid)

        # I(V)
        Ival = current_from_TE(E_grid, T_E, muL, muR, kT=kT)
        I_list.append(Ival)

        # 窗內 ⟨T⟩（導通能力指標）
        _, Tavg = TE_window_stats(E_grid, T_E, muL, muR)
        Tavg_list.append(Tavg)
    return np.array(I_list), np.array(Tavg_list)

if __name__ == "__main__":
    # Baseline（無去相干）
    I_base, Tavg_base = sweep_IV(eta_deph=0.0, deph_region=deph_region)
    RR_pt_base, sat_pt_b = rectification_ratio_point_safe(V_list, I_base, V_probe=V_probe, floor=I_floor)
    RR_ar_base, sat_ar_b = rectification_ratio_area_safe (V_list, I_base, Vmax=Vmax_RR, floor=I_floor)
    Von_base = turn_on_voltage_abs(V_list, I_base, Ith=Ith_abs)

    # Dephasing 對照
    I_deph, Tavg_deph = sweep_IV(eta_deph=eta_deph, deph_region=deph_region)
    RR_pt_deph, sat_pt_d = rectification_ratio_point_safe(V_list, I_deph, V_probe=V_probe, floor=I_floor)
    RR_ar_deph, sat_ar_d = rectification_ratio_area_safe (V_list, I_deph, Vmax=Vmax_RR, floor=I_floor)
    Von_deph = turn_on_voltage_abs(V_list, I_deph, Ith=Ith_abs)

    # ——輸出（若反向<下限 → 用 “>RR” 呈現）——
    tag_b = ">" if sat_pt_b else "="
    tag_d = ">" if sat_pt_d else "="
    tag_ab = ">" if sat_ar_b else "="
    tag_ad = ">" if sat_ar_d else "="
    print(f"[Baseline] RR_point@±{V_probe}V{tag_b}{RR_pt_base:.2e}, RR_area@±{Vmax_RR}V{tag_ab}{RR_ar_base:.2e} (floor={I_floor:g}), "
          f"V_on≈{Von_base:.3f} V, <T>={np.mean(Tavg_base):.3f}")
    print(f"[Dephase ] RR_point@±{V_probe}V{tag_d}{RR_pt_deph:.2e}, RR_area@±{Vmax_RR}V{tag_ad}{RR_ar_deph:.2e} (floor={I_floor:g}), "
          f"V_on≈{Von_deph:.3f} V, <T>={np.mean(Tavg_deph):.3f}")

    # 繪 I–V
    plt.figure(figsize=(7,4))
    plt.plot(V_list, I_base, label=f"baseline (η=0)")
    if eta_deph > 0:
        plt.plot(V_list, I_deph, label=f"dephasing (η={eta_deph} eV, {deph_region})")
    plt.xlabel("Bias V (eV)"); plt.ylabel("Current (arb.)")
    plt.title("I–V of CNT-like Schottky diode @300K")
    plt.grid(True, ls="--", alpha=0.6); plt.legend(); plt.tight_layout(); plt.show()
