# cnt_diode_dephase.py
# CNT-like heterojunction diode with bias-dependent barrier and dephasing (NEGF, 1D TB)

import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# 1) Analytic surface Green's function for semi-infinite 1D chain
# -----------------------------
def surface_gf_1d(E, eps0, t, eta=1e-3):
    """
    Retarded surface Green's function for semi-infinite 1D chain
    onsite = eps0, NN hopping = t (typically negative).
    eta: small positive imaginary part for numerical stability.
    """
    z = E + 1j*eta
    delta = z - eps0
    band_edge = 2.0 * abs(t)

    if np.abs(delta.real) <= band_edge:
        # inside band
        root = np.sqrt((2*t)**2 - delta**2)      # sqrt(4 t^2 - (E-eps0)^2)
        g = (delta - 1j*root) / (2 * (t**2))
    else:
        # outside band
        root = np.sqrt(delta**2 - (2*t)**2)
        sgn = np.sign(delta.real) if delta.real != 0 else 1.0
        g = (delta - sgn*root) / (2 * (t**2))
    return g

def self_energy_scalar(gsurf, t_couple):
    # Σ = τ g_s τ*
    return t_couple * gsurf * np.conjugate(t_couple)

def gamma_from_sigma(Sigma):
    return 1j * (Sigma - np.conjugate(Sigma))

# -----------------------------
# 2) Device H(V): 左右異質 + 偏壓讓右段位能上/下移；可加去相干 -i*eta_dephase
# -----------------------------
def build_hetero_device_biased(epsL, tL,
                               epsR0, tR,
                               NL, NR,
                               barrier_shift_extra,
                               eta_dephase=0.0,
                               apply_region="all"):
    """
    產生 (NL+NR)x(NL+NR) 的裝置哈密頓量 Hc：
      左段 on-site = epsL, hopping = tL
      右段 on-site = epsR0 + barrier_shift_extra, hopping = tR
      交界面 hopping = (tL+tR)/2
    並可在指定區域加入 -i*eta_dephase 模擬去相干/散射。

    apply_region ∈ {"all","left","right","center"}：
      - "all": 全裝置加去相干
      - "left": 只加左段
      - "right": 只加右段（較貼近障壁區散射）
      - "center": 只在交界 1~2 格加（示意）
    """
    N = NL + NR
    Hc = np.zeros((N, N), dtype=complex)

    # 左段 on-site
    for i in range(NL):
        Hc[i, i] = epsL

    # 右段 on-site (含偏壓位移)
    for j in range(NR):
        Hc[NL + j, NL + j] = epsR0 + barrier_shift_extra

    # 左段 hopping
    for i in range(NL - 1):
        Hc[i, i+1] = tL
        Hc[i+1, i] = tL

    # 右段 hopping
    for j in range(NR - 1):
        a = NL + j
        b = NL + j + 1
        Hc[a, b] = tR
        Hc[b, a] = tR

    # 界面耦合
    Hc[NL-1, NL] = 0.5*(tL + tR)
    Hc[NL, NL-1] = 0.5*(tL + tR)

    # 加入 -i*eta_dephase
    if eta_dephase > 0.0:
        diag = np.zeros(N, dtype=complex)
        if apply_region == "all":
            diag[:] = -1j*eta_dephase
        elif apply_region == "left":
            diag[:NL] = -1j*eta_dephase
        elif apply_region == "right":
            diag[NL:] = -1j*eta_dephase
        elif apply_region == "center":
            diag[NL-1] = -1j*eta_dephase
            diag[NL]   = -1j*eta_dephase
        Hc += np.diag(diag)

    t_left_contact  = tL
    t_right_contact = tR
    return Hc, t_left_contact, t_right_contact

# -----------------------------
# 3) NEGF core
# -----------------------------
def retarded_green(E, Hc, SigmaL, SigmaR, eta=1e-4):
    N = Hc.shape[0]
    I = np.eye(N, dtype=complex)
    M = (E + 1j*eta)*I - Hc - SigmaL - SigmaR
    return np.linalg.inv(M)

def transmission_at_E(E, Hc,
                      epsL, tL,
                      epsR_eff, tR,
                      tLc, tRc):
    """
    單能量下的傳輸係數 T(E)
    epsR_eff = epsR0 + barrier_shift_extra (與偏壓一致)
    """
    N = Hc.shape[0]
    gL = surface_gf_1d(E, epsL, tL)
    gR = surface_gf_1d(E, epsR_eff, tR)

    SigmaL = np.zeros((N, N), dtype=complex)
    SigmaR = np.zeros((N, N), dtype=complex)
    SigmaL_edge = self_energy_scalar(gL, tLc)
    SigmaR_edge = self_energy_scalar(gR, tRc)
    SigmaL[0,0]   = SigmaL_edge
    SigmaR[-1,-1] = SigmaR_edge

    G = retarded_green(E, Hc, SigmaL, SigmaR, eta=1e-4)

    GammaL = np.zeros((N, N), dtype=complex)
    GammaR = np.zeros((N, N), dtype=complex)
    GammaL[0,0]   = gamma_from_sigma(SigmaL_edge)
    GammaR[-1,-1] = gamma_from_sigma(SigmaR_edge)

    T = np.trace(GammaL @ G @ GammaR @ G.conjugate().T).real
    return max(T, 0.0)

def compute_TE_for_bias(E_grid,
                        epsL, tL,
                        epsR0, tR,
                        NL, NR,
                        barrier_shift_extra,
                        eta_dephase=0.0,
                        apply_region="right"):
    """
    對「某個偏壓 V ⇒ barrier_shift_extra」計算 T(E)
    """
    epsR_eff = epsR0 + barrier_shift_extra
    Hc, tLc, tRc = build_hetero_device_biased(
        epsL, tL, epsR0, tR, NL, NR,
        barrier_shift_extra,
        eta_dephase=eta_dephase,
        apply_region=apply_region
    )
    Tvals = []
    for E in E_grid:
        Tvals.append(
            transmission_at_E(E, Hc, epsL, tL, epsR_eff, tR, tLc, tRc)
        )
    return np.array(Tvals)

# -----------------------------
# 4) I–V (0 K Landauer)
# -----------------------------
def fermi_step(E, mu):  # 0 K
    return (E < mu).astype(float)

def current_from_TE(E_grid, T_E, muL, muR):
    fL = fermi_step(E_grid, muL)
    fR = fermi_step(E_grid, muR)
    return np.trapz(T_E * (fL - fR), E_grid)

# -----------------------------
# 5) Main: 掃描去相干，畫 T(E) 與 I–V
# -----------------------------
if __name__ == "__main__":
    # 幾何與材料參數（這組在你先前就能看到強整流）
    NL, NR = 8, 16
    epsL, tL = 0.0, -3.0
    epsR0, tR = +2.0, -1.0

    # 偏壓如何改變障壁： barrier_shift_extra = alpha * V
    # alpha > 0: 正偏抬高障壁（壞），負偏壓低障壁（好）
    # 若想讓「正偏是順向」，把 alpha 改成負值。
    alpha = -0.5

    # 電極本底化學勢（功函數差）
    muL0, muR0 = 0.0, -0.8

    # 能量/偏壓掃描網格
    E_grid = np.linspace(-3, 3, 800)
    V_list = np.linspace(-2.0, 2.0, 81)

    # 去相干掃描（eV）
    dephase_list = [0.0, 1e-3, 1e-2, 5e-2]

    # 只在右段（障壁區）加去相干，更貼近「半導體段散射較多」的直覺
    dephase_region = "right"

    # ---- 圖 1：T(E) vs dephasing（取代表性偏壓 V=0 比較） ----
    plt.figure(figsize=(7,4))
    for eta_d in dephase_list:
        barrier_shift_extra = alpha * 0.0  # at V=0
        Tvals = compute_TE_for_bias(
            E_grid, epsL, tL, epsR0, tR, NL, NR,
            barrier_shift_extra,
            eta_dephase=eta_d,
            apply_region=dephase_region
        )
        plt.plot(E_grid, Tvals, label=f"η_deph={eta_d:g} eV")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Transmission T(E)")
    plt.title(f"T(E) vs dephasing (V=0, dephase in {dephase_region})")
    plt.grid(True, ls="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # ---- 圖 2：I–V vs dephasing（每個 V 都重建 H(V) 再積分）----
    plt.figure(figsize=(7,4))
    # --- 掃去相干參數，畫 I–V 曲線並存成圖 + CSV ---
for eta_d in dephase_list:
    I_vals = []
    for V in V_list:
        barrier_shift = alpha * V

        # 對每個偏壓計算 T(E)
        T_E = compute_TE_for_bias(
            E_grid,
            epsL, tL,
            epsR0, tR,
            NL, NR,
            barrier_shift,
            eta_dephase=eta_d,
            apply_region=dephase_region
        )

        # 設置 μ_L, μ_R 並積分得到 I(V)
        muL = muL0 + 0.5*V
        muR = muR0 - 0.5*V
        I_vals.append(current_from_TE(E_grid, T_E, muL, muR))

    I_vals = np.array(I_vals)

    # 在同一張圖上畫不同 η_deph 的 I–V
    plt.plot(V_list, I_vals, label=f"η_deph={eta_d:g} eV")

    # ★ 把每一條曲線存成一個 CSV 檔
    csv_name = f"IV_dephasing_{dephase_region}_eta_{eta_d:g}.csv"
    np.savetxt(
        csv_name,
        np.column_stack([V_list, I_vals]),
        delimiter=",",
        header="V(eV),Current(arb.)",
        comments=""
    )

# === 整張圖的外觀與存檔 ===
plt.xlabel("Bias V (eV)")
plt.ylabel("Relative current (arb.)")
plt.title(f"I–V vs dephasing (dephase in {dephase_region})")
plt.grid(True, ls="--", alpha=0.5)
plt.legend()
plt.tight_layout()

# ★ 存成圖片
fig_name_png = f"IV_dephasing_{dephase_region}.png"
fig_name_pdf = f"IV_dephasing_{dephase_region}.pdf"
plt.savefig(fig_name_png, dpi=300)
plt.savefig(fig_name_pdf)
plt.close()

