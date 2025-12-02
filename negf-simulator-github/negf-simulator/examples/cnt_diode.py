import numpy as np
import matplotlib.pyplot as plt

########################################
# 1. 解析式 semi-infinite 1D lead surface Green's function
########################################
def surface_gf_1d(E, eps0, t, eta=1e-3):
    """
    Analytic surface Green's function for a semi-infinite 1D tight-binding chain.
    Chain onsite = eps0, nearest-neighbor hopping = t (typically negative).
    eta is a small positive imaginary part to ensure causality / numerical stability.
    """
    z = E + 1j*eta
    delta = z - eps0

    # 帶寬條件: band is roughly eps0 ± 2|t|
    band_edge = 2.0 * abs(t)

    if np.abs(delta.real) <= band_edge:
        # Inside band: g has an imaginary part ~ propagating states
        root = np.sqrt((2*t)**2 - delta**2)  # sqrt(4 t^2 - (E-eps0)^2)
        g = (delta - 1j*root) / (2 * (t**2))
    else:
        # Outside band: purely evanescent (imag part ~0 or tiny)
        root = np.sqrt(delta**2 - (2*t)**2)
        # choose sign so Im(g) <= 0 (retarded)
        sgn = np.sign(delta.real) if delta.real != 0 else 1.0
        g = (delta - sgn*root) / (2 * (t**2))

    return g


########################################
# 2. Self-energy and Gamma
########################################
def self_energy_scalar(gsurf, t_couple):
    # Σ = τ g_s τ*
    return t_couple * gsurf * np.conjugate(t_couple)

def gamma_from_sigma(Sigma):
    return 1j * (Sigma - np.conjugate(Sigma))


########################################
# 3. 建構異質結構(左段+右段)
########################################
def build_hetero_device(epsL, tL, epsR, tR,
                        NL, NR,
                        barrier_shift=0.0):
    """
    我們做一條 1D chain:
    - 左段: NL sites, on-site epsL
    - 右段: NR sites, on-site epsR (+ barrier_shift 整段往上抬)
    - 左區 hopping = tL
    - 右區 hopping = tR
    - 界面 hopping = 平均 (tL+tR)/2
    """
    N = NL + NR
    Hc = np.zeros((N, N), dtype=complex)

    # 左半段 on-site
    for i in range(NL):
        Hc[i, i] = epsL

    # 右半段 on-site
    for j in range(NR):
        Hc[NL + j, NL + j] = epsR + barrier_shift

    # 左半段內部 hopping
    for i in range(NL - 1):
        Hc[i, i+1] = tL
        Hc[i+1, i] = tL

    # 右半段內部 hopping
    for j in range(NR - 1):
        a = NL + j
        b = NL + j + 1
        Hc[a, b] = tR
        Hc[b, a] = tR

    # 交界面 hopping
    Hc[NL-1, NL] = 0.5*(tL + tR)
    Hc[NL, NL-1] = 0.5*(tL + tR)

    # leads 假設是各自區段的延伸
    t_left_contact  = tL
    t_right_contact = tR

    return Hc, t_left_contact, t_right_contact, {
        "N":  N,
        "NL": NL,
        "NR": NR
    }


########################################
# 4. Green's function, T(E)
########################################
def retarded_green(E, Hc, SigmaL, SigmaR, eta=1e-9):
    N = Hc.shape[0]
    I = np.eye(N, dtype=complex)
    M = (E + 1j*eta)*I - Hc - SigmaL - SigmaR
    return np.linalg.inv(M)

def transmission_at_E(E,
                      Hc,
                      epsL, tL,
                      epsR, tR,
                      t_left_contact,
                      t_right_contact):
    """
    用 NEGF 計算單能量下的傳輸 T(E)
    """
    N = Hc.shape[0]

    # 1. 左右 leads 的 surface GF
    gL = surface_gf_1d(E, epsL, tL)
    gR = surface_gf_1d(E, epsR, tR)

    # 2. Self-energies
    SigmaL = np.zeros((N, N), dtype=complex)
    SigmaR = np.zeros((N, N), dtype=complex)

    SigmaL_edge = self_energy_scalar(gL, t_left_contact)
    SigmaR_edge = self_energy_scalar(gR, t_right_contact)

    SigmaL[0,0]   = SigmaL_edge
    SigmaR[-1,-1] = SigmaR_edge

    # 3. G^r
    G = retarded_green(E, Hc, SigmaL, SigmaR)

    # 4. Γ
    GammaL = np.zeros((N, N), dtype=complex)
    GammaR = np.zeros((N, N), dtype=complex)
    GammaL[0,0]   = gamma_from_sigma(SigmaL_edge)
    GammaR[-1,-1] = gamma_from_sigma(SigmaR_edge)

    # 5. Landauer 傳輸: T = Tr[ ΓL G ΓR G^\dagger ]
    Tval = np.trace(GammaL @ G @ GammaR @ G.conjugate().T).real
    return max(Tval, 0.0)


def compute_TE_spectrum(Hc,
                        epsL, tL,
                        epsR, tR,
                        t_left_contact,
                        t_right_contact,
                        E_grid):
    Tvals = []
    for E in E_grid:
        Tvals.append(
            transmission_at_E(E, Hc,
                              epsL, tL,
                              epsR, tR,
                              t_left_contact,
                              t_right_contact)
        )
    return np.array(Tvals)


########################################
# 5. I-V 計算，這次加「不同電極功函數」→ Schottky
########################################
def fermi_step(E, mu):
    # T=0K 近似: E < mu → 1, E > mu → 0
    return (E < mu).astype(float)

def current_vs_bias_schottky(V_list,
                             E_grid,
                             T_of_E,
                             muL0,
                             muR0):
    """
    muL0, muR0 是「偏壓=0」時左右接觸的費米能。
    I(V) ∝ ∫ dE T(E) [ f_L(E) - f_R(E) ]
    """
    I_list = []
    for V in V_list:
        muL = muL0 + 0.5*V
        muR = muR0 - 0.5*V

        fL = fermi_step(E_grid, muL)
        fR = fermi_step(E_grid, muR)

        diff = fL - fR
        integrand = T_of_E * diff
        I_val = np.trapz(integrand, E_grid)
        I_list.append(I_val)

    return np.array(I_list)


########################################
# 6. 主程式：建立 CNT-like Schottky 二極體，畫 T(E) 和 I-V
########################################
if __name__ == "__main__":
    # --- 幾何長度設定 ---
    NL = 12   # 左半段長度
    NR = 16   # 右半段長度

    # --- 左邊 (近似金屬 / n-type CNT contact) ---
    epsL = 0.0     # on-site (eV)
    tL   = -3.0    # hopping  (eV)

    # --- 右邊 (近似半導體 / barrier CNT side) ---
    epsR = +2.1    # 把整個能帶往上抬，像有能障/價帶已遠離
    tR   = -1.0    # 稍微不同的有效質量 (改帶寬)

    barrier_shift = 0.0
    Hc, tLc, tRc, info = build_hetero_device(
        epsL, tL,
        epsR, tR,
        NL, NR,
        barrier_shift=barrier_shift
    )

    # --- 能量掃描範圍 ---
    E_grid = np.linspace(-3, 3, 800)

    # --- 計算 T(E) ---
    T_E = compute_TE_spectrum(
        Hc,
        epsL, tL,
        epsR, tR,
        tLc, tRc,
        E_grid
    )

    # === 畫傳輸譜並儲存 ===
    plt.figure(figsize=(7,4))
    plt.plot(E_grid, T_E, lw=1.5)
    plt.xlabel("Energy E (eV)")
    plt.ylabel("Transmission T(E)")
    plt.title("Transmission of CNT-like heterojunction (metal ↔ semiconducting)")
    plt.grid(True, ls="--", alpha=0.5)
    plt.tight_layout()

    plt.savefig("TE_spectrum.png", dpi=300)
    plt.savefig("TE_spectrum.pdf")
    plt.close()

    np.savetxt(
        "TE_spectrum.csv",
        np.column_stack([E_grid, T_E]),
        delimiter=",",
        header="E_grid(eV),T(E)",
        comments=""
    )

    # --- I-V 計算 ---
    muL0 = 0.2   # 左接觸的自然化學勢 (像金屬 work function)
    muR0 = -0.8  # 右接觸比較深，像是能帶下沉 → 注入比較難

    V_list = np.linspace(-2.0, 2.0, 81)
    I_list = current_vs_bias_schottky(
        V_list,
        E_grid,
        T_E,
        muL0,
        muR0
    )

    # === 畫 I-V 並儲存 ===
    plt.figure(figsize=(7,4))
    plt.plot(V_list, I_list, lw=1.8, color="tomato")
    plt.xlabel("Bias V (eV)")
    plt.ylabel("Relative current (arb. units)")
    plt.title("I-V of CNT heterojunction (Schottky-like rectification)")
    plt.grid(True, ls="--", alpha=0.5)
    plt.tight_layout()

    plt.savefig("IV_curve.png", dpi=300)
    plt.savefig("IV_curve.pdf")
    plt.close()

    np.savetxt(
        "IV_curve.csv",
        np.column_stack([V_list, I_list]),
        delimiter=",",
        header="V_list(eV),Current(arbitrary units)",
        comments=""
    )
