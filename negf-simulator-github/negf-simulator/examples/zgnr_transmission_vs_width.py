import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm  # for automatic color map

########################
# tight-binding for ZGNR (t1 + t2)
########################
def build_zgnr_cell_with_t2(Ny, t1, t2):
    dim = 2 * Ny
    H0 = np.zeros((dim, dim), dtype=complex)
    V  = np.zeros((dim, dim), dtype=complex)

    for y in range(Ny):
        A = 2*y
        B = 2*y + 1

        # 最近鄰 A<->B (同 slice)
        H0[A,B] = t1
        H0[B,A] = t1

        # slice 間 B_y -> A_y (右一個 slice)
        V[B,A] = t1

        # slice 間 B_y -> A_(y-1) (斜鍵)
        if y > 0:
            A_prev = 2*(y-1)
            V[B, A_prev] = t1

        # 次近鄰 t2: 同子格 (A_y <-> A_{y+1}, B_y <-> B_{y+1})
        if y < Ny-1:
            A_up = 2*(y+1)
            B_up = 2*(y+1) + 1
            H0[A, A_up] = t2
            H0[A_up, A] = t2
            H0[B, B_up] = t2
            H0[B_up, B] = t2

    return H0, V

########################
# surface GF with safety
########################
def surface_gf_block(E, H0, V, eta=1e-6, max_iter=200, tol=1e-10):
    n = H0.shape[0]
    I = np.eye(n, dtype=complex)

    alpha = V.copy()
    beta  = V.conj().T.copy()
    eps   = H0.copy()

    z = E + 1j*eta

    for _ in range(max_iter):
        M = z*I - eps
        # 防止不可逆或數值爆掉
        try:
            M_inv = np.linalg.inv(M)
        except np.linalg.LinAlgError:
            # 如果發散，回傳零矩陣 -> 傳輸那點可能變小，但程式不會死
            return np.zeros((n,n), dtype=complex)

        eps_new   = eps + alpha @ M_inv @ beta
        alpha_new = alpha @ M_inv @ alpha
        beta_new  = beta  @ M_inv @ beta

        # 收斂就收手
        if (np.linalg.norm(eps_new - eps) < tol and
            np.linalg.norm(alpha_new - alpha) < tol and
            np.linalg.norm(beta_new  - beta)  < tol):
            eps = eps_new
            break

        eps   = eps_new
        alpha = alpha_new
        beta  = beta_new

    # 最終 surface GF
    try:
        gsurf = np.linalg.inv(z*I - eps)
    except np.linalg.LinAlgError:
        gsurf = np.zeros((n,n), dtype=complex)

    return gsurf

########################
# NEGF helpers
########################
def self_energy_block(gsurf, V_couple):
    return V_couple @ gsurf @ V_couple.conj().T

def embed_self_energy(Sigma_edge, N_blocks, block_size, which):
    dim_full = N_blocks * block_size
    Sigma_full = np.zeros((dim_full, dim_full), dtype=complex)
    if which == "left":
        Sigma_full[:block_size, :block_size] = Sigma_edge
    else:
        Sigma_full[-block_size:, -block_size:] = Sigma_edge
    return Sigma_full

def retarded_green(E, Hc, SigmaL_full, SigmaR_full, eta=1e-6):
    I = np.eye(Hc.shape[0], dtype=complex)
    return np.linalg.inv((E + 1j*eta)*I - Hc - SigmaL_full - SigmaR_full)

def gamma_matrix(Sigma):
    return 1j * (Sigma - Sigma.conj().T)

def transmission(E, Hc, VL_ifc, VR_ifc, gL, gR, N_blocks, block_size):
    SigmaL_edge = self_energy_block(gL, VL_ifc)
    SigmaR_edge = self_energy_block(gR, VR_ifc)

    SigmaL_full = embed_self_energy(SigmaL_edge, N_blocks, block_size, "left")
    SigmaR_full = embed_self_energy(SigmaR_edge, N_blocks, block_size, "right")

    G = retarded_green(E, Hc, SigmaL_full, SigmaR_full)

    GammaL_full = gamma_matrix(SigmaL_full)
    GammaR_full = gamma_matrix(SigmaR_full)

    Tval = np.trace(GammaL_full @ G @ GammaR_full @ G.conj().T).real
    return Tval

########################
# 掃 Ny，畫一張圖
########################

# 掃寬度：不是每一個 Ny=4,5,...25 全畫，而是幾個代表
Ny_list = [4, 6, 8, 10, 14, 18, 22, 25]

Nx_c = 20          # channel 長度 (slice 數)
t1   = -3.0       # 最近鄰 hopping
t2   = -0.22      # 次近鄰 hopping
energies = np.linspace(-6, 6, 150)   # 比原來 400 點輕很多

colors = cm.viridis(np.linspace(0, 1, len(Ny_list)))

plt.figure(figsize=(8,5))

for Ny, color in zip(Ny_list, colors):
    print(f"running Ny={Ny} ...")  # <- 讓你在 console 確認它有在跑

    # 1. 建 ribbon slice H0, V
    H0, V = build_zgnr_cell_with_t2(Ny, t1, t2)

    block_size = H0.shape[0]        # = 2*Ny
    N_blocks   = Nx_c               # 幾個 slice
    dim_c      = block_size * N_blocks

    # 2. 組中央 Hc (block-tridiagonal)
    HC = np.zeros((dim_c, dim_c), dtype=complex)
    for n in range(N_blocks):
        i0, i1 = n*block_size, (n+1)*block_size
        HC[i0:i1, i0:i1] = H0
        if n < N_blocks-1:
            j0, j1 = (n+1)*block_size, (n+2)*block_size
            HC[i0:i1, j0:j1] = V
            HC[j0:j1, i0:i1] = V.conj().T

    # 3. 接面耦合
    VL_ifc = V
    VR_ifc = V.conj().T

    # 4. 掃能量
    Tvals = []
    for E in energies:
        gL = surface_gf_block(E, H0, V)
        gR = surface_gf_block(E, H0, V)
        Tvals.append(
            transmission(E, HC, VL_ifc, VR_ifc, gL, gR, N_blocks, block_size)
        )

    # 5. 畫線
    plt.plot(energies, Tvals, color=color, lw=1.5, label=f"Ny={Ny}")

# ====== 美化圖表 ======
plt.xlabel("Energy (eV)", fontsize=12)
plt.ylabel("Transmission T(E)", fontsize=12)
plt.title("ZGNR Transmission vs Width Ny (t2 = -0.22 eV)", fontsize=13)
plt.legend(
    title="Ribbon width (Ny)",
    fontsize=9,
    ncol=2,
    loc="upper left",
    frameon=False
)
plt.grid(True, ls="--", alpha=0.6)
plt.tight_layout()
plt.show()
