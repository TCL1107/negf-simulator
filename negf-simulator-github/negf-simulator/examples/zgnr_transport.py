import numpy as np
import matplotlib.pyplot as plt

########################################
# 1. 建立 ZGNR slice 的 TB
########################################

def build_zgnr_cell(Ny, t1):
    """
    回傳:
    H0 : (2Ny x 2Ny) 單一 slice (column) 的哈密頓量
    V  : (2Ny x 2Ny) 這個 slice -> 下一個 slice 之間的耦合

    拿的是最簡單的最近鄰 π-band 模型。
    Zigzag ribbon 的幾何可以這樣想：
    slice n 跟 slice n+1 之間只有特定 A-B 配對有鍵結。
    """
    dim = 2 * Ny
    H0 = np.zeros((dim, dim), dtype=complex)
    V  = np.zeros((dim, dim), dtype=complex)

    # 我們用一個常見的 ZGNR indexing：
    # 對於每個 "dimer line" y = 0..Ny-1:
    #   site A_y = 2*y
    #   site B_y = 2*y + 1
    #
    # 內部 (同一 slice) 有 A_y 連到 B_y (垂直鍵)
    # slice 間有 B_y 連到 A_y (往右)，以及 B_y 連到 A_{y-1} (斜向)
    #
    # 這是標準 zigzag TB 拓撲。

    for y in range(Ny):
        A = 2*y
        B = 2*y + 1

        # 同一 slice 裡 A_y <-> B_y
        H0[A, B] = t1
        H0[B, A] = t1

        # slice 間耦合: B_y (這 slice) -> A_y (下一 slice)
        V[B, A] = t1

        # slice 間耦合: B_y (這 slice) -> A_{y-1} (下一 slice), if y>0
        if y > 0:
            A_prev = 2*(y-1)
            V[B, A_prev] = t1

    return H0, V

########################################
# 2. surface Green's function for block lead
########################################

def surface_gf_block(E, H0, V, eta=1e-6, max_iter=200, tol=1e-12):
    """
    Sancho-like iterative decimation for multi-orbital semi-infinite lead
    """
    n = H0.shape[0]
    I = np.eye(n, dtype=complex)

    alpha = V.copy()            # forward
    beta  = V.conj().T.copy()   # backward
    eps   = H0.copy()

    z = E + 1j*eta

    for _ in range(max_iter):
        M = z*I - eps
        M_inv = np.linalg.inv(M)

        eps_new   = eps + alpha @ M_inv @ beta
        alpha_new = alpha @ M_inv @ alpha
        beta_new  = beta  @ M_inv @ beta

        if (np.linalg.norm(eps_new - eps) < tol and
            np.linalg.norm(alpha_new - alpha) < tol and
            np.linalg.norm(beta_new  - beta)  < tol):
            eps = eps_new
            break

        eps   = eps_new
        alpha = alpha_new
        beta  = beta_new

    gsurf = np.linalg.inv(z*I - eps)
    return gsurf

########################################
# 3. NEGF core
########################################

def self_energy_block(gsurf, V_couple):
    # Σ = V * g * V^\dagger
    return V_couple @ gsurf @ V_couple.conj().T

def embed_self_energy(Sigma_edge, N_blocks, block_size, which):
    dim_full = N_blocks * block_size
    Sigma_full = np.zeros((dim_full, dim_full), dtype=complex)
    if which == "left":
        Sigma_full[0:block_size, 0:block_size] = Sigma_edge
    elif which == "right":
        start = dim_full - block_size
        Sigma_full[start:dim_full, start:dim_full] = Sigma_edge
    else:
        raise ValueError
    return Sigma_full

def retarded_green(E, Hc, SigmaL_full, SigmaR_full, eta=1e-6):
    dim = Hc.shape[0]
    I = np.eye(dim, dtype=complex)
    M = (E + 1j*eta)*I - Hc - SigmaL_full - SigmaR_full
    return np.linalg.inv(M)

def gamma_matrix(Sigma):
    return 1j * (Sigma - Sigma.conj().T)

def transmission(E, Hc,
                 VL_ifc, VR_ifc,
                 gL, gR,
                 N_blocks, block_size):
    SigmaL_edge = self_energy_block(gL, VL_ifc)
    SigmaR_edge = self_energy_block(gR, VR_ifc)

    SigmaL_full = embed_self_energy(SigmaL_edge, N_blocks, block_size, "left")
    SigmaR_full = embed_self_energy(SigmaR_edge, N_blocks, block_size, "right")

    G = retarded_green(E, Hc, SigmaL_full, SigmaR_full)

    GammaL_full = gamma_matrix(SigmaL_full)
    GammaR_full = gamma_matrix(SigmaR_full)

    Tval = np.trace(GammaL_full @ G @ GammaR_full @ G.conj().T).real
    return Tval

########################################
# 4. 把所有東西串起來跑
########################################

# 參數：帶寬 Ny、中心長度 Nx_c
Ny   = 4          # ribbon 寬度 = 幾條dimer lines
Nx_c = 4          # 中央散射區的 slice 數
t1   = -3.0       # 來自你貼的表的第一個 hopping (最近鄰)

# 建立單一 slice 的 H0, V (lead單元)
H0, V = build_zgnr_cell(Ny, t1)

block_size = H0.shape[0]        # = 2*Ny
N_blocks   = Nx_c               # 中央有多少 slice
dim_c      = block_size * N_blocks  # HC 的總維度

# 組中央散射區 HC：把 Nx_c 個 slice 串起來成一條鏈
HC = np.zeros((dim_c, dim_c), dtype=complex)
for n in range(N_blocks):
    # 放 H0 在對角 block
    i0 = n*block_size
    i1 = i0+block_size
    HC[i0:i1, i0:i1] = H0
    # 與下一個 slice 的耦合 V
    if n < N_blocks-1:
        j0 = (n+1)*block_size
        j1 = j0+block_size
        HC[i0:i1, j0:j1] = V
        HC[j0:j1, i0:i1] = V.conj().T

# interface 耦合矩陣：
# 左電極的 surface slice ↔ 中央散射區最左 slice
VL_ifc = V              # (block_size x block_size)
# 中央散射區最右 slice ↔ 右電極的 surface slice
VR_ifc = V.conj().T     # 用對稱方向

energies = np.linspace(-10, 10, 400)
Tvals = []

debug_once = True
for E in energies:
    gL = surface_gf_block(E, H0, V)
    gR = surface_gf_block(E, H0, V)

    if debug_once:
        SigmaL_edge = self_energy_block(gL, VL_ifc)
        print("||SigmaL_edge|| =", np.linalg.norm(SigmaL_edge))
        debug_once = False

    T_E = transmission(E, HC,
                       VL_ifc, VR_ifc,
                       gL, gR,
                       N_blocks, block_size)
    Tvals.append(T_E)

plt.figure()
plt.plot(energies, Tvals)
plt.xlabel("Energy (eV)")
plt.ylabel("Transmission T(E)")
plt.title("ZGNR NEGF with first-neighbor t1")
plt.grid(True)
plt.ylim(bottom=0)
plt.show()
