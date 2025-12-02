# analyze_sweeps.py (robust variant)
# 讀取 sweep_* 結果，做穩健欄位偵測、變異分解，輸出熱圖與比例摘要。
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ============ 通用工具 ============

def _norm_name(s: str) -> str:
    return str(s).strip().lower().replace(" ", "").replace("-", "").replace("±", "+").replace("__", "_")

def find_col(df, predicates):
    """
    在 df.columns 中以寬鬆規則尋找目標欄位。
    predicates: list[callable(normed_name)->bool]
    """
    for c in df.columns:
        cn = _norm_name(c)
        for p in predicates:
            try:
                if p(cn):
                    return c
            except Exception:
                pass
    return None

def load_csv(fname):
    p = Path(fname)
    if not p.exists():
        raise FileNotFoundError(f"找不到檔案：{p.resolve()}")
    df = pd.read_csv(p)
    return df

def to_numeric_df(df, cols):
    out = df.copy()
    for c in cols:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")
    return out

def pivot_to_grid(df, row_col, col_col, val_col, sort=True, agg="mean"):
    """
    把 (row, col) → value 的寬表變成 numpy 網格；自動平均重複點。
    回傳 (row_values, col_values, Z)
    """
    dfg = df[[row_col, col_col, val_col]].copy()
    dfg = to_numeric_df(dfg, [row_col, col_col, val_col]).dropna()
    if dfg.empty:
        raise ValueError("pivot 資料全為空：請檢查欄位或數值。")
    dfg = dfg.groupby([row_col, col_col], as_index=False)[val_col].mean()
    piv = dfg.pivot_table(index=row_col, columns=col_col, values=val_col, aggfunc=agg)
    if sort:
        piv = piv.sort_index(axis=0).sort_index(axis=1)
    return piv.index.values, piv.columns.values, piv.to_numpy()

def smooth_along_rows(Z, k=3):
    """沿列做移動平均平滑（奇數 k）。"""
    k = int(k); k = k if k % 2 == 1 and k >= 1 else 3
    pad = k // 2
    Zpad = np.pad(Z, ((pad, pad), (0, 0)), mode="edge")
    out = np.zeros_like(Z, dtype=float)
    for i in range(Z.shape[0]):
        out[i] = Zpad[i:i+k].mean(axis=0)
    return out

def decompose_grid(Y, smooth_k=3):
    """
    把 2D 響應 Y(row×col) 分解：
      - row 主效應（SS_row）
      - col 主效應（SS_col）
      - 低頻交互（SS_int）
      - 高頻紋理（SS_qtz, 量子條紋）
    皆以總平方和 SS_tot 正規化成百分比。
    對常數陣列（SS_tot=0）回報 0%，避免 NaN。
    """
    Y = np.asarray(Y, dtype=float)
    m = np.nanmean(Y)
    # 若全 NaN
    if np.all(~np.isfinite(Y)):
        return {
            "tot": 0.0, "row_main_pct": 0.0, "col_main_pct": 0.0,
            "interaction_pct": 0.0, "quantization_pct": 0.0,
            "remainder_pct": 0.0, "parts": (np.zeros_like(Y),)*4, "constant": True
        }

    # 用有限值填補整體平均（避免極少數 NaN 破壞分解）
    Y = np.where(np.isfinite(Y), Y, m)

    r_mean = Y.mean(axis=1, keepdims=True)
    c_mean = Y.mean(axis=0, keepdims=True)
    SS_tot = np.sum((Y - m)**2)

    if SS_tot == 0.0:
        return {
            "tot": 0.0, "row_main_pct": 0.0, "col_main_pct": 0.0,
            "interaction_pct": 0.0, "quantization_pct": 0.0,
            "remainder_pct": 0.0, "parts": (r_mean*0, c_mean*0, Y*0, Y*0), "constant": True
        }

    SS_row = Y.shape[1] * np.sum((r_mean - m)**2)
    SS_col = Y.shape[0] * np.sum((c_mean - m)**2)

    Interact = Y - (r_mean + c_mean - m)
    Inter_s  = smooth_along_rows(Interact, k=smooth_k)
    HighFreq = Interact - Inter_s

    SS_int = np.sum(Inter_s**2)
    SS_qtz = np.sum(HighFreq**2)

    pct = lambda s: 100.0 * s / SS_tot
    rem = max(0.0, 100.0 - sum([pct(SS_row), pct(SS_col), pct(SS_int), pct(SS_qtz)]))
    return {
        "tot": SS_tot,
        "row_main_pct": pct(SS_row),
        "col_main_pct": pct(SS_col),
        "interaction_pct": pct(SS_int),
        "quantization_pct": pct(SS_qtz),
        "remainder_pct": rem,
        "parts": (r_mean, c_mean, Inter_s, HighFreq),
        "constant": False
    }

def _extent_for_axes(xvals, yvals):
    # 支援非等距座標
    x = np.array(xvals, dtype=float); y = np.array(yvals, dtype=float)
    if len(x) > 1:
        dx = np.diff(x).mean()
    else:
        dx = 1.0
    if len(y) > 1:
        dy = np.diff(y).mean()
    else:
        dy = 1.0
    return [x[0]-dx/2, x[-1]+dx/2, y[0]-dy/2, y[-1]+dy/2]

def heatmap(Z, xvals, yvals, title, cbar_label, xlabel, ylabel, fname=None,
            log10=False, vmin=None, vmax=None):
    Zp = np.array(Z, dtype=float)
    if log10:
        Zp = np.log10(np.maximum(Zp, 1e-12))
        cbar_label = f"log10({cbar_label})"
    extent = _extent_for_axes(xvals, yvals)

    plt.figure(figsize=(6.0, 4.6))
    im = plt.imshow(Zp, origin="lower", aspect="auto",
                    extent=extent, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(im); cbar.set_label(cbar_label)
    plt.xlabel(xlabel); plt.ylabel(ylabel); plt.title(title)
    plt.tight_layout()
    if fname:
        plt.savefig(fname, dpi=180)
    plt.show()

def print_report(name, result):
    print(f"\n=== Variance share for {name} ===")
    if result.get("constant", False):
        print(" （此網格無變化：SS_tot=0 → 百分比皆 0）")
        return
    print(f" row(main): {result['row_main_pct']:5.1f}%")
    print(f" col(main): {result['col_main_pct']:5.1f}%")
    print(f" interact : {result['interaction_pct']:5.1f}%")
    print(f" quantize : {result['quantization_pct']:5.1f}%")
    if result["remainder_pct"] > 0.1:
        print(f" remainder: {result['remainder_pct']:5.1f}%")

def append_variance_csv(rows, path="variance_report.csv"):
    df = pd.DataFrame(rows)
    if Path(path).exists():
        # 追加
        old = pd.read_csv(path)
        df = pd.concat([old, df], axis=0, ignore_index=True)
    df.to_csv(path, index=False)

# ============ A) 分析 (NR, alpha) ============

def analyze_NR_alpha(csv="sweep_NR_alpha.csv", smooth_k=3):
    df = load_csv(csv)
    print("CSV columns:", list(df.columns))

    col_NR    = find_col(df, [lambda s: s in {"nr", "n_r"} or s.startswith("nr")])
    col_alpha = find_col(df, [lambda s: s.startswith("alpha")])
    col_von   = find_col(df, [
        lambda s: s in {"von","v_on","vonabs","v_on_abs","vonarbv","vonrel","v_onrel"},
        lambda s: s.startswith("von")
    ])
    col_rr    = find_col(df, [
        lambda s: s in {"rr","rrarea","rr_point","rrpoint","rr_area","rrarea15v","rrarea+15v"},
        lambda s: s.startswith("rr")
    ])

    if col_NR is None or col_alpha is None:
        raise ValueError("需要欄位：NR 與 alpha")

    rows_for_csv = []

    # 1) V_on
    if col_von is not None:
        r_vals, c_vals, Y = pivot_to_grid(df, col_NR, col_alpha, col_von)
        res = decompose_grid(Y, smooth_k=smooth_k)
        print_report(f"V_on vs ({col_NR}, {col_alpha}) [k={smooth_k}]", res)
        rows_for_csv.append({
            "which": "V_on", "row": col_NR, "col": col_alpha,
            "row_main_pct": res["row_main_pct"], "col_main_pct": res["col_main_pct"],
            "interaction_pct": res["interaction_pct"], "quantization_pct": res["quantization_pct"],
            "remainder_pct": res["remainder_pct"], "SS_tot": res["tot"]
        })
        # 全域色階
        vmin, vmax = np.nanmin(Y), np.nanmax(Y)
        heatmap(Y, c_vals, r_vals,
                f"Turn-on voltage vs ({col_NR}, {col_alpha})",
                "V_on (arb. V)", xlabel=col_alpha, ylabel=col_NR,
                fname="decomp_von_raw.png", vmin=vmin, vmax=vmax)
        (r_mean, c_mean, Inter, High) = res["parts"]
        heatmap(Inter, c_vals, r_vals,
                "Interaction (low-freq) of V_on", "arb.",
                xlabel=col_alpha, ylabel=col_NR,
                fname="decomp_von_inter.png")
        heatmap(High, c_vals, r_vals,
                "Quantization texture (high-freq) of V_on", "arb.",
                xlabel=col_alpha, ylabel=col_NR,
                fname="decomp_von_quant.png")
    else:
        print("（提醒）此 CSV 中找不到 V_on 欄位，略過 V_on 分析。")

    # 2) RR
    if col_rr is not None:
        r_vals, c_vals, Y = pivot_to_grid(df, col_NR, col_alpha, col_rr)
        res = decompose_grid(Y, smooth_k=smooth_k)
        print_report(f"RR vs ({col_NR}, {col_alpha}) [k={smooth_k}]", res)
        rows_for_csv.append({
            "which": "RR", "row": col_NR, "col": col_alpha,
            "row_main_pct": res["row_main_pct"], "col_main_pct": res["col_main_pct"],
            "interaction_pct": res["interaction_pct"], "quantization_pct": res["quantization_pct"],
            "remainder_pct": res["remainder_pct"], "SS_tot": res["tot"]
        })
        vmin, vmax = np.nanmin(np.log10(np.maximum(Y, 1e-12))), np.nanmax(np.log10(np.maximum(Y, 1e-12)))
        heatmap(Y, c_vals, r_vals,
                f"log10(RR) vs ({col_NR}, {col_alpha})", "RR",
                xlabel=col_alpha, ylabel=col_NR,
                fname="decomp_rr_raw.png", log10=True, vmin=vmin, vmax=vmax)
        (r_mean, c_mean, Inter, High) = res["parts"]
        heatmap(Inter, c_vals, r_vals,
                "Interaction (low-freq) of RR", "arb.",
                xlabel=col_alpha, ylabel=col_NR,
                fname="decomp_rr_inter.png")
        heatmap(High, c_vals, r_vals,
                "Quantization texture (high-freq) of RR", "arb.",
                xlabel=col_alpha, ylabel=col_NR,
                fname="decomp_rr_quant.png")
    else:
        print("（提醒）此 CSV 中找不到 RR 欄位，略過 RR 分析。")

    if rows_for_csv:
        append_variance_csv(rows_for_csv)

# ============ B) 分析 (epsR, tR) ============

def analyze_epsR_tR(csv="sweep_epsR_tR.csv", smooth_k=3):
    df = load_csv(csv)
    print("CSV columns:", list(df.columns))

    col_epsR = find_col(df, [lambda s: s in {"epsr","epsilonr","eps_r"} or s.startswith("epsr")])
    col_tR   = find_col(df, [lambda s: s in {"tr","t_r"} or s == "tr"])
    col_von  = find_col(df, [
        lambda s: s in {"von","v_on","vonabs","v_on_abs","vonarbv","vonrel","v_onrel"},
        lambda s: s.startswith("von")
    ])
    col_rr   = find_col(df, [
        lambda s: s in {"rr","rrarea","rr_point","rrpoint","rr_area","rrarea15v","rrarea+15v"},
        lambda s: s.startswith("rr")
    ])

    if col_epsR is None or col_tR is None:
        raise ValueError("需要欄位：epsR 與 tR")

    rows_for_csv = []

    if col_von is not None:
        r_vals, c_vals, Y = pivot_to_grid(df, col_epsR, col_tR, col_von)
        res = decompose_grid(Y, smooth_k=smooth_k)
        print_report(f"V_on vs ({col_epsR}, {col_tR}) [k={smooth_k}]", res)
        rows_for_csv.append({
            "which": "V_on", "row": col_epsR, "col": col_tR,
            "row_main_pct": res["row_main_pct"], "col_main_pct": res["col_main_pct"],
            "interaction_pct": res["interaction_pct"], "quantization_pct": res["quantization_pct"],
            "remainder_pct": res["remainder_pct"], "SS_tot": res["tot"]
        })
        vmin, vmax = np.nanmin(Y), np.nanmax(Y)
        heatmap(Y, c_vals, r_vals,
                f"Turn-on voltage vs ({col_epsR}, {col_tR})", "V_on (arb. V)",
                xlabel=col_tR, ylabel=col_epsR,
                fname="decomp_von_epsR_tR_raw.png", vmin=vmin, vmax=vmax)
        (r_mean, c_mean, Inter, High) = res["parts"]
        heatmap(Inter, c_vals, r_vals,
                "Interaction (low-freq) of V_on", "arb.",
                xlabel=col_tR, ylabel=col_epsR,
                fname="decomp_von_epsR_tR_inter.png")
        heatmap(High, c_vals, r_vals,
                "Quantization texture (high-freq) of V_on", "arb.",
                xlabel=col_tR, ylabel=col_epsR,
                fname="decomp_von_epsR_tR_quant.png")
    else:
        print("（提醒）此 CSV 中找不到 V_on 欄位，略過 V_on 分析。")

    if col_rr is not None:
        r_vals, c_vals, Y = pivot_to_grid(df, col_epsR, col_tR, col_rr)
        res = decompose_grid(Y, smooth_k=smooth_k)
        print_report(f"RR vs ({col_epsR}, {col_tR}) [k={smooth_k}]", res)
        rows_for_csv.append({
            "which": "RR", "row": col_epsR, "col": col_tR,
            "row_main_pct": res["row_main_pct"], "col_main_pct": res["col_main_pct"],
            "interaction_pct": res["interaction_pct"], "quantization_pct": res["quantization_pct"],
            "remainder_pct": res["remainder_pct"], "SS_tot": res["tot"]
        })
        vmin, vmax = np.nanmin(np.log10(np.maximum(Y, 1e-12))), np.nanmax(np.log10(np.maximum(Y, 1e-12)))
        heatmap(Y, c_vals, r_vals,
                f"log10(RR) vs ({col_epsR}, {col_tR})", "RR",
                xlabel=col_tR, ylabel=col_epsR,
                fname="decomp_rr_epsR_tR_raw.png", log10=True, vmin=vmin, vmax=vmax)
        (r_mean, c_mean, Inter, High) = res["parts"]
        heatmap(Inter, c_vals, r_vals,
                "Interaction (low-freq) of RR", "arb.",
                xlabel=col_tR, ylabel=col_epsR,
                fname="decomp_rr_epsR_tR_inter.png")
        heatmap(High, c_vals, r_vals,
                "Quantization texture (high-freq) of RR", "arb.",
                xlabel=col_tR, ylabel=col_epsR,
                fname="decomp_rr_epsR_tR_quant.png")
    else:
        print("（提醒）此 CSV 中找不到 RR 欄位，略過 RR 分析。")

    if rows_for_csv:
        append_variance_csv(rows_for_csv)

# ============ 入口 ============

if __name__ == "__main__":
    # k 越大 → 當作更強的「去相干/平均化」，高頻紋理比例會下降
    SMOOTH_K = 3

    try:
        analyze_NR_alpha("sweep_NR_alpha.csv", smooth_k=SMOOTH_K)
    except Exception as e:
        print("[NR, alpha] 分析略過：", e)

    try:
        analyze_epsR_tR("sweep_epsR_tR.csv", smooth_k=SMOOTH_K)
    except Exception as e:
        print("[epsR, tR] 分析略過：", e)

    print("\n完成。輸出 PNG 會存到目前資料夾；比例摘要彙整於 variance_report.csv。")
