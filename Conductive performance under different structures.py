# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 17:38:16 2025

@author: USER
"""

# ============================================
# 第 0 步：載入必要的套件
# ============================================

# pandas：處理表格資料
# scikit-learn：提供ML模型與資料前處理工具
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score
import matplotlib.pyplot as plt

# ============================================
# 第 1 步：建立訓練資料 (模擬 NanoDCAL 輸出)
# ============================================

# 我們先用一小筆模擬資料取代真實模擬結果。
# 這些數值代表不同結構下的電子穿隧行為。

np.random.seed(42)  # 固定亂數讓結果可重現

# -----------------------------
# 1️⃣ 產生結構參數 (features)
# -----------------------------
n_samples = 100

barrier_thickness_nm = np.random.uniform(0.8, 1.6, n_samples)  # 厚度 0.8~1.6 nm
Ef_T = np.random.uniform(0.3, 0.95, n_samples)                  # 傳輸係數 0.3~0.95
materials = np.random.choice(["MgO", "Al2O3"], n_samples)        # 兩種材料
interfaces = np.random.choice(["Fe-O-Mg", "Fe-Mg-O", "Fe-O-Al", "Fe-Al-O"], n_samples)

# -----------------------------
# 2️⃣ 根據物理邏輯模擬電流 (target)
# -----------------------------
# 基本模型：I ~ T(E) * exp(-α * thickness)
alpha = 2.5  # 代表穿隧阻礙強度
base_current = 20 * Ef_T * np.exp(-alpha * (barrier_thickness_nm - 0.8))

# 根據材料調整比例
for i in range(n_samples):
    if materials[i] == "MgO":
        base_current[i] *= 1.0  # MgO 為基準
    else:
        base_current[i] *= 0.8  # Al2O3 導電性稍弱

# 加入少許隨機雜訊（模擬計算或量測誤差）
noise = np.random.normal(0, 0.8, n_samples)
I_at_0_1V = np.clip(base_current + noise, 0, None)  # μA, 不允許負電流

# -----------------------------
# 3️⃣ 組成 DataFrame
# -----------------------------
df = pd.DataFrame({
    "barrier_thickness_nm": barrier_thickness_nm,
    "barrier_material": materials,
    "interface_atom_order": interfaces,
    "Ef_T(E)": Ef_T,
    "I_at_0.1V_uA": I_at_0_1V
})

print(df.head(10))
print(f"\n✅ 資料筆數: {len(df)}")

# ============================================
# 🧮 第 2 步：定義特徵與目標
# ============================================

# 我們想預測在 0.1 V 下的電流 (I_at_0.1V_uA)
# 特徵（inputs）包含厚度、材料、界面排列、T(E)
X = df[["barrier_thickness_nm", "barrier_material", "interface_atom_order", "Ef_T(E)"]]
y = df["I_at_0.1V_uA"]

# ============================================
# 第 3 步：對文字特徵做 One-Hot 編碼
# ============================================

# 這是把文字轉成機器可理解的數值特徵。
# 例如：MgO → [1,0]；Al2O3 → [0,1]
categorical_features = ["barrier_material", "interface_atom_order"]
preprocessor = ColumnTransformer(
    transformers=[
        ("cat", OneHotEncoder(), categorical_features)
    ],
    remainder="passthrough"  # 其他數值欄位保持原樣
)

# ============================================
# 第 4 步：建立隨機森林迴歸模型
# ============================================

# RandomForest 是一種穩定且解釋性強的非線性模型
model = Pipeline(steps=[
    ("preprocessor", preprocessor),
    ("regressor", RandomForestRegressor(
        n_estimators=200,  # 樹的數量
        random_state=42,
        max_depth=None     # 讓模型自動決定深度
    ))
])

# ============================================
# 第 5 步：切分訓練與測試資料
# ============================================

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# ============================================
# 第 6 步：訓練模型
# ============================================

model.fit(X_train, y_train)

# ============================================
# 第 7 步：模型評估
# ============================================

y_pred = model.predict(X_test)

mae = mean_absolute_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"\n📊 模型評估結果：")
print(f"平均絕對誤差 (MAE): {mae:.3f}")
print(f"決定係數 R²: {r2:.3f}")

# 畫出預測 vs 實際 的散佈圖
plt.figure(figsize=(6,6))
plt.scatter(y_test, y_pred, color='royalblue', s=70)
plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--')
plt.xlabel("real current (μA)")
plt.ylabel("predic current (μA)")
plt.title("predic vs real")
plt.text(min(y_test)+0.1, max(y_test)-0.5, "Slope = 1.0", color='red')
plt.grid(True)
plt.show()

# ============================================
# 第 8 步：用模型預測新結構
# ============================================

# 假設有兩個新結構，想看 AI 覺得導電表現如何
new_samples = pd.DataFrame({
    "barrier_thickness_nm": [1.3, 0.9],
    "barrier_material": ["MgO", "Al2O3"],
    "interface_atom_order": ["Fe-O-Mg", "Fe-O-Al"],
    "Ef_T(E)": [0.45, 0.60]
})

predicted_currents = model.predict(new_samples)
print("\n🔮 新結構預測結果：")
for i, val in enumerate(predicted_currents):
    print(f"結構 {i+1} 的預測電流 ≈ {val:.2f} μA")

# ============================================
# 第 9 步（進階）：查看特徵重要性
# ============================================

# Random Forest 可以告訴我們哪個特徵對模型影響最大
feature_names = model.named_steps["preprocessor"].get_feature_names_out()
importances = model.named_steps["regressor"].feature_importances_

plt.figure(figsize=(8,4))
plt.barh(feature_names, importances, color='teal')
plt.xlabel("Feature Importance")
plt.title("The contribution of different features to conductivity prediction")
plt.show()
