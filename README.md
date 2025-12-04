# 📘 NEGF-Simulator  
*A Materials-Aware Quantum Transport Toolkit for Nanoelectronic Devices*

---

# 🌄 Pipeline Overview (Materials → Hamiltonian → NEGF → Device)

<img src="CNT_pipeline_composite.png" width="100%">

---

## 🚀 Overview

**NEGF-Simulator** is a clean, research-oriented Python toolkit connecting:

- **materials-level physics** (CNT bandstructure, zone-folding)  
- **effective tight-binding modeling** (Δ = Eg/2 chain model)  
- **quantum transport** (NEGF Green’s functions)  
- **device-level signatures** (gap, rectification, resonance)

The core methodology:

```
CNT (n,m) → Eg(n,m)
      → Δ = Eg/2
      → effective Hamiltonian
      → NEGF transport
      → T(E) → device behavior
```

This end-to-end flow mirrors how ECE device-modeling groups perform physics-based simulation.

---

## 🧭 Model Hierarchy

```
┌─────────────────────────────┐
│  Materials (CNT)            │
│  - zone-folding Eg(n,m)     │
└─────────────────────────────┘
               ↓ mapping
┌─────────────────────────────┐
│  Effective Hamiltonian      │
│  - 1D dimerized chain       │
│  - Eg ≈ 2Δ                  │
└─────────────────────────────┘
               ↓ NEGF
┌─────────────────────────────┐
│  Device Behavior            │
│  - T(E)                     │
│  - gap extraction           │
│  - transport resonance      │
└─────────────────────────────┘
```

---

## 1️⃣ CNT Bandstructure via Zone Folding

Scripts:
- `src/negf/cnt_zonefold.py`  
- `examples/cnt_zonefold_bandstructure.py`  
- `examples/cnt_zonefold_fullTB.py`

Outputs include:
- CNT diameter  
- bandgap Eg(n,m)  
- subbands  
- metallic/semiconducting classification  

Example:

```
CNT (17,0)
  d = 1.331 nm
  Eg ≈ 0.619 eV
```

*(see CNT_bandstructure.png for reference)*

---

## 2️⃣ Mapping CNT Bandgap → Effective Chain Hamiltonian

Using CNT’s zone-folding gap:

Eg_CNT → Δ = Eg/2

Construct a 1D dimerized chain:

H_ii = ±Δ,  H_{i,i+1} = t

This preserves the CNT’s band-edge physics while greatly simplifying computation.

---

## 3️⃣ NEGF Transport: T(E), G(E), Σ_L/R

Module:
- `src/negf/negf.py`

Features:
- retarded Green’s function  
- analytic 1D surface self-energies  
- T(E)  
- gap extraction  
- interface effects  
- resonance transport  

Your true T(E) figure is embedded in the composite banner above.

---

## 4️⃣ Device-Level Examples

Folder: `examples/`

Includes:
- CNT bandstructure  
- CNT→Δ→NEGF mapping  
- uniform chain transport  
- heterojunction rectifiers  
- dephasing α study  
- parameter sweeps  

These show how electronic structure translates to device-level ON/OFF behavior,  
turn-on voltage, and rectification ratio.

---

## 🧱 Repository Structure

```
negf-simulator/
│
├── src/negf/
│     ├── negf.py
│     ├── surface_gf.py
│     ├── cnt_zonefold.py
│     └── __init__.py
│
├── examples/
│     ├── cnt_zonefold_bandstructure.py
│     ├── cnt_zonefold_fullTB.py
│     ├── cnt_gap_to_chain.py
│     ├── uniform_chain.py
│     └── rectifier_demo.py
│
└── README.md
```

---

## 🧠 Why This Toolkit Matters for ECE Device Modeling

This project demonstrates essential research skills:

✔ Physics-based modeling  
✔ Effective Hamiltonian construction  
✔ NEGF implementation from scratch  
✔ Model validation (Eg_CNT ≈ Eg_chain)  
✔ Modular, reproducible architecture  

These abilities are exactly what ECE device-physics advisors expect.

---

## 🔧 Future Extensions

- Poisson-NEGF self-consistent I–V  
- Electron–phonon scattering  
- Multi-orbital CNT TB model  
- Graphene nanoribbon transport  
- Metal/CNT Schottky barrier extraction  
- ML-based parameter tuning  
- Pareto optimization (RR, V_on, energy cost)

---

## 👤 Maintainer

**Te-Chang Liu**  
Quantum Transport • Nanoelectronics • Device Modeling  
GitHub: https://github.com/TCL1107/negf-simulator
