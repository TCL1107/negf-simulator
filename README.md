# NEGF Quantum Transport Simulator

This repository provides a transparent, research-oriented Python toolkit for simulating quantum electron transport using the Non-Equilibrium Green’s Function (NEGF) formalism.  
It supports studies on **1D tight-binding chains, graphene-like ribbons, and CNT-like heterojunctions**, focusing on physical interpretability and device-level insight.

## Motivation & Purpose

Modern nanoelectronic devices—graphene nanoribbons, CNT diodes, and engineered heterostructures—
are all governed by quantum-coherent electron transport.  
Yet most NEGF tools are:

- too complex for quick prototyping,
- too opaque to modify,
- or too heavy for conceptual studies.

**This project aims to bridge that gap** by providing a clean, minimal, and fully transparent NEGF pipeline.  
Every step—from Hamiltonian construction, surface Green’s function, self-energies, and transmission—is visible and editable, so users can directly connect *mathematical expressions → numerical behavior*.

---

## Key Features

- **Analytic surface Green’s function** for semi-infinite 1D leads  
- **Custom heterojunction construction** (metal/semiconductor segments, potential steps)
- **Landauer–Büttiker transport calculation**:
  - Transmission spectrum **T(E)**
  - Current–voltage (I–V) curves with asymmetric contacts
- **Parameter sweeps & variance decomposition** for device design rules
- **Exploration tools** for:
  - Fabry–Pérot interference  
  - Kronig–Penney–like structures  
  - Geometry-dependent conduction in ribbon-like systems

---

## Repository Structure

```
src/negf/         # Core NEGF modules (GF, self-energy, device assembly)
examples/         # Toy chains, CNT diodes, GNR transport, Fabry–Pérot
figures/          # Generated plots (T(E), IV curves, sweep maps)
data/             # CSV/numpy outputs
notebooks/        # (optional) workflow demos
```

---

## Getting Started

### Install dependencies

```bash
pip install -r requirements.txt
pip install -e .
```

### Run a CNT-like diode example

```bash
python examples/cnt_diode.py
```

This produces:

- `TE_spectrum.png`, `TE_spectrum.csv`
- `IV_curve.png`, `IV_curve.csv`

in the working directory.

---

## Minimal Working Example (MWE)

```python
from negf.negf import transmission
from negf.surface_gf import surface_gf_1d
import numpy as np

# Simple 1D chain
N = 5
t = -3.0
eps0 = 0.0

Hc = np.zeros((N, N), dtype=complex)
for i in range(N):
    Hc[i, i] = eps0
    if i < N - 1:
        Hc[i, i+1] = Hc[i+1, i] = t

E = np.linspace(-6, 6, 400)
SigmaL = surface_gf_1d(E, eps0, t)
SigmaR = surface_gf_1d(E, eps0, t)

T = transmission(E, Hc, SigmaL, SigmaR)
```

---

## What You Can Learn

With small code modifications, users can explore:

- band-edge shifts from onsite potentials  
- Fabry–Pérot oscillations from multi-cell interference  
- effect of coupling strength on absolute current  
- rectification mechanisms in heterojunction structures  
- dephasing-induced broadening of resonances  

These tools provide intuition essential for nanoelectronics research.

---

## Background

This toolkit originated from an undergraduate research project on quantum transport and rectification in CNT-like systems.  
It has since been reorganized into a reusable, transparent, and research-ready form for future device studies.
