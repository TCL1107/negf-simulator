# NEGF Quantum Transport Simulator

*Python‑based quantum transport toolkit using the Non-Equilibrium Green’s Function (NEGF) formalism.*  
Supports 1D chains, graphene‑like ribbons, and CNT‑heterojunction diodes, enabling transparent, customizable, research‑ready transport simulations.

---

## Motivation & Purpose  

Nano‑electronic devices — graphene ribbons, carbon nanotubes (CNTs), and engineered heterostructures — rely critically on quantum‑coherent electron transport.  
However, existing NEGF tools are often:

- computationally heavy or tightly coupled to DFT packages,  
- opaque black‑boxes that are hard to modify,  
- or too complex for rapid prototyping and conceptual studies.

**This project fills that gap** by providing a clean, minimal, and fully transparent NEGF pipeline.  
Every step — from Hamiltonian construction, surface Green’s function, contact self‑energy, to transmission — is visible and editable, making it easy to connect *mathematical expressions → numerical behavior → device‑level intuition*.

Use this toolkit to:

- prototype device concepts quickly,  
- test physics hypotheses,  
- generate transport datasets,  
- and build intuition before heavier simulations.

---

## Key Features  

- **Analytic surface Green’s function** for semi‑infinite 1D tight‑binding leads  
- **Custom heterojunction construction** (metal/semiconductor segments, tunable onsite energies, tunable hoppings)  
- **Landauer–Büttiker transport calculations**  
  - Transmission spectrum **T(E)**  
  - Current–voltage (I–V) curves with asymmetric Fermi levels  
- **Parameter sweeps & variance analysis** for extracting device design rules  
- **Exploration tools** for:  
  - Fabry–Pérot oscillations  
  - Kronig–Penney periodic potentials  
  - Geometry‑dependent conduction in ribbon‑like systems  

---

## Repository Structure  

```
src/negf/         # Core NEGF modules (Green’s function, self‑energy, device builder)
examples/         # Toy chains, CNT diode, GNR transport, Fabry–Pérot scripts
figures/          # Generated plots (T(E), I–V curves, sweep maps)
data/             # Exported CSV/numpy datasets
notebooks/        # (Optional) workflow demos and tutorials
```

---

## Getting Started  

### 1. Install dependencies  

```bash
pip install -r requirements.txt
pip install -e .
```

### 2. Run a demonstration example  

```bash
python examples/cnt_diode.py
```

This script computes quantum transport through a CNT‑like heterojunction and outputs:

- `TE_spectrum.png`, `TE_spectrum.csv`  
- `IV_curve.png`, `IV_curve.csv`  

---

## Minimal Working Example (MWE)

```python
from negf.negf import transmission
from negf.surface_gf import surface_gf_1d
import numpy as np

# Simple 1‑D chain model
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

This demonstrates the complete workflow:  
**Hamiltonian → Contact self‑energies → Retarded Green’s function → Transmission.**

---

## What You Can Learn / Explore  

This simulator allows rapid exploration of:

- Band‑edge behavior from onsite shifts  
- Fabry–Pérot resonance structures in multi‑cell devices  
- Contact‑coupling effects on absolute current scale  
- Schottky‑like rectification in heterojunction devices  
- Dephasing‑induced resonance broadening  

These insights are essential for understanding nanoscale quantum transport.

---

## Background & Research Context  

This toolkit originated from an undergraduate research project on quantum transport and rectification in CNT‑like systems.  
It has been reorganized into a clean and reusable form to support future device studies, reproducible simulations, and educational demonstrations.

---

## License  

Released under the **MIT License**. Contributions and feature requests are welcome.
