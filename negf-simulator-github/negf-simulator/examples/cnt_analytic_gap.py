"""
cnt_analytic_gap.py

Analytic CNT bandgap vs diameter using a simple
zone-folding model:

    Eg ≈ 2 * gamma0 * d_CC / d

This script is independent of NEGF; it just uses geometry.
"""

import sys
import os
import numpy as np

# Make sure we can import the local "negf" package
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(ROOT, "src"))

from negf.cnt_geometry import cnt_diameter, cnt_bandgap_nm, GAMMA0, D_CC


def main():
    # A few example chiralities (you可以自己改成你論文中常見的那幾顆)
    chiralities = [
        (17, 0),
        (19, 0),
        (13, 6),
        (10, 9),
        (10, 10),  # metallic (Eg ≈ 0) just to show the formula
    ]

    print("Analytic CNT bandgaps (zone-folding)")
    print("gamma0 = {:.3f} eV, d_CC = {:.3f} nm".format(GAMMA0, D_CC))
    print("-" * 60)
    print("{:>8s} {:>10s} {:>10s}".format("(n,m)", "d (nm)", "Eg (eV)"))
    print("-" * 60)

    for (n, m) in chiralities:
        d = cnt_diameter(n, m)
        Eg = cnt_bandgap_nm(n, m)
        print("({:2d},{:2d}) {:10.3f} {:10.3f}".format(n, m, d, Eg))

    print("-" * 60)
    print("Note: metallic tubes should have Eg ≈ 0; this model mainly")
    print("captures the 1/d scaling trend for semiconducting tubes.")


if __name__ == "__main__":
    main()
