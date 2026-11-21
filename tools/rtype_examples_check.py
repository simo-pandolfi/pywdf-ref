# pywdf-ref: checks of the corrections to pywdf (E1, E2, E3), see tools/README.md
# Copyright (c) 2026 Simone Pandolfi
# SPDX-License-Identifier: MIT

# Checks the pywdf R-type examples against the E1/E2 corrections (rswdf VALIDATION.md).
# Usage: python -I rtype_examples_check.py <path to the pywdf package dir that contains core/>
#   1. |S[up][up]| of every non-root RTypeAdaptor in the examples (E2 requires 0).
#   2. BaxandallEQ (RTypeAdaptor under IdealVoltageSource) vs UnadaptedBaxandallEQ
#      (RootRTypeAdaptor): same circuit, outputs must agree. Expected: ~1e-1 on
#      pywdf 3a1dac9, ~1e-14 with E1/E2.

import os, sys, importlib, numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")   # examples import matplotlib
sys.path.insert(0, sys.argv[1])               # pywdf package dir (contains core/)
from core.wdf import PolarityInverter
from core.rtype import RTypeAdaptor, RootRTypeAdaptor
fs = 48000
cases = [
    ("examples.bassmantonestack", "BassmanToneStack", (fs, .5, .5, .5)),
    ("examples.baxandalleq", "BaxandallEQ", (fs, .5, .5)),
    ("examples.sallenkeyfilter", "SallenKeyFilter", (fs, 1000, .707)),
    ("examples.sallenkeyfilter_hpf", "SallenKeyFilter", (fs, 1000, .707)),
    ("examples.passive_apf", "PassiveAPF", (fs, 1000)),
]
def rtypes(obj):
    for name, v in vars(obj).items():
        if isinstance(v, PolarityInverter):   # passive_apf wraps its adaptor
            v = v.p1
        if isinstance(v, RTypeAdaptor) and not isinstance(v, RootRTypeAdaptor):
            yield name, v
print("1. |S[up][up]| of non-root RTypeAdaptor")
for mod, cls, args in cases:
    c = getattr(importlib.import_module(mod), cls)(*args)
    found = False
    for name, ad in rtypes(c):
        u = ad.up_port_idx; found = True
        print(f"   {mod:30s} {name:10s} up={u} n={ad.n_ports} |S[up][up]| = {abs(ad.S_matrix[u][u]):.3e}")
    if not found:
        print(f"   {mod:30s} no RTypeAdaptor found as attribute")
n = 4800
x = np.zeros(n); x[0] = 1.0
x[n // 2:] = 0.5 * np.sin(2 * np.pi * 440 * np.arange(n - n // 2) / fs)
bax = importlib.import_module("examples.baxandalleq")
ya = np.array([s for s in map(bax.BaxandallEQ(fs, .5, .5).process_sample, x)])
yr = np.array([s for s in map(bax.UnadaptedBaxandallEQ(fs, .5, .5).process_sample, x)])
print("2. BaxandallEQ vs UnadaptedBaxandallEQ (impulse + 440 Hz tone)")
print(f"   peak {np.max(np.abs(ya)):.3e}   max abs diff {np.max(np.abs(ya - yr)):.2e}   RMS diff {np.sqrt(np.mean((ya - yr)**2)):.2e}")
