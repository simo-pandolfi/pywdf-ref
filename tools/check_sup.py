# pywdf-ref: checks of the corrections to pywdf (E1, E2, E3), see tools/README.md
# Copyright (c) 2026 Simone Pandolfi
# SPDX-License-Identifier: MIT

# Prints |S[up][up]| for every RTypeAdaptor (non-root) in the pywdf examples
# (the condition of E2); also part 1 of rtype_examples_check.py.
# Usage: python -I check_sup.py <pywdf package dir>
import sys, importlib, numpy as np
sys.path.insert(0, sys.argv[1])
from core.rtype import RTypeAdaptor, RootRTypeAdaptor

cases = [
    ("examples.bassmantonestack", "BassmanToneStack", (48000, .5, .5, .5)),
    ("examples.baxandalleq", "BaxandallEQ", (48000, .5, .5)),
    ("examples.sallenkeyfilter", "SallenKeyFilter", (48000, 1000, .707)),
    ("examples.sallenkeyfilter_hpf", "SallenKeyFilter", (48000, 1000, .707)),
    ("examples.passive_apf", "PassiveAPF", (48000, 1000)),
]
for mod, cls, args in cases:
    m = importlib.import_module(mod)
    c = getattr(m, cls)(*args)
    for name, v in vars(c).items():
        if isinstance(v, RTypeAdaptor) and not isinstance(v, RootRTypeAdaptor):
            u = v.up_port_idx
            print(f"{mod:32s} {name:10s} up={u} n={v.n_ports} |S[up][up]|={abs(v.S_matrix[u][u]):.3e}")
