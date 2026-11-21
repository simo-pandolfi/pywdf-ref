# pywdf-ref: checks of the corrections to pywdf (E1, E2, E3), see tools/README.md
# Copyright (c) 2026 Simone Pandolfi
# SPDX-License-Identifier: MIT

# Runs the R-type pywdf examples on a test signal and saves the outputs.
# Usage: python -I smoke_examples.py <pywdf package dir> <out.npz>
import sys, importlib, numpy as np
sys.path.insert(0, sys.argv[1])
fs, n = 48000, 4800
t = np.arange(n) / fs
x = np.zeros(n); x[0] = 1.0
x[n // 2:] = 0.5 * np.sin(2 * np.pi * 440 * t[n // 2:])  # impulse, then a 440 Hz tone
cases = [
    ("bassman", "examples.bassmantonestack", "BassmanToneStack", (fs, .5, .5, .5)),
    ("baxandall", "examples.baxandalleq", "BaxandallEQ", (fs, .5, .5)),
    ("baxandall_root", "examples.baxandalleq", "UnadaptedBaxandallEQ", (fs, .5, .5)),
    ("sallenkey", "examples.sallenkeyfilter", "SallenKeyFilter", (fs, 1000, .707)),
    ("sallenkey_hpf", "examples.sallenkeyfilter_hpf", "SallenKeyFilter", (fs, 1000, .707)),
    ("passive_apf", "examples.passive_apf", "PassiveAPF", (fs, 1000)),
    ("tr808_hat", "examples.tr_808_hatresonator", "TR_808_HatResonator", (fs, 1000, .5)),
]
out = {}
for key, mod, cls, args in cases:
    try:
        c = getattr(importlib.import_module(mod), cls)(*args)
        out[key] = np.array([c.process_sample(v) for v in x])
    except Exception as e:
        print(f"{key}: {type(e).__name__}: {e}")
np.savez(sys.argv[2], **out)
print("ok:", ", ".join(out))
