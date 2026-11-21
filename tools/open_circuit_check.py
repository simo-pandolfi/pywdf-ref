# pywdf-ref: checks of the corrections to pywdf (E1, E2, E3), see tools/README.md
# Copyright (c) 2026 Simone Pandolfi
# SPDX-License-Identifier: MIT

# Checks pywdf OpenCircuit (rswdf VALIDATION.md, D3) in parallel and in series.
# Usage: python -I open_circuit_check.py <path to the pywdf package dir that contains core/>
# Expected: pywdf 3a1dac9 (Rp = 0, G = 1e16) fails both cases; the corrected
# OpenCircuit (Rp = 1e15) gives errors ~1e-14 (parallel) and ~1e-12 (series).
#   parallel: RC low-pass with OpenCircuit || C must equal the plain RC low-pass
#   series:   Vin - R - C - OpenCircuit carries no current: V_C = 0, V_OC = Vin

import sys, numpy as np
sys.path.insert(0, sys.argv[1])          # pywdf package dir (contains core/)
from core.wdf import (Resistor, Capacitor, OpenCircuit, SeriesAdaptor,
                      ParallelAdaptor, IdealVoltageSource)
fs, R, C, n = 48000.0, 1e3, 1e-6, 2000
x = np.zeros(n); x[0] = 1.0                                       # impulse...
x[n // 2:] = np.sin(2 * np.pi * 440 * np.arange(n - n // 2) / fs)  # ...then a 440 Hz tone
def drive(s, probes):
    v = IdealVoltageSource(s)
    y = np.zeros((n, len(probes)))
    for k in range(n):
        v.set_voltage(x[k])
        v.accept_incident_wave(s.propagate_reflected_wave())
        s.accept_incident_wave(v.propagate_reflected_wave())
        y[k] = [p.wave_to_voltage() for p in probes]
    return y
def rc(with_oc):
    c = Capacitor(C, fs)
    node = ParallelAdaptor(c, OpenCircuit()) if with_oc else c
    return drive(SeriesAdaptor(Resistor(R), node), [c])[:, 0]
c, oc = Capacitor(C, fs), OpenCircuit()
ys = drive(SeriesAdaptor(Resistor(R), SeriesAdaptor(c, oc)), [c, oc])
print("OpenCircuit Rp = %g" % oc.Rp)
print("parallel  max|V_C - V_C(plain RC)| %.2e" % np.max(np.abs(rc(True) - rc(False))))
print("series    max|V_C| %.2e   max|V_OC - Vin| %.2e" % (np.max(np.abs(ys[:, 0])), np.max(np.abs(ys[:, 1] - x))))
