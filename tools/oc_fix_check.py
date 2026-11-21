# pywdf-ref: checks of the corrections to pywdf (E1, E2, E3), see tools/README.md
# Copyright (c) 2026 Simone Pandolfi
# SPDX-License-Identifier: MIT

# Candidate OpenCircuit fix (Rp = 1e15, G = 1/Rp), patched in memory only: the
# check used to choose E3 before 7cdf940; open_circuit_check.py tests the fix.
# Parallel: OC || C in an RC low-pass must equal the plain RC.
# Series: Vin - R - C - OC: no current, so V_C must stay 0 and V_OC must follow -Vin (series sign convention).
# Usage: python -I oc_fix_check.py <pywdf package dir> <Rp>
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
from core.wdf import *
RP = float(sys.argv[2])
def ci(self): self.Rp = RP; self.G = 1.0 / RP
OpenCircuit.calc_impedance = ci
fs, R, C, n = 48000.0, 1e3, 1e-6, 2000
x = np.zeros(n); x[0] = 1.0; x[n//2:] = np.sin(2*np.pi*440*np.arange(n - n//2)/fs)
def drive(v, s, probes):
    y = np.zeros((n, len(probes)))
    for k in range(n):
        v.set_voltage(x[k]); v.accept_incident_wave(s.propagate_reflected_wave())
        s.accept_incident_wave(v.propagate_reflected_wave())
        y[k] = [p.wave_to_voltage() for p in probes]
    return y
def rc(with_oc):
    c = Capacitor(C, fs); node = ParallelAdaptor(c, OpenCircuit()) if with_oc else c
    s = SeriesAdaptor(Resistor(R), node); return drive(IdealVoltageSource(s), s, [c])[:, 0]
a, b = rc(False), rc(True)
print(f"Rp={RP:g} parallel: maxdiff vs plain RC {np.max(abs(a-b)):.2e}")
c, oc = Capacitor(C, fs), OpenCircuit()
s = SeriesAdaptor(Resistor(R), SeriesAdaptor(c, oc))
y = drive(IdealVoltageSource(s), s, [c, oc])
print(f"Rp={RP:g} series:   max|V_C| {np.max(abs(y[:,0])):.2e}   max|V_OC + Vin| {np.max(abs(y[:,1] + x)):.2e}  max|V_OC - Vin| {np.max(abs(y[:,1] - x)):.2e}")
