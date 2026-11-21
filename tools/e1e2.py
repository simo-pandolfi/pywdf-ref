# pywdf-ref: checks of the corrections to pywdf (E1, E2, E3), see tools/README.md
# Copyright (c) 2026 Simone Pandolfi
# SPDX-License-Identifier: MIT

# Reproduces the R-type evaluation-order errors (E1, E2) of pywdf on an RC low-pass.
# Usage: python e1e2.py <path to the pywdf package dir that contains core/>
# Reference: analytic bilinear response of the RC low-pass (scipy.signal.bilinear).

import sys, numpy as np, scipy.signal as sg
sys.path.insert(0, sys.argv[1])          # pywdf package dir (contains core/)
from core.wdf import Resistor, Capacitor, SeriesAdaptor, IdealVoltageSource, ResistiveVoltageSource
from core.rtype import RTypeAdaptor, RootRTypeAdaptor
fs, R, C, n = 48000.0, 1e3, 1e-6, 2000
x = np.zeros(n); x[0] = 1.0
b, a = sg.bilinear([1.0], [R*C, 1.0], fs); ref = sg.lfilter(b, a, x)
def run(root, src, out):
    y = np.zeros(n)
    for k in range(n):
        src.set_voltage(x[k]); root()
        y[k] = out.wave_to_voltage()
    return y
# A: binary series
r, c = Resistor(R), Capacitor(C, fs); s = SeriesAdaptor(r, c); v = IdealVoltageSource(s)
yA = run(lambda: (v.accept_incident_wave(s.propagate_reflected_wave()), s.accept_incident_wave(v.propagate_reflected_wave())), v, c)
# B: RTypeAdaptor series junction, up port 0
r, c = Resistor(R), Capacitor(C, fs)
def imp(ad):
    rb, rc = ad.get_port_impedances(); ra = rb + rc; rv = np.array([ra, rb, rc])
    ad.set_S_matrix(np.eye(3) - 2*np.outer(rv, np.ones(3))/rv.sum()); return ra
rt = RTypeAdaptor([r, c], imp, 0); v = IdealVoltageSource(rt)
yB = run(lambda: (v.accept_incident_wave(rt.propagate_reflected_wave()), rt.accept_incident_wave(v.propagate_reflected_wave())), v, c)
# C: RootRTypeAdaptor parallel: resistive source || C
vs, c = ResistiveVoltageSource(R), Capacitor(C, fs)
def impr(ad):
    g = 1/np.array(ad.get_port_impedances()); ad.set_S_matrix(2*np.outer(np.ones(2), g)/g.sum() - np.eye(2))
root = RootRTypeAdaptor([vs, c], impr)
yC = run(root.compute, vs, c)
e = lambda y: np.max(np.abs(y - ref))
# The series junction with an ideal source at port 0 gives -Vc (sign convention)
yA, yB = -yA, -yB
print("A binary series     maxerr %.2e" % e(yA))
print("B RTypeAdaptor      maxerr %.2e" % e(yB))
print("C RootRTypeAdaptor  maxerr %.2e   shifted by 1: %.2e" % (e(yC), np.max(np.abs(yC[1:] - ref[:-1]))))
