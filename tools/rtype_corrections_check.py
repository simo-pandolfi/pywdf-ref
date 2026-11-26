# pywdf-ref: checks of the corrections to pywdf (E8-E11), see tools/README.md
# Copyright (c) 2026 Simone Pandolfi
# SPDX-License-Identifier: MIT

# Checks the example corrections E8 (Bassman output) and E9-E11 (Baxandall)
# against a nodal analysis of the circuits, discretised with the bilinear map
# s = 2 fs (z - 1)/(z + 1) (the pywdf capacitor discretisation, so pywdf and
# the nodal analysis agree to rounding once E1/E2 are in).
# Usage: python -I rtype_corrections_check.py <path to the pywdf package dir that contains core/>
#   1. E8: the original Bassman output, the sum of the port voltages
#      v(R3_minus) + v(R1_minus) + v(R3_plus) + v(R2), against
#      2 V(tb) - V(w) - 2 V(x) and against the treble wiper voltage V(w)
#      (nodes: w treble wiper, tb between C2/R1_minus and R2, x above
#      R3_minus), at the resistor values of the checkout's instance.
#   2. Baxandall: BaxandallEQ and UnadaptedBaxandallEQ of the checkout against
#      the nodal analysis, 7 settings, component values read from the instance.
#   3. Gains (dB) of four versions of the Baxandall network: original
#      (set_treble sets Pt_plus = Pt (1 - t) and leaves Pt_minus at 50 kOhm,
#      Ca = 1e6 F, Cc = 22 nF), E9, E9 + E10, E9-E11; computed by pywdf
#      (BaxandallEQ with the components changed in memory) and by the nodal
#      analysis.
#   4. E10: effect of Ca = 1 uF instead of 1e6 F, gain at 1, 5 and 20 Hz and
#      max abs change of the impulse response, with Cc = 22 nF and 220 nF.
# Expected with E1/E2 (branch fixes): parts 1-3 agree to ~1e-13; on the
# original library (3a1dac9) the pywdf columns carry the E1/E2 errors.

import os, sys, importlib, numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")   # examples import matplotlib
sys.path.insert(0, sys.argv[1])               # pywdf package dir (contains core/)
bassman = importlib.import_module("examples.bassmantonestack")
bax = importlib.import_module("examples.baxandalleq")

FS = 44100
GND, VIN = -1, -2


def nodal(n, branches):
    """Node voltages for VIN = 1 V; branches: (node a, node b, admittance)."""
    y = np.zeros((n, n), complex)
    rhs = np.zeros(n, complex)
    for a, b, g in branches:
        for p, q in ((a, b), (b, a)):
            if p >= 0:
                y[p, p] += g
                if q >= 0:
                    y[p, q] -= g
                elif q == VIN:
                    rhs[p] += g
    return np.linalg.solve(y, rhs)


def s_at(f, fs=FS):
    z = np.exp(2j * np.pi * f / fs)
    return 2 * fs * (z - 1) / (z + 1)


def dft(x, f, fs=FS):
    return np.sum(x * np.exp(-2j * np.pi * f / fs * np.arange(len(x))))


def db(h):
    return 20 * np.log10(abs(h))


# 1. E8 ----------------------------------------------------------------------
print("1. E8: original Bassman output (sum of four port voltages) vs nodal analysis")
freqs = [20, 100, 440, 1000, 5000, 15000]
for b, m, t in [(0.5, 0.5, 0.5), (0.01, 0.99, 0.25), (0.99, 0.01, 0.75)]:
    bts = bassman.BassmanToneStack(FS, b, m, t)
    n = 2 * FS
    x = np.zeros(n); x[0] = 1.0
    bts.reset()
    y = np.empty(n)
    for k in range(n):
        bts.process_sample(x[k])
        y[k] = (bts.R3_minus.wave_to_voltage() + bts.R1_minus.wave_to_voltage()
                + bts.R3_plus.wave_to_voltage() + bts.R2.wave_to_voltage())
    e_formula = e_wiper = 0.0
    for f in freqs:
        s = s_at(f)
        # nodes: c 0, w 1, tb 2, n1 3, mt 4, x 5
        v = nodal(6, [
            (VIN, 0, s * bts.C1.C), (0, 1, 1 / bts.R1_plus.Rp), (1, 2, 1 / bts.R1_minus.Rp),
            (VIN, 3, 1 / bts.R4.Rp), (3, 2, s * bts.C2.C), (3, 5, s * bts.C3.C),
            (2, 4, 1 / bts.R2.Rp), (4, 5, 1 / bts.R3_plus.Rp), (5, GND, 1 / bts.R3_minus.Rp),
        ])
        h = dft(y, f)
        e_formula = max(e_formula, abs(h - (2 * v[2] - v[1] - 2 * v[5])))
        e_wiper = max(e_wiper, abs(h - v[1]))
    print(f"   bass {b} mid {m} treble {t}: vs 2V(tb)-V(w)-2V(x) {e_formula:.2e}, vs V(w) {e_wiper:.2e}")


# 2-4. Baxandall --------------------------------------------------------------
def bax_values(eq):
    return dict(ptp=eq.Pt_plus.Rp, ptm=eq.Pt_minus.Rp, pbp=eq.Pb_plus.Rp, pbm=eq.Pb_minus.Rp,
                ca=eq.Ca.C, cb=eq.Cb.C, cc=eq.Cc.C, cd=eq.Cd.C, ce=eq.Ce.C, rl=eq.Rl.Rp,
                resa=eq.Resa.Rp, resb=eq.Resb.Rp, resc=eq.Resc.Rp, resd=eq.Resd.Rp, rese=eq.Rese.Rp)


def bax_nodal(s, v, rs=0.0):
    """Voltage across Rl; nodes: in 0, Resa/Pb+ 1, m 2, Pb-/Resb 3, out 4, Cd/Pt+ 5, Pt-/Ce 6, source/Ca 7."""
    if rs > 0:
        n, src = 8, [(VIN, 7, 1 / rs), (7, 0, s * v["ca"])]
    else:
        n, src = 7, [(VIN, 0, s * v["ca"])]
    return nodal(n, src + [
        (0, 1, 1 / v["resa"]), (1, 2, 1 / v["pbp"]), (1, 2, s * v["cb"]),
        (2, 3, 1 / v["pbm"]), (2, 3, s * v["cc"]), (3, GND, 1 / v["resb"]),
        (2, 4, 1 / v["resc"]),
        (0, 5, s * v["cd"]), (5, 4, 1 / v["ptp"]), (5, 4, 1 / v["resd"]),
        (4, 6, 1 / v["ptm"]), (4, 6, 1 / v["rese"]), (6, GND, s * v["ce"]),
        (4, GND, 1 / v["rl"]),
    ])[4]


N4 = 4 * FS   # Ca = 1 uF: time constant ~0.1 s
print("2. Baxandall: pywdf vs nodal analysis (4 s impulse responses)")
freqs = [5, 20, 50, 200, 1000, 5000, 15000, 20000]
for cls, rs_attr in [("BaxandallEQ", None), ("UnadaptedBaxandallEQ", "Vin")]:
    worst = 0.0
    for b, t in [(.5, .5), (.01, .5), (.99, .5), (.5, .01), (.5, .99), (.2, .8), (.8, .2)]:
        eq = getattr(bax, cls)(FS, b, t)
        rs = getattr(eq, rs_attr).Rp if rs_attr else 0.0
        v = bax_values(eq)
        ir = eq.get_impulse_response(N4 / FS)
        worst = max(worst, max(abs(dft(ir, f) - bax_nodal(s_at(f), v, rs)) for f in freqs))
    print(f"   {cls:22s} 7 settings x 8 frequencies: max abs err {worst:.2e}")

VERSIONS = [  # (label, E9 applied, Ca, Cc)
    ("original", False, 1e6, 22e-9),
    ("E9", True, 1e6, 22e-9),
    ("E9 + E10", True, 1e-6, 22e-9),
    ("E9-E11", True, 1e-6, 220e-9),
]


def make_version(b, t, e9, ca, cc):
    """BaxandallEQ of the checkout with the components of a version."""
    eq = bax.BaxandallEQ(FS, b, t)
    pt = eq.Pt
    eq.Pt_plus.set_resistance(pt * t if e9 else pt * (1 - t))
    eq.Pt_minus.set_resistance(pt * (1 - t) if e9 else pt * .5)
    eq.Ca.set_capacitance(ca)
    eq.Cc.set_capacitance(cc)
    return eq


print("3. Baxandall gains (dB), pywdf / nodal; max |H_pywdf - H_nodal| per version")
cases = [("b=t=.5, 50 Hz", .5, .5, 50), ("b=t=.5, 200 Hz", .5, .5, 200),
         ("b=t=.5, 20 kHz", .5, .5, 20000), ("t=.01, 20 kHz", .5, .01, 20000),
         ("t=.99, 20 kHz", .5, .99, 20000), ("b=.01, 20 Hz", .01, .5, 20),
         ("b=.99, 20 Hz", .99, .5, 20), ("b=t=.5, 20 Hz", .5, .5, 20)]
irs = {}
for label, e9, ca, cc in VERSIONS:
    worst = 0.0
    cells = []
    for name, b, t, f in cases:
        key = (label, b, t)
        if key not in irs:
            eq = make_version(b, t, e9, ca, cc)
            irs[key] = (eq.get_impulse_response(N4 / FS), bax_values(eq))
        ir, v = irs[key]
        hp, hn = dft(ir, f), bax_nodal(s_at(f), v)
        worst = max(worst, abs(hp - hn))
        cells.append(f"{name}: {db(hp):.2f} / {db(hn):.2f}")
    print(f"   {label:9s} max err {worst:.1e}")
    for c in cells:
        print(f"      {c}")

print("4. E10: Ca 1 uF vs 1e6 F (gain change in dB, nodal; impulse response change, pywdf)")
for e9, cc in [(False, 22e-9), (True, 220e-9)]:
    for b in (.5, .01):
        with_ca = make_version(b, .5, e9, 1e-6, cc)
        without = make_version(b, .5, e9, 1e6, cc)
        d = [db(bax_nodal(s_at(f), bax_values(with_ca))) - db(bax_nodal(s_at(f), bax_values(without)))
             for f in (1, 5, 20)]
        ir_change = np.max(np.abs(with_ca.get_impulse_response(N4 / FS) - without.get_impulse_response(N4 / FS)))
        circuit = "corrected (E9-E11)" if e9 else "original, Cc 22 nF"
        print(f"   {circuit:19s} bass {b}: 1 / 5 / 20 Hz {d[0]:.2f} / {d[1]:.2f} / {d[2]:.2f} dB, "
              f"impulse response {ir_change:.1e}")
