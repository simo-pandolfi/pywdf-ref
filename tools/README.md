# tools: checks of the corrections

Scripts that reproduce the errors of the original pywdf (tag `orig-3a1dac9`) and
check the corrections of the `fixes` branch (E2 `7911766`, E1 `412bb44`, E3
`7cdf940`; `rtype_corrections_check.py` also the example corrections E8-E11).
Each script takes the pywdf package folder (the one that contains `core/`) as
its first argument, so the same script can be run on both versions:

```
python -I tools/e1e2.py <checkout>/pywdf
```

Requirements: numpy, scipy (and matplotlib, imported by the examples).

| Script | What it checks | `orig-3a1dac9` | `fixes` (`7cdf940`) |
|---|---|---|---|
| `e1e2.py` | RC low-pass (R = 1 kOhm, C = 1 uF, fs = 48 kHz), 2000-sample impulse response against the analytic bilinear one: binary series adaptor (A), `RTypeAdaptor` (B, E2), `RootRTypeAdaptor` (C, E1) | A 6.9e-17, B 2.1e-4, C 1.0e-2 | A 6.9e-17, B 7.5e-17, C 9.0e-17 |
| `rtype_examples_check.py` | 1. \|S[up][up]\| of every non-root `RTypeAdaptor` in the examples (E2 requires 0); 2. `BaxandallEQ` against `UnadaptedBaxandallEQ`, the same circuit with the R-type adaptor under the source or at the root | 1. all 0; 2. max diff 1.05e-1 (peak 1.05e-1) | 1. all 0; 2. max diff 1.8e-14 |
| `open_circuit_check.py` | `OpenCircuit` in parallel (OC \|\| C in an RC low-pass equals the plain RC) and in series (R - C - OC: no current, V_OC = Vin), impulse then 440 Hz | parallel 4.6e-1; series \|V_C\| 4.6e-1, \|V_OC - Vin\| 1.0 | parallel 4.9e-15; series 1.0e-14, 2.0e-12 |
| `rtype_corrections_check.py` | Example corrections E8-E11 against a nodal analysis (bilinear transform): 1. original Bassman output (sum of four port voltages) = 2 V(tb) - V(w) - 2 V(x), not the treble wiper V(w); 2. `BaxandallEQ` and `UnadaptedBaxandallEQ` against the nodal analysis; 3. gains of the Baxandall network in four versions (original, E9, E9 + E10, E9-E11), pywdf and nodal; 4. effect of `Ca` = 1 uF (E10) | pywdf columns carry the E1/E2 errors: 1. 5.6e-3 to 2.9e-2; 2. 2.5e-1 / 9.8e-1 | with E1-E13 (`e4f5d25`): 1. 1.1e-14; 2. 4.0e-14 / 1.0e-13; 3. treble at 20 kHz, t = 0.01 / 0.99: -20.28 / -6.11 dB (original), -6.07 / -25.38 dB (E9); b = t = 0.5 at 50 / 200 Hz: -9.73 / -14.06 dB (E9 + E10), -16.01 / -22.63 dB (E9-E11); 4. -0.55 / -0.22 dB at 5 / 20 Hz, impulse response 3.6e-4 (Cc 22 nF) |
| `smoke_examples.py` | runs the R-type examples (and the TR-808 hat resonator) on an impulse and a 440 Hz tone, saves the outputs to a `.npz` | runs | runs |
| `check_sup.py` | part 1 of `rtype_examples_check.py` alone | all 0 | all 0 |
| `oc_fix_check.py` | the candidate E3 (`Rp` given on the command line) patched in memory, used to choose E3; `open_circuit_check.py` checks the committed fix | with `1e15`: as `fixes` | same |

Change of the example outputs between the two versions (`smoke_examples.py`,
max abs difference, peak in brackets): bassman 7.8e-3 (0.84), baxandall
1.3e-2 (0.10), baxandall_root 1.0e-1 (0.10), sallenkey 9.1e-2 (0.025),
sallenkey_hpf 3.1e-3 (0.082), passive_apf 1.0e-2 (0.98), tr808_hat 0.99 (2.6).

Measured on 2026-11-21 with Python 3.13, numpy 2.5, scipy 1.18.
