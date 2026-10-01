# OPAL 1994 data / MC validation

`CompareDataMC` is a compiled ROOT program that compares the converted OPAL 1994
collision data with the reconstructed JETSET 7.4 + GOPAL Monte Carlo
(`jt74mh/r2792`). It is written to be read by students: the source file starts
with a plain-language description of the physics and of the input files, and
is divided into eight commented parts (constants, options, histograms, track
selection, event loop, phi fit, drawing, `main`). Read `CompareDataMC.cpp`
from the top; every ROOT class used is named in the include list with a
one-line explanation.

## What it compares

1. **Charged-particle pT spectra** — `(1/N_ev) dN/dpT` with linear and
   logarithmic binning, data vs reconstructed MC with a data/MC ratio panel,
   and the generator-level MC (`tgen`) overlaid to show the size of the
   detector effect. The good-track multiplicity is drawn as well because it
   sets the per-event normalisation.
2. **phi(1020) -> K+K-** — the invariant mass of opposite-sign good-track
   pairs under the kaon-mass hypothesis. The converted trees carry no dE/dx, so
   there is no kaon identification; the phi sits on a large combinatorial
   background. Each spectrum is fitted with a Voigtian (Breit-Wigner with the
   PDG width, convolved with a Gaussian resolution) on top of a background
   built from the like-sign pairs of the same sample (same combinatorics, no
   phi) multiplied by a quadratic polynomial. The fit returns the phi yield per
   event, the fitted mass and the mass resolution for data and MC, and the
   summary prints their ratios. By default only pairs with momentum above
   6 GeV are used (`--phiPMin`), which raises the signal-to-background ratio
   near the peak from about 1% to about 7%; `--phiPMin 0` gives the inclusive
   spectrum.

## Inputs

Converted files on grendel01:

```
/raid5/data/yjlee/OPAL/converted/1994/merged/OPAL_1994_data.root        tree t
/raid5/data/yjlee/OPAL/converted/1994/merged/OPAL_1994_mc_jt74mh.root   trees t, tgen, tgenBefore
```

## Build and run

```bash
make
bash run_grendel01.sh output                     # full statistics
bash run_grendel01.sh test --maxEvents 100000    # quick look
./CompareDataMC --help                           # all options
```

Outputs in the chosen directory:

| File | Content |
|---|---|
| `pt_spectrum_linear.pdf`, `pt_spectrum_log.pdf` | pT spectra with data/MC ratio |
| `ngood_tracks.pdf` | good-track multiplicity (fraction of events) |
| `phi_kk_mass.pdf` | K+K- mass with fits, data (left) and MC (right) |
| `phi_kk_subtracted.pdf` | background-subtracted phi peak per event, data vs MC |
| `summary.txt` | event counts, fit results, MC/data ratios |
| `OPAL_1994_DataMC_validation.root` | all histograms, ratios and fit functions |

## Selection

Event: `mtError == 0`, `passesOPALHIMSTrack` (OPAL's standard "at least five
good tracks" hadronic-event criterion — already built into the MC trees by the
converter, applied here to data as well), and at least `--nGoodTracksMin` (5)
good charged tracks.

Good charged track: `pwflag == 0` (charged central track) and `highPurity`
(the converter's HIMS definition: pT >= 0.1 GeV, |cos theta| <= 0.966,
|d0| <= 2.5 cm, |z0| <= 50 cm, >= 20 jet-chamber hits). The cuts can be
tightened with `--trackPtMin`, `--trackPMin`, `--cosThetaMax`, `--d0Max`,
`--z0Max`. At generator level `highPurity` is the kinematic acceptance
pT >= 0.1 GeV, |cos theta| <= 0.966, and the d0/z0 cuts do not apply.

phi candidates: every opposite-sign pair of good tracks with pair momentum
above `--phiPMin` (6 GeV); `--kaonPMin` additionally cuts on each track.
The fit range starts at `--fitLow` (0.995 GeV), a few MeV above the K+K-
threshold where like-sign and opposite-sign pairs differ most.

## Results with the full 1994 sample

`bash run_grendel01.sh output` on 1,890,842 data and 651,015 MC events
(1,861,295 and 649,881 selected):

| | OPAL 1994 data | JETSET 7.4 + GOPAL | MC / data |
|---|---|---|---|
| good charged tracks per event | 17.74 | 18.74 | 1.056 |
| pT spectrum shape | — | agrees within ±5% for 0.2 < pT < 20 GeV | — |
| phi yield per event (pair p > 6 GeV) | (10.98 ± 0.69) × 10⁻³ | (9.09 ± 1.13) × 10⁻³ | 0.83 ± 0.12 |
| fitted phi mass | 1020.14 ± 0.25 MeV | 1019.59 ± 0.51 MeV | shift −0.5 MeV |
| mass resolution sigma | 3.20 ± 0.34 MeV | 3.06 ± 0.61 MeV | 0.96 |
| chi2 / ndf | 116 / 89 | 114 / 89 | |

What to look at: the data/MC ratio of the pT spectrum is flat at about 0.93
between 0.3 and 8 GeV (the 5–7% multiplicity difference), rises above one
below 0.2 GeV, and drops to about 0.6 above 25 GeV; the data have an excess of
events with only 5–8 good tracks (tau pairs and two-photon events, which the
quark-pair MC does not contain); the phi mass in data is 0.7 ± 0.25 MeV above
the PDG value, i.e. the momentum scale agrees with the simulation to better
than one per mille.

## Things to keep in mind

- The MC uses one nominal 1994 condition set (`EXPT 1005`, fixed
  sqrt(s) = 91.208 GeV); the data span runs 5014–5712 with sqrt(s) drifting
  between 91.19 and 91.45 GeV.
- The generator-level overlay uses the `tgen` tree, i.e. the generator record
  of the events that passed the reconstruction-level selection; it illustrates
  the detector effect, it is not an unfolding.
- Expect an offset in the per-event pT normalisation: the JETSET 7.4 sample has
  about 7% more charged particles per event than the data. Use
  `--normalize area` for a shape-only comparison.
- A fit "status" of 0 in `summary.txt` means the minimiser converged; chi2/ndf
  close to 1 means the model describes the histogram within statistics.
