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
2. **phi(1020) -> K+K-** — the invariant mass of all opposite-sign good-track
   pairs under the kaon-mass hypothesis. The converted trees carry no dE/dx, so
   there is no kaon identification; the phi sits on a combinatorial
   background. Each spectrum is fitted with a Voigtian (Breit-Wigner with the
   PDG width, convolved with a Gaussian resolution) on top of a threshold
   background `b0 (m - 2mK)^a exp(c (m - 2mK))`. Like-sign pairs are shown for
   reference. The fit returns the phi yield per event, the fitted mass and the
   mass resolution for data and MC, and the summary prints their ratios.

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

phi candidates: every opposite-sign pair of good tracks; `--kaonPMin` and
`--phiPMin` can be used to suppress the combinatorial background.

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
