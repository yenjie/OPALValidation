# The OPAL 1994 ntuples: every branch explained

This page documents the ROOT files used by `CompareDataMC`:

```
/raid5/data/yjlee/OPAL/converted/1994/merged/OPAL_1994_data.root        1,890,842 events, tree t
/raid5/data/yjlee/OPAL/converted/1994/merged/OPAL_1994_mc_jt74mh.root     651,015 events, trees t, tgen; 653,749 events, tree tgenBefore
```

Each tree has one entry per e+e- collision event. An entry holds a set of
**scalar** branches (one number per event: run number, flags, missing
momentum, ...) and a set of **per-particle** branches (one array per event
with one element per particle: `px`, `py`, `pz`, `charge`, ...). The three
trees have exactly the same 84 branches, so the same analysis code can run
on any of them:

| Tree | Present in | Contents |
|---|---|---|
| `t` | data and MC | the reconstructed event: charged tracks and neutral calorimeter energy as the detector saw them |
| `tgen` | MC only | the generator's final-state particles for the *same events* as `t`, before any detector effect |
| `tgenBefore` | MC only | generator particles for *all* simulated events, including the 0.4% that fail the reconstruction-level selection applied to `t`/`tgen` |

The branch names follow the ParticleTree interface of the MIT ALEPH analysis
code, so that the same readers work for both experiments. Some names are
therefore ALEPH detector words (`ntpc`, `nitc`, `nvdet`, `bFlag`, ...) that
either carry the closest OPAL quantity or a placeholder; this is spelled out
below. The numbers quoted ("observed ...") were measured on the first 50,000
events of each tree with `inspect_ntuple.py`.

## Where the numbers come from

The chain is: OPAL DST on CERN EOS → the ROPE user processor `opal2csv`
(reads the OPAL `OD` banks, runs OPAL's standard **MT** package to build
energy-flow particles) → CSV tables (`events`, `particles`, `tracks`,
`ecal_clusters`, `hcal_clusters`, `mc_particles`) → `opal-convert`
(`src/opal_open_data/aleph.py`) → these ROOT files.

**MT** ("matching") is the OPAL algorithm that combines tracks and calorimeter
clusters into a list of particles without double counting: every good charged
track becomes one particle, the calorimeter energy that the track is expected
to deposit is subtracted from the clusters it points to, and the remaining
cluster energy becomes neutral particles. MT assigns the **pion mass to every
charged particle and zero mass to every neutral one**.

## Reading the trees

Per-particle branches are variable-length arrays. In ROOT, `TTreeReader`
handles them; in Python, `uproot` returns awkward arrays:

```cpp
TTreeReader reader("t", file);
TTreeReaderValue<Int_t>   nParticle(reader, "nParticle");
TTreeReaderArray<Float_t> pt(reader, "pt");
TTreeReaderArray<Short_t> charge(reader, "charge");
while (reader.Next())
  for (unsigned i = 0; i < pt.GetSize(); ++i)
    if (charge[i] != 0) h->Fill(pt[i]);
```

```python
import uproot
t = uproot.open("OPAL_1994_data.root")["t"]
arrays = t.arrays(["pt", "charge"], entry_stop=10000)
```

Every per-particle branch `X` comes with a counter branch `nX` (`npx`, `npy`,
...). They are an artefact of how uproot writes variable-length arrays; all of
them equal `nParticle` and can be ignored.

## Per-event (scalar) branches

Types: `int32`/`int64` integers, `float32` real numbers, `bool` true/false.
"ALEPH placeholder" means the branch exists only so that ALEPH readers work;
its value carries no OPAL information.

| Branch | Type | Meaning | Data | MC |
|---|---|---|---|---|
| `EventNo` | int32 | OPAL event number within the run | 25 ... 182354 (first chunk) | 245001 ... per simulated partition |
| `RunNo` | int32 | OPAL run number. Data runs 5014–5712 cover periods p59–p68 of 1994. The MC carries the assigned production run number 2792 for every event; it is a bookkeeping number, not a 1991 data run | 5014–5712 | 2792 |
| `year` | int32 | data-taking year | 1994 | 1994 |
| `subDir` | int32 | ALEPH placeholder (ALEPH sub-directory code) | −999 | −1 |
| `process` | int32 | ALEPH placeholder (ALEPH process code) | −1 | −1 |
| `isMC` | bool | true for simulation | false | true |
| `uniqueID` | uint64 | one number that identifies the event across all files: `(sourceID << 48) | (RunNo << 24) | EventNo` | | |
| `Energy` | float32 | centre-of-mass energy √s in GeV, twice the LEP beam energy recorded for the run. Varies run by run in data (91.19–91.45 GeV, mean 91.22) because LEP's energy drifted and period p68 included off-peak scan points | 91.212, 91.214, ... | 91.208 or 91.210 |
| `bFlag` | int32 | ALEPH placeholder (b-tag flag) | −1 | −999 |
| `particleWeight` | float32 | ALEPH placeholder (event weight); always 1 | 1 | 1 |
| `bx`, `by`, `ebx`, `eby` | float32 | ALEPH placeholder (beam-spot position and its error). The real OPAL beam spot is in the `SUBEAM` conditions on CVMFS, not here | 0 | −999 |
| `nParticle` | int32 | number of entries in the per-particle arrays (charged + neutral). `t`: 0–95 observed (a handful of data events have none); `tgen`: 9–122 | | |
| `passesNTupleAfterCut` | bool | true if the converter kept the event (see "Event selection" below). Always true in `t` and `tgen`; in `tgenBefore` false marks the events that `t` does not contain | true | true / mixed in `tgenBefore` |
| `passesNTrkMin` | bool | same as `passesOPALHIMSTrack` | mixed | true |
| `passesLEP1TwoPC` | bool | `nParticle >= 2` (ALEPH two-particle requirement) | mixed | true |
| `passesAll` | bool | same as `passesNTupleAfterCut` | true | true |
| `passesOPALHIMSTrack` | bool | the "track arm" of OPAL's standard hadronic-event selection (HIMS): at least five charged tracks passing the good-track cuts listed under `highPurity`. The published HIMS lists also accept events through a calorimeter-based arm (GPMH bit) that is not stored in the ntuple; 0.9% of the data events fail this track arm but are in the official list. False for those 0.9% of data; true for every MC event in `t`/`tgen` because the converter applied this cut to the MC | 99.1% true | true |
| `missP` | float32 | magnitude of the missing momentum vector (GeV), from the OPAL event summary: minus the vector sum of all reconstructed particle momenta | mean 15 GeV; rare unphysical outliers (> 1000 GeV) exist | |
| `missPt` | float32 | transverse component of the missing momentum (GeV) | | |
| `missTheta`, `missPhi` | float32 | polar and azimuthal angle of the missing momentum vector (rad) | | |
| `nChargedHadrons` | int32 | number of charged particles in the arrays (`charge != 0`) | 0–53 | 2–56 |
| `nChargedHadronsHP` | int32 | number of charged particles with `highPurity` true; equal to `nChargedHadrons` in `t` (see `highPurity`) | | |
| `nChargedHadronsHP_Corrected` | float32 | ALEPH placeholder (efficiency-corrected count); a copy of `nChargedHadronsHP` | | |
| `sourceID` | int32 | which extraction job produced the event. Data: periods converted in one job carry the period number (60, 66, 68 for p60, p66, p68); periods converted in chunks carry `period × 100 + chunk` (5900 = p59 chunk 0, ..., 6711 = p67 chunk 11); 68 values in total. MC: the DDST partition number, 197–800, 523 values | 60, 66, 68, 5900–6711 | 197–800 |
| `eventID` | int64 | running index of the event inside its extraction job, starting at 0; together with `sourceID` it is the key that links all converter tables | | |
| `mtError` | int32 | error code returned by MT; 0 = the particle list was built without problems. Always 0 in `t`, `tgen` and `tgenBefore` (events with an MT error are dropped by the converter) | 0 | 0 |

## Per-particle branches

In `t` a "particle" is an MT energy-flow object: a charged track or a neutral
calorimeter deposit. In `tgen`/`tgenBefore` it is a generator-level
final-state particle: a JETSET particle that was produced as a primary or in
a generator decay and did not itself decay inside the generator. Charged
objects are stored first in every event, then the neutral ones.

### Kinematics (GeV, rad, cm)

| Branch | Type | `t` (reconstructed) | `tgen`, `tgenBefore` (generator) |
|---|---|---|---|
| `px`, `py`, `pz` | float32 | momentum components. z is the beam axis (electron direction), x points to the centre of the LEP ring, y upwards | same, true generator momentum |
| `pt` | float32 | transverse momentum √(px²+py²). Charged tracks: ≥ 0.12 GeV (an MT quality cut); neutrals can be arbitrarily soft | ≥ 0 |
| `pmag` | float32 | magnitude of the momentum | |
| `theta` | float32 | polar angle to the +z beam axis, 0–π. Tracks: 0.17–2.97 (the chamber acceptance) | 0–π, particles can go straight down the beam pipe |
| `phi` | float32 | azimuthal angle, −π to π | |
| `eta` | float32 | pseudorapidity −ln tan(θ/2) | up to ±14 for particles along the beam |
| `rap` | float32 | rapidity ½ ln((E+pz)/(E−pz)) using `particleEnergy`, which depends on the assigned mass | with the true mass |
| `mass` | float32 | the mass MT assigned: 0.1396 GeV (pion) for every charged object, 0 for every neutral one. This is **not** a measurement; use the kaon or proton mass yourself when you need another hypothesis | the true mass: 0 (γ, ν), 0.1396 (π), 0.4937 (K±), 0.4976 (K⁰), 0.938 (p), ... up to 1.672 (Ω) |
| `particleEnergy` | float32 | energy √(p²+m²) with the assigned mass | the true energy |
| `charge` | int16 | electric charge −1, 0, +1 | same |

### Identity and provenance

| Branch | Type | `t` (reconstructed) | `tgen`, `tgenBefore` (generator) |
|---|---|---|---|
| `pwflag` | int16 | object class, from the ALEPH convention: **0** = charged track, **4** = neutral energy in the electromagnetic calorimeter or the forward detectors, **5** = neutral energy in the hadron calorimeter | **0** = charged particle, **4** = neutral particle, **−11** = neutrino (ν_e, ν_μ, ν_τ; invisible in the detector) |
| `pid` | int32 | not available: −999 (the ntuple carries no particle identification) | PDG particle code with the sign of the antiparticle: 22 γ, ±211 π, ±321 K, 130 K⁰_L, 310 K⁰_S, ±2212 p, ±2112 n, ±11 e, ±13 μ, ±12/14/16 ν, ±3122 Λ, ±3112/3222 Σ, ±3312/3322 Ξ, ±3334 Ω. Note that K⁰_S and the hyperons are generator-stable here: their decays were simulated inside GOPAL, so `tgen` contains the K⁰_S, not its two pions |
| `opalOriginType` | int32 | MT provenance code: **1** = charged central track, **3** = unassociated ECAL cluster, **4** = residual energy of an ECAL cluster that a track pointed to, **5** = unassociated HCAL cluster, **6** = residual energy of a track-associated HCAL cluster (codes 7, 9, 11 for the forward detectors do not occur in these samples) | −1 |
| `opalOriginIndex` | int32 | one-based index of the object in the corresponding OPAL bank (`CTRK` for tracks, `ECAL`/`HCAL` for clusters), as in the original DST | one-based index of the particle in the OPAL `TREE` generator record |

### Track quality (reconstructed tracks only; −999 / −127 elsewhere)

| Branch | Type | Meaning | Values for tracks in `t` |
|---|---|---|---|
| `d0` | float32 | signed distance of closest approach of the track to the beam axis in the transverse (x–y) plane, cm. Large values indicate particles not produced at the collision point (decays in flight, conversions, cosmic rays) | −2.5 ... +2.5 (MT quality cut) |
| `z0` | float32 | z coordinate of the track at that point of closest approach, cm | −30 ... +30 (MT quality cut) |
| `ntpc` | int16 | ALEPH name (TPC hits); here the number of hits on the track in the OPAL **jet chamber CJ** (`CTRK` word `JCNHCJ`). The jet chamber is OPAL's main tracking detector, a 4 m long cylindrical drift chamber with up to 159 sense wires along a track | 20 ... 210, mean 126 |
| `nitc` | int16 | ALEPH name (ITC hits); here the number of hits in the OPAL **vertex chamber CV** (axial + stereo wires, unpacked from `opalCVHitsPacked`) | 0 ... 18, mean 10 |
| `nvdet` | int16 | ALEPH name (silicon vertex detector hits). OPAL's silicon microvertex hit counts are not in the staging tables, so this is always −127 | −127 |
| `opalCVHitsPacked` | int32 | the raw OPAL encoding of the CV hits: `axial + 100 × stereo`, e.g. 612 = 12 axial and 6 stereo hits | 0 ... 612 |
| `highPurity` | bool | the converter's **HIMS good-track** definition: pT ≥ 0.1 GeV, \|cos θ\| ≤ 0.966, \|d0\| ≤ 2.5 cm, \|z0\| ≤ 50 cm, ≥ 20 CJ hits. Always true for neutrals. In practice also true for every track in `t`, because MT's own quality selection is already at least as tight — this is why `nChargedHadronsHP == nChargedHadrons`. At generator level it means "inside the tracking acceptance": pT ≥ 0.1 GeV and \|cos θ\| ≤ 0.966 (93% of charged generator particles); in `tgenBefore` it is set to true for all | true |

### Generator-level vertex and weights

| Branch | Type | `t` | `tgen`, `tgenBefore` |
|---|---|---|---|
| `vx`, `vy`, `vz` | float32 | −999 | production point of the particle, cm. Primaries sit at the simulated 1994 beam spot, (−0.021, 0.043, 0.46) cm on average with the beam spread; decay products of K⁰_S etc. are displaced by up to a few cm |
| `weight` | float32 | ALEPH placeholder, always 1 | 1 |

## Event selection built into the files

Data: every valid event of the official pass-7 HIMS multihadron lists for
periods p59–p68 is kept (`selection: "none"` in the metadata); the only
events dropped are those where MT reported an error. Therefore the data `t`
tree still contains the 0.9% of events that fail the track arm
(`passesOPALHIMSTrack == false`).

MC: `t` and `tgen` contain only events that pass `mtError == 0` **and** the
HIMS track arm (`selection: "hims-track"`), 651,015 of 653,749 generated
events. `tgenBefore` contains all 653,749, flagged by `passesAll`.

To compare data with MC on equal footing, apply `passesOPALHIMSTrack` to the
data yourself — `CompareDataMC` does this.

## The `opal_metadata` objects

Besides the trees, each file holds `TNamed` objects called `opal_metadata`
whose title is a JSON string, one per input file that went into the merge
(68 cycles in the data file, 523 in the MC file, one per partition). Each
records the converter version, the staging prefix, the `source_ids`, the
event counts per tree, the selection, and that truth codes are PDG.

## What is *not* in these files

- **No particle identification**: no dE/dx, no electron/muon/kaon/proton
  weights, no time of flight. The extractor's `tracks` table does contain the
  jet-chamber dE/dx and OPAL's five dE/dx hypothesis weights, but the ntuple
  writer does not copy them. Adding them requires re-running the CERN
  extraction.
- **No calorimeter detail**: cluster energies, shapes and track–cluster
  associations are reduced to the single MT particle list.
- **No silicon microvertex hit counts**, no secondary-vertex finding.
- **No run conditions**: beam spot, detector status and the exact beam energy
  per run live in the OPAL conditions databases (`opalcal.rzexp5/6`,
  `sudb.rzdata` on `/cvmfs/opal.cern.ch`), not in the ntuple.
- **No generator history**: `tgen` has final-state particles only, without
  mother–daughter links, so a φ or a ρ cannot be identified at generator level.
- The fraction of events with an MT error, and the GPMH calorimeter arm of
  the hadronic selection, are not recoverable from the ntuple.

## Known quirks

- `missP`/`missPt` have rare unphysical outliers of several TeV; guard against
  them before using the missing momentum.
- `mass` is an assignment, not a measurement (see above).
- `Energy` in data is the recorded LEP energy of the run, so it takes
  slightly different values from run to run; the MC uses 91.208 or
  91.210 GeV throughout.
- `eventID` restarts at 0 in every extraction job; use `uniqueID` to identify
  an event globally.
- The `n…` counter branches duplicate `nParticle`.

## Reproducing the numbers on this page

```bash
python3 inspect_ntuple.py 50000      # branch types, ranges and unique values per tree
```
