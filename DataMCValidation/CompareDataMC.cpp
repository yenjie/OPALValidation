// =============================================================================
//  CompareDataMC.cpp
//
//  Compare OPAL 1994 collision data with the OPAL Monte Carlo simulation.
//
//  This program was written to be read by undergraduate students who are new
//  to ROOT and to e+e- physics.  The comments explain both the physics and the
//  programming choices.  Every ROOT class used here is documented in the ROOT
//  reference guide: https://root.cern/doc/master/
// =============================================================================
//
//  PHYSICS BACKGROUND
//
//  In 1994 the LEP collider at CERN collided electrons and positrons at a
//  centre-of-mass energy of about 91.2 GeV, the mass of the Z boson.  The Z
//  decays into a quark-antiquark pair, and the quarks turn into sprays ("jets")
//  of hadrons -- pions, kaons, protons, ... -- that the OPAL detector measured.
//
//  To interpret such data physicists rely on a Monte Carlo (MC) simulation:
//  the JETSET 7.4 event generator produces the hadrons, and GOPAL (a GEANT 3
//  model of the OPAL detector) simulates how the detector responds to them.
//  Before using the simulation in a measurement we must check that it
//  reproduces the data.  This program performs two such checks:
//
//   1. The transverse-momentum (pT) spectrum of charged particles.
//      pT is the component of a particle's momentum perpendicular to the beam
//      axis.  Its distribution tests both the physics model (how quarks turn
//      into hadrons) and the detector simulation (tracking efficiency and
//      momentum resolution).
//
//   2. The phi(1020) meson.  The phi decays to a K+ K- pair about half of the
//      time.  If we compute the invariant mass of every pair of oppositely
//      charged tracks, pretending that both tracks are kaons, the pairs that
//      really come from a phi pile up in a narrow peak at 1019.5 MeV, while
//      all other pairs form a smooth "combinatorial" background.  Comparing the
//      position, width and size of the peak in data and MC tests the momentum
//      scale, the momentum resolution and the simulated rate of strange
//      particles.
//
//  THE INPUT FILES
//
//  The original OPAL data were converted into plain ROOT files with one entry
//  per collision event.  Each file contains a TTree called "t" with the
//  reconstructed particles.  The MC file has two more trees:
//
//     t           reconstructed particles (what the detector saw)
//     tgen        the generator's final-state particles for the same events
//                 (what was actually produced, before the detector)
//     tgenBefore  the same, but for all generated events (not used here)
//
//  Each event stores a few scalars (run number, flags, ...) and, for each
//  particle, a set of variable-length arrays: px, py, pz, pt, charge, ...
//  The quantities we use are:
//
//     mtError              0 if the event was reconstructed without problems
//     passesOPALHIMSTrack  true if OPAL's standard hadronic-event track
//                          selection found at least five good tracks
//     px, py, pz           momentum components [GeV]
//     pt, pmag             transverse momentum and |p| [GeV]
//     charge               electric charge (+1, -1, or 0 for neutral objects)
//     pwflag               what kind of object: 0 = charged track,
//                          4/5 = neutral calorimeter energy
//     highPurity           true for tracks that pass OPAL's good-track cuts
//                          (pT >= 0.1 GeV, |cos theta| <= 0.966, |d0| <= 2.5 cm,
//                          |z0| <= 50 cm, >= 20 hits in the jet chamber)
//     d0, z0               distance of closest approach of the track to the
//                          beam line, in the transverse plane and along the
//                          beam [cm]; large values indicate particles that did
//                          not come from the collision point
//
//  HOW THE CODE IS ORGANISED
//
//     Part 1   Physics constants
//     Part 2   Command-line options
//     Part 3   Histograms for one sample (data, MC, generator level)
//     Part 4   Reading events and selecting good tracks
//     Part 5   The event loop that fills the histograms
//     Part 6   Fitting the phi -> K+K- peak
//     Part 7   Drawing the comparison plots
//     Part 8   main(): putting everything together
//
//  HOW TO BUILD AND RUN
//
//     make
//     ./CompareDataMC --data OPAL_1994_data.root --mc OPAL_1994_mc_jt74mh.root
//     ./CompareDataMC --help          (lists all options)
//
//  or simply ./run_grendel01.sh, which fills in the file locations.
// =============================================================================

#include <algorithm>   // std::max, std::min
#include <cmath>       // std::sqrt, std::fabs, std::pow, std::exp, std::log
#include <fstream>     // writing summary.txt
#include <iostream>    // printing to the terminal
#include <map>         // used by the command-line parser
#include <memory>      // std::unique_ptr
#include <sstream>     // building the summary text
#include <string>
#include <vector>

#include "TCanvas.h"           // a drawing area that can be saved as PDF
#include "TError.h"            // gErrorIgnoreLevel: silence ROOT's chatter
#include "TF1.h"               // one-dimensional functions, used for fitting
#include "TFile.h"             // reading and writing ROOT files
#include "TFitResult.h"        // the result of a fit (status, parameters)
#include "TFitResultPtr.h"
#include "TH1D.h"              // one-dimensional histograms with double precision
#include "TLatex.h"            // text with LaTeX-like formatting on a plot
#include "TLegend.h"           // the legend box of a plot
#include "TLine.h"             // a straight line on a plot
#include "TMath.h"             // TMath::Voigt
#include "TNamed.h"            // a named text object, used to save the summary
#include "TPad.h"              // sub-areas of a canvas
#include "TROOT.h"             // gROOT: global ROOT settings
#include "TStyle.h"            // gStyle: global drawing style
#include "TSystem.h"           // gSystem->mkdir
#include "TTree.h"
#include "TTreeReader.h"       // convenient, safe way to loop over a TTree
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

// Everything except main() lives in an anonymous namespace.  This is C++ for
// "these names are private to this file".
namespace {

// =============================================================================
//  Part 1: physics constants (all masses and widths in GeV)
// =============================================================================

const double kKaonMass    = 0.493677;              // charged kaon mass
const double kKKThreshold = 2.0 * kKaonMass;       // smallest possible K+K- mass
const double kPhiMassPDG  = 1.019461;              // phi(1020) mass
const double kPhiWidthPDG = 0.004249;              // phi natural width (FWHM)

// When we quote "signal over background" we count pairs within this distance
// of the fitted peak position.
const double kSignalWindow = 0.010;                // +- 10 MeV

// =============================================================================
//  Part 2: command-line options
//
//  The program is steered from the command line with "--name value" pairs,
//  for example:   ./CompareDataMC --data d.root --mc m.root --maxEvents 1000
//  A flag without a value (e.g. "--help") is stored as "true".
// =============================================================================

class CommandLine {
public:
  CommandLine(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
      std::string argument = argv[i];
      if (argument.size() < 3 || argument.compare(0, 2, "--") != 0) {
        std::cerr << "Ignoring unexpected argument: " << argument << std::endl;
        continue;
      }
      std::string name = argument.substr(2);    // drop the leading "--"
      std::string value = "true";
      // If the next word does not start with "--" it is the value of this option.
      if (i + 1 < argc) {
        std::string next = argv[i + 1];
        bool nextIsOption = (next.size() >= 2 && next.compare(0, 2, "--") == 0);
        if (!nextIsOption) {
          value = next;
          ++i;
        }
      }
      fOptions[name] = value;
    }
  }

  bool Has(const std::string& name) const { return fOptions.count(name) > 0; }

  std::string GetString(const std::string& name, const std::string& defaultValue) const {
    if (!Has(name)) return defaultValue;
    return fOptions.find(name)->second;
  }
  double GetDouble(const std::string& name, double defaultValue) const {
    if (!Has(name)) return defaultValue;
    return std::stod(GetString(name, ""));
  }
  int GetInt(const std::string& name, int defaultValue) const {
    if (!Has(name)) return defaultValue;
    return std::stoi(GetString(name, ""));
  }
  long long GetLong(const std::string& name, long long defaultValue) const {
    if (!Has(name)) return defaultValue;
    return std::stoll(GetString(name, ""));
  }
  bool GetBool(const std::string& name, bool defaultValue) const {
    if (!Has(name)) return defaultValue;
    std::string value = GetString(name, "");
    return !(value == "false" || value == "0" || value == "no" || value == "off");
  }

private:
  std::map<std::string, std::string> fOptions;
};

// All settings of one run of the program, collected in one place so that the
// functions below do not need long argument lists.
struct Config {
  // Input and output
  std::string dataFile;
  std::string mcFile;
  std::string outputFile;
  std::string plotDir;
  long long maxEvents = -1;         // -1 means "all events"
  bool useGen = true;               // overlay the generator-level MC on the pT plots

  // Event selection
  bool requireHIMS = true;          // require passesOPALHIMSTrack
  int nGoodTracksMin = 5;           // at least this many good charged tracks

  // Track selection.  The defaults are exactly the cuts already contained in
  // the highPurity flag, so by default no track is removed beyond that.  Use
  // the command line to tighten them.
  double trackPtMin = 0.10;         // GeV
  double trackPMin = 0.0;           // GeV
  double cosThetaMax = 0.966;       // |cos theta| of the track direction
  double d0Max = 2.5;               // cm
  double z0Max = 50.0;              // cm

  // phi -> K+K- candidate selection
  double kaonPMin = 0.0;            // minimum momentum of each kaon candidate [GeV]
  double phiPMin = 0.0;             // minimum momentum of the pair [GeV]

  // Histogram binning and fit range for the K+K- mass
  double massLow = 0.98;            // GeV
  double massHigh = 1.10;           // GeV
  double massBinWidth = 0.001;      // GeV, i.e. 1 MeV bins
  double fitLow = 0.990;            // GeV
  double fitHigh = 1.090;           // GeV

  // How to normalise the pT spectra: "event" or "area" (see NormalizedCopy)
  std::string normalize = "event";
};

// =============================================================================
//  Part 3: the histograms of one sample
//
//  We fill the same set of histograms three times: for the data, for the
//  reconstructed MC and (optionally) for the generator-level MC.  Keeping
//  them in a struct lets us write the filling code once.
// =============================================================================

struct Sample {
  std::string tag;                  // short name used in histogram names, e.g. "data"
  std::string label;                // text shown in plot legends

  long long nEventsRead = 0;        // how many events we looked at
  long long nEventsSelected = 0;    // how many passed the event selection
  double sumGoodTracks = 0.0;       // to compute the mean number of good tracks

  TH1D* hNGoodTracks = nullptr;     // number of good charged tracks per event
  TH1D* hPt = nullptr;              // pT of good tracks, uniform bins
  TH1D* hPtLog = nullptr;           // pT of good tracks, logarithmic bins
  TH1D* hMassOS = nullptr;          // K+K- mass of opposite-sign pairs (contains the phi)
  TH1D* hMassLS = nullptr;          // K+K- mass of like-sign pairs (background only)
};

// Bin edges that are equally spaced in log(x).  A steeply falling spectrum such
// as dN/dpT is much easier to look at with such bins on a logarithmic axis.
std::vector<double> LogarithmicBinEdges(int nBins, double low, double high) {
  std::vector<double> edges(nBins + 1);
  double step = std::log(high / low) / nBins;
  for (int i = 0; i <= nBins; ++i) edges[i] = low * std::exp(step * i);
  return edges;
}

Sample MakeSample(const std::string& tag, const std::string& label, const Config& cfg) {
  Sample s;
  s.tag = tag;
  s.label = label;

  // The histogram title string ";x title;y title" sets the axis titles.
  s.hNGoodTracks = new TH1D(("hNGoodTracks_" + tag).c_str(),
                            ";N_{good charged tracks};Events", 71, -0.5, 70.5);

  s.hPt = new TH1D(("hPt_" + tag).c_str(),
                   ";p_{T} [GeV];(1/N_{ev}) dN/dp_{T} [GeV^{-1}]", 100, 0.0, 10.0);

  std::vector<double> edges = LogarithmicBinEdges(50, 0.1, 50.0);
  s.hPtLog = new TH1D(("hPtLog_" + tag).c_str(),
                      ";p_{T} [GeV];(1/N_{ev}) dN/dp_{T} [GeV^{-1}]",
                      static_cast<int>(edges.size()) - 1, edges.data());

  // Form() is ROOT's printf that returns a string: handy for titles.
  int nMassBins = static_cast<int>(std::lround((cfg.massHigh - cfg.massLow) / cfg.massBinWidth));
  std::string massTitle = Form(";m(K^{+}K^{-}) [GeV];Pairs / %.1f MeV", cfg.massBinWidth * 1000.0);
  s.hMassOS = new TH1D(("hMassOS_" + tag).c_str(), massTitle.c_str(), nMassBins, cfg.massLow, cfg.massHigh);
  s.hMassLS = new TH1D(("hMassLS_" + tag).c_str(), massTitle.c_str(), nMassBins, cfg.massLow, cfg.massHigh);
  return s;
}

// =============================================================================
//  Part 4: reading events and selecting good tracks
//
//  TTreeReader is ROOT's modern way to read a TTree.  You declare one
//  TTreeReaderValue per scalar branch and one TTreeReaderArray per
//  per-particle branch, and then call Next() to move from event to event.
//  A TTreeReaderArray behaves like a std::vector that is refilled with the
//  particles of the current event; its size is the number of particles.
// =============================================================================

struct EventReader {
  TTreeReader reader;
  TTreeReaderValue<Int_t> mtError;
  TTreeReaderValue<Bool_t> passesHIMS;
  TTreeReaderArray<Float_t> px;
  TTreeReaderArray<Float_t> py;
  TTreeReaderArray<Float_t> pz;
  TTreeReaderArray<Float_t> pt;
  TTreeReaderArray<Float_t> pmag;
  TTreeReaderArray<Float_t> d0;
  TTreeReaderArray<Float_t> z0;
  TTreeReaderArray<Short_t> charge;
  TTreeReaderArray<Short_t> pwflag;
  TTreeReaderArray<Bool_t> highPurity;

  // The strings are the branch names inside the ROOT file.
  explicit EventReader(TTree* tree)
      : reader(tree),
        mtError(reader, "mtError"),
        passesHIMS(reader, "passesOPALHIMSTrack"),
        px(reader, "px"),
        py(reader, "py"),
        pz(reader, "pz"),
        pt(reader, "pt"),
        pmag(reader, "pmag"),
        d0(reader, "d0"),
        z0(reader, "z0"),
        charge(reader, "charge"),
        pwflag(reader, "pwflag"),
        highPurity(reader, "highPurity") {}

  bool Next() { return reader.Next(); }
  int NumberOfParticles() { return static_cast<int>(px.GetSize()); }
};

// A copy of the few numbers we need about one good track.  Copying them out of
// the TTreeReaderArrays makes the rest of the code easier to read.
struct Track {
  double px, py, pz;    // momentum components [GeV]
  double pt;            // transverse momentum [GeV]
  double p;             // magnitude of the momentum [GeV]
  int charge;           // +1 or -1
};

// Does the event as a whole qualify?  We require a clean reconstruction and
// OPAL's standard "at least five good tracks" hadronic-event criterion.  The
// MC trees were produced with these two requirements already applied, so
// applying them to the data puts both samples on the same footing.
bool PassesEventSelection(EventReader& ev, const Config& cfg) {
  if (*ev.mtError != 0) return false;
  if (cfg.requireHIMS && !*ev.passesHIMS) return false;
  return true;
}

// Is particle number i a good charged track?
//
// At generator level (isGeneratorLevel == true) the "tracks" are the true
// charged final-state particles.  They have no d0/z0, so those cuts are
// skipped; highPurity then simply means "inside the detector acceptance".
bool IsGoodChargedTrack(EventReader& ev, int i, const Config& cfg, bool isGeneratorLevel) {
  if (ev.charge[i] == 0) return false;          // neutral: not a track
  if (ev.pwflag[i] != 0) return false;          // not a charged-track object
  if (!ev.highPurity[i]) return false;          // fails OPAL's good-track cuts

  // Optional tighter kinematic cuts (by default identical to highPurity).
  if (ev.pt[i] < cfg.trackPtMin) return false;
  if (ev.pmag[i] < cfg.trackPMin) return false;
  if (ev.pmag[i] <= 0.0) return false;
  double cosTheta = ev.pz[i] / ev.pmag[i];      // theta = angle to the beam axis
  if (std::fabs(cosTheta) > cfg.cosThetaMax) return false;

  if (!isGeneratorLevel) {
    if (std::fabs(ev.d0[i]) > cfg.d0Max) return false;
    if (std::fabs(ev.z0[i]) > cfg.z0Max) return false;
  }
  return true;
}

// Collect the good charged tracks of the current event.
std::vector<Track> SelectGoodTracks(EventReader& ev, const Config& cfg, bool isGeneratorLevel) {
  std::vector<Track> tracks;
  for (int i = 0; i < ev.NumberOfParticles(); ++i) {
    if (!IsGoodChargedTrack(ev, i, cfg, isGeneratorLevel)) continue;
    Track t;
    t.px = ev.px[i];
    t.py = ev.py[i];
    t.pz = ev.pz[i];
    t.pt = ev.pt[i];
    t.p = ev.pmag[i];
    t.charge = ev.charge[i];
    tracks.push_back(t);
  }
  return tracks;
}

// Invariant mass of two tracks under the hypothesis that both are kaons.
//
// Special relativity: for a system of particles the quantity
//     m^2 = (sum of energies)^2 - |sum of momenta|^2
// is the same in every reference frame.  For the two daughters of a phi decay
// it equals the phi mass.  The energy of each track is E = sqrt(p^2 + m_K^2);
// we do not know which tracks are kaons, so we simply assume all of them are.
// The pair momentum |p1 + p2| is also returned because it can be used as a cut.
double InvariantMassKK(const Track& a, const Track& b, double& pairMomentum) {
  double energyA = std::sqrt(a.p * a.p + kKaonMass * kKaonMass);
  double energyB = std::sqrt(b.p * b.p + kKaonMass * kKaonMass);
  double sumPx = a.px + b.px;
  double sumPy = a.py + b.py;
  double sumPz = a.pz + b.pz;
  double sumP2 = sumPx * sumPx + sumPy * sumPy + sumPz * sumPz;
  pairMomentum = std::sqrt(sumP2);
  double mass2 = (energyA + energyB) * (energyA + energyB) - sumP2;
  return mass2 > 0.0 ? std::sqrt(mass2) : 0.0;
}

// =============================================================================
//  Part 5: the event loop
// =============================================================================

void FillHistogramsFromTree(TTree* tree, const Config& cfg, bool isGeneratorLevel, Sample& s) {
  EventReader ev(tree);

  long long nEntries = tree->GetEntries();
  long long nToRead = nEntries;
  if (cfg.maxEvents > 0 && cfg.maxEvents < nEntries) nToRead = cfg.maxEvents;
  std::cout << "[" << s.tag << "] reading " << nToRead << " of " << nEntries
            << " events from tree '" << tree->GetName() << "'" << std::endl;

  while (ev.Next()) {
    if (s.nEventsRead >= nToRead) break;
    ++s.nEventsRead;
    if (s.nEventsRead % 200000 == 0) std::cout << "   ... " << s.nEventsRead << " events" << std::endl;

    // ---- event selection -------------------------------------------------
    if (!PassesEventSelection(ev, cfg)) continue;
    std::vector<Track> tracks = SelectGoodTracks(ev, cfg, isGeneratorLevel);
    if (static_cast<int>(tracks.size()) < cfg.nGoodTracksMin) continue;

    ++s.nEventsSelected;
    s.sumGoodTracks += tracks.size();
    s.hNGoodTracks->Fill(tracks.size());

    // ---- pT spectrum: one entry per good track ---------------------------
    for (const Track& t : tracks) {
      s.hPt->Fill(t.pt);
      s.hPtLog->Fill(t.pt);
    }

    // ---- K+K- invariant mass: one entry per pair of good tracks ----------
    // The double loop visits every pair exactly once (j > i).  Pairs with
    // opposite charges may contain a phi; pairs with equal charges cannot and
    // show what pure combinatorial background looks like.
    for (size_t i = 0; i < tracks.size(); ++i) {
      if (tracks[i].p < cfg.kaonPMin) continue;
      for (size_t j = i + 1; j < tracks.size(); ++j) {
        if (tracks[j].p < cfg.kaonPMin) continue;

        double pairMomentum = 0.0;
        double mass = InvariantMassKK(tracks[i], tracks[j], pairMomentum);
        if (pairMomentum < cfg.phiPMin) continue;
        if (mass >= cfg.massHigh) continue;       // outside the histogram anyway

        bool oppositeSign = (tracks[i].charge * tracks[j].charge < 0);
        if (oppositeSign)
          s.hMassOS->Fill(mass);
        else
          s.hMassLS->Fill(mass);
      }
    }
  }

  std::cout << "[" << s.tag << "] selected " << s.nEventsSelected << " of " << s.nEventsRead
            << " events" << std::endl;
}

// =============================================================================
//  Part 6: fitting the phi -> K+K- peak
//
//  The opposite-sign mass spectrum is modelled as  signal + background.
//
//  Signal.  An unstable particle does not have one sharp mass: its mass
//  distribution is a Breit-Wigner (Lorentzian) curve whose full width at half
//  maximum is the natural width Gamma (4.25 MeV for the phi).  On top of that
//  the detector measures each momentum with a finite resolution, which smears
//  the peak with a Gaussian of width sigma.  The convolution of a Breit-Wigner
//  with a Gaussian is called a Voigt profile; ROOT provides it as
//  TMath::Voigt(x, sigma, Gamma), normalised to unit area.  We fix Gamma to
//  its known value and let the fit determine sigma, i.e. the mass resolution.
//
//  Background.  Random pairs of tracks give a smooth spectrum that starts at
//  the K+K- threshold 2*m_K = 0.987 GeV and rises from there.  We use the
//  empirical shape  b0 * (m - 2m_K)^a * exp(c * (m - 2m_K)),  which is zero
//  at threshold and flexible enough to follow the data in the fit range.
//
//  Fit parameters:  [0] N      number of phi mesons (the signal yield)
//                   [1] m0     peak position
//                   [2] sigma  Gaussian resolution
//                   [3] Gamma  natural width (fixed)
//                   [4] b0, [5] a, [6] c   background shape
// =============================================================================

// TF1 evaluates user functions through this C-style signature:
// x[0] is the mass, p[] are the parameters.  The signal is multiplied by the
// bin width so that N counts entries, not entries per GeV.
double gMassBinWidth = 0.001;

double BackgroundFunction(double* x, double* p) {
  double dm = x[0] - kKKThreshold;
  if (dm <= 0.0) return 0.0;
  return p[0] * std::pow(dm, p[1]) * std::exp(p[2] * dm);
}

double SignalFunction(double* x, double* p) {
  return p[0] * gMassBinWidth * TMath::Voigt(x[0] - p[1], p[2], p[3]);
}

double SignalPlusBackground(double* x, double* p) {
  return SignalFunction(x, p) + BackgroundFunction(x, p + 4);
}

// Everything we want to know about one fitted peak.
struct PhiFitResult {
  int status = -1;                  // 0 means the fit converged
  double yield = 0, yieldError = 0; // number of phi mesons
  double mass = 0, massError = 0;   // GeV
  double sigma = 0, sigmaError = 0; // GeV
  double chi2 = 0;                  // goodness of fit ...
  int ndf = 0;                      // ... per degree of freedom (chi2/ndf ~ 1 is good)
  double signalInWindow = 0;        // fitted signal within +- kSignalWindow of the peak
  double backgroundInWindow = 0;    // fitted background in the same window
  TF1* total = nullptr;             // signal + background, for drawing
  TF1* signal = nullptr;
  TF1* background = nullptr;
};

PhiFitResult FitPhiPeak(TH1D* h, const Config& cfg, const std::string& tag) {
  PhiFitResult result;
  gMassBinWidth = h->GetBinWidth(1);

  TF1* total = new TF1(("fit_total_" + tag).c_str(), SignalPlusBackground, cfg.fitLow, cfg.fitHigh, 7);
  total->SetParNames("N_{#phi}", "m_{#phi}", "#sigma", "#Gamma_{#phi}", "b_{0}", "a", "c");

  // A fit needs reasonable starting values.  We estimate the background
  // normalisation from the region above the peak (1.05-1.085 GeV), assuming
  // the shape (m - 2m_K)^0.5, and the signal from the excess over that estimate
  // within +- 10 MeV of the known phi mass.
  const double aStart = 0.5;
  const double cStart = 0.0;
  double b0Start = 0.0;
  int nReferenceBins = 0;
  for (int bin = h->FindBin(1.050); bin <= h->FindBin(std::min(1.085, cfg.fitHigh)); ++bin) {
    double dm = h->GetBinCenter(bin) - kKKThreshold;
    if (dm <= 0.0) continue;
    b0Start += h->GetBinContent(bin) / std::pow(dm, aStart);
    ++nReferenceBins;
  }
  b0Start = nReferenceBins > 0 ? b0Start / nReferenceBins : 1.0;

  double countsInWindow = 0.0;
  double backgroundInWindow = 0.0;
  for (int bin = h->FindBin(kPhiMassPDG - kSignalWindow); bin <= h->FindBin(kPhiMassPDG + kSignalWindow); ++bin) {
    countsInWindow += h->GetBinContent(bin);
    double dm = h->GetBinCenter(bin) - kKKThreshold;
    if (dm > 0.0) backgroundInWindow += b0Start * std::pow(dm, aStart);
  }
  double yieldStart = std::max(countsInWindow - backgroundInWindow, 0.05 * countsInWindow);
  if (yieldStart <= 0.0) yieldStart = 1.0;

  total->SetParameters(yieldStart, kPhiMassPDG, 0.0025, kPhiWidthPDG, b0Start, aStart, cStart);
  total->SetParLimits(0, 0.0, 1.0e10);        // the yield cannot be negative
  total->SetParLimits(1, 1.010, 1.030);       // keep the peak near the phi mass
  total->SetParLimits(2, 0.0005, 0.0100);     // resolution between 0.5 and 10 MeV
  total->FixParameter(3, kPhiWidthPDG);       // natural width is known
  total->SetParLimits(4, 0.0, 1.0e12);
  total->SetParLimits(5, 0.05, 4.0);
  total->SetParLimits(6, -50.0, 50.0);

  // Fit options:  L = use a Poisson likelihood (correct for counting data),
  // R = only fit inside the function's range, Q = quiet, S = return the
  // result object, 0 = do not draw.  The first pass only improves the starting
  // values for the second pass.
  h->Fit(total, "LRQ0");
  TFitResultPtr fitResult = h->Fit(total, "LRSQ");

  result.status = fitResult.Get() ? fitResult->Status() : -1;
  result.yield = total->GetParameter(0);
  result.yieldError = total->GetParError(0);
  result.mass = total->GetParameter(1);
  result.massError = total->GetParError(1);
  result.sigma = total->GetParameter(2);
  result.sigmaError = total->GetParError(2);
  result.chi2 = total->GetChisquare();
  result.ndf = total->GetNDF();
  result.total = total;

  // Separate copies of the two components, for drawing and for integrals.
  result.signal = new TF1(("fit_signal_" + tag).c_str(), SignalFunction, cfg.fitLow, cfg.fitHigh, 4);
  result.signal->SetParameters(total->GetParameter(0), total->GetParameter(1),
                               total->GetParameter(2), total->GetParameter(3));
  result.background = new TF1(("fit_background_" + tag).c_str(), BackgroundFunction, cfg.fitLow, cfg.fitHigh, 3);
  result.background->SetParameters(total->GetParameter(4), total->GetParameter(5), total->GetParameter(6));

  // TF1::Integral gives the area under the curve; dividing by the bin width
  // converts it back into a number of entries.
  double windowLow = result.mass - kSignalWindow;
  double windowHigh = result.mass + kSignalWindow;
  result.signalInWindow = result.signal->Integral(windowLow, windowHigh) / gMassBinWidth;
  result.backgroundInWindow = result.background->Integral(windowLow, windowHigh) / gMassBinWidth;
  return result;
}

// Subtract the fitted background from the mass histogram and divide by the
// number of events, so that data and MC peaks can be overlaid directly.
TH1D* BackgroundSubtractedCopy(const TH1D* h, const PhiFitResult& fit, double nEvents,
                               const Config& cfg, const std::string& name) {
  TH1D* copy = static_cast<TH1D*>(h->Clone(name.c_str()));
  for (int bin = 1; bin <= copy->GetNbinsX(); ++bin) {
    double center = copy->GetBinCenter(bin);
    bool insideFitRange = (center >= cfg.fitLow && center <= cfg.fitHigh);
    if (!insideFitRange || fit.background == nullptr) {
      copy->SetBinContent(bin, 0.0);
      copy->SetBinError(bin, 0.0);
    } else {
      copy->SetBinContent(bin, h->GetBinContent(bin) - fit.background->Eval(center));
      copy->SetBinError(bin, h->GetBinError(bin));   // the subtraction adds no statistical error
    }
  }
  if (nEvents > 0.0) copy->Scale(1.0 / nEvents);
  copy->GetYaxis()->SetTitle(Form("(1/N_{ev}) pairs / %.1f MeV", gMassBinWidth * 1000.0));
  return copy;
}

// =============================================================================
//  Part 7: drawing
// =============================================================================

// The data contain 1.9 million events and the MC 0.65 million, so the raw
// histograms cannot be compared directly.  Two normalisations are offered:
//
//   "event": divide by the number of selected events and by the bin width.
//            The result is the number of tracks per event per GeV, so a
//            difference in the charged multiplicity shows up as a different
//            overall height.
//   "area":  divide by the total number of entries instead, which compares
//            only the shapes.
TH1D* NormalizedCopy(const TH1D* h, double nEvents, const std::string& mode, const std::string& suffix) {
  TH1D* copy = static_cast<TH1D*>(h->Clone((std::string(h->GetName()) + suffix).c_str()));
  if (mode == "area") {
    double integral = copy->Integral();
    if (integral > 0.0) copy->Scale(1.0 / integral, "width");   // "width" also divides by the bin width
    copy->GetYaxis()->SetTitle("(1/N_{trk}) dN/dp_{T} [GeV^{-1}]");
  } else {
    if (nEvents > 0.0) copy->Scale(1.0 / nEvents, "width");
  }
  return copy;
}

// Two lines of text in the upper left corner of the current pad.
// Coordinates are "NDC": fractions of the pad size, from 0 to 1.
void DrawHeader(const std::string& line1, const std::string& line2) {
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.045);
  latex.DrawLatex(0.16, 0.86, line1.c_str());
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.16, 0.81, line2.c_str());
}

std::string SelectionText(const Config& cfg) {
  return Form("p_{T} > %.2f GeV, |cos#theta| < %.3f, N_{good} #geq %d",
              cfg.trackPtMin, cfg.cosThetaMax, cfg.nGoodTracksMin);
}

// Draw data (points) and MC (lines) in the upper pad and the ratio data/MC in
// the lower pad.  The ratio makes small differences visible that would be
// hidden on a logarithmic scale.  Returns the ratio histogram so that it can
// be saved.
TH1D* DrawSpectrumComparison(TH1D* hData, TH1D* hMC, TH1D* hGen, const Config& cfg,
                             bool logX, bool logY, const std::string& pdfName,
                             const std::string& ratioName) {
  TCanvas canvas("canvas", "", 800, 900);
  TPad upper("upper", "", 0.0, 0.32, 1.0, 1.0);
  TPad lower("lower", "", 0.0, 0.0, 1.0, 0.32);
  upper.SetBottomMargin(0.025);
  upper.SetTopMargin(0.06);
  upper.SetLeftMargin(0.14);
  upper.SetRightMargin(0.04);
  lower.SetTopMargin(0.03);
  lower.SetBottomMargin(0.36);
  lower.SetLeftMargin(0.14);
  lower.SetRightMargin(0.04);
  if (logX) {
    upper.SetLogx();
    lower.SetLogx();
  }
  if (logY) upper.SetLogy();
  upper.Draw();
  lower.Draw();

  // Styles: black points for data, red line for reconstructed MC, dashed blue
  // line for generator level.
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(0.8);
  hData->SetMarkerColor(kBlack);
  hData->SetLineColor(kBlack);
  hMC->SetLineColor(kRed + 1);
  hMC->SetLineWidth(2);
  if (hGen) {
    hGen->SetLineColor(kBlue + 1);
    hGen->SetLineStyle(2);
    hGen->SetLineWidth(2);
  }

  // ---- upper pad ---------------------------------------------------------
  upper.cd();
  double yMax = std::max(hData->GetMaximum(), hMC->GetMaximum());
  if (hGen) yMax = std::max(yMax, hGen->GetMaximum());
  if (logY) {
    // On a log scale the lower edge must be positive: find the smallest
    // non-zero bin content among the histograms.
    double yMin = 1.0e300;
    for (TH1D* h : {hData, hMC, hGen}) {
      if (!h) continue;
      for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        double value = h->GetBinContent(bin);
        if (value > 0.0) yMin = std::min(yMin, value);
      }
    }
    if (yMin >= 1.0e299) yMin = 1.0e-6;
    hMC->SetMaximum(yMax * 8.0);
    hMC->SetMinimum(yMin * 0.3);
  } else {
    hMC->SetMaximum(yMax * 1.35);
    hMC->SetMinimum(0.0);
  }
  hMC->GetXaxis()->SetLabelSize(0);    // the x axis is labelled in the lower pad
  hMC->GetXaxis()->SetTitleSize(0);
  hMC->GetYaxis()->SetTitleOffset(1.3);
  hMC->Draw("HIST");                   // "HIST": draw as a line without error bars
  if (hGen) hGen->Draw("HIST SAME");
  hData->Draw("E1 P SAME");            // "E1 P": points with error bars

  TLegend legend(0.55, 0.66, 0.93, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(hData, "OPAL 1994 data", "lp");
  legend.AddEntry(hMC, "JETSET 7.4 + GOPAL (reconstructed)", "l");
  if (hGen) legend.AddEntry(hGen, "JETSET 7.4 (generator level)", "l");
  legend.Draw();
  DrawHeader("OPAL 1994, e^{+}e^{-} #rightarrow hadrons, #sqrt{s} #approx 91.2 GeV", SelectionText(cfg));

  // ---- lower pad: ratio ----------------------------------------------------
  lower.cd();
  TH1D* ratio = static_cast<TH1D*>(hData->Clone(ratioName.c_str()));
  ratio->Divide(hMC);                  // bin-by-bin division, errors propagated
  ratio->SetTitle("");
  ratio->GetYaxis()->SetTitle("Data / MC");
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetTitleSize(0.11);
  ratio->GetYaxis()->SetTitleOffset(0.55);
  ratio->GetYaxis()->SetLabelSize(0.10);
  ratio->GetXaxis()->SetTitleSize(0.13);
  ratio->GetXaxis()->SetTitleOffset(1.1);
  ratio->GetXaxis()->SetLabelSize(0.10);
  if (logX) ratio->GetXaxis()->SetMoreLogLabels();
  ratio->SetMinimum(0.4);
  ratio->SetMaximum(1.6);
  ratio->Draw("E1 P");
  TLine unity(ratio->GetXaxis()->GetXmin(), 1.0, ratio->GetXaxis()->GetXmax(), 1.0);
  unity.SetLineStyle(2);
  unity.Draw();

  canvas.SaveAs(pdfName.c_str());
  return ratio;
}

// Draw one K+K- mass spectrum with its fit into the given pad.
void DrawPhiFit(TPad* pad, Sample& s, const PhiFitResult& fit) {
  pad->cd();
  pad->SetLeftMargin(0.14);
  pad->SetRightMargin(0.04);
  pad->SetTopMargin(0.06);
  pad->SetBottomMargin(0.13);

  TH1D* oppositeSign = s.hMassOS;
  TH1D* likeSign = s.hMassLS;
  oppositeSign->SetMarkerStyle(20);
  oppositeSign->SetMarkerSize(0.6);
  oppositeSign->SetMarkerColor(kBlack);
  oppositeSign->SetLineColor(kBlack);
  likeSign->SetFillColor(kGray);
  likeSign->SetLineColor(kGray + 1);
  oppositeSign->SetMinimum(0.0);
  oppositeSign->SetMaximum(1.3 * std::max(oppositeSign->GetMaximum(), likeSign->GetMaximum()));
  oppositeSign->GetYaxis()->SetTitleOffset(1.5);

  oppositeSign->Draw("AXIS");          // first only the axes, to fix the frame
  likeSign->Draw("HIST SAME");         // grey background reference
  oppositeSign->Draw("E1 P SAME");     // the spectrum with the phi
  if (fit.total) {
    fit.total->SetLineColor(kRed + 1);
    fit.total->SetLineWidth(2);
    fit.total->SetNpx(1000);           // draw the curve smoothly
    fit.total->Draw("SAME");
  }
  if (fit.background) {
    fit.background->SetLineColor(kBlue + 1);
    fit.background->SetLineStyle(2);
    fit.background->SetLineWidth(2);
    fit.background->SetNpx(1000);
    fit.background->Draw("SAME");
  }
  oppositeSign->Draw("AXIS SAME");     // redraw the axes on top

  // The legend is created with "new" because the canvas is only saved after
  // this function returns; a local object would already have been destroyed.
  TLegend* legend = new TLegend(0.50, 0.70, 0.94, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(oppositeSign, (s.label + ", K^{+}K^{-} pairs").c_str(), "lp");
  legend->AddEntry(likeSign, "like-sign pairs", "f");
  if (fit.total) legend->AddEntry(fit.total, "Voigtian + background", "l");
  if (fit.background) legend->AddEntry(fit.background, "background", "l");
  legend->Draw();

  double nEvents = std::max(1LL, s.nEventsSelected);
  double signalOverBackground = fit.backgroundInWindow > 0 ? fit.signalInWindow / fit.backgroundInWindow : 0.0;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.032);
  double y = 0.62;
  const double lineStep = 0.045;
  latex.DrawLatex(0.50, y, Form("N_{#phi} = %.0f #pm %.0f", fit.yield, fit.yieldError));
  y -= lineStep;
  latex.DrawLatex(0.50, y, Form("N_{#phi}/event = (%.3f #pm %.3f) #times 10^{-3}",
                                1000.0 * fit.yield / nEvents, 1000.0 * fit.yieldError / nEvents));
  y -= lineStep;
  latex.DrawLatex(0.50, y, Form("m_{#phi} = %.2f #pm %.2f MeV", 1000.0 * fit.mass, 1000.0 * fit.massError));
  y -= lineStep;
  latex.DrawLatex(0.50, y, Form("#sigma_{res} = %.2f #pm %.2f MeV", 1000.0 * fit.sigma, 1000.0 * fit.sigmaError));
  y -= lineStep;
  latex.DrawLatex(0.50, y, Form("#Gamma_{#phi} = %.2f MeV (fixed)", 1000.0 * kPhiWidthPDG));
  y -= lineStep;
  latex.DrawLatex(0.50, y, Form("S/B (#pm%.0f MeV) = %.3f", 1000.0 * kSignalWindow, signalOverBackground));
  y -= lineStep;
  latex.DrawLatex(0.50, y, Form("#chi^{2}/ndf = %.1f / %d", fit.chi2, fit.ndf));
  DrawHeader("OPAL 1994, #sqrt{s} #approx 91.2 GeV", "");
}

// Overlay the background-subtracted phi peaks of data and MC.
void DrawSubtractedPeaks(TH1D* subData, TH1D* subMC, const Config& cfg, const std::string& pdfName) {
  TCanvas canvas("canvasSubtracted", "", 800, 650);
  canvas.SetLeftMargin(0.14);
  canvas.SetRightMargin(0.04);

  subData->SetMarkerStyle(20);
  subData->SetMarkerSize(0.7);
  subData->SetMarkerColor(kBlack);
  subData->SetLineColor(kBlack);
  subMC->SetMarkerStyle(24);
  subMC->SetMarkerSize(0.7);
  subMC->SetMarkerColor(kRed + 1);
  subMC->SetLineColor(kRed + 1);

  double yMax = std::max(subData->GetMaximum(), subMC->GetMaximum());
  double yMin = std::min(subData->GetMinimum(), subMC->GetMinimum());
  subData->SetMaximum(yMax * 1.4);
  subData->SetMinimum(std::min(0.0, yMin * 1.2));
  subData->GetYaxis()->SetTitleOffset(1.5);
  subData->Draw("E1 P");
  subMC->Draw("E1 P SAME");

  TLine zero(cfg.massLow, 0.0, cfg.massHigh, 0.0);
  zero.SetLineStyle(2);
  zero.Draw();

  TLegend legend(0.50, 0.72, 0.94, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(subData, "OPAL 1994 data", "lp");
  legend.AddEntry(subMC, "JETSET 7.4 + GOPAL", "lp");
  legend.Draw();
  DrawHeader("OPAL 1994, #phi #rightarrow K^{+}K^{-}", "fitted background subtracted, per selected event");

  canvas.SaveAs(pdfName.c_str());
}

// =============================================================================
//  Text summary
// =============================================================================

void AddSampleSummary(std::ostringstream& out, const Sample& s) {
  double meanGoodTracks = s.nEventsSelected > 0 ? s.sumGoodTracks / s.nEventsSelected : 0.0;
  out << Form("  %-28s read %10lld   selected %10lld   <N_good> = %6.3f   <pT> = %6.4f GeV\n",
              s.label.c_str(), s.nEventsRead, s.nEventsSelected, meanGoodTracks, s.hPt->GetMean());
}

void AddPhiSummary(std::ostringstream& out, const Sample& s, const PhiFitResult& fit) {
  double nEvents = std::max(1LL, s.nEventsSelected);
  double signalOverBackground = fit.backgroundInWindow > 0 ? fit.signalInWindow / fit.backgroundInWindow : 0.0;
  out << Form("  %-28s fit status %d\n", s.label.c_str(), fit.status);
  out << Form("      N_phi            = %10.0f +- %7.0f   (%.4f +- %.4f) x 10^-3 per event\n",
              fit.yield, fit.yieldError, 1000.0 * fit.yield / nEvents, 1000.0 * fit.yieldError / nEvents);
  out << Form("      mass             = %8.3f +- %.3f MeV   (PDG %.3f MeV)\n",
              1000.0 * fit.mass, 1000.0 * fit.massError, 1000.0 * kPhiMassPDG);
  out << Form("      resolution sigma = %8.3f +- %.3f MeV\n", 1000.0 * fit.sigma, 1000.0 * fit.sigmaError);
  out << Form("      chi2/ndf         = %.1f / %d,   S/B within +-%.0f MeV = %.3f\n",
              fit.chi2, fit.ndf, 1000.0 * kSignalWindow, signalOverBackground);
}

std::string BuildSummary(const Config& cfg, const Sample& data, const Sample& mc, const Sample& gen,
                         const PhiFitResult& fitData, const PhiFitResult& fitMC) {
  std::ostringstream out;
  out << "OPAL 1994 data / MC validation\n"
      << "  data : " << cfg.dataFile << "\n"
      << "  MC   : " << cfg.mcFile << "\n"
      << "  event selection : mtError == 0" << (cfg.requireHIMS ? ", passesOPALHIMSTrack" : "")
      << ", N_good >= " << cfg.nGoodTracksMin << "\n"
      << "  track selection : charged, pwflag == 0, highPurity, pT > " << cfg.trackPtMin
      << " GeV, p > " << cfg.trackPMin << " GeV, |cos theta| < " << cfg.cosThetaMax
      << ", |d0| < " << cfg.d0Max << " cm, |z0| < " << cfg.z0Max << " cm\n"
      << "  phi candidates  : kaon p > " << cfg.kaonPMin << " GeV, pair p > " << cfg.phiPMin
      << " GeV, fit range " << cfg.fitLow << "-" << cfg.fitHigh << " GeV, Gamma fixed to PDG\n\n";

  out << "Events and tracks\n";
  AddSampleSummary(out, data);
  AddSampleSummary(out, mc);
  if (cfg.useGen) AddSampleSummary(out, gen);
  if (data.nEventsSelected > 0 && mc.nEventsSelected > 0) {
    double tracksPerEventData = data.sumGoodTracks / data.nEventsSelected;
    double tracksPerEventMC = mc.sumGoodTracks / mc.nEventsSelected;
    out << Form("  MC / data ratio of good tracks per event : %.4f\n", tracksPerEventMC / tracksPerEventData);
  }

  out << "\nphi -> K+K- (opposite-sign pairs, kaon mass hypothesis)\n";
  AddPhiSummary(out, data, fitData);
  AddPhiSummary(out, mc, fitMC);
  if (data.nEventsSelected > 0 && mc.nEventsSelected > 0 && fitData.yield > 0 && fitMC.yield > 0) {
    double perEventData = fitData.yield / data.nEventsSelected;
    double perEventMC = fitMC.yield / mc.nEventsSelected;
    double ratio = perEventMC / perEventData;
    // Relative errors add in quadrature for a ratio of independent numbers.
    double relativeError = std::sqrt(std::pow(fitData.yieldError / fitData.yield, 2) +
                                     std::pow(fitMC.yieldError / fitMC.yield, 2));
    out << Form("  MC / data ratio of phi per event : %.4f +- %.4f\n", ratio, ratio * relativeError);
    out << Form("  mass shift MC - data             : %.3f MeV\n", 1000.0 * (fitMC.mass - fitData.mass));
    out << Form("  resolution ratio MC / data       : %.3f\n", fitData.sigma > 0 ? fitMC.sigma / fitData.sigma : 0.0);
  }
  return out.str();
}

}  // namespace

// =============================================================================
//  Part 8: main()
// =============================================================================

void PrintUsage() {
  std::cout
      << "Usage: CompareDataMC --data <data.root> --mc <mc.root> [options]\n"
         "\n"
         "Required:\n"
         "  --data FILE          converted OPAL 1994 data file (tree t)\n"
         "  --mc FILE            converted OPAL 1994 MC file (trees t, tgen)\n"
         "Output:\n"
         "  --output FILE        output ROOT file [OPAL_1994_DataMC_validation.root]\n"
         "  --plotDir DIR        directory for the PDFs and summary.txt [.]\n"
         "  --maxEvents N        read at most N events per tree [-1 = all]\n"
         "  --gen BOOL           overlay generator-level MC (tgen) on the pT plots [true]\n"
         "Event selection:\n"
         "  --requireHIMS BOOL   require passesOPALHIMSTrack [true]\n"
         "  --nGoodTracksMin N   minimum number of good charged tracks [5]\n"
         "Track selection (defaults = the cuts already inside highPurity):\n"
         "  --trackPtMin X       [0.10 GeV]     --trackPMin X   [0 GeV]\n"
         "  --cosThetaMax X      [0.966]        --d0Max X [2.5 cm]   --z0Max X [50 cm]\n"
         "phi -> K+K- candidates:\n"
         "  --kaonPMin X         minimum momentum of each kaon candidate [0 GeV]\n"
         "  --phiPMin X          minimum momentum of the pair [0 GeV]\n"
         "  --massLow/--massHigh/--massBinWidth   mass histogram [0.98 / 1.10 / 0.001 GeV]\n"
         "  --fitLow/--fitHigh   fit range [0.990 / 1.090 GeV]\n"
         "pT normalisation:\n"
         "  --normalize MODE     event = (1/N_ev) dN/dpT,  area = unit area [event]\n"
      << std::endl;
}

int main(int argc, char** argv) {
  // ---- 1. read the command line ---------------------------------------------
  CommandLine cl(argc, argv);
  if (cl.Has("help") || !cl.Has("data") || !cl.Has("mc")) {
    PrintUsage();
    return cl.Has("help") ? 0 : 1;
  }

  Config cfg;
  cfg.dataFile = cl.GetString("data", "");
  cfg.mcFile = cl.GetString("mc", "");
  cfg.outputFile = cl.GetString("output", "OPAL_1994_DataMC_validation.root");
  cfg.plotDir = cl.GetString("plotDir", ".");
  cfg.maxEvents = cl.GetLong("maxEvents", -1);
  cfg.useGen = cl.GetBool("gen", true);
  cfg.requireHIMS = cl.GetBool("requireHIMS", true);
  cfg.nGoodTracksMin = cl.GetInt("nGoodTracksMin", 5);
  cfg.trackPtMin = cl.GetDouble("trackPtMin", 0.10);
  cfg.trackPMin = cl.GetDouble("trackPMin", 0.0);
  cfg.cosThetaMax = cl.GetDouble("cosThetaMax", 0.966);
  cfg.d0Max = cl.GetDouble("d0Max", 2.5);
  cfg.z0Max = cl.GetDouble("z0Max", 50.0);
  cfg.kaonPMin = cl.GetDouble("kaonPMin", 0.0);
  cfg.phiPMin = cl.GetDouble("phiPMin", 0.0);
  cfg.massLow = cl.GetDouble("massLow", 0.98);
  cfg.massHigh = cl.GetDouble("massHigh", 1.10);
  cfg.massBinWidth = cl.GetDouble("massBinWidth", 0.001);
  cfg.fitLow = cl.GetDouble("fitLow", 0.990);
  cfg.fitHigh = cl.GetDouble("fitHigh", 1.090);
  cfg.normalize = cl.GetString("normalize", "event");
  if (cfg.normalize != "event" && cfg.normalize != "area") {
    std::cerr << "--normalize must be 'event' or 'area'" << std::endl;
    return 1;
  }

  // ---- 2. global ROOT settings ----------------------------------------------
  gROOT->SetBatch(true);                // no graphics windows, just files
  gErrorIgnoreLevel = kWarning;         // hide ROOT's "Info" messages
  gStyle->SetOptStat(0);                // no statistics box on histograms
  gStyle->SetOptTitle(0);               // no histogram title above the plot
  gStyle->SetPadTickX(1);               // tick marks on all four sides
  gStyle->SetPadTickY(1);
  TH1::SetDefaultSumw2(true);           // keep proper statistical errors when scaling
  TH1::AddDirectory(false);             // histograms are not tied to the input files

  // ---- 3. open the input files ----------------------------------------------
  std::unique_ptr<TFile> dataFile(TFile::Open(cfg.dataFile.c_str()));
  if (!dataFile || dataFile->IsZombie()) {
    std::cerr << "Cannot open data file " << cfg.dataFile << std::endl;
    return 1;
  }
  std::unique_ptr<TFile> mcFile(TFile::Open(cfg.mcFile.c_str()));
  if (!mcFile || mcFile->IsZombie()) {
    std::cerr << "Cannot open MC file " << cfg.mcFile << std::endl;
    return 1;
  }
  TTree* dataTree = dynamic_cast<TTree*>(dataFile->Get("t"));
  TTree* mcTree = dynamic_cast<TTree*>(mcFile->Get("t"));
  TTree* genTree = cfg.useGen ? dynamic_cast<TTree*>(mcFile->Get("tgen")) : nullptr;
  if (!dataTree || !mcTree) {
    std::cerr << "Tree 't' is missing in the data or MC file" << std::endl;
    return 1;
  }
  if (cfg.useGen && !genTree) {
    std::cerr << "Tree 'tgen' not found in the MC file; continuing without generator level" << std::endl;
    cfg.useGen = false;
  }

  // ---- 4. fill the histograms -----------------------------------------------
  Sample data = MakeSample("data", "OPAL 1994 data", cfg);
  Sample mc = MakeSample("mc", "JETSET 7.4 + GOPAL", cfg);
  Sample gen = MakeSample("gen", "JETSET 7.4 generator level", cfg);

  FillHistogramsFromTree(dataTree, cfg, false, data);
  FillHistogramsFromTree(mcTree, cfg, false, mc);
  if (cfg.useGen) FillHistogramsFromTree(genTree, cfg, true, gen);

  // ---- 5. fit the phi peaks -------------------------------------------------
  PhiFitResult fitData = FitPhiPeak(data.hMassOS, cfg, "data");
  PhiFitResult fitMC = FitPhiPeak(mc.hMassOS, cfg, "mc");

  // ---- 6. normalised copies for the comparison plots -------------------------
  TH1D* ptData = NormalizedCopy(data.hPt, data.nEventsSelected, cfg.normalize, "_norm");
  TH1D* ptMC = NormalizedCopy(mc.hPt, mc.nEventsSelected, cfg.normalize, "_norm");
  TH1D* ptGen = cfg.useGen ? NormalizedCopy(gen.hPt, gen.nEventsSelected, cfg.normalize, "_norm") : nullptr;

  TH1D* ptLogData = NormalizedCopy(data.hPtLog, data.nEventsSelected, cfg.normalize, "_norm");
  TH1D* ptLogMC = NormalizedCopy(mc.hPtLog, mc.nEventsSelected, cfg.normalize, "_norm");
  TH1D* ptLogGen = cfg.useGen ? NormalizedCopy(gen.hPtLog, gen.nEventsSelected, cfg.normalize, "_norm") : nullptr;

  // The multiplicity is always compared as a shape (fraction of events).
  TH1D* nGoodData = NormalizedCopy(data.hNGoodTracks, data.nEventsSelected, "area", "_norm");
  TH1D* nGoodMC = NormalizedCopy(mc.hNGoodTracks, mc.nEventsSelected, "area", "_norm");
  TH1D* nGoodGen = cfg.useGen ? NormalizedCopy(gen.hNGoodTracks, gen.nEventsSelected, "area", "_norm") : nullptr;
  for (TH1D* h : {nGoodData, nGoodMC, nGoodGen}) {
    if (h) h->GetYaxis()->SetTitle("Fraction of events");
  }

  TH1D* subtractedData = BackgroundSubtractedCopy(data.hMassOS, fitData, data.nEventsSelected, cfg, "hMassSubtracted_data");
  TH1D* subtractedMC = BackgroundSubtractedCopy(mc.hMassOS, fitMC, mc.nEventsSelected, cfg, "hMassSubtracted_mc");

  // ---- 7. draw ---------------------------------------------------------------
  std::string dir = cfg.plotDir.empty() ? "." : cfg.plotDir;
  gSystem->mkdir(dir.c_str(), true);

  TH1D* ratioPt = DrawSpectrumComparison(ptData, ptMC, ptGen, cfg, false, true,
                                         dir + "/pt_spectrum_linear.pdf", "ratioPt_dataOverMC");
  TH1D* ratioPtLog = DrawSpectrumComparison(ptLogData, ptLogMC, ptLogGen, cfg, true, true,
                                            dir + "/pt_spectrum_log.pdf", "ratioPtLog_dataOverMC");
  TH1D* ratioNGood = DrawSpectrumComparison(nGoodData, nGoodMC, nGoodGen, cfg, false, false,
                                            dir + "/ngood_tracks.pdf", "ratioNGoodTracks_dataOverMC");
  {
    TCanvas canvas("canvasPhi", "", 1400, 600);
    canvas.Divide(2, 1);                                 // data on the left, MC on the right
    DrawPhiFit(static_cast<TPad*>(canvas.cd(1)), data, fitData);
    DrawPhiFit(static_cast<TPad*>(canvas.cd(2)), mc, fitMC);
    canvas.SaveAs((dir + "/phi_kk_mass.pdf").c_str());
  }
  DrawSubtractedPeaks(subtractedData, subtractedMC, cfg, dir + "/phi_kk_subtracted.pdf");

  // ---- 8. summary text -------------------------------------------------------
  std::string summary = BuildSummary(cfg, data, mc, gen, fitData, fitMC);
  std::cout << "\n" << summary << std::endl;
  std::ofstream summaryFile(dir + "/summary.txt");
  summaryFile << summary;
  summaryFile.close();

  // ---- 9. save everything in a ROOT file ------------------------------------
  // Each object is written under its own name, so the plots can be remade
  // later (for example interactively in ROOT) without rerunning the event loop.
  std::unique_ptr<TFile> outputFile(TFile::Open(cfg.outputFile.c_str(), "RECREATE"));
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Cannot create output file " << cfg.outputFile << std::endl;
    return 1;
  }
  outputFile->cd();
  for (Sample* s : {&data, &mc, &gen}) {
    if (s == &gen && !cfg.useGen) continue;
    for (TH1D* h : {s->hNGoodTracks, s->hPt, s->hPtLog, s->hMassOS, s->hMassLS}) h->Write();
  }
  for (TH1D* h : {ptData, ptMC, ptGen, ptLogData, ptLogMC, ptLogGen, nGoodData, nGoodMC, nGoodGen,
                  ratioPt, ratioPtLog, ratioNGood, subtractedData, subtractedMC}) {
    if (h) h->Write();
  }
  for (const PhiFitResult* fit : {&fitData, &fitMC}) {
    if (fit->total) fit->total->Write();
    if (fit->signal) fit->signal->Write();
    if (fit->background) fit->background->Write();
  }
  TNamed summaryObject("summary", summary.c_str());
  summaryObject.Write();
  outputFile->Close();

  std::cout << "Wrote " << cfg.outputFile << " and the plots in " << dir << std::endl;
  return 0;
}
