# TopKinFit

![C++ version](https://img.shields.io/badge/C++-11-blue.svg)

Kinematic reconstruction of processes with top quark decays includes
the following signal hypotheses:

- Top quark pair production in semileptonic channel
- Top quark pair production in dilepton channel
- Single top quark pair production in single-lepton channel

The kinematic fit can also account for the photon radiation from any
charged hard-process particle in the aforementioned processes.
Implementation is realized as a C++ shared library. A Python wrapper
is also provided.

The full API documentation is available at [this page](io).

Installation:

```
git clone https://github.com/kskovpen/TopKinFit
cd TopKinFit
make
./buildPyLib.sh
```

The latter command is used to build the Python wrapper (if needed).

How to use (in C++):

```c++
#include "include/kinfit.h"

KINFIT::kfit *kf = new KINFIT::kfit();

// Specify ROOT file with transfer function histograms
std::string pdfFileName = "pdf.root";

// Specify transfer function histogram names
kf->SetPDF("TopMass", pdfFileName, "top_mass_tf");
kf->SetPDF("TopWMass", pdfFileName, "w_mass_tf");
kf->SetPDF("TopHadMass", pdfFileName, "top_mass_tf");
kf->SetPDF("TopWHadMass", pdfFileName, "w_mass_tf");
kf->SetPDF("MetPx", pdfFileName, "emiss_px_tf");
kf->SetPDF("MetPy", pdfFileName, "emiss_py_tf");
kf->SetPDF("BJetPx", pdfFileName, "bjet_px_tf");
kf->SetPDF("BJetPy", pdfFileName, "bjet_py_tf");
kf->SetPDF("BJetPz", pdfFileName, "bjet_pz_tf");
kf->SetPDF("NonBJetPx", pdfFileName, "ljet_px_tf");
kf->SetPDF("NonBJetPy", pdfFileName, "ljet_py_tf");
kf->SetPDF("NonBJetPz", pdfFileName, "ljet_pz_tf");
kf->SetPDF("ElecPx", pdfFileName, "elec_px_tf");
kf->SetPDF("ElecPy", pdfFileName, "elec_py_tf");
kf->SetPDF("ElecPz", pdfFileName, "elec_pz_tf");
kf->SetPDF("MuonPx", pdfFileName, "muon_px_tf");
kf->SetPDF("MuonPy", pdfFileName, "muon_py_tf");
kf->SetPDF("MuonPz", pdfFileName, "muon_pz_tf");

// Set signal hypothesis (ttbar semileptonic)
kf->Init(HYPO::TOPTOPLEPHAD);

// Pass the four-momenta of reconstructed objects in event
kf->SetBJet(BJetPt, BJetEta, BJetPhi, BJetE);
kf->SetNonBJet(NonBJetPt, NonBJetEta, NonBJetPhi, NonBJetE);
kf->SetElectron(ElectronPt, ElectronEta, ElectronPhi, ElectronE, ElectronCharge);
kf->SetMuon(MuonPt, MuonEta, MuonPhi, MuonE, MuonCharge);
kf->SetMet(MetPx, MetPy);

// One can specify reconstructed photons (not required)
kf->SetPhoton(PhotonPt, PhotonEta, PhotonPhi, PhotonE);

// Run the tool
kf->Run();

// Get the minimized log-likelihood value corresponding to the
best permutation
float disc = kf->GetDisc(0);

// Get the reconstructed top quark (with leptonic W boson decay) pT
float top1Pt = kf->GetTopPt(0, 0);

// Get the reconstructed top quark (with hadronic W boson decay) pT
float top2Pt = kf->GetTopPt(0, 1);

// Get the minimized log-likelihood value for the second-best permutation
float disc = kf->GetDisc(1);

// Total number of permutations
int nperm = kf->GetNPerm();

// Get photon origin in the best permutation
int origin = kf->GetPhotonOrigin(0);

```
