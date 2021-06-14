## About

The TopPhit implementation is a kinematic fit that reconstructs the neutrino momentum in top quark decays. The algorithm can be applied to various event topologies, including top quark pair production, as well as singly-produced top quarks. If specified, the photon radiation from the hard-process particles is included in the kinematic fit.

## Installation

```
git clone https://github.com/kskovpen/TopKinFit
cd TopKinFit
make
./buildPyLib.sh
```
## Main tool

An instance of the main class can be added with:
```c++
#include "include/kinfit.h"
KINFIT::kfit *kf = new KINFIT::kfit();
```

## Signal hypothesis

The reconstruction algorithm is implemented for:
- **TOPTOPLEPHAD**: semileptonic ttbar
- **TOPTOPLEPLEP**: dilepton ttbar
- **TOPLEP**: single-lepton top quark decay

and is specified with:
```c++
kf->Init(HYPO::TOPTOPLEPHAD);
```

The decaying top quarks and W bosons can also be defined for non standard model particles, provided the decay chain is the same.

## Transfer functions

The transfer functions (**TFs**) and probability distribution functions (**PDFs**) are defined as the relative difference between reconstructed and generated kinematic variables for leptons, jets, transverse missing energy and photons. The input PDFs correspond to the top quark and W boson mass distributions. The exact definition is (gen-rec)/gen. The TFs and PDFs can be provided as **ROOT TH1** histograms in a single input file.

```c++
kf->SetPDF([internal TF name], [input file name], [input file histogram name]);
```
The first argument corresponds to one of the internal definitions of TFs. The possible choices are: 
- **TopMass**: top quark mass distribution (for leptonic W boson decay)
- **TopWMass**: W boson mass distribution (for leptonic decay)
- **TopHadMass**: top quark mass distribution (for hadronic W boson decay)
- **TopWHadMass**: W boson mass distribution (for hadronic decay)
- **MetPx**: x component of transverse missing energy
- **MetPy**: y component of transverse missing energy
- **BJetPx**: x component of b jet momentum
- **BJetPy**: y component of b jet momentum
- **BJetPz**: z component of b jet momentum
- **BJetPt**: pT of b jet
- **BJetEta**: eta of b jet
- **BJetPhi**: phi of b jet
- **NonBJetPx**: x component of non-b jet momentum
- **NonBJetPy**: y component of non-b jet momentum
- **NonBJetPz**: z component of non-b jet momentum
- **NonBJetPt**: pT of non-b jet
- **NonBJetEta**: eta of non-b jet
- **NonBJetPhi**: phi of non-b jet
- **ElecPx**: x component of electron momentum
- **ElecPy**: y component of electron momentum
- **ElecPz**: z component of electron momentum
- **ElecPt**: pT of electron
- **ElecEta**: eta of electron
- **ElecPhi**: phi of electon
- **MuonPx**: x component of muon momentum
- **MuonPy**: y component of muon momentum
- **MuonPz**: z component of muon momentum
- **MuonPt**: pT of muon
- **MuonEta**: eta of muon
- **MuonPhi**: phi of muon
- **PhotonPx**: x component of photon momentum
- **PhotonPy**: y component of photon momentum
- **PhotonPz**: z component of photon momentum
- **PhotonPt**: pT of photon
- **PhotonEta**: eta of photon
- **PhotonPhi**: phi of photon

In the case of lepton, jet and photon objects, the TFs can be specified either for (Px, Py, Pz) or (Pt, Eta, Phi).
