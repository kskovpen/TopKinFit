## About

The TopPhit implementation is a kinematic fit that reconstructs the neutrino momentum in top quark decays. The algorithm can be applied to various event topologies, including top quark pair production, as well as singly-produced top quarks. If specified, the photon radiation from the hard-process particles is included in the kinematic fit.

## Installation

```
git clone https://github.com/kskovpen/TopPhit
cd TopPhit
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

The transfer functions (**TFs**) are defined as the relative difference between reconstructed and generated kinematic variables for leptons, jets, transverse missing energy and photons. The probability distribution functions (**PDFs**) have to be defined for the top quark and W boson mass distributions. The exact definition is (gen-rec)/gen. The TFs and PDFs can be provided as **ROOT TH1** histograms in a single input file.

```c++
kf->SetPDF([internal PDF/TF name], [input file name], [input file histogram name]);
```
The first argument corresponds to one of the internal definitions of TFs. The possible choices are: 
- **TopMass**: PDF, top quark mass distribution (for leptonic W boson decay)
- **TopWMass**: PDF, W boson mass distribution (for leptonic decay)
- **TopHadMass**: PDF, top quark mass distribution (for hadronic W boson decay)
- **TopWHadMass**: PDF, W boson mass distribution (for hadronic decay)
- **MetPx**: TF, x component of transverse missing energy
- **MetPy**: TF, y component of transverse missing energy
- **BJetPx**: TF, x component of b jet momentum
- **BJetPy**: TF, y component of b jet momentum
- **BJetPz**: TF, z component of b jet momentum
- **BJetPt**: TF, pT of b jet
- **BJetEta**: TF, eta of b jet
- **BJetPhi**: TF, phi of b jet
- **NonBJetPx**: TF, x component of non-b jet momentum
- **NonBJetPy**: TF, y component of non-b jet momentum
- **NonBJetPz**: TF, z component of non-b jet momentum
- **NonBJetPt**: TF, pT of non-b jet
- **NonBJetEta**: TF, eta of non-b jet
- **NonBJetPhi**: TF, phi of non-b jet
- **ElecPx**: TF, x component of electron momentum
- **ElecPy**: TF, y component of electron momentum
- **ElecPz**: TF, z component of electron momentum
- **ElecPt**: TF, pT of electron
- **ElecEta**: TF, eta of electron
- **ElecPhi**: TF, phi of electon
- **MuonPx**: TF, x component of muon momentum
- **MuonPy**: TF, y component of muon momentum
- **MuonPz**: TF, z component of muon momentum
- **MuonPt**: TF, pT of muon
- **MuonEta**: TF, eta of muon
- **MuonPhi**: TF, phi of muon
- **PhotonPx**: TF, x component of photon momentum
- **PhotonPy**: TF, y component of photon momentum
- **PhotonPz**: TF, z component of photon momentum
- **PhotonPt**: TF, pT of photon
- **PhotonEta**: TF, eta of photon
- **PhotonPhi**: TF, phi of photon

In the case of lepton, jet and photon objects, the TFs can be specified either for (Px, Py, Pz) or (Pt, Eta, Phi).

## Input data

Specify input collections. Each input variable corresponds to std::vector<float>. Floats are used in case of missing ET. If SetPhoton() method is called, photons are included in the kinematic fit with corresponding permutations.

```c++
kf->SetBJet(Pt, Eta, Phi, E);
kf->SetNonBJet(Pt, Eta, Phi, E);
kf->SetElectron(Pt, Eta, Phi, E);
kf->SetMuon(Pt, Eta, Phi, E);
kf->SetMet(Pt, Eta, Phi, E);
kf->SetPhoton(Pt, Eta, Phi, E);
```

## Run options

Run the kinematic fit:
  
```c++
kf->Run();
```
  
Set the number of toys in the generic minimization algorithm:

```c++
kf->SetNToy(10000);
```

Set the maximum number of toys to retain from generic minimization:

```c++
kf->SetNGrid(50);
```

Run MINUIT fits (default = true for TOPTOPLEPLEP, default = false for TOPLEP and TOPTOPLEPHAD).

```c++
kf->DoFit(true or false);
```

Set the maximum number of MINUIT fits. By default, this option is only relevant for TOPTOPLEPLEP. If DoFit() is set to True, the MINUIT fits can also be included in the TOPLEP and TOPTOPLEPHAD hypotheses evaluation.

```c++
kf->SetNFitMax(50);
```

Only retain the variations that correspond to the minimized log-likelihood value below the specified cut-off in generic minimization and fits.

```c++
kf->SetLHMaxGeneric(20.);
kf->SetLHMaxMinuit(0.01);
```

## Output results
  
Get the minimized log-likelihood value for a specific permutation. Permutations are sorted according to the minimized likelihood value, and the first permutation (N = 0) in the list therefore corresponds to the best fit result.
  
```c++
kf->GetDisc(int N);
```

Get the number of permutations in event:

```c++
kf->GetNPerm();
```
  
Get object's momenta for the Nth permutation. The second argument corresponds to the object's index (0 for TOPLEP and TOPTOPLEPHAD, 0 or 1 for TOPTOPLEPLEP).

```c++
kf->GetNuPx(int N, int i);
kf->GetNuPy(int N, int i);
kf->GetNuPz(int N, int i);
```
```c++
kf->GetTopPt(int N, int i);
kf->GetTopEta(int N, int i);
kf->GetTopPhi(int N, int i);
kf->GetTopPx(int N, int i);
kf->GetTopPy(int N, int i);
kf->GetTopPz(int N, int i);
kf->GetTopP(int N, int i);
kf->GetTopE(int N, int i);
kf->GetTopMass(int N, int i);
```
```c++
kf->GetWPt(int N, int i);
kf->GetWEta(int N, int i);
kf->GetWPhi(int N, int i);
kf->GetWPx(int N, int i);
kf->GetWPy(int N, int i);
kf->GetWPz(int N, int i);
kf->GetWP(int N, int i);
kf->GetWE(int N, int i);
kf->GetWMass(int N, int i);
```
```c++
kf->GetDrTopTop(int N);
kf->GetMTopTop(int N);
kf->GetPtTopTop(int N);
kf->GetPTopTop(int N);
kf->GetEtaTopTop(int N);
kf->GetPhiTopTop(int N);
```

Get index of an object in the input collection that is matched to a specific hard-process particle. The possible options are listed below:

- **ELECTRON1_TOPTOPLEPLEP**: Electron from the first top quark's W boson decay
- **MUON1_TOPTOPLEPLEP**: Muon from the first top quark's W boson decay
- **ELECTRON2_TOPTOPLEPLEP**: Electron from the second top quark's W boson decay
- **MUON2_TOPTOPLEPLEP**: Muon from the second top quark's W boson decay
- **BJET1_TOPTOPLEPLEP**: b jet from the first top quark decay
- **BJET2_TOPTOPLEPLEP**: b jet from the second top quark decay
- **PHOTON_TOPTOPLEPLEP**: radiated photon
  
- **ELECTRON_TOPTOPLEPHAD**: Electron from the leptonic top quark's W boson decay
- **MUON_TOPTOPLEPHAD**: Muon from the leptonic top quark's W boson decay
- **BJETLEP_TOPTOPLEPHAD**: b jet from the leptonic top quark decay
- **BJETHAD_TOPTOPLEPHAD**: b jet from the hadronic top quark decay
- **NONBJET1_TOPTOPLEPHAD**: first non b jet from the hadronic top quark decay
- **NONBJET2_TOPTOPLEPHAD**: second non b jet from the hadronic top quark decay
- **PHOTON_TOPTOPLEPHAD**: radiated photon
  
- **ELECTRON_TOPLEP**: Electron from the top quark's W boson decay
- **MUON_TOPLEP**: Muon from the top quark's W boson decay
- **BJETLEP_TOPLEP**: b jet from the top quark decay
- **PHOTON_TOPLEP**: radiated photon
  
```c++
kf->GetIndex(objType, int N);
```
  
## Photon origin
  
Possible returned values are mentioned below.
  
- **PHOTON_FROM_TOP1_COMB_TOPTOPLEPLEP**: first top quark decay (from top or b quark)
- **PHOTON_FROM_W1_COMB_TOPTOPLEPLEP**: W boson from the first top quark decay (from W boson or lepton)
- **PHOTON_FROM_TOP2_COMB_TOPTOPLEPLEP**: second top quark decay (from top or b quark)
- **PHOTON_FROM_W2_COMB_TOPTOPLEPLEP**: W boson from the second top quark decay (from W boson or lepton)

- **PHOTON_FROM_TOPLEP_COMB_TOPTOPLEPHAD**: leptonic top quark decay (from top or b quark)
- **PHOTON_FROM_WLEP_COMB_TOPTOPLEPHAD**: W boson from the leptonic top quark decay (from W boson or lepton)
- **PHOTON_FROM_TOPHAD_COMB_TOPTOPLEPHAD**: hadronic top quark decay (from top or b quark)
- **PHOTON_FROM_WHAD_COMB_TOPTOPLEPHAD**: W boson from the hadronic top quark decay (from W boson or quarks)
  
- **PHOTON_FROM_TOPLEP_COMB_TOPLEP**: top quark decay (from top or b quark)
- **PHOTON_FROM_WLEP_COMB_TOPLEP**: W boson from the top quark decay (from W boson or lepton)
  
```c++
kf->GetPhotonOrigin(int i);
```
