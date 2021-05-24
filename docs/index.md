## Installation

```
git clone https://github.com/kskovpen/TopKinFit
cd TopKinFit
make
./buildPyLib.sh
```
## Main tool

An instance of the main class can be added with:
```
#include "include/kinfit.h"
KINFIT::kfit *kf = new KINFIT::kfit();
```

## Transfer functions

The transfer functions (TFs) and probability distribution functions (PDFs) are defined as the relative difference between reconstructed and generated kinematic variables for leptons, jets, transverse missing energy and photons. The input PDFs corresponds to the top quark and W boson mass distributions. The exact definition is (gen-rec)/gen. The TFs and PDFs can be provided as ROOT TH1 histograms in a single input file.

```c++
kf->SetPDF([internal TF name], [input file name], [input file histogram name]);
```
The first argument corresponds to one of the internal definitions of TFs. Possible choices are: 
- TopMass, PDF, top quark mass distribution (for leptonic W boson decay)
- TopWMass, PDF, W boson mass distribution (for leptonic decay)
- TopHadMass, PDF, top quark mass distribution (for hadronic W boson decay)
- TopWHadMass, PDF, W boson mass distribution (for hadronic decay)
- MetPx, TF, x component of transverse missing energy
- MetPy, TF, y component of transverse missing energy
- BJetPx, TF, x component of b jet momentum
- BJetPy, TF, y component of b jet momentum
- BJetPz, TF, z component of b jet momentum
- NonBJetPx, TF, x component of non-b jet momentum
- NonBJetPy, TF, y component of non-b jet momentum
- NonBJetPz, TF, z component of non-b jet momentum
- ElecPx, TF, x component of electron momentum
- ElecPy, TF, y component of electron momentum
- ElecPz, TF, z component of electron momentum
- MuonPx, TF, x component of muon momentum
- MuonPy, TF, y component of muon momentum
- MuonPz, TF, z component of muon momentum
- PhotonPx, TF, x component of photon momentum
- PhotonPy, TF, y component of photon momentum
- PhotonPz, TF, z component of photon momentum
