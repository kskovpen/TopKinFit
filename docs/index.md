### Installation

```
git clone https://github.com/kskovpen/TopKinFit
cd TopKinFit
make
./buildPyLib.sh
```

### Transfer functions

The transfer functions (TFs) are defined as a relative difference between reconstructed and generated kinematic variables for leptons, jets, transverse missing energy and photons. The exact definition is (gen-rec)/gen. The TFs can be provided as ROOT TH1 histograms in a single input file.
