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
git clone
cd TopKinFit
make
./buildPyLib.sh
```

The latter command is used to build the Python wrapper (if needed).
