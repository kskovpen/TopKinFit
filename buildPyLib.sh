#/bin/env bash

if [[ $(type -P python) == "" ]]; then
  echo "Python not found"
  exit 1
fi

PYTHONINC=$(python -c "from distutils.sysconfig import get_python_inc;print get_python_inc()")
BDIR=$(pwd)

swig -c++ -python python/kfit.i

g++ -fPIC -Wno-register -c -std=c++11 -Wno-deprecated-declarations -m64 python/kfit_wrap.cxx -I${PYTHONINC} -I$(root-config --cflags) -I${BDIR}

g++ kfit_wrap.o -shared -fPIC -O3 -pthread -m64 -std=c++11 -o _kfit.so -L${BDIR} -lKinFit

rm -f kfit_wrap.o kinfitDict.cxx_*
mv _kfit.so python/
