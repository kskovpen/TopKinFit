ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lMinuit

CXX           = gcc -fpic -O3
SOFLAGS       = -shared

CXXFLAGS       = $(ROOTCFLAGS)
INCLUDE_FLAGS  = 
LDLIBS         = $(ROOTLIBS)

EXE           = kinfit

INC 	      = include/kinfit.h include/TopTopLepLep.h include/TopTopLepHad.h include/TopLep.h

SRC	      = src/kinfit.cxx src/TopTopLepLep.cxx src/TopTopLepHad.cxx src/TopLep.cxx

OBJS          = kinfit.o TopTopLepLep.o TopTopLepHad.o TopLep.o

LIB           = libKinFit.so

all: 	      $(LIB)

$(LIB):	      $(INC) $(SRC)
	      @echo "####### Generating dictionary"
	      @rootcint -f kinfitDict.cxx -c -p $(CXXFLAGS) \
	      $(INCLUDE_FLAGS) -I. $(INC) include/LinkDef.h

	      @echo "####### Building library $(LIB)"
	      @$(CXX) $(SOFLAGS) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDE_FLAGS) -I. $(SRC) \
	      kinfitDict.cxx -o $(LIB) $(ROOTLIBS)
	      
	      @echo  "####### Removing generated dictionary"
	      @rm -f kinfitDict.cxx kinfitDict.h
	      @rm -f *.o

clean:
	      @rm -f $(OBJS) $(EXE) kinfitDict.cxx kinfitDict.h $(LIB)
