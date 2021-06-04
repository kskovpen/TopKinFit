/// KINFIT class implementation

#include "../include/kinfit.h"

#include "../include/TopTopLepLep.h"
#include "../include/TopTopLepHad.h"
#include "../include/TopLep.h"

#include <boost/bind.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

ClassImp(KINFIT::kfit)

const int NNUMAX = 2;    ///< Maximum number of neutrinos
int FPARAM_N;

KINFIT::TopTopLepLep *hypTopTopLepLep;
KINFIT::TopTopLepHad *hypTopTopLepHad;
KINFIT::TopLep *hypTopLep;

//bool KINFIT::kfit::AddTopTopMassToNLL_;
bool KINFIT::kfit::SortByTopTopMass_;
bool KINFIT::kfit::SortByTopTopPz_;
bool KINFIT::kfit::SortByTopTopPzAbs_;
bool KINFIT::kfit::SortByNeutrinoPz_;
bool KINFIT::kfit::SortByNeutrinoPt_;
bool KINFIT::kfit::SortByNeutrinoP_;

bool KINFIT::kfit::IncludePhotons_;
bool KINFIT::kfit::CheckAllPhotonOrigins_;

std::shared_ptr<double> KINFIT::kfit::CHISQ;

std::shared_ptr<double> KINFIT::kfit::EtMissX;
std::shared_ptr<double> KINFIT::kfit::EtMissY;

std::shared_ptr<double> KINFIT::kfit::PxLepton1;
std::shared_ptr<double> KINFIT::kfit::PyLepton1;
std::shared_ptr<double> KINFIT::kfit::PzLepton1;
std::shared_ptr<double> KINFIT::kfit::ELepton1;
std::shared_ptr<double> KINFIT::kfit::PtLepton1;
std::shared_ptr<double> KINFIT::kfit::EtaLepton1;
std::shared_ptr<double> KINFIT::kfit::PhiLepton1;
std::shared_ptr<double> KINFIT::kfit::MassLepton1;
std::shared_ptr<int> KINFIT::kfit::LabelLepton1;

std::shared_ptr<double> KINFIT::kfit::PxLepton2;
std::shared_ptr<double> KINFIT::kfit::PyLepton2;
std::shared_ptr<double> KINFIT::kfit::PzLepton2;
std::shared_ptr<double> KINFIT::kfit::ELepton2;
std::shared_ptr<double> KINFIT::kfit::PtLepton2;
std::shared_ptr<double> KINFIT::kfit::EtaLepton2;
std::shared_ptr<double> KINFIT::kfit::PhiLepton2;
std::shared_ptr<double> KINFIT::kfit::MassLepton2;
std::shared_ptr<int> KINFIT::kfit::LabelLepton2;

std::shared_ptr<double> KINFIT::kfit::PxBJet1;
std::shared_ptr<double> KINFIT::kfit::PyBJet1;
std::shared_ptr<double> KINFIT::kfit::PzBJet1;
std::shared_ptr<double> KINFIT::kfit::EBJet1;
std::shared_ptr<double> KINFIT::kfit::PtBJet1;
std::shared_ptr<double> KINFIT::kfit::EtaBJet1;
std::shared_ptr<double> KINFIT::kfit::PhiBJet1;
std::shared_ptr<double> KINFIT::kfit::MassBJet1;

std::shared_ptr<double> KINFIT::kfit::PxBJet2;
std::shared_ptr<double> KINFIT::kfit::PyBJet2;
std::shared_ptr<double> KINFIT::kfit::PzBJet2;
std::shared_ptr<double> KINFIT::kfit::EBJet2;
std::shared_ptr<double> KINFIT::kfit::PtBJet2;
std::shared_ptr<double> KINFIT::kfit::EtaBJet2;
std::shared_ptr<double> KINFIT::kfit::PhiBJet2;
std::shared_ptr<double> KINFIT::kfit::MassBJet2;

std::shared_ptr<double> KINFIT::kfit::PxNonBJet1;
std::shared_ptr<double> KINFIT::kfit::PyNonBJet1;
std::shared_ptr<double> KINFIT::kfit::PzNonBJet1;
std::shared_ptr<double> KINFIT::kfit::ENonBJet1;
std::shared_ptr<double> KINFIT::kfit::PtNonBJet1;
std::shared_ptr<double> KINFIT::kfit::EtaNonBJet1;
std::shared_ptr<double> KINFIT::kfit::PhiNonBJet1;
std::shared_ptr<double> KINFIT::kfit::MassNonBJet1;

std::shared_ptr<double> KINFIT::kfit::PxNonBJet2;
std::shared_ptr<double> KINFIT::kfit::PyNonBJet2;
std::shared_ptr<double> KINFIT::kfit::PzNonBJet2;
std::shared_ptr<double> KINFIT::kfit::ENonBJet2;
std::shared_ptr<double> KINFIT::kfit::PtNonBJet2;
std::shared_ptr<double> KINFIT::kfit::EtaNonBJet2;
std::shared_ptr<double> KINFIT::kfit::PhiNonBJet2;
std::shared_ptr<double> KINFIT::kfit::MassNonBJet2;

std::shared_ptr<double> KINFIT::kfit::PxPhoton;
std::shared_ptr<double> KINFIT::kfit::PyPhoton;
std::shared_ptr<double> KINFIT::kfit::PzPhoton;
std::shared_ptr<double> KINFIT::kfit::EPhoton;
std::shared_ptr<double> KINFIT::kfit::PtPhoton;
std::shared_ptr<double> KINFIT::kfit::EtaPhoton;
std::shared_ptr<double> KINFIT::kfit::PhiPhoton;
std::shared_ptr<int> KINFIT::kfit::PhotonOrigin;

std::shared_ptr<double> KINFIT::kfit::PxNu1;
std::shared_ptr<double> KINFIT::kfit::PyNu1;
std::shared_ptr<double> KINFIT::kfit::PzNu1;

std::shared_ptr<double> KINFIT::kfit::PxNu2;
std::shared_ptr<double> KINFIT::kfit::PyNu2;
std::shared_ptr<double> KINFIT::kfit::PzNu2;

std::shared_ptr<double> KINFIT::kfit::TopTopMass;
std::shared_ptr<double> KINFIT::kfit::TopTopPz;
std::shared_ptr<double> KINFIT::kfit::TopTopPzAbs;

std::shared_ptr<std::vector<double> > KINFIT::kfit::FitParam;
std::shared_ptr<std::vector<double> > KINFIT::kfit::ChiTerm;
std::shared_ptr<std::vector<std::string> > KINFIT::kfit::ChiTermName;
std::shared_ptr<std::vector<bool> > KINFIT::kfit::IsParFixed;

std::shared_ptr<std::vector<float> > KINFIT::kfit::ParMin;
std::shared_ptr<std::vector<float> > KINFIT::kfit::ParMax;

std::shared_ptr<std::map<int, bool> > KINFIT::kfit::foundPDF;

std::shared_ptr<bool> KINFIT::kfit::usePDFBJetPxPyPz;
std::shared_ptr<bool> KINFIT::kfit::usePDFBJetPtEtaPhi;
std::shared_ptr<bool> KINFIT::kfit::usePDFNonBJetPxPyPz;
std::shared_ptr<bool> KINFIT::kfit::usePDFNonBJetPtEtaPhi;
std::shared_ptr<bool> KINFIT::kfit::usePDFLeptonPxPyPz;
std::shared_ptr<bool> KINFIT::kfit::usePDFLeptonPtEtaPhi;
std::shared_ptr<bool> KINFIT::kfit::usePDFPhotonPxPyPz;
std::shared_ptr<bool> KINFIT::kfit::usePDFPhotonPtEtaPhi;

std::shared_ptr<TF1> KINFIT::kfit::hPDFTopWMass;
std::shared_ptr<TF1> KINFIT::kfit::hPDFTopMass;
std::shared_ptr<TF1> KINFIT::kfit::hPDFTopWHadMass;
std::shared_ptr<TF1> KINFIT::kfit::hPDFTopHadMass;
//std::shared_ptr<TF1> KINFIT::kfit::hPDFTopTopMass;

std::shared_ptr<TF1> KINFIT::kfit::hPDFMetPx;
std::shared_ptr<TF1> KINFIT::kfit::hPDFMetPy;

std::shared_ptr<TF1> KINFIT::kfit::hPDFBJetPx;
std::shared_ptr<TF1> KINFIT::kfit::hPDFBJetPy;
std::shared_ptr<TF1> KINFIT::kfit::hPDFBJetPz;
std::shared_ptr<TF1> KINFIT::kfit::hPDFBJetPt;
std::shared_ptr<TF1> KINFIT::kfit::hPDFBJetEta;
std::shared_ptr<TF1> KINFIT::kfit::hPDFBJetPhi;

std::shared_ptr<TF1> KINFIT::kfit::hPDFElecPx;
std::shared_ptr<TF1> KINFIT::kfit::hPDFElecPy;
std::shared_ptr<TF1> KINFIT::kfit::hPDFElecPz;
std::shared_ptr<TF1> KINFIT::kfit::hPDFElecPt;
std::shared_ptr<TF1> KINFIT::kfit::hPDFElecEta;
std::shared_ptr<TF1> KINFIT::kfit::hPDFElecPhi;

std::shared_ptr<TF1> KINFIT::kfit::hPDFMuonPx;
std::shared_ptr<TF1> KINFIT::kfit::hPDFMuonPy;
std::shared_ptr<TF1> KINFIT::kfit::hPDFMuonPz;
std::shared_ptr<TF1> KINFIT::kfit::hPDFMuonPt;
std::shared_ptr<TF1> KINFIT::kfit::hPDFMuonEta;
std::shared_ptr<TF1> KINFIT::kfit::hPDFMuonPhi;

std::shared_ptr<TF1> KINFIT::kfit::hPDFNonBJetPx;
std::shared_ptr<TF1> KINFIT::kfit::hPDFNonBJetPy;
std::shared_ptr<TF1> KINFIT::kfit::hPDFNonBJetPz;
std::shared_ptr<TF1> KINFIT::kfit::hPDFNonBJetPt;
std::shared_ptr<TF1> KINFIT::kfit::hPDFNonBJetEta;
std::shared_ptr<TF1> KINFIT::kfit::hPDFNonBJetPhi;

std::shared_ptr<TF1> KINFIT::kfit::hPDFPhotonPx;
std::shared_ptr<TF1> KINFIT::kfit::hPDFPhotonPy;
std::shared_ptr<TF1> KINFIT::kfit::hPDFPhotonPz;
std::shared_ptr<TF1> KINFIT::kfit::hPDFPhotonPt;
std::shared_ptr<TF1> KINFIT::kfit::hPDFPhotonEta;
std::shared_ptr<TF1> KINFIT::kfit::hPDFPhotonPhi;

std::string KINFIT::kfit::FPARAM_NAME[NPARMAX];
std::string KINFIT::kfit::PDF_NAME[PDF_N];

bool KINFIT::kfit::isDeltaFuncPDFTopWMass;
bool KINFIT::kfit::isDeltaFuncPDFTopMass;
bool KINFIT::kfit::isDeltaFuncPDFTopWHadMass;
bool KINFIT::kfit::isDeltaFuncPDFTopHadMass;
bool KINFIT::kfit::isDeltaFuncPDFMetPx;
bool KINFIT::kfit::isDeltaFuncPDFMetPy;
bool KINFIT::kfit::isDeltaFuncPDFBJetPx;
bool KINFIT::kfit::isDeltaFuncPDFBJetPy;
bool KINFIT::kfit::isDeltaFuncPDFBJetPz;
bool KINFIT::kfit::isDeltaFuncPDFBJetPt;
bool KINFIT::kfit::isDeltaFuncPDFBJetEta;
bool KINFIT::kfit::isDeltaFuncPDFBJetPhi;
bool KINFIT::kfit::isDeltaFuncPDFNonBJetPx;
bool KINFIT::kfit::isDeltaFuncPDFNonBJetPy;
bool KINFIT::kfit::isDeltaFuncPDFNonBJetPz;
bool KINFIT::kfit::isDeltaFuncPDFNonBJetPt;
bool KINFIT::kfit::isDeltaFuncPDFNonBJetEta;
bool KINFIT::kfit::isDeltaFuncPDFNonBJetPhi;
bool KINFIT::kfit::isDeltaFuncPDFElecPx;
bool KINFIT::kfit::isDeltaFuncPDFElecPy;
bool KINFIT::kfit::isDeltaFuncPDFElecPz;
bool KINFIT::kfit::isDeltaFuncPDFElecPt;
bool KINFIT::kfit::isDeltaFuncPDFElecEta;
bool KINFIT::kfit::isDeltaFuncPDFElecPhi;
bool KINFIT::kfit::isDeltaFuncPDFMuonPx;
bool KINFIT::kfit::isDeltaFuncPDFMuonPy;
bool KINFIT::kfit::isDeltaFuncPDFMuonPz;
bool KINFIT::kfit::isDeltaFuncPDFMuonPt;
bool KINFIT::kfit::isDeltaFuncPDFMuonEta;
bool KINFIT::kfit::isDeltaFuncPDFMuonPhi;
bool KINFIT::kfit::isDeltaFuncPDFPhotonPx;
bool KINFIT::kfit::isDeltaFuncPDFPhotonPy;
bool KINFIT::kfit::isDeltaFuncPDFPhotonPz;
bool KINFIT::kfit::isDeltaFuncPDFPhotonPt;
bool KINFIT::kfit::isDeltaFuncPDFPhotonEta;
bool KINFIT::kfit::isDeltaFuncPDFPhotonPhi;

float KINFIT::kfit::maxPDFTopWMass;
float KINFIT::kfit::meanPDFTopWMass;
float KINFIT::kfit::sigmaPDFTopWMass;
double KINFIT::kfit::xminPDFTopWMass;
double KINFIT::kfit::xmaxPDFTopWMass;

float KINFIT::kfit::maxPDFTopMass;
float KINFIT::kfit::meanPDFTopMass;
float KINFIT::kfit::sigmaPDFTopMass;
double KINFIT::kfit::xminPDFTopMass;
double KINFIT::kfit::xmaxPDFTopMass;

float KINFIT::kfit::maxPDFTopWHadMass;
float KINFIT::kfit::meanPDFTopWHadMass;
float KINFIT::kfit::sigmaPDFTopWHadMass;
double KINFIT::kfit::xminPDFTopWHadMass;
double KINFIT::kfit::xmaxPDFTopWHadMass;

float KINFIT::kfit::maxPDFTopHadMass;
float KINFIT::kfit::meanPDFTopHadMass;
float KINFIT::kfit::sigmaPDFTopHadMass;
double KINFIT::kfit::xminPDFTopHadMass;
double KINFIT::kfit::xmaxPDFTopHadMass;

//float KINFIT::kfit::maxPDFTopTopMass;
//float KINFIT::kfit::meanPDFTopTopMass;
//float KINFIT::kfit::sigmaPDFTopTopMass;
//double KINFIT::kfit::xminPDFTopTopMass;
//double KINFIT::kfit::xmaxPDFTopTopMass;

float KINFIT::kfit::maxPDFBJetPx;
float KINFIT::kfit::meanPDFBJetPx;
float KINFIT::kfit::sigmaPDFBJetPx;
double KINFIT::kfit::xminPDFBJetPx;
double KINFIT::kfit::xmaxPDFBJetPx;

float KINFIT::kfit::maxPDFBJetPy;
float KINFIT::kfit::meanPDFBJetPy;
float KINFIT::kfit::sigmaPDFBJetPy;
double KINFIT::kfit::xminPDFBJetPy;
double KINFIT::kfit::xmaxPDFBJetPy;

float KINFIT::kfit::maxPDFBJetPz;
float KINFIT::kfit::meanPDFBJetPz;
float KINFIT::kfit::sigmaPDFBJetPz;
double KINFIT::kfit::xminPDFBJetPz;
double KINFIT::kfit::xmaxPDFBJetPz;

float KINFIT::kfit::maxPDFBJetPt;
float KINFIT::kfit::meanPDFBJetPt;
float KINFIT::kfit::sigmaPDFBJetPt;
double KINFIT::kfit::xminPDFBJetPt;
double KINFIT::kfit::xmaxPDFBJetPt;

float KINFIT::kfit::maxPDFBJetEta;
float KINFIT::kfit::meanPDFBJetEta;
float KINFIT::kfit::sigmaPDFBJetEta;
double KINFIT::kfit::xminPDFBJetEta;
double KINFIT::kfit::xmaxPDFBJetEta;

float KINFIT::kfit::maxPDFBJetPhi;
float KINFIT::kfit::meanPDFBJetPhi;
float KINFIT::kfit::sigmaPDFBJetPhi;
double KINFIT::kfit::xminPDFBJetPhi;
double KINFIT::kfit::xmaxPDFBJetPhi;

float KINFIT::kfit::maxPDFMetPx;
float KINFIT::kfit::meanPDFMetPx;
float KINFIT::kfit::sigmaPDFMetPx;
double KINFIT::kfit::xminPDFMetPx;
double KINFIT::kfit::xmaxPDFMetPx;

float KINFIT::kfit::maxPDFMetPy;
float KINFIT::kfit::meanPDFMetPy;
float KINFIT::kfit::sigmaPDFMetPy;
double KINFIT::kfit::xminPDFMetPy;
double KINFIT::kfit::xmaxPDFMetPy;

float KINFIT::kfit::maxPDFElecPx;
float KINFIT::kfit::meanPDFElecPx;
float KINFIT::kfit::sigmaPDFElecPx;
double KINFIT::kfit::xminPDFElecPx;
double KINFIT::kfit::xmaxPDFElecPx;

float KINFIT::kfit::maxPDFElecPy;
float KINFIT::kfit::meanPDFElecPy;
float KINFIT::kfit::sigmaPDFElecPy;
double KINFIT::kfit::xminPDFElecPy;
double KINFIT::kfit::xmaxPDFElecPy;

float KINFIT::kfit::maxPDFElecPz;
float KINFIT::kfit::meanPDFElecPz;
float KINFIT::kfit::sigmaPDFElecPz;
double KINFIT::kfit::xminPDFElecPz;
double KINFIT::kfit::xmaxPDFElecPz;

float KINFIT::kfit::maxPDFElecPt;
float KINFIT::kfit::meanPDFElecPt;
float KINFIT::kfit::sigmaPDFElecPt;
double KINFIT::kfit::xminPDFElecPt;
double KINFIT::kfit::xmaxPDFElecPt;

float KINFIT::kfit::maxPDFElecEta;
float KINFIT::kfit::meanPDFElecEta;
float KINFIT::kfit::sigmaPDFElecEta;
double KINFIT::kfit::xminPDFElecEta;
double KINFIT::kfit::xmaxPDFElecEta;

float KINFIT::kfit::maxPDFElecPhi;
float KINFIT::kfit::meanPDFElecPhi;
float KINFIT::kfit::sigmaPDFElecPhi;
double KINFIT::kfit::xminPDFElecPhi;
double KINFIT::kfit::xmaxPDFElecPhi;

float KINFIT::kfit::maxPDFMuonPx;
float KINFIT::kfit::meanPDFMuonPx;
float KINFIT::kfit::sigmaPDFMuonPx;
double KINFIT::kfit::xminPDFMuonPx;
double KINFIT::kfit::xmaxPDFMuonPx;

float KINFIT::kfit::maxPDFMuonPy;
float KINFIT::kfit::meanPDFMuonPy;
float KINFIT::kfit::sigmaPDFMuonPy;
double KINFIT::kfit::xminPDFMuonPy;
double KINFIT::kfit::xmaxPDFMuonPy;

float KINFIT::kfit::maxPDFMuonPz;
float KINFIT::kfit::meanPDFMuonPz;
float KINFIT::kfit::sigmaPDFMuonPz;
double KINFIT::kfit::xminPDFMuonPz;
double KINFIT::kfit::xmaxPDFMuonPz;

float KINFIT::kfit::maxPDFMuonPt;
float KINFIT::kfit::meanPDFMuonPt;
float KINFIT::kfit::sigmaPDFMuonPt;
double KINFIT::kfit::xminPDFMuonPt;
double KINFIT::kfit::xmaxPDFMuonPt;

float KINFIT::kfit::maxPDFMuonEta;
float KINFIT::kfit::meanPDFMuonEta;
float KINFIT::kfit::sigmaPDFMuonEta;
double KINFIT::kfit::xminPDFMuonEta;
double KINFIT::kfit::xmaxPDFMuonEta;

float KINFIT::kfit::maxPDFMuonPhi;
float KINFIT::kfit::meanPDFMuonPhi;
float KINFIT::kfit::sigmaPDFMuonPhi;
double KINFIT::kfit::xminPDFMuonPhi;
double KINFIT::kfit::xmaxPDFMuonPhi;

float KINFIT::kfit::maxPDFNonBJetPx;
float KINFIT::kfit::meanPDFNonBJetPx;
float KINFIT::kfit::sigmaPDFNonBJetPx;
double KINFIT::kfit::xminPDFNonBJetPx;
double KINFIT::kfit::xmaxPDFNonBJetPx;

float KINFIT::kfit::maxPDFNonBJetPy;
float KINFIT::kfit::meanPDFNonBJetPy;
float KINFIT::kfit::sigmaPDFNonBJetPy;
double KINFIT::kfit::xminPDFNonBJetPy;
double KINFIT::kfit::xmaxPDFNonBJetPy;

float KINFIT::kfit::maxPDFNonBJetPz;
float KINFIT::kfit::meanPDFNonBJetPz;
float KINFIT::kfit::sigmaPDFNonBJetPz;
double KINFIT::kfit::xminPDFNonBJetPz;
double KINFIT::kfit::xmaxPDFNonBJetPz;

float KINFIT::kfit::maxPDFNonBJetPt;
float KINFIT::kfit::meanPDFNonBJetPt;
float KINFIT::kfit::sigmaPDFNonBJetPt;
double KINFIT::kfit::xminPDFNonBJetPt;
double KINFIT::kfit::xmaxPDFNonBJetPt;

float KINFIT::kfit::maxPDFNonBJetEta;
float KINFIT::kfit::meanPDFNonBJetEta;
float KINFIT::kfit::sigmaPDFNonBJetEta;
double KINFIT::kfit::xminPDFNonBJetEta;
double KINFIT::kfit::xmaxPDFNonBJetEta;

float KINFIT::kfit::maxPDFNonBJetPhi;
float KINFIT::kfit::meanPDFNonBJetPhi;
float KINFIT::kfit::sigmaPDFNonBJetPhi;
double KINFIT::kfit::xminPDFNonBJetPhi;
double KINFIT::kfit::xmaxPDFNonBJetPhi;

float KINFIT::kfit::maxPDFPhotonPx;
float KINFIT::kfit::meanPDFPhotonPx;
float KINFIT::kfit::sigmaPDFPhotonPx;
double KINFIT::kfit::xminPDFPhotonPx;
double KINFIT::kfit::xmaxPDFPhotonPx;

float KINFIT::kfit::maxPDFPhotonPy;
float KINFIT::kfit::meanPDFPhotonPy;
float KINFIT::kfit::sigmaPDFPhotonPy;
double KINFIT::kfit::xminPDFPhotonPy;
double KINFIT::kfit::xmaxPDFPhotonPy;

float KINFIT::kfit::maxPDFPhotonPz;
float KINFIT::kfit::meanPDFPhotonPz;
float KINFIT::kfit::sigmaPDFPhotonPz;
double KINFIT::kfit::xminPDFPhotonPz;
double KINFIT::kfit::xmaxPDFPhotonPz;

float KINFIT::kfit::maxPDFPhotonPt;
float KINFIT::kfit::meanPDFPhotonPt;
float KINFIT::kfit::sigmaPDFPhotonPt;
double KINFIT::kfit::xminPDFPhotonPt;
double KINFIT::kfit::xmaxPDFPhotonPt;

float KINFIT::kfit::maxPDFPhotonEta;
float KINFIT::kfit::meanPDFPhotonEta;
float KINFIT::kfit::sigmaPDFPhotonEta;
double KINFIT::kfit::xminPDFPhotonEta;
double KINFIT::kfit::xmaxPDFPhotonEta;

float KINFIT::kfit::maxPDFPhotonPhi;
float KINFIT::kfit::meanPDFPhotonPhi;
float KINFIT::kfit::sigmaPDFPhotonPhi;
double KINFIT::kfit::xminPDFPhotonPhi;
double KINFIT::kfit::xmaxPDFPhotonPhi;

/// Main constructor
KINFIT::kfit::kfit()
{
   /// Suppress MINUIT output
   gErrorIgnoreLevel = 2000;
   gPrintViaErrorHandler = true;
   
   /// Not initialized yet
   ISINITIALIZED_ = 0;

   /// Create the holder of flags indicating if a specific parameter is fixed
   IsParFixed.reset(); IsParFixed = std::make_shared<std::vector<bool> >(std::vector<bool>(NPARMAX, false));
   
   /// Use default NLL definition
//   AddTopTopMassToNLL_ = false;

   /// Default sorting option for best fits
   SortByTopTopMass_ = false;
   SortByTopTopPz_ = false;
   SortByTopTopPzAbs_ = false;
   SortByNeutrinoPz_ = false;
   SortByNeutrinoPt_ = false;
   SortByNeutrinoP_  = false;
   
   /// Do not measure runtime by default
   DoTiming_ = false;
   
   DoFit_ = false;
   
   /// Default settings
   NToy_ = 10000;
   LHMaxGeneric_ = 20.;
   LHMaxMinuit_ = 0.01;
   NFitMax_ = 50;
   NGrid_ = 50;
   
   /// Default transfer function variation ranges
   NMetRMS_ = 5;

   NBJetPxRMS_ = 3;
   NBJetPyRMS_ = 3;
   NBJetPzRMS_ = 3;
   NBJetPtRMS_ = 3;
   NBJetEtaRMS_ = 3;
   NBJetPhiRMS_ = 3;

   NNonBJetPxRMS_ = 3;
   NNonBJetPyRMS_ = 3;
   NNonBJetPzRMS_ = 3;
   NNonBJetPtRMS_ = 3;
   NNonBJetEtaRMS_ = 3;
   NNonBJetPhiRMS_ = 3;
   
   NElecPxRMS_ = 3;
   NElecPyRMS_ = 3;
   NElecPzRMS_ = 3;
   NElecPtRMS_ = 3;
   NElecEtaRMS_ = 3;
   NElecPhiRMS_ = 3;

   NMuonPxRMS_ = 3;
   NMuonPyRMS_ = 3;
   NMuonPzRMS_ = 3;
   NMuonPtRMS_ = 3;
   NMuonEtaRMS_ = 3;
   NMuonPhiRMS_ = 3;

   NPhotonPxRMS_ = 3;
   NPhotonPyRMS_ = 3;
   NPhotonPzRMS_ = 3;
   NPhotonPtRMS_ = 3;
   NPhotonEtaRMS_ = 3;
   NPhotonPhiRMS_ = 3;
   
   LimNRMS_ = 3;
   
   EtxMin_ = -500.;
   EtxMax_ = 500.;
   
   Etx = 1E+10;
   Ety = 1E+10;

   WMass1 = -1;
   WMass2 = -1;
   
   foundPDF.reset(); foundPDF = std::make_shared<std::map<int, bool> >();
   
   usePDFBJetPxPyPz.reset(); usePDFBJetPxPyPz = std::make_shared<bool>();
   usePDFBJetPtEtaPhi.reset(); usePDFBJetPtEtaPhi = std::make_shared<bool>();
   usePDFNonBJetPxPyPz.reset(); usePDFNonBJetPxPyPz = std::make_shared<bool>();
   usePDFNonBJetPtEtaPhi.reset(); usePDFNonBJetPtEtaPhi = std::make_shared<bool>();
   usePDFLeptonPxPyPz.reset(); usePDFLeptonPxPyPz = std::make_shared<bool>();
   usePDFLeptonPtEtaPhi.reset(); usePDFLeptonPtEtaPhi = std::make_shared<bool>();
   usePDFPhotonPxPyPz.reset(); usePDFPhotonPxPyPz = std::make_shared<bool>();
   usePDFPhotonPtEtaPhi.reset(); usePDFPhotonPtEtaPhi = std::make_shared<bool>();
   
   hPDFTopWMass.reset(); hPDFTopWMass = std::make_shared<TF1>();
   hPDFTopMass.reset(); hPDFTopMass = std::make_shared<TF1>();
   hPDFTopWHadMass.reset(); hPDFTopWHadMass = std::make_shared<TF1>();
   hPDFTopHadMass.reset(); hPDFTopHadMass = std::make_shared<TF1>();
   
   hPDFMetPx.reset(); hPDFMetPx = std::make_shared<TF1>();
   hPDFMetPy.reset(); hPDFMetPy = std::make_shared<TF1>();

//   hPDFTopTopMass.reset(); hPDFTopTopMass = std::make_shared<TF1>();
   
   hPDFBJetPx.reset(); hPDFBJetPx = std::make_shared<TF1>();
   hPDFBJetPy.reset(); hPDFBJetPy = std::make_shared<TF1>();
   hPDFBJetPz.reset(); hPDFBJetPz = std::make_shared<TF1>();
   hPDFBJetPt.reset(); hPDFBJetPt = std::make_shared<TF1>();
   hPDFBJetEta.reset(); hPDFBJetEta = std::make_shared<TF1>();
   hPDFBJetPhi.reset(); hPDFBJetPhi = std::make_shared<TF1>();
   
   hPDFElecPx.reset(); hPDFElecPx = std::make_shared<TF1>();
   hPDFElecPy.reset(); hPDFElecPy = std::make_shared<TF1>();
   hPDFElecPz.reset(); hPDFElecPz = std::make_shared<TF1>();
   hPDFElecPt.reset(); hPDFElecPt = std::make_shared<TF1>();
   hPDFElecEta.reset(); hPDFElecEta = std::make_shared<TF1>();
   hPDFElecPhi.reset(); hPDFElecPhi = std::make_shared<TF1>();
   
   hPDFMuonPx.reset(); hPDFMuonPx = std::make_shared<TF1>();
   hPDFMuonPy.reset(); hPDFMuonPy = std::make_shared<TF1>();
   hPDFMuonPz.reset(); hPDFMuonPz = std::make_shared<TF1>();
   hPDFMuonPt.reset(); hPDFMuonPt = std::make_shared<TF1>();
   hPDFMuonEta.reset(); hPDFMuonEta = std::make_shared<TF1>();
   hPDFMuonPhi.reset(); hPDFMuonPhi = std::make_shared<TF1>();
   
   hPDFNonBJetPx.reset(); hPDFNonBJetPx = std::make_shared<TF1>();
   hPDFNonBJetPy.reset(); hPDFNonBJetPy = std::make_shared<TF1>();
   hPDFNonBJetPz.reset(); hPDFNonBJetPz = std::make_shared<TF1>();
   hPDFNonBJetPt.reset(); hPDFNonBJetPt = std::make_shared<TF1>();
   hPDFNonBJetEta.reset(); hPDFNonBJetEta = std::make_shared<TF1>();
   hPDFNonBJetPhi.reset(); hPDFNonBJetPhi = std::make_shared<TF1>();

   hPDFPhotonPx.reset(); hPDFPhotonPx = std::make_shared<TF1>();
   hPDFPhotonPy.reset(); hPDFPhotonPy = std::make_shared<TF1>();
   hPDFPhotonPz.reset(); hPDFPhotonPz = std::make_shared<TF1>();
   hPDFPhotonPt.reset(); hPDFPhotonPt = std::make_shared<TF1>();
   hPDFPhotonEta.reset(); hPDFPhotonEta = std::make_shared<TF1>();
   hPDFPhotonPhi.reset(); hPDFPhotonPhi = std::make_shared<TF1>();
   
   MetPx_ = NULL;
   MetPy_ = NULL;
   
   PhotonOrigin_ = NULL;
   
   WMass_ = NULL;
   WP_ = NULL;
   WPt_ = NULL;
   WPx_ = NULL;
   WPy_ = NULL;
   WPz_ = NULL;
   WEta_ = NULL;
   WPhi_ = NULL;
   WE_ = NULL;
   
   TopMass_ = NULL;
   TopPt_ = NULL;
   TopPx_ = NULL;
   TopPy_ = NULL;
   TopPz_ = NULL;
   TopP_ = NULL;
   TopEta_ = NULL;
   TopPhi_ = NULL;
   TopE_ = NULL;

   chi_ = NULL;
   chiTerm_ = NULL;
   chiTermName_ = NULL;
   NGeneric_ = NULL;
   par_ = NULL;
   nuPx_ = NULL;
   nuPy_ = NULL;
   nuPz_ = NULL;
   idxMin_ = NULL;
   
   timerWallGeneric_ = NULL;
   timerCPUGeneric_ = NULL;
   timerWallMinuit_ = NULL;
   timerCPUMinuit_ = NULL;
   
   TopTopLepLep_Electron1Idx = NULL;
   TopTopLepLep_Muon1Idx = NULL;
   TopTopLepLep_Electron2Idx = NULL;
   TopTopLepLep_Muon2Idx = NULL;
   TopTopLepLep_BJet1Idx = NULL;
   TopTopLepLep_BJet2Idx = NULL;
   TopTopLepLep_PhotonIdx = NULL;

   TopTopLepHad_ElectronIdx = NULL;
   TopTopLepHad_MuonIdx = NULL;
   TopTopLepHad_BJetLepIdx = NULL;
   TopTopLepHad_BJetHadIdx = NULL;
   TopTopLepHad_NonBJet1Idx = NULL;
   TopTopLepHad_NonBJet2Idx = NULL;
   TopTopLepHad_PhotonIdx = NULL;

   TopLep_ElectronIdx = NULL;
   TopLep_MuonIdx = NULL;
   TopLep_BJetLepIdx = NULL;
   TopLep_PhotonIdx = NULL;
   
   drTopTop_ = NULL;
   detaTopTop_ = NULL;
   dphiTopTop_ = NULL;
   mTopTop_ = NULL;
   ptTopTop_ = NULL;
   pTopTop_ = NULL;
   etaTopTop_ = NULL;
   phiTopTop_ = NULL;
}

/// Main destructor
KINFIT::kfit::~kfit()
{
   if( ISINITIALIZED_ ) delete rnd;
}

/// Initialize process for kinematic reconstruction
void KINFIT::kfit::Init(HYPO hypo)
{
   Reset();
   
   bool hypoFound = 0;   

   for( int i=0;i<HYPO_N;i++ )
     {
	HYPO cur = static_cast<HYPO>(i);
	if( hypo == cur )
	  {
	     hypoFound = 1;
	     break;
	  }
     }   
   if( !hypoFound )
     {
	std::cout << "Process " << hypo << " is not defined" << std::endl;
	exit(1);
     }   
   else
     hypoMode = hypo;
   
   if( hypoMode == TOPTOPLEPLEP )
     {	
	hypTopTopLepLep = (KINFIT::TopTopLepLep*)(this);

	/// Define fit parameters
	
	FPARAM_NAME[FPARAM_Etx_TOPTOPLEPLEP]        = "Etx";
	FPARAM_NAME[FPARAM_Ety_TOPTOPLEPLEP]        = "Ety";
	FPARAM_NAME[FPARAM_Sign1_TOPTOPLEPLEP]      = "Sign1";
	FPARAM_NAME[FPARAM_Sign2_TOPTOPLEPLEP]      = "Sign2";
	FPARAM_NAME[FPARAM_EtRealX_TOPTOPLEPLEP]    = "EtRealX";
	FPARAM_NAME[FPARAM_EtRealY_TOPTOPLEPLEP]    = "EtRealY";
	FPARAM_NAME[FPARAM_mW1_TOPTOPLEPLEP]        = "mW1";
	FPARAM_NAME[FPARAM_mW2_TOPTOPLEPLEP]        = "mW2";
	FPARAM_NAME[FPARAM_BJet1Px_TOPTOPLEPLEP]    = "BJet1Px";
	FPARAM_NAME[FPARAM_BJet1Py_TOPTOPLEPLEP]    = "BJet1Py";
	FPARAM_NAME[FPARAM_BJet1Pz_TOPTOPLEPLEP]    = "BJet1Pz";
	FPARAM_NAME[FPARAM_BJet1Pt_TOPTOPLEPLEP]    = "BJet1Pt";
	FPARAM_NAME[FPARAM_BJet1Eta_TOPTOPLEPLEP]   = "BJet1Eta";
	FPARAM_NAME[FPARAM_BJet1Phi_TOPTOPLEPLEP]   = "BJet1Phi";
	FPARAM_NAME[FPARAM_BJet1E_TOPTOPLEPLEP]     = "BJet1E";
	FPARAM_NAME[FPARAM_BJet2Px_TOPTOPLEPLEP]    = "BJet2Px";
	FPARAM_NAME[FPARAM_BJet2Py_TOPTOPLEPLEP]    = "BJet2Py";
	FPARAM_NAME[FPARAM_BJet2Pz_TOPTOPLEPLEP]    = "BJet2Pz";
	FPARAM_NAME[FPARAM_BJet2Pt_TOPTOPLEPLEP]    = "BJet2Pt";
	FPARAM_NAME[FPARAM_BJet2Eta_TOPTOPLEPLEP]   = "BJet2Eta";
	FPARAM_NAME[FPARAM_BJet2Phi_TOPTOPLEPLEP]   = "BJet2Phi";
	FPARAM_NAME[FPARAM_BJet2E_TOPTOPLEPLEP]     = "BJet2E";
	FPARAM_NAME[FPARAM_Lepton1Px_TOPTOPLEPLEP]  = "Lepton1Px";
	FPARAM_NAME[FPARAM_Lepton1Py_TOPTOPLEPLEP]  = "Lepton1Py";
	FPARAM_NAME[FPARAM_Lepton1Pz_TOPTOPLEPLEP]  = "Lepton1Pz";
	FPARAM_NAME[FPARAM_Lepton1Pt_TOPTOPLEPLEP]  = "Lepton1Pt";
	FPARAM_NAME[FPARAM_Lepton1Eta_TOPTOPLEPLEP] = "Lepton1Eta";
	FPARAM_NAME[FPARAM_Lepton1Phi_TOPTOPLEPLEP] = "Lepton1Phi";
	FPARAM_NAME[FPARAM_Lepton1E_TOPTOPLEPLEP]   = "Lepton1E";
	FPARAM_NAME[FPARAM_Lepton2Px_TOPTOPLEPLEP]  = "Lepton2Px";
	FPARAM_NAME[FPARAM_Lepton2Py_TOPTOPLEPLEP]  = "Lepton2Py";
	FPARAM_NAME[FPARAM_Lepton2Pz_TOPTOPLEPLEP]  = "Lepton2Pz";
	FPARAM_NAME[FPARAM_Lepton2Pt_TOPTOPLEPLEP]  = "Lepton2Pt";
	FPARAM_NAME[FPARAM_Lepton2Eta_TOPTOPLEPLEP] = "Lepton2Eta";
	FPARAM_NAME[FPARAM_Lepton2Phi_TOPTOPLEPLEP] = "Lepton2Phi";
	FPARAM_NAME[FPARAM_Lepton2E_TOPTOPLEPLEP]   = "Lepton2E";
	FPARAM_NAME[FPARAM_PhotonPx_TOPTOPLEPLEP]   = "PhotonPx";
	FPARAM_NAME[FPARAM_PhotonPy_TOPTOPLEPLEP]   = "PhotonPy";
	FPARAM_NAME[FPARAM_PhotonPz_TOPTOPLEPLEP]   = "PhotonPz";
	FPARAM_NAME[FPARAM_PhotonPt_TOPTOPLEPLEP]   = "PhotonPt";
	FPARAM_NAME[FPARAM_PhotonEta_TOPTOPLEPLEP]  = "PhotonEta";
	FPARAM_NAME[FPARAM_PhotonPhi_TOPTOPLEPLEP]  = "PhotonPhi";
	FPARAM_NAME[FPARAM_PhotonE_TOPTOPLEPLEP]    = "PhotonE";
	FPARAM_N = FPARAM_N_TOPTOPLEPLEP;

	/// Transfer functions
	
	PDF_NAME[PDF_TopWMass]     = "PDFTopWMass";
	PDF_NAME[PDF_TopMass]      = "PDFTopMass";
//	PDF_NAME[PDF_TopTopMass]   = "PDFTopTopMass";
	PDF_NAME[PDF_MetPx]        = "PDFMetPx";
	PDF_NAME[PDF_MetPy]        = "PDFMetPy";
	PDF_NAME[PDF_BJetPx]       = "PDFBJetPx";
	PDF_NAME[PDF_BJetPy]       = "PDFBJetPy";
	PDF_NAME[PDF_BJetPz]       = "PDFBJetPz";
	PDF_NAME[PDF_BJetPt]       = "PDFBJetPt";
	PDF_NAME[PDF_BJetEta]      = "PDFBJetEta";
	PDF_NAME[PDF_BJetPhi]      = "PDFBJetPhi";
	PDF_NAME[PDF_ElecPx]       = "PDFElecPx";
	PDF_NAME[PDF_ElecPy]       = "PDFElecPy";
	PDF_NAME[PDF_ElecPz]       = "PDFElecPz";
	PDF_NAME[PDF_ElecPt]       = "PDFElecPt";
	PDF_NAME[PDF_ElecEta]      = "PDFElecEta";
	PDF_NAME[PDF_ElecPhi]      = "PDFElecPhi";
	PDF_NAME[PDF_MuonPx]       = "PDFMuonPx";
	PDF_NAME[PDF_MuonPy]       = "PDFMuonPy";
	PDF_NAME[PDF_MuonPz]       = "PDFMuonPz";
	PDF_NAME[PDF_MuonPt]       = "PDFMuonPt";
	PDF_NAME[PDF_MuonEta]      = "PDFMuonEta";
	PDF_NAME[PDF_MuonPhi]      = "PDFMuonPhi";
	PDF_NAME[PDF_PhotonPx]     = "PDFPhotonPx";
	PDF_NAME[PDF_PhotonPy]     = "PDFPhotonPy";
	PDF_NAME[PDF_PhotonPz]     = "PDFPhotonPz";
	PDF_NAME[PDF_PhotonPt]     = "PDFPhotonPt";
	PDF_NAME[PDF_PhotonEta]    = "PDFPhotonEta";
	PDF_NAME[PDF_PhotonPhi]    = "PDFPhotonPhi";
	
	/// Fix all parameters except one neutrino (px, py) and two W boson masses
	
	(*IsParFixed)[FPARAM_Etx_TOPTOPLEPLEP] = false;
	(*IsParFixed)[FPARAM_Ety_TOPTOPLEPLEP] = false;
	
	(*IsParFixed)[FPARAM_Sign1_TOPTOPLEPLEP] = true;
	(*IsParFixed)[FPARAM_Sign2_TOPTOPLEPLEP] = true;
	
	(*IsParFixed)[FPARAM_EtRealX_TOPTOPLEPLEP] = true;
	(*IsParFixed)[FPARAM_EtRealY_TOPTOPLEPLEP] = true;

	(*IsParFixed)[FPARAM_mW1_TOPTOPLEPLEP] = false;
	(*IsParFixed)[FPARAM_mW2_TOPTOPLEPLEP] = false;
	
	(*IsParFixed)[FPARAM_BJet1Px_TOPTOPLEPLEP]    = true;
	(*IsParFixed)[FPARAM_BJet1Py_TOPTOPLEPLEP]    = true;
	(*IsParFixed)[FPARAM_BJet1Pz_TOPTOPLEPLEP]    = true;
	(*IsParFixed)[FPARAM_BJet1Pt_TOPTOPLEPLEP]    = true;
	(*IsParFixed)[FPARAM_BJet1Eta_TOPTOPLEPLEP]    = true;
	(*IsParFixed)[FPARAM_BJet1Phi_TOPTOPLEPLEP]    = true;
	(*IsParFixed)[FPARAM_BJet1E_TOPTOPLEPLEP]     = true; /// Must always be fixed
	(*IsParFixed)[FPARAM_BJet2Px_TOPTOPLEPLEP]    = true;
	(*IsParFixed)[FPARAM_BJet2Py_TOPTOPLEPLEP]    = true;
	(*IsParFixed)[FPARAM_BJet2Pz_TOPTOPLEPLEP]    = true;
	(*IsParFixed)[FPARAM_BJet2Pt_TOPTOPLEPLEP]    = true;
	(*IsParFixed)[FPARAM_BJet2Eta_TOPTOPLEPLEP]    = true;
	(*IsParFixed)[FPARAM_BJet2Phi_TOPTOPLEPLEP]    = true;
	(*IsParFixed)[FPARAM_BJet2E_TOPTOPLEPLEP]     = true; /// Must always be fixed
	(*IsParFixed)[FPARAM_Lepton1Px_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_Lepton1Py_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_Lepton1Pz_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_Lepton1Pt_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_Lepton1Eta_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_Lepton1Phi_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_Lepton1E_TOPTOPLEPLEP]   = true; /// Must always be fixed
	(*IsParFixed)[FPARAM_Lepton2Px_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_Lepton2Py_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_Lepton2Pz_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_Lepton2Pt_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_Lepton2Eta_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_Lepton2Phi_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_Lepton2E_TOPTOPLEPLEP]   = true; /// Must always be fixed
	(*IsParFixed)[FPARAM_PhotonPx_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_PhotonPy_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_PhotonPz_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_PhotonPt_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_PhotonEta_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_PhotonPhi_TOPTOPLEPLEP]  = true;
	(*IsParFixed)[FPARAM_PhotonE_TOPTOPLEPLEP]   = true; /// Must always be fixed
	
	/// Check that transfer functions are present in the input file
	
	bool foundPDFTopWMass = (! isDeltaFuncPDFTopWMass) ? checkPDF(hPDFTopWMass.get(), PDF_NAME[PDF_TopWMass]) : true; foundPDF.get()->insert({PDF_TopWMass, foundPDFTopWMass});
	bool foundPDFTopMass = (! isDeltaFuncPDFTopMass) ? checkPDF(hPDFTopMass.get(), PDF_NAME[PDF_TopMass]) : true; foundPDF.get()->insert({PDF_TopMass, foundPDFTopMass});
//	checkPDF(hPDFTopTopMass.get(), PDF_NAME[PDF_TopTopMass]);
	
	bool foundPDFMetPx = (! isDeltaFuncPDFMetPx) ? checkPDF(hPDFMetPx.get(), PDF_NAME[PDF_MetPx]) : true; foundPDF.get()->insert({PDF_MetPx, foundPDFMetPx});
	bool foundPDFMetPy = (! isDeltaFuncPDFMetPy) ? checkPDF(hPDFMetPy.get(), PDF_NAME[PDF_MetPy]) : true; foundPDF.get()->insert({PDF_MetPy, foundPDFMetPy});

	bool foundPDFBJetPx = (! isDeltaFuncPDFBJetPx) ? checkPDF(hPDFBJetPx.get(), PDF_NAME[PDF_BJetPx]) : true;
	bool foundPDFBJetPy = (! isDeltaFuncPDFBJetPy) ? checkPDF(hPDFBJetPy.get(), PDF_NAME[PDF_BJetPy]) : true;
	bool foundPDFBJetPz = (! isDeltaFuncPDFBJetPz) ? checkPDF(hPDFBJetPz.get(), PDF_NAME[PDF_BJetPz]) : true;
	bool foundPDFBJetPt = (! isDeltaFuncPDFBJetPt) ? checkPDF(hPDFBJetPt.get(), PDF_NAME[PDF_BJetPt]) : true;
	bool foundPDFBJetEta = (! isDeltaFuncPDFBJetEta) ? checkPDF(hPDFBJetEta.get(), PDF_NAME[PDF_BJetEta]) : true;
	bool foundPDFBJetPhi = (! isDeltaFuncPDFBJetPhi) ? checkPDF(hPDFBJetPhi.get(), PDF_NAME[PDF_BJetPhi]) : true;
	
	bool foundPDFBJetPxPyPz = (foundPDFBJetPx || foundPDFBJetPy || foundPDFBJetPz);
	bool foundPDFBJetPtEtaPhi = (foundPDFBJetPt || foundPDFBJetEta || foundPDFBJetPhi);
	
	if( foundPDFBJetPxPyPz && foundPDFBJetPtEtaPhi ) printPDFWarning();
	else if( foundPDFBJetPxPyPz ) {
	   (*usePDFBJetPxPyPz) = true;
	   (*usePDFBJetPtEtaPhi) = false;
	   foundPDF.get()->insert({PDF_BJetPx, foundPDFBJetPx});
	   foundPDF.get()->insert({PDF_BJetPy, foundPDFBJetPy});
	   foundPDF.get()->insert({PDF_BJetPz, foundPDFBJetPz});
	}
	else {
	   (*usePDFBJetPxPyPz) = false;
	   (*usePDFBJetPtEtaPhi) = true;
	   foundPDF.get()->insert({PDF_BJetPt, foundPDFBJetPt});
	   foundPDF.get()->insert({PDF_BJetEta, foundPDFBJetEta});
	   foundPDF.get()->insert({PDF_BJetPhi, foundPDFBJetPhi});
	}
	
	bool foundPDFElecPx = (! isDeltaFuncPDFElecPx) ? checkPDF(hPDFElecPx.get(), PDF_NAME[PDF_ElecPx]) : true;
	bool foundPDFElecPy = (! isDeltaFuncPDFElecPy) ? checkPDF(hPDFElecPy.get(), PDF_NAME[PDF_ElecPy]) : true;
	bool foundPDFElecPz = (! isDeltaFuncPDFElecPz) ? checkPDF(hPDFElecPz.get(), PDF_NAME[PDF_ElecPz]) : true;
	bool foundPDFElecPt = (! isDeltaFuncPDFElecPt) ? checkPDF(hPDFElecPt.get(), PDF_NAME[PDF_ElecPt]) : true;
	bool foundPDFElecEta = (! isDeltaFuncPDFElecEta) ? checkPDF(hPDFElecEta.get(), PDF_NAME[PDF_ElecEta]) : true;
	bool foundPDFElecPhi = (! isDeltaFuncPDFElecPhi) ? checkPDF(hPDFElecPhi.get(), PDF_NAME[PDF_ElecPhi]) : true;

	bool foundPDFElecPxPyPz = (foundPDFElecPx || foundPDFElecPy || foundPDFElecPz);
	bool foundPDFElecPtEtaPhi = (foundPDFElecPt || foundPDFElecEta || foundPDFElecPhi);
	
	if( foundPDFElecPxPyPz && foundPDFElecPtEtaPhi ) printPDFWarning();
	else if( foundPDFElecPxPyPz ) {
	   foundPDF.get()->insert({PDF_ElecPx, foundPDFElecPx});
	   foundPDF.get()->insert({PDF_ElecPy, foundPDFElecPy});
	   foundPDF.get()->insert({PDF_ElecPz, foundPDFElecPz});
	}
	else {
	   foundPDF.get()->insert({PDF_ElecPt, foundPDFElecPt});
	   foundPDF.get()->insert({PDF_ElecEta, foundPDFElecEta});
	   foundPDF.get()->insert({PDF_ElecPhi, foundPDFElecPhi});
	}
	
	bool foundPDFMuonPx = (! isDeltaFuncPDFMuonPx) ? checkPDF(hPDFMuonPx.get(), PDF_NAME[PDF_MuonPx]) : true;
	bool foundPDFMuonPy = (! isDeltaFuncPDFMuonPy) ? checkPDF(hPDFMuonPy.get(), PDF_NAME[PDF_MuonPy]) : true;
	bool foundPDFMuonPz = (! isDeltaFuncPDFMuonPz) ? checkPDF(hPDFMuonPz.get(), PDF_NAME[PDF_MuonPz]) : true;
	bool foundPDFMuonPt = (! isDeltaFuncPDFMuonPt) ? checkPDF(hPDFMuonPt.get(), PDF_NAME[PDF_MuonPt]) : true;
	bool foundPDFMuonEta = (! isDeltaFuncPDFMuonEta) ? checkPDF(hPDFMuonEta.get(), PDF_NAME[PDF_MuonEta]) : true;
	bool foundPDFMuonPhi = (! isDeltaFuncPDFMuonPhi) ? checkPDF(hPDFMuonPhi.get(), PDF_NAME[PDF_MuonPhi]) : true;

	bool foundPDFMuonPxPyPz = (foundPDFMuonPx || foundPDFMuonPy || foundPDFMuonPz);
	bool foundPDFMuonPtEtaPhi = (foundPDFMuonPt || foundPDFMuonEta || foundPDFMuonPhi);
	
	if( foundPDFMuonPxPyPz && foundPDFMuonPtEtaPhi ) printPDFWarning();
	else if( foundPDFMuonPxPyPz ) {
	   foundPDF.get()->insert({PDF_MuonPx, foundPDFMuonPx});
	   foundPDF.get()->insert({PDF_MuonPy, foundPDFMuonPy});
	   foundPDF.get()->insert({PDF_MuonPz, foundPDFMuonPz});
	}
	else {
	   foundPDF.get()->insert({PDF_MuonPt, foundPDFMuonPt});
	   foundPDF.get()->insert({PDF_MuonEta, foundPDFMuonEta});
	   foundPDF.get()->insert({PDF_MuonPhi, foundPDFMuonPhi});
	}
	
	if( foundPDFElecPxPyPz && foundPDFMuonPxPyPz ) {
	   (*usePDFLeptonPxPyPz) = true;
	   (*usePDFLeptonPtEtaPhi) = false;
	}
	else if( foundPDFElecPtEtaPhi && foundPDFMuonPtEtaPhi ) {
	   (*usePDFLeptonPxPyPz) = false;
	   (*usePDFLeptonPtEtaPhi) = true;
	}
	else {
	   std::cout << "Please define the same type PDFs for leptons" << std::endl;
	   exit(1);
	}
	
	DoFit_ = true;
	
	NNu_ = 2; ///< Number of neutrinos in the final state
     }
   else if( hypoMode == TOPTOPLEPHAD )
     {	
	hypTopTopLepHad = (KINFIT::TopTopLepHad*)(this);

	/// Define fit parameters
	
	FPARAM_NAME[FPARAM_Sign_TOPTOPLEPHAD]         = "Sign";
	FPARAM_NAME[FPARAM_EtRealX_TOPTOPLEPHAD]      = "EtRealX";
	FPARAM_NAME[FPARAM_EtRealY_TOPTOPLEPHAD]      = "EtRealY";
	FPARAM_NAME[FPARAM_mWLep_TOPTOPLEPHAD]        = "mWLep";
	FPARAM_NAME[FPARAM_mWHad_TOPTOPLEPHAD]        = "mWHad";
	FPARAM_NAME[FPARAM_BJetLepPx_TOPTOPLEPHAD]    = "BJetLepPx";
	FPARAM_NAME[FPARAM_BJetLepPy_TOPTOPLEPHAD]    = "BJetLepPy";
	FPARAM_NAME[FPARAM_BJetLepPz_TOPTOPLEPHAD]    = "BJetLepPz";
	FPARAM_NAME[FPARAM_BJetLepPt_TOPTOPLEPHAD]    = "BJetLepPt";
	FPARAM_NAME[FPARAM_BJetLepEta_TOPTOPLEPHAD]   = "BJetLepEta";
	FPARAM_NAME[FPARAM_BJetLepPhi_TOPTOPLEPHAD]   = "BJetLepPhi";
	FPARAM_NAME[FPARAM_BJetLepE_TOPTOPLEPHAD]     = "BJetLepE";
	FPARAM_NAME[FPARAM_BJetHadPx_TOPTOPLEPHAD]    = "BJetHadPx";
	FPARAM_NAME[FPARAM_BJetHadPy_TOPTOPLEPHAD]    = "BJetHadPy";
	FPARAM_NAME[FPARAM_BJetHadPz_TOPTOPLEPHAD]    = "BJetHadPz";
	FPARAM_NAME[FPARAM_BJetHadPt_TOPTOPLEPHAD]    = "BJetHadPt";
	FPARAM_NAME[FPARAM_BJetHadEta_TOPTOPLEPHAD]   = "BJetHadEta";
	FPARAM_NAME[FPARAM_BJetHadPhi_TOPTOPLEPHAD]   = "BJetHadPhi";
	FPARAM_NAME[FPARAM_BJetHadE_TOPTOPLEPHAD]     = "BJetHadE";
	FPARAM_NAME[FPARAM_NonBJet1Px_TOPTOPLEPHAD]   = "NonBJet1Px";
	FPARAM_NAME[FPARAM_NonBJet1Py_TOPTOPLEPHAD]   = "NonBJet1Py";
	FPARAM_NAME[FPARAM_NonBJet1Pz_TOPTOPLEPHAD]   = "NonBJet1Pz";
	FPARAM_NAME[FPARAM_NonBJet1Pt_TOPTOPLEPHAD]   = "NonBJet1Pt";
	FPARAM_NAME[FPARAM_NonBJet1Eta_TOPTOPLEPHAD]  = "NonBJet1Eta";
	FPARAM_NAME[FPARAM_NonBJet1Phi_TOPTOPLEPHAD]  = "NonBJet1Phi";
	FPARAM_NAME[FPARAM_NonBJet1E_TOPTOPLEPHAD]    = "NonBJet1E";
	FPARAM_NAME[FPARAM_NonBJet2Px_TOPTOPLEPHAD]   = "NonBJet2Px";
	FPARAM_NAME[FPARAM_NonBJet2Py_TOPTOPLEPHAD]   = "NonBJet2Py";
	FPARAM_NAME[FPARAM_NonBJet2Pz_TOPTOPLEPHAD]   = "NonBJet2Pz";
	FPARAM_NAME[FPARAM_NonBJet2Pt_TOPTOPLEPHAD]   = "NonBJet2Pt";
	FPARAM_NAME[FPARAM_NonBJet2Eta_TOPTOPLEPHAD]  = "NonBJet2Eta";
	FPARAM_NAME[FPARAM_NonBJet2Phi_TOPTOPLEPHAD]  = "NonBJet2Phi";
	FPARAM_NAME[FPARAM_NonBJet2E_TOPTOPLEPHAD]    = "NonBJet2E";
	FPARAM_NAME[FPARAM_LeptonPx_TOPTOPLEPHAD]     = "LeptonPx";
	FPARAM_NAME[FPARAM_LeptonPy_TOPTOPLEPHAD]     = "LeptonPy";
	FPARAM_NAME[FPARAM_LeptonPz_TOPTOPLEPHAD]     = "LeptonPz";
	FPARAM_NAME[FPARAM_LeptonPt_TOPTOPLEPHAD]     = "LeptonPt";
	FPARAM_NAME[FPARAM_LeptonEta_TOPTOPLEPHAD]    = "LeptonEta";
	FPARAM_NAME[FPARAM_LeptonPhi_TOPTOPLEPHAD]    = "LeptonPhi";
	FPARAM_NAME[FPARAM_LeptonE_TOPTOPLEPHAD]      = "LeptonE";
	FPARAM_NAME[FPARAM_PhotonPx_TOPTOPLEPHAD]     = "PhotonPx";
	FPARAM_NAME[FPARAM_PhotonPy_TOPTOPLEPHAD]     = "PhotonPy";
	FPARAM_NAME[FPARAM_PhotonPz_TOPTOPLEPHAD]     = "PhotonPz";
	FPARAM_NAME[FPARAM_PhotonPt_TOPTOPLEPHAD]     = "PhotonPt";
	FPARAM_NAME[FPARAM_PhotonEta_TOPTOPLEPHAD]    = "PhotonEta";
	FPARAM_NAME[FPARAM_PhotonPhi_TOPTOPLEPHAD]    = "PhotonPhi";
	FPARAM_NAME[FPARAM_PhotonE_TOPTOPLEPHAD]      = "PhotonE";
	FPARAM_N = FPARAM_N_TOPTOPLEPHAD;

	/// Transfer functions
	
	PDF_NAME[PDF_TopWMass]        = "PDFTopWMass";
	PDF_NAME[PDF_TopMass]         = "PDFTopMass";
	PDF_NAME[PDF_TopWHadMass]     = "PDFTopWHadMass";
	PDF_NAME[PDF_TopHadMass]      = "PDFTopHadMass";
//	PDF_NAME[PDF_TopTopMass]      = "PDFTopTopMass";
	PDF_NAME[PDF_MetPx]           = "PDFMetPx";
	PDF_NAME[PDF_MetPy]           = "PDFMetPy";
	PDF_NAME[PDF_BJetPx]          = "PDFBJetPx";
	PDF_NAME[PDF_BJetPy]          = "PDFBJetPy";
	PDF_NAME[PDF_BJetPz]          = "PDFBJetPz";
	PDF_NAME[PDF_BJetPt]          = "PDFBJetPt";
	PDF_NAME[PDF_BJetEta]         = "PDFBJetEta";
	PDF_NAME[PDF_BJetPhi]         = "PDFBJetPhi";
	PDF_NAME[PDF_NonBJetPx]       = "PDFNonBJetPx";
	PDF_NAME[PDF_NonBJetPy]       = "PDFNonBJetPy";
	PDF_NAME[PDF_NonBJetPz]       = "PDFNonBJetPz";
	PDF_NAME[PDF_NonBJetPt]       = "PDFNonBJetPt";
	PDF_NAME[PDF_NonBJetEta]      = "PDFNonBJetEta";
	PDF_NAME[PDF_NonBJetPhi]      = "PDFNonBJetPhi";
	PDF_NAME[PDF_ElecPx]          = "PDFElecPx";
	PDF_NAME[PDF_ElecPy]          = "PDFElecPy";
	PDF_NAME[PDF_ElecPz]          = "PDFElecPz";
	PDF_NAME[PDF_ElecPt]          = "PDFElecPt";
	PDF_NAME[PDF_ElecEta]         = "PDFElecEta";
	PDF_NAME[PDF_ElecPhi]         = "PDFElecPhi";
	PDF_NAME[PDF_MuonPx]          = "PDFMuonPx";
	PDF_NAME[PDF_MuonPy]          = "PDFMuonPy";
	PDF_NAME[PDF_MuonPz]          = "PDFMuonPz";
	PDF_NAME[PDF_MuonPt]          = "PDFMuonPt";
	PDF_NAME[PDF_MuonEta]         = "PDFMuonEta";
	PDF_NAME[PDF_MuonPhi]         = "PDFMuonPhi";
	PDF_NAME[PDF_PhotonPx]        = "PDFPhotonPx";
	PDF_NAME[PDF_PhotonPy]        = "PDFPhotonPy";
	PDF_NAME[PDF_PhotonPz]        = "PDFPhotonPz";
	PDF_NAME[PDF_PhotonPt]        = "PDFPhotonPt";
	PDF_NAME[PDF_PhotonEta]       = "PDFPhotonEta";
	PDF_NAME[PDF_PhotonPhi]       = "PDFPhotonPhi";
	
	/// Fix all parameters except leptonic W boson mass
	
	(*IsParFixed)[FPARAM_Sign_TOPTOPLEPHAD] = true;
	
	(*IsParFixed)[FPARAM_EtRealX_TOPTOPLEPHAD] = true;
	(*IsParFixed)[FPARAM_EtRealY_TOPTOPLEPHAD] = true;

	(*IsParFixed)[FPARAM_mWLep_TOPTOPLEPHAD] = false;
	(*IsParFixed)[FPARAM_mWHad_TOPTOPLEPHAD] = true;
	
	(*IsParFixed)[FPARAM_BJetLepPx_TOPTOPLEPHAD]    = true;
	(*IsParFixed)[FPARAM_BJetLepPy_TOPTOPLEPHAD]    = true;
	(*IsParFixed)[FPARAM_BJetLepPz_TOPTOPLEPHAD]    = true;
	(*IsParFixed)[FPARAM_BJetLepPt_TOPTOPLEPHAD]    = true;
	(*IsParFixed)[FPARAM_BJetLepEta_TOPTOPLEPHAD]   = true;
	(*IsParFixed)[FPARAM_BJetLepPhi_TOPTOPLEPHAD]   = true;
	(*IsParFixed)[FPARAM_BJetLepE_TOPTOPLEPHAD]     = true; /// Must always be fixed
	(*IsParFixed)[FPARAM_BJetHadPx_TOPTOPLEPHAD]    = true;
	(*IsParFixed)[FPARAM_BJetHadPy_TOPTOPLEPHAD]    = true;
	(*IsParFixed)[FPARAM_BJetHadPz_TOPTOPLEPHAD]    = true;
	(*IsParFixed)[FPARAM_BJetHadPt_TOPTOPLEPHAD]    = true;
	(*IsParFixed)[FPARAM_BJetHadEta_TOPTOPLEPHAD]   = true;
	(*IsParFixed)[FPARAM_BJetHadPhi_TOPTOPLEPHAD]   = true;
	(*IsParFixed)[FPARAM_BJetHadE_TOPTOPLEPHAD]     = true; /// Must always be fixed
	(*IsParFixed)[FPARAM_NonBJet1Px_TOPTOPLEPHAD]   = true;
	(*IsParFixed)[FPARAM_NonBJet1Py_TOPTOPLEPHAD]   = true;
	(*IsParFixed)[FPARAM_NonBJet1Pz_TOPTOPLEPHAD]   = true;
	(*IsParFixed)[FPARAM_NonBJet1Pt_TOPTOPLEPHAD]   = true;
	(*IsParFixed)[FPARAM_NonBJet1Eta_TOPTOPLEPHAD]  = true;
	(*IsParFixed)[FPARAM_NonBJet1Phi_TOPTOPLEPHAD]  = true;
	(*IsParFixed)[FPARAM_NonBJet1E_TOPTOPLEPHAD]    = true; /// Must always be fixed
	(*IsParFixed)[FPARAM_NonBJet2Px_TOPTOPLEPHAD]   = true;
	(*IsParFixed)[FPARAM_NonBJet2Py_TOPTOPLEPHAD]   = true;
	(*IsParFixed)[FPARAM_NonBJet2Pz_TOPTOPLEPHAD]   = true;
	(*IsParFixed)[FPARAM_NonBJet2Pt_TOPTOPLEPHAD]   = true;
	(*IsParFixed)[FPARAM_NonBJet2Eta_TOPTOPLEPHAD]  = true;
	(*IsParFixed)[FPARAM_NonBJet2Phi_TOPTOPLEPHAD]  = true;
	(*IsParFixed)[FPARAM_NonBJet2E_TOPTOPLEPHAD]    = true; /// Must always be fixed
	(*IsParFixed)[FPARAM_LeptonPx_TOPTOPLEPHAD]     = true;
	(*IsParFixed)[FPARAM_LeptonPy_TOPTOPLEPHAD]     = true;
	(*IsParFixed)[FPARAM_LeptonPz_TOPTOPLEPHAD]     = true;
	(*IsParFixed)[FPARAM_LeptonPt_TOPTOPLEPHAD]     = true;
	(*IsParFixed)[FPARAM_LeptonEta_TOPTOPLEPHAD]    = true;
	(*IsParFixed)[FPARAM_LeptonPhi_TOPTOPLEPHAD]    = true;
	(*IsParFixed)[FPARAM_LeptonE_TOPTOPLEPHAD]      = true; /// Must always be fixed
	(*IsParFixed)[FPARAM_PhotonPx_TOPTOPLEPHAD]     = true;
	(*IsParFixed)[FPARAM_PhotonPy_TOPTOPLEPHAD]     = true;
	(*IsParFixed)[FPARAM_PhotonPz_TOPTOPLEPHAD]     = true;
	(*IsParFixed)[FPARAM_PhotonPt_TOPTOPLEPHAD]     = true;
	(*IsParFixed)[FPARAM_PhotonEta_TOPTOPLEPHAD]    = true;
	(*IsParFixed)[FPARAM_PhotonPhi_TOPTOPLEPHAD]    = true;
	(*IsParFixed)[FPARAM_PhotonE_TOPTOPLEPHAD]      = true; /// Must always be fixed

	/// Check that transfer functions are present in the input file
	
	bool foundPDFTopWMass = (! isDeltaFuncPDFTopWMass) ? checkPDF(hPDFTopWMass.get(), PDF_NAME[PDF_TopWMass]) : true; foundPDF.get()->insert({PDF_TopWMass, foundPDFTopWMass});
	bool foundPDFTopMass = (! isDeltaFuncPDFTopMass) ? checkPDF(hPDFTopMass.get(), PDF_NAME[PDF_TopMass]) : true; foundPDF.get()->insert({PDF_TopMass, foundPDFTopMass});
	bool foundPDFTopWHadMass = (! isDeltaFuncPDFTopWHadMass) ? checkPDF(hPDFTopWHadMass.get(), PDF_NAME[PDF_TopWHadMass]) : true; foundPDF.get()->insert({PDF_TopWHadMass, foundPDFTopWHadMass});
	bool foundPDFTopHadMass = (! isDeltaFuncPDFTopHadMass) ? checkPDF(hPDFTopHadMass.get(), PDF_NAME[PDF_TopHadMass]) : true; foundPDF.get()->insert({PDF_TopHadMass, foundPDFTopHadMass});
//	checkPDF(hPDFTopTopMass.get(), PDF_NAME[PDF_TopTopMass]);
	
	bool foundPDFMetPx = (! isDeltaFuncPDFMetPx) ? checkPDF(hPDFMetPx.get(), PDF_NAME[PDF_MetPx]) : true; foundPDF.get()->insert({PDF_MetPx, foundPDFMetPx});
	bool foundPDFMetPy = (! isDeltaFuncPDFMetPy) ? checkPDF(hPDFMetPy.get(), PDF_NAME[PDF_MetPy]) : true; foundPDF.get()->insert({PDF_MetPy, foundPDFMetPy});
	
	bool foundPDFBJetPx = (! isDeltaFuncPDFBJetPx) ? checkPDF(hPDFBJetPx.get(), PDF_NAME[PDF_BJetPx]) : true;
	bool foundPDFBJetPy = (! isDeltaFuncPDFBJetPy) ? checkPDF(hPDFBJetPy.get(), PDF_NAME[PDF_BJetPy]) : true;
	bool foundPDFBJetPz = (! isDeltaFuncPDFBJetPz) ? checkPDF(hPDFBJetPz.get(), PDF_NAME[PDF_BJetPz]) : true;
	bool foundPDFBJetPt = (! isDeltaFuncPDFBJetPt) ? checkPDF(hPDFBJetPt.get(), PDF_NAME[PDF_BJetPt]) : true;
	bool foundPDFBJetEta = (! isDeltaFuncPDFBJetEta) ? checkPDF(hPDFBJetEta.get(), PDF_NAME[PDF_BJetEta]) : true;
	bool foundPDFBJetPhi = (! isDeltaFuncPDFBJetPhi) ? checkPDF(hPDFBJetPhi.get(), PDF_NAME[PDF_BJetPhi]) : true;

	bool foundPDFBJetPxPyPz = (foundPDFBJetPx || foundPDFBJetPy || foundPDFBJetPz);
	bool foundPDFBJetPtEtaPhi = (foundPDFBJetPt || foundPDFBJetEta || foundPDFBJetPhi);
	
	if( foundPDFBJetPxPyPz && foundPDFBJetPtEtaPhi ) printPDFWarning();
	else if( foundPDFBJetPxPyPz ) {
	   (*usePDFBJetPxPyPz) = true;
	   (*usePDFBJetPtEtaPhi) = false;
	   foundPDF.get()->insert({PDF_BJetPx, foundPDFBJetPx});
	   foundPDF.get()->insert({PDF_BJetPy, foundPDFBJetPy});
	   foundPDF.get()->insert({PDF_BJetPz, foundPDFBJetPz});
	}
	else {
	   (*usePDFBJetPxPyPz) = false;
	   (*usePDFBJetPtEtaPhi) = true;
	   foundPDF.get()->insert({PDF_BJetPt, foundPDFBJetPt});
	   foundPDF.get()->insert({PDF_BJetEta, foundPDFBJetEta});
	   foundPDF.get()->insert({PDF_BJetPhi, foundPDFBJetPhi});
	}
	
	bool foundPDFNonBJetPx = (! isDeltaFuncPDFNonBJetPx) ? checkPDF(hPDFNonBJetPx.get(), PDF_NAME[PDF_NonBJetPx]) : true;
	bool foundPDFNonBJetPy = (! isDeltaFuncPDFNonBJetPy) ? checkPDF(hPDFNonBJetPy.get(), PDF_NAME[PDF_NonBJetPy]) : true;
	bool foundPDFNonBJetPz = (! isDeltaFuncPDFNonBJetPz) ? checkPDF(hPDFNonBJetPz.get(), PDF_NAME[PDF_NonBJetPz]) : true;
	bool foundPDFNonBJetPt = (! isDeltaFuncPDFNonBJetPt) ? checkPDF(hPDFNonBJetPt.get(), PDF_NAME[PDF_NonBJetPt]) : true;
	bool foundPDFNonBJetEta = (! isDeltaFuncPDFNonBJetEta) ? checkPDF(hPDFNonBJetEta.get(), PDF_NAME[PDF_NonBJetEta]) : true;
	bool foundPDFNonBJetPhi = (! isDeltaFuncPDFNonBJetPhi) ? checkPDF(hPDFNonBJetPhi.get(), PDF_NAME[PDF_NonBJetPhi]) : true;

	bool foundPDFNonBJetPxPyPz = (foundPDFNonBJetPx || foundPDFNonBJetPy || foundPDFNonBJetPz);
	bool foundPDFNonBJetPtEtaPhi = (foundPDFNonBJetPt || foundPDFNonBJetEta || foundPDFNonBJetPhi);
	
	if( foundPDFNonBJetPxPyPz && foundPDFNonBJetPtEtaPhi ) printPDFWarning();
	else if( foundPDFNonBJetPxPyPz ) {
	   (*usePDFNonBJetPxPyPz) = true;
	   (*usePDFNonBJetPtEtaPhi) = false;
	   foundPDF.get()->insert({PDF_NonBJetPx, foundPDFNonBJetPx});
	   foundPDF.get()->insert({PDF_NonBJetPy, foundPDFNonBJetPy});
	   foundPDF.get()->insert({PDF_NonBJetPz, foundPDFNonBJetPz});
	}
	else {
	   (*usePDFNonBJetPxPyPz) = false;
	   (*usePDFNonBJetPtEtaPhi) = true;
	   foundPDF.get()->insert({PDF_NonBJetPt, foundPDFNonBJetPt});
	   foundPDF.get()->insert({PDF_NonBJetEta, foundPDFNonBJetEta});
	   foundPDF.get()->insert({PDF_NonBJetPhi, foundPDFNonBJetPhi});
	}
	
	bool foundPDFElecPx = (! isDeltaFuncPDFElecPx) ? checkPDF(hPDFElecPx.get(), PDF_NAME[PDF_ElecPx]) : true;
	bool foundPDFElecPy = (! isDeltaFuncPDFElecPy) ? checkPDF(hPDFElecPy.get(), PDF_NAME[PDF_ElecPy]) : true;
	bool foundPDFElecPz = (! isDeltaFuncPDFElecPz) ? checkPDF(hPDFElecPz.get(), PDF_NAME[PDF_ElecPz]) : true;
	bool foundPDFElecPt = (! isDeltaFuncPDFElecPt) ? checkPDF(hPDFElecPt.get(), PDF_NAME[PDF_ElecPt]) : true;
	bool foundPDFElecEta = (! isDeltaFuncPDFElecEta) ? checkPDF(hPDFElecEta.get(), PDF_NAME[PDF_ElecEta]) : true;
	bool foundPDFElecPhi = (! isDeltaFuncPDFElecPhi) ? checkPDF(hPDFElecPhi.get(), PDF_NAME[PDF_ElecPhi]) : true;

	bool foundPDFElecPxPyPz = (foundPDFElecPx || foundPDFElecPy || foundPDFElecPz);
	bool foundPDFElecPtEtaPhi = (foundPDFElecPt || foundPDFElecEta || foundPDFElecPhi);
	
	if( foundPDFElecPxPyPz && foundPDFElecPtEtaPhi ) printPDFWarning();
	else if( foundPDFElecPxPyPz ) {
	   foundPDF.get()->insert({PDF_ElecPx, foundPDFElecPx});
	   foundPDF.get()->insert({PDF_ElecPy, foundPDFElecPy});
	   foundPDF.get()->insert({PDF_ElecPz, foundPDFElecPz});
	}
	else {
	   foundPDF.get()->insert({PDF_ElecPt, foundPDFElecPt});
	   foundPDF.get()->insert({PDF_ElecEta, foundPDFElecEta});
	   foundPDF.get()->insert({PDF_ElecPhi, foundPDFElecPhi});
	}
	
	bool foundPDFMuonPx = (! isDeltaFuncPDFMuonPx) ? checkPDF(hPDFMuonPx.get(), PDF_NAME[PDF_MuonPx]) : true;
	bool foundPDFMuonPy = (! isDeltaFuncPDFMuonPy) ? checkPDF(hPDFMuonPy.get(), PDF_NAME[PDF_MuonPy]) : true;
	bool foundPDFMuonPz = (! isDeltaFuncPDFMuonPz) ? checkPDF(hPDFMuonPz.get(), PDF_NAME[PDF_MuonPz]) : true;
	bool foundPDFMuonPt = (! isDeltaFuncPDFMuonPt) ? checkPDF(hPDFMuonPt.get(), PDF_NAME[PDF_MuonPt]) : true;
	bool foundPDFMuonEta = (! isDeltaFuncPDFMuonEta) ? checkPDF(hPDFMuonEta.get(), PDF_NAME[PDF_MuonEta]) : true;
	bool foundPDFMuonPhi = (! isDeltaFuncPDFMuonPhi) ? checkPDF(hPDFMuonPhi.get(), PDF_NAME[PDF_MuonPhi]) : true;

	bool foundPDFMuonPxPyPz = (foundPDFMuonPx || foundPDFMuonPy || foundPDFMuonPz);
	bool foundPDFMuonPtEtaPhi = (foundPDFMuonPt || foundPDFMuonEta || foundPDFMuonPhi);
	
	if( foundPDFMuonPxPyPz && foundPDFMuonPtEtaPhi ) printPDFWarning();
	else if( foundPDFMuonPxPyPz ) {
	   foundPDF.get()->insert({PDF_MuonPx, foundPDFMuonPx});
	   foundPDF.get()->insert({PDF_MuonPy, foundPDFMuonPy});
	   foundPDF.get()->insert({PDF_MuonPz, foundPDFMuonPz});
	}
	else {
	   foundPDF.get()->insert({PDF_MuonPt, foundPDFMuonPt});
	   foundPDF.get()->insert({PDF_MuonEta, foundPDFMuonEta});
	   foundPDF.get()->insert({PDF_MuonPhi, foundPDFMuonPhi});
	}

	if( foundPDFElecPxPyPz && foundPDFMuonPxPyPz ) {
	   (*usePDFLeptonPxPyPz) = true;
	   (*usePDFLeptonPtEtaPhi) = false;
	}
	else if( foundPDFElecPtEtaPhi && foundPDFMuonPtEtaPhi ) {
	   (*usePDFLeptonPxPyPz) = false;
	   (*usePDFLeptonPtEtaPhi) = true;
	}
	else {
	   std::cout << "Please define the same type PDFs for leptons" << std::endl;
	   exit(1);
	}
	
	DoFit_ = false;
	
	NNu_ = 1; ///< Number of neutrinos in the final state
     }
   else if( hypoMode == TOPLEP )
     {
	hypTopLep = (KINFIT::TopLep*)(this);

	/// Define fit parameters
	
	FPARAM_NAME[FPARAM_Sign_TOPLEP]         = "Sign";
	FPARAM_NAME[FPARAM_EtRealX_TOPLEP]      = "EtRealX";
	FPARAM_NAME[FPARAM_EtRealY_TOPLEP]      = "EtRealY";
	FPARAM_NAME[FPARAM_mWLep_TOPLEP]        = "mWLep";
	FPARAM_NAME[FPARAM_BJetLepPx_TOPLEP]    = "BJetLepPx";
	FPARAM_NAME[FPARAM_BJetLepPy_TOPLEP]    = "BJetLepPy";
	FPARAM_NAME[FPARAM_BJetLepPz_TOPLEP]    = "BJetLepPz";
	FPARAM_NAME[FPARAM_BJetLepE_TOPLEP]     = "BJetLepE";
	FPARAM_NAME[FPARAM_BJetLepPt_TOPLEP]    = "BJetLepPt";
	FPARAM_NAME[FPARAM_BJetLepEta_TOPLEP]   = "BJetLepEta";
	FPARAM_NAME[FPARAM_BJetLepPhi_TOPLEP]   = "BJetLepPhi";
	FPARAM_NAME[FPARAM_LeptonPx_TOPLEP]     = "LeptonPx";
	FPARAM_NAME[FPARAM_LeptonPy_TOPLEP]     = "LeptonPy";
	FPARAM_NAME[FPARAM_LeptonPz_TOPLEP]     = "LeptonPz";
	FPARAM_NAME[FPARAM_LeptonE_TOPLEP]      = "LeptonE";
	FPARAM_NAME[FPARAM_LeptonPt_TOPLEP]     = "LeptonPt";
	FPARAM_NAME[FPARAM_LeptonEta_TOPLEP]    = "LeptonEta";
	FPARAM_NAME[FPARAM_LeptonPhi_TOPLEP]    = "LeptonPhi";
	FPARAM_NAME[FPARAM_PhotonPx_TOPLEP]     = "PhotonPx";
	FPARAM_NAME[FPARAM_PhotonPy_TOPLEP]     = "PhotonPy";
	FPARAM_NAME[FPARAM_PhotonPz_TOPLEP]     = "PhotonPz";
	FPARAM_NAME[FPARAM_PhotonE_TOPLEP]      = "PhotonE";
	FPARAM_NAME[FPARAM_PhotonPt_TOPLEP]     = "PhotonPt";
	FPARAM_NAME[FPARAM_PhotonEta_TOPLEP]    = "PhotonEta";
	FPARAM_NAME[FPARAM_PhotonPhi_TOPLEP]    = "PhotonPhi";
	FPARAM_N = FPARAM_N_TOPLEP;

	/// Transfer functions
	
	PDF_NAME[PDF_TopWMass]       = "PDFTopWMass";
	PDF_NAME[PDF_TopMass]        = "PDFTopMass";
	PDF_NAME[PDF_MetPx]          = "PDFMetPx";
	PDF_NAME[PDF_MetPy]          = "PDFMetPy";
	PDF_NAME[PDF_BJetPx]         = "PDFBJetPx";
	PDF_NAME[PDF_BJetPy]         = "PDFBJetPy";
	PDF_NAME[PDF_BJetPz]         = "PDFBJetPz";
	PDF_NAME[PDF_BJetPt]         = "PDFBJetPt";
	PDF_NAME[PDF_BJetEta]        = "PDFBJetEta";
	PDF_NAME[PDF_BJetPhi]        = "PDFBJetPhi";
	PDF_NAME[PDF_ElecPx]         = "PDFElecPx";
	PDF_NAME[PDF_ElecPy]         = "PDFElecPy";
	PDF_NAME[PDF_ElecPz]         = "PDFElecPz";
	PDF_NAME[PDF_ElecPt]         = "PDFElecPt";
	PDF_NAME[PDF_ElecEta]        = "PDFElecEta";
	PDF_NAME[PDF_ElecPhi]        = "PDFElecPhi";
	PDF_NAME[PDF_MuonPx]         = "PDFMuonPx";
	PDF_NAME[PDF_MuonPy]         = "PDFMuonPy";
	PDF_NAME[PDF_MuonPz]         = "PDFMuonPz";
	PDF_NAME[PDF_MuonPt]         = "PDFMuonPt";
	PDF_NAME[PDF_MuonEta]        = "PDFMuonEta";
	PDF_NAME[PDF_MuonPhi]        = "PDFMuonPhi";
	PDF_NAME[PDF_PhotonPx]       = "PDFPhotonPx";
	PDF_NAME[PDF_PhotonPy]       = "PDFPhotonPy";
	PDF_NAME[PDF_PhotonPz]       = "PDFPhotonPz";
	PDF_NAME[PDF_PhotonPt]       = "PDFPhotonPt";
	PDF_NAME[PDF_PhotonEta]      = "PDFPhotonEta";
	PDF_NAME[PDF_PhotonPhi]      = "PDFPhotonPhi";
	
	/// Fix all parameters except leptonic W boson mass
	
	(*IsParFixed)[FPARAM_Sign_TOPLEP] = true;
	
	(*IsParFixed)[FPARAM_EtRealX_TOPLEP] = true;
	(*IsParFixed)[FPARAM_EtRealY_TOPLEP] = true;

	(*IsParFixed)[FPARAM_mWLep_TOPLEP] = false;
	
	(*IsParFixed)[FPARAM_BJetLepPx_TOPLEP]    = true;
	(*IsParFixed)[FPARAM_BJetLepPy_TOPLEP]    = true;
	(*IsParFixed)[FPARAM_BJetLepPz_TOPLEP]    = true;
	(*IsParFixed)[FPARAM_BJetLepE_TOPLEP]     = true; /// Must always be fixed
	(*IsParFixed)[FPARAM_BJetLepPt_TOPLEP]    = true;
	(*IsParFixed)[FPARAM_BJetLepEta_TOPLEP]   = true;
	(*IsParFixed)[FPARAM_BJetLepPhi_TOPLEP]   = true;
	(*IsParFixed)[FPARAM_LeptonPx_TOPLEP]     = true;
	(*IsParFixed)[FPARAM_LeptonPy_TOPLEP]     = true;
	(*IsParFixed)[FPARAM_LeptonPz_TOPLEP]     = true;
	(*IsParFixed)[FPARAM_LeptonE_TOPLEP]      = true; /// Must always be fixed
	(*IsParFixed)[FPARAM_LeptonPt_TOPLEP]     = true;
	(*IsParFixed)[FPARAM_LeptonEta_TOPLEP]    = true;
	(*IsParFixed)[FPARAM_LeptonPhi_TOPLEP]    = true;
	(*IsParFixed)[FPARAM_PhotonPx_TOPLEP]     = true;
	(*IsParFixed)[FPARAM_PhotonPy_TOPLEP]     = true;
	(*IsParFixed)[FPARAM_PhotonPz_TOPLEP]     = true;
	(*IsParFixed)[FPARAM_PhotonE_TOPLEP]      = true; /// Must always be fixed
	(*IsParFixed)[FPARAM_PhotonPt_TOPLEP]     = true;
	(*IsParFixed)[FPARAM_PhotonEta_TOPLEP]    = true;
	(*IsParFixed)[FPARAM_PhotonPhi_TOPLEP]    = true;

	/// Check that transfer functions are present in the input file
	
	bool foundPDFTopWMass = (! isDeltaFuncPDFTopWMass) ? checkPDF(hPDFTopWMass.get(), PDF_NAME[PDF_TopWMass]) : true; foundPDF.get()->insert({PDF_TopWMass, foundPDFTopWMass});
	bool foundPDFTopMass = (! isDeltaFuncPDFTopMass) ? checkPDF(hPDFTopMass.get(), PDF_NAME[PDF_TopMass]) : true; foundPDF.get()->insert({PDF_TopMass, foundPDFTopMass});
	
	bool foundPDFMetPx = (! isDeltaFuncPDFMetPx) ? checkPDF(hPDFMetPx.get(), PDF_NAME[PDF_MetPx]) : true; foundPDF.get()->insert({PDF_MetPx, foundPDFMetPx});
	bool foundPDFMetPy = (! isDeltaFuncPDFMetPy) ? checkPDF(hPDFMetPy.get(), PDF_NAME[PDF_MetPy]) : true; foundPDF.get()->insert({PDF_MetPy, foundPDFMetPy});
	
	bool foundPDFBJetPx = (! isDeltaFuncPDFBJetPx) ? checkPDF(hPDFBJetPx.get(), PDF_NAME[PDF_BJetPx]) : true;
	bool foundPDFBJetPy = (! isDeltaFuncPDFBJetPy) ? checkPDF(hPDFBJetPy.get(), PDF_NAME[PDF_BJetPy]) : true;
	bool foundPDFBJetPz = (! isDeltaFuncPDFBJetPz) ? checkPDF(hPDFBJetPz.get(), PDF_NAME[PDF_BJetPz]) : true;
	bool foundPDFBJetPt = (! isDeltaFuncPDFBJetPt) ? checkPDF(hPDFBJetPt.get(), PDF_NAME[PDF_BJetPt]) : true;
	bool foundPDFBJetEta = (! isDeltaFuncPDFBJetEta) ? checkPDF(hPDFBJetEta.get(), PDF_NAME[PDF_BJetEta]) : true;
	bool foundPDFBJetPhi = (! isDeltaFuncPDFBJetPhi) ? checkPDF(hPDFBJetPhi.get(), PDF_NAME[PDF_BJetPhi]) : true;

	bool foundPDFBJetPxPyPz = (foundPDFBJetPx || foundPDFBJetPy || foundPDFBJetPz);
	bool foundPDFBJetPtEtaPhi = (foundPDFBJetPt || foundPDFBJetEta || foundPDFBJetPhi);
	
	if( foundPDFBJetPxPyPz && foundPDFBJetPtEtaPhi ) printPDFWarning();
	else if( foundPDFBJetPxPyPz ) {
	   (*usePDFBJetPxPyPz) = true;
	   (*usePDFBJetPtEtaPhi) = false;
	   foundPDF.get()->insert({PDF_BJetPx, foundPDFBJetPx});
	   foundPDF.get()->insert({PDF_BJetPy, foundPDFBJetPy});
	   foundPDF.get()->insert({PDF_BJetPz, foundPDFBJetPz});
	}
	else {
	   (*usePDFBJetPxPyPz) = false;
	   (*usePDFBJetPtEtaPhi) = true;
	   foundPDF.get()->insert({PDF_BJetPt, foundPDFBJetPt});
	   foundPDF.get()->insert({PDF_BJetEta, foundPDFBJetEta});
	   foundPDF.get()->insert({PDF_BJetPhi, foundPDFBJetPhi});
	}
	
	bool foundPDFElecPx = (! isDeltaFuncPDFElecPx) ? checkPDF(hPDFElecPx.get(), PDF_NAME[PDF_ElecPx]) : true;
	bool foundPDFElecPy = (! isDeltaFuncPDFElecPy) ? checkPDF(hPDFElecPy.get(), PDF_NAME[PDF_ElecPy]) : true;
	bool foundPDFElecPz = (! isDeltaFuncPDFElecPz) ? checkPDF(hPDFElecPz.get(), PDF_NAME[PDF_ElecPz]) : true;
	bool foundPDFElecPt = (! isDeltaFuncPDFElecPt) ? checkPDF(hPDFElecPt.get(), PDF_NAME[PDF_ElecPt]) : true;
	bool foundPDFElecEta = (! isDeltaFuncPDFElecEta) ? checkPDF(hPDFElecEta.get(), PDF_NAME[PDF_ElecEta]) : true;
	bool foundPDFElecPhi = (! isDeltaFuncPDFElecPhi) ? checkPDF(hPDFElecPhi.get(), PDF_NAME[PDF_ElecPhi]) : true;

	bool foundPDFElecPxPyPz = (foundPDFElecPx || foundPDFElecPy || foundPDFElecPz);
	bool foundPDFElecPtEtaPhi = (foundPDFElecPt || foundPDFElecEta || foundPDFElecPhi);
	
	if( foundPDFElecPxPyPz && foundPDFElecPtEtaPhi ) printPDFWarning();
	else if( foundPDFElecPxPyPz ) {
	   foundPDF.get()->insert({PDF_ElecPx, foundPDFElecPx});
	   foundPDF.get()->insert({PDF_ElecPy, foundPDFElecPy});
	   foundPDF.get()->insert({PDF_ElecPz, foundPDFElecPz});
	}
	else {
	   foundPDF.get()->insert({PDF_ElecPt, foundPDFElecPt});
	   foundPDF.get()->insert({PDF_ElecEta, foundPDFElecEta});
	   foundPDF.get()->insert({PDF_ElecPhi, foundPDFElecPhi});
	}
	
	bool foundPDFMuonPx = (! isDeltaFuncPDFMuonPx) ? checkPDF(hPDFMuonPx.get(), PDF_NAME[PDF_MuonPx]) : true;
	bool foundPDFMuonPy = (! isDeltaFuncPDFMuonPy) ? checkPDF(hPDFMuonPy.get(), PDF_NAME[PDF_MuonPy]) : true;
	bool foundPDFMuonPz = (! isDeltaFuncPDFMuonPz) ? checkPDF(hPDFMuonPz.get(), PDF_NAME[PDF_MuonPz]) : true;
	bool foundPDFMuonPt = (! isDeltaFuncPDFMuonPt) ? checkPDF(hPDFMuonPt.get(), PDF_NAME[PDF_MuonPt]) : true;
	bool foundPDFMuonEta = (! isDeltaFuncPDFMuonEta) ? checkPDF(hPDFMuonEta.get(), PDF_NAME[PDF_MuonEta]) : true;
	bool foundPDFMuonPhi = (! isDeltaFuncPDFMuonPhi) ? checkPDF(hPDFMuonPhi.get(), PDF_NAME[PDF_MuonPhi]) : true;

	bool foundPDFMuonPxPyPz = (foundPDFMuonPx || foundPDFMuonPy || foundPDFMuonPz);
	bool foundPDFMuonPtEtaPhi = (foundPDFMuonPt || foundPDFMuonEta || foundPDFMuonPhi);
	
	if( foundPDFMuonPxPyPz && foundPDFMuonPtEtaPhi ) printPDFWarning();
	else if( foundPDFMuonPxPyPz ) {
	   foundPDF.get()->insert({PDF_MuonPx, foundPDFMuonPx});
	   foundPDF.get()->insert({PDF_MuonPy, foundPDFMuonPy});
	   foundPDF.get()->insert({PDF_MuonPz, foundPDFMuonPz});
	}
	else {
	   foundPDF.get()->insert({PDF_MuonPt, foundPDFMuonPt});
	   foundPDF.get()->insert({PDF_MuonEta, foundPDFMuonEta});
	   foundPDF.get()->insert({PDF_MuonPhi, foundPDFMuonPhi});
	}

	if( foundPDFElecPxPyPz && foundPDFMuonPxPyPz ) {
	   (*usePDFLeptonPxPyPz) = true;
	   (*usePDFLeptonPtEtaPhi) = false;
	}
	else if( foundPDFElecPtEtaPhi && foundPDFMuonPtEtaPhi ) {
	   (*usePDFLeptonPxPyPz) = false;
	   (*usePDFLeptonPtEtaPhi) = true;
	}
	else {
	   std::cout << "Please define the same type PDFs for leptons" << std::endl;
	   exit(1);
	}
	
	DoFit_ = false;
	
	NNu_ = 1; ///< Number of neutrinos in the final state
     }
   
   rnd = new TRandom3(777);
   
   ISINITIALIZED_ = 1;
}

/// Read transfer functions from file
void KINFIT::kfit::SetPDF(std::string obj, std::string fileName, std::string hName, bool isDeltaFunc)
{
   if( isDeltaFunc )
     {
	if( strcmp(obj.c_str(), "TopWMass") == 0 ) isDeltaFuncPDFTopWMass = true;
	else if( strcmp(obj.c_str(), "TopMass") == 0 ) isDeltaFuncPDFTopMass = true;
	else if( strcmp(obj.c_str(), "TopWHadMass") == 0 ) isDeltaFuncPDFTopWHadMass = true;
	else if( strcmp(obj.c_str(), "TopHadMass") == 0 ) isDeltaFuncPDFTopHadMass = true;
	else if( strcmp(obj.c_str(), "MetPx") == 0 ) isDeltaFuncPDFMetPx = true;
	else if( strcmp(obj.c_str(), "MetPy") == 0 ) isDeltaFuncPDFMetPy = true;
	else if( strcmp(obj.c_str(), "BJetPx") == 0 ) isDeltaFuncPDFBJetPx = true;
	else if( strcmp(obj.c_str(), "BJetPy") == 0 ) isDeltaFuncPDFBJetPy = true;
	else if( strcmp(obj.c_str(), "BJetPz") == 0 ) isDeltaFuncPDFBJetPz = true;
	else if( strcmp(obj.c_str(), "BJetPt") == 0 ) isDeltaFuncPDFBJetPt = true;
	else if( strcmp(obj.c_str(), "BJetEta") == 0 ) isDeltaFuncPDFBJetEta = true;
	else if( strcmp(obj.c_str(), "BJetPhi") == 0 ) isDeltaFuncPDFBJetPhi = true;
	else if( strcmp(obj.c_str(), "NonBJetPx") == 0 ) isDeltaFuncPDFNonBJetPx = true;
	else if( strcmp(obj.c_str(), "NonBJetPy") == 0 ) isDeltaFuncPDFNonBJetPy = true;
	else if( strcmp(obj.c_str(), "NonBJetPz") == 0 ) isDeltaFuncPDFNonBJetPz = true;
	else if( strcmp(obj.c_str(), "NonBJetPt") == 0 ) isDeltaFuncPDFNonBJetPt = true;
	else if( strcmp(obj.c_str(), "NonBJetEta") == 0 ) isDeltaFuncPDFNonBJetEta = true;
	else if( strcmp(obj.c_str(), "NonBJetPhi") == 0 ) isDeltaFuncPDFNonBJetPhi = true;
	else if( strcmp(obj.c_str(), "ElecPx") == 0 ) isDeltaFuncPDFElecPx = true;
	else if( strcmp(obj.c_str(), "ElecPy") == 0 ) isDeltaFuncPDFElecPy = true;
	else if( strcmp(obj.c_str(), "ElecPz") == 0 ) isDeltaFuncPDFElecPz = true;
	else if( strcmp(obj.c_str(), "ElecPt") == 0 ) isDeltaFuncPDFElecPt = true;
	else if( strcmp(obj.c_str(), "ElecEta") == 0 ) isDeltaFuncPDFElecEta = true;
	else if( strcmp(obj.c_str(), "ElecPhi") == 0 ) isDeltaFuncPDFElecPhi = true;
	else if( strcmp(obj.c_str(), "MuonPx") == 0 ) isDeltaFuncPDFMuonPx = true;
	else if( strcmp(obj.c_str(), "MuonPy") == 0 ) isDeltaFuncPDFMuonPy = true;
	else if( strcmp(obj.c_str(), "MuonPz") == 0 ) isDeltaFuncPDFMuonPz = true;
	else if( strcmp(obj.c_str(), "MuonPt") == 0 ) isDeltaFuncPDFMuonPt = true;
	else if( strcmp(obj.c_str(), "MuonEta") == 0 ) isDeltaFuncPDFMuonEta = true;
	else if( strcmp(obj.c_str(), "MuonPhi") == 0 ) isDeltaFuncPDFMuonPhi = true;
	else if( strcmp(obj.c_str(), "PhotonPx") == 0 ) isDeltaFuncPDFPhotonPx = true;
	else if( strcmp(obj.c_str(), "PhotonPy") == 0 ) isDeltaFuncPDFPhotonPy = true;
	else if( strcmp(obj.c_str(), "PhotonPz") == 0 ) isDeltaFuncPDFPhotonPz = true;
	else if( strcmp(obj.c_str(), "PhotonPt") == 0 ) isDeltaFuncPDFPhotonPt = true;
	else if( strcmp(obj.c_str(), "PhotonEta") == 0 ) isDeltaFuncPDFPhotonEta = true;
	else if( strcmp(obj.c_str(), "PhotonPhi") == 0 ) isDeltaFuncPDFPhotonPhi = true;
	else
	  {
	     std::cout << "Provided transfer function " << obj << " is not defined for this process" << std::endl;
	     exit(1);
	  }
	
	return;
     }   
   
   TFile *fPDF = TFile::Open(fileName.c_str());
   if( ! fPDF->IsOpen() )
     {
	std::cout << "File " << fileName << " does not exist" << std::endl;
	exit(1);
     }   
   if( ! fPDF->GetListOfKeys()->Contains(hName.c_str()) )
     {
	std::cout << "Histogram " << hName << " does not exist" << std::endl;
	exit(1);
     }
   
   TF1 *hPDF = (TF1*)fPDF->Get(hName.c_str());

   if( strcmp(obj.c_str(), "TopWMass") == 0 )
     {
	hPDFTopWMass = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_TopWMass].c_str()))) );

	maxPDFTopWMass = hPDFTopWMass->GetMaximum();
	meanPDFTopWMass = hPDFTopWMass->GetMaximumX();
	sigmaPDFTopWMass = fabs(hPDFTopWMass->GetX(maxPDFTopWMass/2.));
	hPDFTopWMass->GetRange(xminPDFTopWMass, xmaxPDFTopWMass);
     }
   else if( strcmp(obj.c_str(), "TopMass") == 0 )
     {
	hPDFTopMass = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_TopMass].c_str()))) );

	maxPDFTopMass = hPDFTopMass->GetMaximum();
	meanPDFTopMass = hPDFTopMass->GetMaximumX();
	sigmaPDFTopMass = fabs(hPDFTopMass->GetX(maxPDFTopMass/2.));
	hPDFTopMass->GetRange(xminPDFTopMass, xmaxPDFTopMass);
     }
   else if( strcmp(obj.c_str(), "TopWHadMass") == 0 )
     {
	hPDFTopWHadMass = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_TopWHadMass].c_str()))) );

	maxPDFTopWHadMass = hPDFTopWHadMass->GetMaximum();
	meanPDFTopWHadMass = hPDFTopWHadMass->GetMaximumX();
	sigmaPDFTopWHadMass = fabs(hPDFTopWHadMass->GetX(maxPDFTopWHadMass/2.));
	hPDFTopWHadMass->GetRange(xminPDFTopWHadMass, xmaxPDFTopWHadMass);
     }
   else if( strcmp(obj.c_str(), "TopHadMass") == 0 )
     {
	hPDFTopHadMass = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_TopHadMass].c_str()))) );

	maxPDFTopHadMass = hPDFTopHadMass->GetMaximum();
	meanPDFTopHadMass = hPDFTopHadMass->GetMaximumX();
	sigmaPDFTopHadMass = fabs(hPDFTopHadMass->GetX(maxPDFTopHadMass/2.));
	hPDFTopHadMass->GetRange(xminPDFTopHadMass, xmaxPDFTopHadMass);
     }
//   else if( strcmp(obj.c_str(), "TopTopMass") == 0 )
//     {
//	hPDFTopTopMass = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_TopTopMass].c_str()))) );
//
//	maxPDFTopTopMass = hPDFTopTopMass->GetMaximum();
//	meanPDFTopTopMass = hPDFTopTopMass->GetMaximumX();
//	sigmaPDFTopTopMass = fabs(hPDFTopTopMass->GetX(maxPDFTopTopMass/2.));
//	hPDFTopTopMass->GetRange(xminPDFTopTopMass, xmaxPDFTopTopMass);
//     }
   else if( strcmp(obj.c_str(), "MetPx") == 0 )
     {
	hPDFMetPx = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_MetPx].c_str()))) );

	maxPDFMetPx = hPDFMetPx->GetMaximum();
	meanPDFMetPx = hPDFMetPx->GetMaximumX();
	sigmaPDFMetPx = fabs(hPDFMetPx->GetX(maxPDFMetPx/2.));
	hPDFMetPx->GetRange(xminPDFMetPx, xmaxPDFMetPx);
     }
   else if( strcmp(obj.c_str(), "MetPy") == 0 )
     {
	hPDFMetPy = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_MetPy].c_str()))) );

	maxPDFMetPy = hPDFMetPy->GetMaximum();
	meanPDFMetPy = hPDFMetPy->GetMaximumX();
	sigmaPDFMetPy = fabs(hPDFMetPy->GetX(maxPDFMetPy/2.));
	hPDFMetPy->GetRange(xminPDFMetPy, xmaxPDFMetPy);
     }   
   else if( strcmp(obj.c_str(), "BJetPx") == 0 )
     {
	hPDFBJetPx = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_BJetPx].c_str()))) );

	maxPDFBJetPx = hPDFBJetPx->GetMaximum();
	meanPDFBJetPx = hPDFBJetPx->GetMaximumX();
	sigmaPDFBJetPx = fabs(hPDFBJetPx->GetX(maxPDFBJetPx/2.));
	hPDFBJetPx->GetRange(xminPDFBJetPx, xmaxPDFBJetPx);
     }   
   else if( strcmp(obj.c_str(), "BJetPy") == 0 )
     {
	hPDFBJetPy = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_BJetPy].c_str()))) );

	maxPDFBJetPy = hPDFBJetPy->GetMaximum();
	meanPDFBJetPy = hPDFBJetPy->GetMaximumX();
	sigmaPDFBJetPy = fabs(hPDFBJetPy->GetX(maxPDFBJetPy/2.));
	hPDFBJetPy->GetRange(xminPDFBJetPy, xmaxPDFBJetPy);
     }   
   else if( strcmp(obj.c_str(), "BJetPz") == 0 )
     {
	hPDFBJetPz = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_BJetPz].c_str()))) );

	maxPDFBJetPz = hPDFBJetPz->GetMaximum();
	meanPDFBJetPz = hPDFBJetPz->GetMaximumX();
	sigmaPDFBJetPz = fabs(hPDFBJetPz->GetX(maxPDFBJetPz/2.));
	hPDFBJetPz->GetRange(xminPDFBJetPz, xmaxPDFBJetPz);
     }
   else if( strcmp(obj.c_str(), "BJetPt") == 0 )
     {
	hPDFBJetPt = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_BJetPt].c_str()))) );

	maxPDFBJetPt = hPDFBJetPt->GetMaximum();
	meanPDFBJetPt = hPDFBJetPt->GetMaximumX();
	sigmaPDFBJetPt = fabs(hPDFBJetPt->GetX(maxPDFBJetPt/2.));
	hPDFBJetPt->GetRange(xminPDFBJetPt, xmaxPDFBJetPt);
     }
   else if( strcmp(obj.c_str(), "BJetEta") == 0 )
     {
	hPDFBJetEta = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_BJetEta].c_str()))) );

	maxPDFBJetEta = hPDFBJetEta->GetMaximum();
	meanPDFBJetEta = hPDFBJetEta->GetMaximumX();
	sigmaPDFBJetEta = fabs(hPDFBJetEta->GetX(maxPDFBJetEta/2.));
	hPDFBJetEta->GetRange(xminPDFBJetEta, xmaxPDFBJetEta);
     }
   else if( strcmp(obj.c_str(), "BJetPhi") == 0 )
     {
	hPDFBJetPhi = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_BJetPhi].c_str()))) );

	maxPDFBJetPhi = hPDFBJetPhi->GetMaximum();
	meanPDFBJetPhi = hPDFBJetPhi->GetMaximumX();
	sigmaPDFBJetPhi = fabs(hPDFBJetPhi->GetX(maxPDFBJetPhi/2.));
	hPDFBJetPhi->GetRange(xminPDFBJetPhi, xmaxPDFBJetPhi);
     }
   else if( strcmp(obj.c_str(), "NonBJetPx") == 0 )
     {
	hPDFNonBJetPx = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_NonBJetPx].c_str()))) );

	maxPDFNonBJetPx = hPDFNonBJetPx->GetMaximum();
	meanPDFNonBJetPx = hPDFNonBJetPx->GetMaximumX();
	sigmaPDFNonBJetPx = fabs(hPDFNonBJetPx->GetX(maxPDFNonBJetPx/2.));
	hPDFNonBJetPx->GetRange(xminPDFNonBJetPx, xmaxPDFNonBJetPx);
     }   
   else if( strcmp(obj.c_str(), "NonBJetPy") == 0 )
     {
	hPDFNonBJetPy = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_NonBJetPy].c_str()))) );

	maxPDFNonBJetPy = hPDFNonBJetPy->GetMaximum();
	meanPDFNonBJetPy = hPDFNonBJetPy->GetMaximumX();
	sigmaPDFNonBJetPy = fabs(hPDFNonBJetPy->GetX(maxPDFNonBJetPy/2.));
	hPDFNonBJetPy->GetRange(xminPDFNonBJetPy, xmaxPDFNonBJetPy);
     }   
   else if( strcmp(obj.c_str(), "NonBJetPz") == 0 )
     {
	hPDFNonBJetPz = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_NonBJetPz].c_str()))) );

	maxPDFNonBJetPz = hPDFNonBJetPz->GetMaximum();
	meanPDFNonBJetPz = hPDFNonBJetPz->GetMaximumX();
	sigmaPDFNonBJetPz = fabs(hPDFNonBJetPz->GetX(maxPDFNonBJetPz/2.));
	hPDFNonBJetPz->GetRange(xminPDFNonBJetPz, xmaxPDFNonBJetPz);
     }
   else if( strcmp(obj.c_str(), "NonBJetPt") == 0 )
     {
	hPDFNonBJetPt = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_NonBJetPt].c_str()))) );

	maxPDFNonBJetPt = hPDFNonBJetPt->GetMaximum();
	meanPDFNonBJetPt = hPDFNonBJetPt->GetMaximumX();
	sigmaPDFNonBJetPt = fabs(hPDFNonBJetPt->GetX(maxPDFNonBJetPt/2.));
	hPDFNonBJetPt->GetRange(xminPDFNonBJetPt, xmaxPDFNonBJetPt);
     }
   else if( strcmp(obj.c_str(), "NonBJetEta") == 0 )
     {
	hPDFNonBJetEta = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_NonBJetEta].c_str()))) );

	maxPDFNonBJetEta = hPDFNonBJetEta->GetMaximum();
	meanPDFNonBJetEta = hPDFNonBJetEta->GetMaximumX();
	sigmaPDFNonBJetEta = fabs(hPDFNonBJetEta->GetX(maxPDFNonBJetEta/2.));
	hPDFNonBJetEta->GetRange(xminPDFNonBJetEta, xmaxPDFNonBJetEta);
     }
   else if( strcmp(obj.c_str(), "NonBJetPhi") == 0 )
     {
	hPDFNonBJetPhi = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_NonBJetPhi].c_str()))) );

	maxPDFNonBJetPhi = hPDFNonBJetPhi->GetMaximum();
	meanPDFNonBJetPhi = hPDFNonBJetPhi->GetMaximumX();
	sigmaPDFNonBJetPhi = fabs(hPDFNonBJetPhi->GetX(maxPDFNonBJetPhi/2.));
	hPDFNonBJetPhi->GetRange(xminPDFNonBJetPhi, xmaxPDFNonBJetPhi);
     }
   else if( strcmp(obj.c_str(), "ElecPx") == 0 )
     {
	hPDFElecPx = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_ElecPx].c_str()))) );

	maxPDFElecPx = hPDFElecPx->GetMaximum();
	meanPDFElecPx = hPDFElecPx->GetMaximumX();
	sigmaPDFElecPx = fabs(hPDFElecPx->GetX(maxPDFElecPx/2.));
	hPDFElecPx->GetRange(xminPDFElecPx, xmaxPDFElecPx);
     }   
   else if( strcmp(obj.c_str(), "ElecPy") == 0 )
     {
	hPDFElecPy = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_ElecPy].c_str()))) );

	maxPDFElecPy = hPDFElecPy->GetMaximum();
	meanPDFElecPy = hPDFElecPy->GetMaximumX();
	sigmaPDFElecPy = fabs(hPDFElecPy->GetX(maxPDFElecPy/2.));
	hPDFElecPy->GetRange(xminPDFElecPy, xmaxPDFElecPy);
     }	    
   else if( strcmp(obj.c_str(), "ElecPz") == 0 )
     {
	hPDFElecPz = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_ElecPz].c_str()))) );
	
	maxPDFElecPz = hPDFElecPz->GetMaximum();
	meanPDFElecPz = hPDFElecPz->GetMaximumX();
	sigmaPDFElecPz = fabs(hPDFElecPz->GetX(maxPDFElecPz/2.));
	hPDFElecPz->GetRange(xminPDFElecPz, xmaxPDFElecPz);
     }
   else if( strcmp(obj.c_str(), "ElecPt") == 0 )
     {
	hPDFElecPt = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_ElecPt].c_str()))) );
	
	maxPDFElecPt = hPDFElecPt->GetMaximum();
	meanPDFElecPt = hPDFElecPt->GetMaximumX();
	sigmaPDFElecPt = fabs(hPDFElecPt->GetX(maxPDFElecPt/2.));
	hPDFElecPt->GetRange(xminPDFElecPt, xmaxPDFElecPt);
     }
   else if( strcmp(obj.c_str(), "ElecEta") == 0 )
     {
	hPDFElecEta = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_ElecEta].c_str()))) );
	
	maxPDFElecEta = hPDFElecEta->GetMaximum();
	meanPDFElecEta = hPDFElecEta->GetMaximumX();
	sigmaPDFElecEta = fabs(hPDFElecEta->GetX(maxPDFElecEta/2.));
	hPDFElecEta->GetRange(xminPDFElecEta, xmaxPDFElecEta);
     }
   else if( strcmp(obj.c_str(), "ElecPhi") == 0 )
     {
	hPDFElecPhi = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_ElecPhi].c_str()))) );
	
	maxPDFElecPhi = hPDFElecPhi->GetMaximum();
	meanPDFElecPhi = hPDFElecPhi->GetMaximumX();
	sigmaPDFElecPhi = fabs(hPDFElecPhi->GetX(maxPDFElecPhi/2.));
	hPDFElecPhi->GetRange(xminPDFElecPhi, xmaxPDFElecPhi);
     }
   else if( strcmp(obj.c_str(), "MuonPx") == 0 )
     {
	hPDFMuonPx = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_MuonPx].c_str()))) );

	maxPDFMuonPx = hPDFMuonPx->GetMaximum();
	meanPDFMuonPx = hPDFMuonPx->GetMaximumX();
	sigmaPDFMuonPx = fabs(hPDFMuonPx->GetX(maxPDFMuonPx/2.));
	hPDFMuonPx->GetRange(xminPDFMuonPx, xmaxPDFMuonPx);
     }   
   else if( strcmp(obj.c_str(), "MuonPy") == 0 )
     {
	hPDFMuonPy = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_MuonPy].c_str()))) );

	maxPDFMuonPy = hPDFMuonPy->GetMaximum();
	meanPDFMuonPy = hPDFMuonPy->GetMaximumX();
	sigmaPDFMuonPy = fabs(hPDFMuonPy->GetX(maxPDFMuonPy/2.));
	hPDFMuonPy->GetRange(xminPDFMuonPy, xmaxPDFMuonPy);
     }   
   else if( strcmp(obj.c_str(), "MuonPz") == 0 )
     {
	hPDFMuonPz = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_MuonPz].c_str()))) );

	maxPDFMuonPz = hPDFMuonPz->GetMaximum();
	meanPDFMuonPz = hPDFMuonPz->GetMaximumX();
	sigmaPDFMuonPz = fabs(hPDFMuonPz->GetX(maxPDFMuonPz/2.));
	hPDFMuonPz->GetRange(xminPDFMuonPz, xmaxPDFMuonPz);
     }
   else if( strcmp(obj.c_str(), "MuonPt") == 0 )
     {
	hPDFMuonPt = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_MuonPt].c_str()))) );

	maxPDFMuonPt = hPDFMuonPt->GetMaximum();
	meanPDFMuonPt = hPDFMuonPt->GetMaximumX();
	sigmaPDFMuonPt = fabs(hPDFMuonPt->GetX(maxPDFMuonPt/2.));
	hPDFMuonPt->GetRange(xminPDFMuonPt, xmaxPDFMuonPt);
     }
   else if( strcmp(obj.c_str(), "MuonEta") == 0 )
     {
	hPDFMuonEta = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_MuonEta].c_str()))) );

	maxPDFMuonEta = hPDFMuonEta->GetMaximum();
	meanPDFMuonEta = hPDFMuonEta->GetMaximumX();
	sigmaPDFMuonEta = fabs(hPDFMuonEta->GetX(maxPDFMuonEta/2.));
	hPDFMuonEta->GetRange(xminPDFMuonEta, xmaxPDFMuonEta);
     }
   else if( strcmp(obj.c_str(), "MuonPhi") == 0 )
     {
	hPDFMuonPhi = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_MuonPhi].c_str()))) );

	maxPDFMuonPhi = hPDFMuonPhi->GetMaximum();
	meanPDFMuonPhi = hPDFMuonPhi->GetMaximumX();
	sigmaPDFMuonPhi = fabs(hPDFMuonPhi->GetX(maxPDFMuonPhi/2.));
	hPDFMuonPhi->GetRange(xminPDFMuonPhi, xmaxPDFMuonPhi);
     }
   else if( strcmp(obj.c_str(), "PhotonPx") == 0 )
     {
	hPDFPhotonPx = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_PhotonPx].c_str()))) );

	maxPDFPhotonPx = hPDFPhotonPx->GetMaximum();
	meanPDFPhotonPx = hPDFPhotonPx->GetMaximumX();
	sigmaPDFPhotonPx = fabs(hPDFPhotonPx->GetX(maxPDFPhotonPx/2.));
	hPDFPhotonPx->GetRange(xminPDFPhotonPx, xmaxPDFPhotonPx);
     }   
   else if( strcmp(obj.c_str(), "PhotonPy") == 0 )
     {
	hPDFPhotonPy = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_PhotonPy].c_str()))) );

	maxPDFPhotonPy = hPDFPhotonPy->GetMaximum();
	meanPDFPhotonPy = hPDFPhotonPy->GetMaximumX();
	sigmaPDFPhotonPy = fabs(hPDFPhotonPy->GetX(maxPDFPhotonPy/2.));
	hPDFPhotonPy->GetRange(xminPDFPhotonPy, xmaxPDFPhotonPy);
     }   
   else if( strcmp(obj.c_str(), "PhotonPz") == 0 )
     {
	hPDFPhotonPz = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_PhotonPz].c_str()))) );

	maxPDFPhotonPz = hPDFPhotonPz->GetMaximum();
	meanPDFPhotonPz = hPDFPhotonPz->GetMaximumX();
	sigmaPDFPhotonPz = fabs(hPDFPhotonPz->GetX(maxPDFPhotonPz/2.));
	hPDFPhotonPz->GetRange(xminPDFPhotonPz, xmaxPDFPhotonPz);
     }
   else if( strcmp(obj.c_str(), "PhotonPt") == 0 )
     {
	hPDFPhotonPt = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_PhotonPt].c_str()))) );

	maxPDFPhotonPt = hPDFPhotonPt->GetMaximum();
	meanPDFPhotonPt = hPDFPhotonPt->GetMaximumX();
	sigmaPDFPhotonPt = fabs(hPDFPhotonPt->GetX(maxPDFPhotonPt/2.));
	hPDFPhotonPt->GetRange(xminPDFPhotonPt, xmaxPDFPhotonPt);
     }
   else if( strcmp(obj.c_str(), "PhotonEta") == 0 )
     {
	hPDFPhotonEta = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_PhotonEta].c_str()))) );

	maxPDFPhotonEta = hPDFPhotonEta->GetMaximum();
	meanPDFPhotonEta = hPDFPhotonEta->GetMaximumX();
	sigmaPDFPhotonEta = fabs(hPDFPhotonEta->GetX(maxPDFPhotonEta/2.));
	hPDFPhotonEta->GetRange(xminPDFPhotonEta, xmaxPDFPhotonEta);
     }
   else if( strcmp(obj.c_str(), "PhotonPhi") == 0 )
     {
	hPDFPhotonPhi = std::shared_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone(PDF_NAME[PDF_PhotonPhi].c_str()))) );

	maxPDFPhotonPhi = hPDFPhotonPhi->GetMaximum();
	meanPDFPhotonPhi = hPDFPhotonPhi->GetMaximumX();
	sigmaPDFPhotonPhi = fabs(hPDFPhotonPhi->GetX(maxPDFPhotonPhi/2.));
	hPDFPhotonPhi->GetRange(xminPDFPhotonPhi, xmaxPDFPhotonPhi);
     }
   else
     {
	std::cout << "Provided transfer function " << obj << " is not defined for this process" << std::endl;
	exit(1);
     }

   if( fPDF->IsOpen() ) fPDF->Close();
}

/// Main routine
void KINFIT::kfit::Run()
{
   /// Reset all variables used in the fit
   ResetFitInfo();

   /// Put electrons and muons into one collection
   FillLeptons();
   
   /// Define process
   if( hypoMode == TOPTOPLEPLEP ) TopTopLepLep();
   else if( hypoMode == TOPTOPLEPHAD ) TopTopLepHad();
   else if( hypoMode == TOPLEP ) TopLep();
   
   /// Sort permutations according to the NLL results
   sortPermIndex();
   
   /// Clean start at next event
   Reset();
}

/// Create a single collection for all leptons
void KINFIT::kfit::FillLeptons()
{
   for(int i=0;i<nElectron;i++)
     {	
	LeptonLabel.push_back(0);
	LeptonIdx.push_back(i);
     }   
   
   nLepton += nElectron;
   LeptonPt.insert(LeptonPt.end(),ElectronPt.begin(),ElectronPt.end());
   LeptonEta.insert(LeptonEta.end(),ElectronEta.begin(),ElectronEta.end());
   LeptonPhi.insert(LeptonPhi.end(),ElectronPhi.begin(),ElectronPhi.end());
   LeptonE.insert(LeptonE.end(),ElectronE.begin(),ElectronE.end());   
   LeptonPx.insert(LeptonPx.end(),ElectronPx.begin(),ElectronPx.end());
   LeptonPy.insert(LeptonPy.end(),ElectronPy.begin(),ElectronPy.end());
   LeptonPz.insert(LeptonPz.end(),ElectronPz.begin(),ElectronPz.end());
   LeptonCharge.insert(LeptonCharge.end(),ElectronCharge.begin(),ElectronCharge.end());

   for(int i=0;i<nMuon;i++)
     {	
	LeptonLabel.push_back(1);
	LeptonIdx.push_back(i);
     }   
   
   nLepton += nMuon;
   LeptonPt.insert(LeptonPt.end(),MuonPt.begin(),MuonPt.end());
   LeptonEta.insert(LeptonEta.end(),MuonEta.begin(),MuonEta.end());
   LeptonPhi.insert(LeptonPhi.end(),MuonPhi.begin(),MuonPhi.end());
   LeptonE.insert(LeptonE.end(),MuonE.begin(),MuonE.end());
   LeptonPx.insert(LeptonPx.end(),MuonPx.begin(),MuonPx.end());
   LeptonPy.insert(LeptonPy.end(),MuonPy.begin(),MuonPy.end());
   LeptonPz.insert(LeptonPz.end(),MuonPz.begin(),MuonPz.end());
   LeptonCharge.insert(LeptonCharge.end(),MuonCharge.begin(),MuonCharge.end());
}

/// Restore fit settings
void KINFIT::kfit::ResetFitInfo()
{
   ParMin.reset(); ParMin = std::make_shared<std::vector<float> >(std::vector<float>(FPARAM_N, 1.));
   ParMax.reset(); ParMax = std::make_shared<std::vector<float> >(std::vector<float>(FPARAM_N, 1.));
}

/// Clean interim fit variables
void KINFIT::kfit::Reset()
{
   CHISQ.reset(); CHISQ = std::make_shared<double>(0);
   
   EtMissX.reset(); EtMissX = std::make_shared<double>(0);
   EtMissY.reset(); EtMissY = std::make_shared<double>(0);
   
   PxLepton1.reset(); PxLepton1 = std::make_shared<double>(0);
   PyLepton1.reset(); PyLepton1 = std::make_shared<double>(0);
   PzLepton1.reset(); PzLepton1 = std::make_shared<double>(0);
   ELepton1.reset(); ELepton1 = std::make_shared<double>(0);
   PtLepton1.reset(); PtLepton1 = std::make_shared<double>(0);
   EtaLepton1.reset(); EtaLepton1 = std::make_shared<double>(0);
   PhiLepton1.reset(); PhiLepton1 = std::make_shared<double>(0);
   MassLepton1.reset(); MassLepton1 = std::make_shared<double>(0);
   LabelLepton1.reset(); LabelLepton1 = std::make_shared<int>(0);
   
   PxLepton2.reset(); PxLepton2 = std::make_shared<double>(0);
   PyLepton2.reset(); PyLepton2 = std::make_shared<double>(0);
   PzLepton2.reset(); PzLepton2 = std::make_shared<double>(0);
   ELepton2.reset(); ELepton2 = std::make_shared<double>(0);
   PtLepton2.reset(); PtLepton2 = std::make_shared<double>(0);
   EtaLepton2.reset(); EtaLepton2 = std::make_shared<double>(0);
   PhiLepton2.reset(); PhiLepton2 = std::make_shared<double>(0);
   MassLepton2.reset(); MassLepton2 = std::make_shared<double>(0);
   LabelLepton2.reset(); LabelLepton2 = std::make_shared<int>(0);
   
   PxBJet1.reset(); PxBJet1 = std::make_shared<double>(0);
   PyBJet1.reset(); PyBJet1 = std::make_shared<double>(0);
   PzBJet1.reset(); PzBJet1 = std::make_shared<double>(0);
   EBJet1.reset(); EBJet1 = std::make_shared<double>(0);
   PtBJet1.reset(); PtBJet1 = std::make_shared<double>(0);
   EtaBJet1.reset(); EtaBJet1 = std::make_shared<double>(0);
   PhiBJet1.reset(); PhiBJet1 = std::make_shared<double>(0);
   MassBJet1.reset(); MassBJet1 = std::make_shared<double>(0);
   
   PxBJet2.reset(); PxBJet2 = std::make_shared<double>(0);
   PyBJet2.reset(); PyBJet2 = std::make_shared<double>(0);
   PzBJet2.reset(); PzBJet2 = std::make_shared<double>(0);
   EBJet2.reset(); EBJet2 = std::make_shared<double>(0);
   PtBJet2.reset(); PtBJet2 = std::make_shared<double>(0);
   EtaBJet2.reset(); EtaBJet2 = std::make_shared<double>(0);
   PhiBJet2.reset(); PhiBJet2 = std::make_shared<double>(0);
   MassBJet2.reset(); MassBJet2 = std::make_shared<double>(0);

   PxNonBJet1.reset(); PxNonBJet1 = std::make_shared<double>(0);
   PyNonBJet1.reset(); PyNonBJet1 = std::make_shared<double>(0);
   PzNonBJet1.reset(); PzNonBJet1 = std::make_shared<double>(0);
   ENonBJet1.reset(); ENonBJet1 = std::make_shared<double>(0);
   PtNonBJet1.reset(); PtNonBJet1 = std::make_shared<double>(0);
   EtaNonBJet1.reset(); EtaNonBJet1 = std::make_shared<double>(0);
   PhiNonBJet1.reset(); PhiNonBJet1 = std::make_shared<double>(0);
   MassNonBJet1.reset(); MassNonBJet1 = std::make_shared<double>(0);
   
   PxNonBJet2.reset(); PxNonBJet2 = std::make_shared<double>(0);
   PyNonBJet2.reset(); PyNonBJet2 = std::make_shared<double>(0);
   PzNonBJet2.reset(); PzNonBJet2 = std::make_shared<double>(0);
   ENonBJet2.reset(); ENonBJet2 = std::make_shared<double>(0);
   PtNonBJet2.reset(); PtNonBJet2 = std::make_shared<double>(0);
   EtaNonBJet2.reset(); EtaNonBJet2 = std::make_shared<double>(0);
   PhiNonBJet2.reset(); PhiNonBJet2 = std::make_shared<double>(0);
   MassNonBJet2.reset(); MassNonBJet2 = std::make_shared<double>(0);

   PxPhoton.reset(); PxPhoton = std::make_shared<double>(0);
   PyPhoton.reset(); PyPhoton = std::make_shared<double>(0);
   PzPhoton.reset(); PzPhoton = std::make_shared<double>(0);
   EPhoton.reset(); EPhoton = std::make_shared<double>(0);
   PtPhoton.reset(); PtPhoton = std::make_shared<double>(0);
   EtaPhoton.reset(); EtaPhoton = std::make_shared<double>(0);
   PhiPhoton.reset(); PhiPhoton = std::make_shared<double>(0);
   PhotonOrigin.reset(); PhotonOrigin = std::make_shared<int>(0);
   
   PxNu1.reset(); PxNu1 = std::make_shared<double>(0);
   PyNu1.reset(); PyNu1 = std::make_shared<double>(0);
   PzNu1.reset(); PzNu1 = std::make_shared<double>(0);
   
   PxNu2.reset(); PxNu2 = std::make_shared<double>(0);
   PyNu2.reset(); PyNu2 = std::make_shared<double>(0);
   PzNu2.reset(); PzNu2 = std::make_shared<double>(0);
   
   TopTopMass.reset(); TopTopMass = std::make_shared<double>(0);
   TopTopPz.reset(); TopTopPz = std::make_shared<double>(0);
   TopTopPzAbs.reset(); TopTopPzAbs = std::make_shared<double>(0);
   
   FitParam.reset(); FitParam = std::make_shared<std::vector<double> >();
   ChiTerm.reset(); ChiTerm = std::make_shared<std::vector<double> >();
   ChiTermName.reset(); ChiTermName = std::make_shared<std::vector<std::string> >();
   
   nLepton = 0;   
   LeptonPt.clear();
   LeptonEta.clear();
   LeptonPhi.clear();
   LeptonE.clear();
   LeptonPx.clear();
   LeptonPy.clear();
   LeptonPz.clear();
   LeptonCharge.clear();
   LeptonLabel.clear();
   LeptonIdx.clear();

   nElectron = 0;
   ElectronPt.clear();
   ElectronEta.clear();
   ElectronPhi.clear();
   ElectronE.clear();
   ElectronPx.clear();
   ElectronPy.clear();
   ElectronPz.clear();
   ElectronCharge.clear();

   nMuon = 0;
   MuonPt.clear();
   MuonEta.clear();
   MuonPhi.clear();
   MuonE.clear();
   MuonPx.clear();
   MuonPy.clear();
   MuonPz.clear();
   MuonCharge.clear();

   nBJet = 0;
   BJetPt.clear();
   BJetEta.clear();
   BJetPhi.clear();
   BJetE.clear();
   BJetPx.clear();
   BJetPy.clear();
   BJetPz.clear();

   nNonBJet = 0;
   NonBJetPt.clear();
   NonBJetEta.clear();
   NonBJetPhi.clear();
   NonBJetE.clear();
   NonBJetPx.clear();
   NonBJetPy.clear();
   NonBJetPz.clear();

   nPhoton = 0;
   PhotonPt.clear();
   PhotonEta.clear();
   PhotonPhi.clear();
   PhotonE.clear();
   PhotonPx.clear();
   PhotonPy.clear();
   PhotonPz.clear();
   
   IncludePhotons_ = false;
   CheckAllPhotonOrigins_ = false;
}

/// Call the run method
void KINFIT::kfit::TopTopLepLep()
{
   hypTopTopLepLep->TopTopLepLepRun();
}

/// Call the run method
void KINFIT::kfit::TopTopLepHad()
{
   hypTopTopLepHad->TopTopLepHadRun();
}

/// Call the run method
void KINFIT::kfit::TopLep()
{
   hypTopLep->TopLepRun();
}

/// Match objects in input collections to the fit result
int KINFIT::kfit::GetIndex(OBJ objType, int idx)
{
   // TOPTOPLEPLEP
   if( objType == ELECTRON1_TOPTOPLEPLEP )       return getPermIndex(idx, TopTopLepLep_Electron1Idx);
   else if( objType == MUON1_TOPTOPLEPLEP )      return getPermIndex(idx, TopTopLepLep_Muon1Idx);
   else if( objType == ELECTRON2_TOPTOPLEPLEP )  return getPermIndex(idx, TopTopLepLep_Electron2Idx);
   else if( objType == MUON2_TOPTOPLEPLEP )      return getPermIndex(idx, TopTopLepLep_Muon2Idx);
   else if( objType == BJET1_TOPTOPLEPLEP )      return getPermIndex(idx, TopTopLepLep_BJet1Idx);
   else if( objType == BJET2_TOPTOPLEPLEP )      return getPermIndex(idx, TopTopLepLep_BJet2Idx);
   else if( objType == PHOTON_TOPTOPLEPLEP )     return getPermIndex(idx, TopTopLepLep_PhotonIdx);

   // TOPTOPLEPHAD
   if( objType == ELECTRON_TOPTOPLEPHAD )        return getPermIndex(idx, TopTopLepHad_ElectronIdx);
   else if( objType == MUON_TOPTOPLEPHAD )       return getPermIndex(idx, TopTopLepHad_MuonIdx);
   else if( objType == BJETLEP_TOPTOPLEPHAD )    return getPermIndex(idx, TopTopLepHad_BJetLepIdx);
   else if( objType == BJETHAD_TOPTOPLEPHAD )    return getPermIndex(idx, TopTopLepHad_BJetHadIdx);
   else if( objType == NONBJET1_TOPTOPLEPHAD )   return getPermIndex(idx, TopTopLepHad_NonBJet1Idx);
   else if( objType == NONBJET2_TOPTOPLEPHAD )   return getPermIndex(idx, TopTopLepHad_NonBJet2Idx);
   else if( objType == PHOTON_TOPTOPLEPHAD )     return getPermIndex(idx, TopTopLepHad_PhotonIdx);

   // TOPLEP
   if( objType == ELECTRON_TOPLEP )        return getPermIndex(idx, TopLep_ElectronIdx);
   else if( objType == MUON_TOPLEP )       return getPermIndex(idx, TopLep_MuonIdx);
   else if( objType == BJETLEP_TOPLEP )    return getPermIndex(idx, TopLep_BJetLepIdx);
   else if( objType == PHOTON_TOPLEP )     return getPermIndex(idx, TopLep_PhotonIdx);
}

/// Get permutation index sorted in minimized NLL value
int KINFIT::kfit::getPermIndex(int idx, int *arr)
{
   if( NPerm_ > 0 )
     {	     
	if( idx < 0 ) return arr[0];
	else if( idx < NPerm_ ) return arr[idxMin_[idx]];
	else
	  {
	     std::cout << "Max number of permutations is " << NPerm_ << std::endl;
	     exit(1);
	  }	     
     }
   else return -1;
}

/// Sort all permutations according to the minimized NLL value
void KINFIT::kfit::sortPermIndex()
{
   float arrSort[NPerm_];
   
   for(int i=0;i<NPerm_;i++)
     {
	idxMin_[i] = i;
	arrSort[i] = chi_[i];
     }
   
   for(int i=0;i<NPerm_;i++)
     {
	for(int j=i;j<NPerm_;j++)
	  {
	     if( arrSort[j] < arrSort[i] )
	       {
		  float temp = arrSort[i];
		  arrSort[i] = arrSort[j];
		  arrSort[j] = temp;
		  
		  int tempIdx = idxMin_[i];
		  idxMin_[i] = idxMin_[j];
		  idxMin_[j] = tempIdx;
	       }
	  }	
     }   
}

/// Sort fit parametes according to the NLL value in the first-layer minimization
void KINFIT::kfit::sortPermVector(std::vector<double> vRef, std::vector<double> &vSort)
{
   std::vector<std::pair<double, double> > vMap;
   
   int nRef = vRef.size();
   
   for(int i=0;i<nRef;i++)
     {
	vMap.push_back(std::make_pair(vRef[i], vSort[i]));
     }   
   
   std::sort(vMap.begin(), vMap.end(), 
	     boost::bind(&std::pair<double, double>::first,_1) <
	     boost::bind(&std::pair<double, double>::first,_2));
   
   vSort.clear();
   
   for(int i=0;i<nRef;i++)
     {
	vSort.push_back(vMap[i].second);
     }   
}

/// Calculate normalized probability on a given transfer function
/// Assume the tails of exp(-0.1*|v-xmin(xmax)|) in case the variable is outside the function range
float KINFIT::kfit::getProb(TF1 *hPDF, float var, float max, double xmin, double xmax)
{
   float prob = 0.;
   float lm = 0.1;

   if( var < xmin ) prob = hPDF->Eval(xmin) * exp(-lm*fabs(var-xmin));
   else if( var > xmax ) prob = hPDF->Eval(xmax) * exp(-lm*fabs(var-xmax));
   else prob = hPDF->Eval(var);

   return prob/max;
}

/// Perform a variation on a given transfer function using rejection method
float KINFIT::kfit::getProbGaus(TF1 *hPDF, float max, float mean, float sigma, TRandom3 *rnd, float nSigma)
{
   float x = 0.;
  
   while(1)
     {
	float r1 = rnd->Rndm();
	float r2 = rnd->Rndm();

	x = mean - nSigma*sigma + 2*nSigma*sigma*r1;
	
	if( hPDF->Eval(x) > max*r2 ) break;
     }

   return x;
}

/// Read b jets
void KINFIT::kfit::SetBJet(std::vector<float> pt,
			   std::vector<float> eta,
			   std::vector<float> phi,
			   std::vector<float> E)
{
   static const bool samesize = (pt.size() == eta.size()) && (pt.size() == phi.size()) && (pt.size() == E.size());
   
   if( !samesize )
     {
	std::cout << "Sizes of input collections for b jets do not match" << std::endl;
	exit(1);
     }
   
   nBJet = pt.size();
   BJetPt = pt;
   BJetEta = eta;
   BJetPhi = phi;
   BJetE = E;
   
   BJetPx.clear();
   BJetPy.clear();
   BJetPz.clear();
   
   for(int i=0;i<nBJet;i++)
     {
	float px = BJetPt[i]*cos(BJetPhi[i]);
	float py = BJetPt[i]*sin(BJetPhi[i]);
	float pz = BJetPt[i]*sinh(BJetEta[i]);
	
	BJetPx.push_back(px);
	BJetPy.push_back(py);
	BJetPz.push_back(pz);
     }      
}

/// Read non-b jets
void KINFIT::kfit::SetNonBJet(std::vector<float> pt,
			      std::vector<float> eta,
			      std::vector<float> phi,
			      std::vector<float> E)
{
   static const bool samesize = (pt.size() == eta.size()) && (pt.size() == phi.size()) && (pt.size() == E.size());
   
   if( !samesize )
     {
	std::cout << "Sizes of input collections for non-b jets do not match" << std::endl;
	exit(1);
     }

   nNonBJet = pt.size();
   NonBJetPt = pt;
   NonBJetEta = eta;
   NonBJetPhi = phi;
   NonBJetE = E;

   NonBJetPx.clear();
   NonBJetPy.clear();
   NonBJetPz.clear();
   
   for(int i=0;i<nNonBJet;i++)
     {
	float px = NonBJetPt[i]*cos(NonBJetPhi[i]);
	float py = NonBJetPt[i]*sin(NonBJetPhi[i]);
	float pz = NonBJetPt[i]*sinh(NonBJetEta[i]);
	
	NonBJetPx.push_back(px);
	NonBJetPy.push_back(py);
	NonBJetPz.push_back(pz);
     }      
}

/// Read electrons
void KINFIT::kfit::SetElectron(std::vector<float> pt,
			       std::vector<float> eta,
			       std::vector<float> phi,
			       std::vector<float> E,
			       std::vector<int> charge)
{
   static const bool samesize = (pt.size() == eta.size()) && (pt.size() == phi.size()) && (pt.size() == E.size());
   
   if( !samesize )
     {
	std::cout << "Sizes of input collections for electrons do not match" << std::endl;
	exit(1);
     }

   nElectron = pt.size();
   ElectronPt = pt;
   ElectronEta = eta;
   ElectronPhi = phi;
   ElectronE = E;
   ElectronCharge = charge;

   ElectronPx.clear();
   ElectronPy.clear();
   ElectronPz.clear();
   
   for(int i=0;i<nElectron;i++)
     {
	float px = ElectronPt[i]*cos(ElectronPhi[i]);
	float py = ElectronPt[i]*sin(ElectronPhi[i]);
	float pz = ElectronPt[i]*sinh(ElectronEta[i]);
	
	ElectronPx.push_back(px);
	ElectronPy.push_back(py);
	ElectronPz.push_back(pz);
     }   
}

/// Read muons
void KINFIT::kfit::SetMuon(std::vector<float> pt,
			   std::vector<float> eta,
			   std::vector<float> phi,
			   std::vector<float> E,
			   std::vector<int> charge)
{
   static const bool samesize = (pt.size() == eta.size()) && (pt.size() == phi.size()) && (pt.size() == E.size());
   
   if( !samesize )
     {
	std::cout << "Sizes of input collections for muons do not match" << std::endl;
	exit(1);
     }

   nMuon = pt.size();
   MuonPt = pt;
   MuonEta = eta;
   MuonPhi = phi;
   MuonE = E;
   MuonCharge = charge;

   MuonPx.clear();
   MuonPy.clear();
   MuonPz.clear();
   
   for(int i=0;i<nMuon;i++)
     {
	float px = MuonPt[i]*cos(MuonPhi[i]);
	float py = MuonPt[i]*sin(MuonPhi[i]);
	float pz = MuonPt[i]*sinh(MuonEta[i]);
	
	MuonPx.push_back(px);
	MuonPy.push_back(py);
	MuonPz.push_back(pz);
     }   
}

/// Read photons
void KINFIT::kfit::SetPhoton(std::vector<float> pt,
			     std::vector<float> eta,
			     std::vector<float> phi,
			     std::vector<float> E)
{
   static const bool samesize = (pt.size() == eta.size()) && (pt.size() == phi.size()) && (pt.size() == E.size());
   
   if( !samesize )
     {
	std::cout << "Sizes of input collections for photons do not match" << std::endl;
	exit(1);
     }

   nPhoton = pt.size();
   PhotonPt = pt;
   PhotonEta = eta;
   PhotonPhi = phi;
   PhotonE = E;

   PhotonPx.clear();
   PhotonPy.clear();
   PhotonPz.clear();
   
   for(int i=0;i<nPhoton;i++)
     {
	float px = PhotonPt[i]*cos(PhotonPhi[i]);
	float py = PhotonPt[i]*sin(PhotonPhi[i]);
	float pz = PhotonPt[i]*sinh(PhotonEta[i]);
	
	PhotonPx.push_back(px);
	PhotonPy.push_back(py);
	PhotonPz.push_back(pz);
     }
   
   if( nPhoton > 0 ) 
     {
	IncludePhotons_ = true;
	
	bool foundPDFPhotonPx = checkPDF(hPDFPhotonPx.get(), PDF_NAME[PDF_PhotonPx]);
	bool foundPDFPhotonPy = checkPDF(hPDFPhotonPy.get(), PDF_NAME[PDF_PhotonPy]);
	bool foundPDFPhotonPz = checkPDF(hPDFPhotonPz.get(), PDF_NAME[PDF_PhotonPz]);
	bool foundPDFPhotonPt = checkPDF(hPDFPhotonPt.get(), PDF_NAME[PDF_PhotonPt]);
	bool foundPDFPhotonEta = checkPDF(hPDFPhotonEta.get(), PDF_NAME[PDF_PhotonEta]);
	bool foundPDFPhotonPhi = checkPDF(hPDFPhotonPhi.get(), PDF_NAME[PDF_PhotonPhi]);

	bool foundPDFPhotonPxPyPz = (foundPDFPhotonPx || foundPDFPhotonPy || foundPDFPhotonPz);
	bool foundPDFPhotonPtEtaPhi = (foundPDFPhotonPt || foundPDFPhotonEta || foundPDFPhotonPhi);
	
	if( foundPDFPhotonPxPyPz && foundPDFPhotonPtEtaPhi ) printPDFWarning();
	else if( foundPDFPhotonPxPyPz ) {
	   (*usePDFPhotonPxPyPz) = true;
	   (*usePDFPhotonPtEtaPhi) = false;
	   foundPDF.get()->insert({PDF_PhotonPx, foundPDFPhotonPx});
	   foundPDF.get()->insert({PDF_PhotonPy, foundPDFPhotonPy});
	   foundPDF.get()->insert({PDF_PhotonPz, foundPDFPhotonPz});
	}
	else {
	   (*usePDFPhotonPxPyPz) = false;
	   (*usePDFPhotonPtEtaPhi) = true;
	   foundPDF.get()->insert({PDF_PhotonPt, foundPDFPhotonPt});
	   foundPDF.get()->insert({PDF_PhotonEta, foundPDFPhotonEta});
	   foundPDF.get()->insert({PDF_PhotonPhi, foundPDFPhotonPhi});
	}	
     }
   else IncludePhotons_ = false;
}

/// Read missing transverse energy
void KINFIT::kfit::SetMet(float px,
			  float py)
{
   MetPx = px;
   MetPy = py;
}

/// Manually set missing transverse energy
void KINFIT::kfit::SetEtxEty(float etx,
			     float ety)
{
   Etx = etx;
   Ety = ety;
}

/// Manually set W boson mass
void KINFIT::kfit::SetWMass(float wm1,
			    float wm2)
{
   WMass1 = wm1;
   WMass2 = wm2;
}

/// Calcualte delta phi between two objects
float KINFIT::kfit::getDeltaPhi(float phi1, float phi2)
{
   float deltaPhi = fabs(phi1 - phi2);
   if( deltaPhi > M_PI ) deltaPhi = 2*M_PI - deltaPhi;
   return deltaPhi;
}

/// Calculate delta R between two objects
float KINFIT::kfit::getDeltaR(float eta1, float phi1, float eta2, float phi2)
{
   float deltaPhi = fabs(phi1 - phi2);
   if( deltaPhi > M_PI ) deltaPhi = 2*M_PI - deltaPhi;
   return sqrt( (eta1-eta2)*(eta1-eta2) + deltaPhi*deltaPhi );
}

/// Calculate eta
float KINFIT::kfit::getEta(float pt, float pz)
{   
   return asinh(pz/pt);
}

/// Check if the transfer function is defined
bool KINFIT::kfit::checkPDF(TF1 *tf, std::string tfname, bool isDeltaFunc)
{   
   if( strcmp(tf->GetName(), "") == 0 && ! isDeltaFunc )
     return false;
   else return true;
}

/// Print error message for transfer function check
void KINFIT::kfit::printPDFWarning()
{   
   std::cout << "Please specify transfer functions either for (px, py, pz) or for (pt, eta, phi)" << std::endl;
   exit(1);
}
