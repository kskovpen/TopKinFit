/// KINFIT class definition

#ifndef KINFIT_H
#define KINFIT_H

/// ROOT includes
#include "TRandom3.h"
#include "TMinuit.h"
#include "TFile.h"
#include "TF1.h"

/// Other includes
#include <vector>
#include <time.h>

/// Global variables
const int NTERMMAX = 30; ///< Maximum number of NLL terms
const int NPARMAX = 30;   ///< Maximum number of fit parameters
extern const int NNUMAX;
extern int FPARAM_N;

/// Result of minimization
struct FRESULT
{
   double par[NPARMAX];
   double lh;   
   double chiTerm[NTERMMAX];
   std::string chiTermName[NTERMMAX];
   int NTERM;
   float PzNu1;
   float PzNu2;
   float PzNuSum;
   float PtNuSum;
   float PNuSum;
   float TopTopMass;
   float TopTopPz;
   float TopTopPzAbs;
};

/// Model parameters
enum 
{
   FPARAM_Etx_TOPTOPLEPLEP,        ///< Absolute difference in neutrino px, (nu1_px-nu2_px)/2
   FPARAM_Ety_TOPTOPLEPLEP,        ///< Absolute difference in neutrino py, (nu1_py-nu2_py)/2
   FPARAM_Sign1_TOPTOPLEPLEP,      ///< Sign of the sqrt term for the first neutrino
   FPARAM_Sign2_TOPTOPLEPLEP,      ///< Sign of the sqrt term for the second neutrino
   FPARAM_EtRealX_TOPTOPLEPLEP,    ///< Total missing transverse energy, px
   FPARAM_EtRealY_TOPTOPLEPLEP,    ///< Total missing transverse energy, py
   FPARAM_mW1_TOPTOPLEPLEP,        ///< Mass of the first W boson
   FPARAM_mW2_TOPTOPLEPLEP,        ///< Mass of the second W boson
   FPARAM_BJet1Px_TOPTOPLEPLEP,    ///< First b jet momentum, px
   FPARAM_BJet1Py_TOPTOPLEPLEP,    ///< First b jet momentum, py
   FPARAM_BJet1Pz_TOPTOPLEPLEP,    ///< First b jet momentum, pz
   FPARAM_BJet1E_TOPTOPLEPLEP,     ///< First b jet energy
   FPARAM_BJet2Px_TOPTOPLEPLEP,    ///< Second b jet momentum, px
   FPARAM_BJet2Py_TOPTOPLEPLEP,    ///< Second b jet momentum, py
   FPARAM_BJet2Pz_TOPTOPLEPLEP,    ///< Second b jet momentum, pz
   FPARAM_BJet2E_TOPTOPLEPLEP,     ///< Second b jet energy
   FPARAM_Lepton1Px_TOPTOPLEPLEP,  ///< First lepton momentum, px
   FPARAM_Lepton1Py_TOPTOPLEPLEP,  ///< First lepton momentum, py
   FPARAM_Lepton1Pz_TOPTOPLEPLEP,  ///< First lepton momentum, pz
   FPARAM_Lepton1E_TOPTOPLEPLEP,   ///< First lepton energy
   FPARAM_Lepton2Px_TOPTOPLEPLEP,  ///< Second lepton momentum, px
   FPARAM_Lepton2Py_TOPTOPLEPLEP,  ///< Second lepton momentum, py
   FPARAM_Lepton2Pz_TOPTOPLEPLEP,  ///< Second lepton momentum, pz
   FPARAM_Lepton2E_TOPTOPLEPLEP,   ///< Second lepton energy
   FPARAM_PhotonPx_TOPTOPLEPLEP,   ///< Photon momentum, px
   FPARAM_PhotonPy_TOPTOPLEPLEP,   ///< Photon momentum, py
   FPARAM_PhotonPz_TOPTOPLEPLEP,   ///< Photon momentum, pz
   FPARAM_PhotonE_TOPTOPLEPLEP,    ///< Photon energy
   FPARAM_N_TOPTOPLEPLEP           ///< Total number of parameters
};

enum 
{   
   FPARAM_Sign_TOPTOPLEPHAD,         ///< Sign of the sqrt term for the neutrino
   FPARAM_EtRealX_TOPTOPLEPHAD,      ///< Total missing transverse energy, px
   FPARAM_EtRealY_TOPTOPLEPHAD,      ///< Total missing transverse energy, py
   FPARAM_mWLep_TOPTOPLEPHAD,        ///< Mass of the leptonic W boson
   FPARAM_mWHad_TOPTOPLEPHAD,        ///< Mass of the hadronic W boson
   FPARAM_BJetLepPx_TOPTOPLEPHAD,    ///< Leptonic top b jet momentum, px
   FPARAM_BJetLepPy_TOPTOPLEPHAD,    ///< Leptonic top b jet momentum, py
   FPARAM_BJetLepPz_TOPTOPLEPHAD,    ///< Leptonic top b jet momentum, pz
   FPARAM_BJetLepE_TOPTOPLEPHAD,     ///< Leptonic top b jet energy
   FPARAM_BJetHadPx_TOPTOPLEPHAD,    ///< Hadronic top b jet momentum, px
   FPARAM_BJetHadPy_TOPTOPLEPHAD,    ///< Hadronic top b jet momentum, py
   FPARAM_BJetHadPz_TOPTOPLEPHAD,    ///< Hadronic top b jet momentum, pz
   FPARAM_BJetHadE_TOPTOPLEPHAD,     ///< Hadronic top b jet energy
   FPARAM_NonBJet1Px_TOPTOPLEPHAD,   ///< Hadronic top first non b jet momentum, px
   FPARAM_NonBJet1Py_TOPTOPLEPHAD,   ///< Hadronic top first non b jet momentum, py
   FPARAM_NonBJet1Pz_TOPTOPLEPHAD,   ///< Hadronic top first non b jet momentum, pz
   FPARAM_NonBJet1E_TOPTOPLEPHAD,    ///< Hadronic top first non b jet energy
   FPARAM_NonBJet2Px_TOPTOPLEPHAD,   ///< Hadronic top second non b jet momentum, px
   FPARAM_NonBJet2Py_TOPTOPLEPHAD,   ///< Hadronic top second non b jet momentum, py
   FPARAM_NonBJet2Pz_TOPTOPLEPHAD,   ///< Hadronic top second non b jet momentum, pz
   FPARAM_NonBJet2E_TOPTOPLEPHAD,    ///< Hadronic top second non b jet energy
   FPARAM_LeptonPx_TOPTOPLEPHAD,     ///< Lepton momentum, px
   FPARAM_LeptonPy_TOPTOPLEPHAD,     ///< Lepton momentum, py
   FPARAM_LeptonPz_TOPTOPLEPHAD,     ///< Lepton momentum, pz
   FPARAM_LeptonE_TOPTOPLEPHAD,      ///< Lepton energy
   FPARAM_PhotonPx_TOPTOPLEPHAD,     ///< Photon momentum, px
   FPARAM_PhotonPy_TOPTOPLEPHAD,     ///< Photon momentum, py
   FPARAM_PhotonPz_TOPTOPLEPHAD,     ///< Photon momentum, pz
   FPARAM_PhotonE_TOPTOPLEPHAD,      ///< Photon energy
   FPARAM_N_TOPTOPLEPHAD             ///< Total number of parameters
};

enum 
{   
   FPARAM_Sign_TOPLEP,         ///< Sign of the sqrt term for the neutrino
   FPARAM_EtRealX_TOPLEP,      ///< Total missing transverse energy, px
   FPARAM_EtRealY_TOPLEP,      ///< Total missing transverse energy, py
   FPARAM_mWLep_TOPLEP,        ///< Mass of the leptonic W boson
   FPARAM_BJetLepPx_TOPLEP,    ///< Leptonic top b jet momentum, px
   FPARAM_BJetLepPy_TOPLEP,    ///< Leptonic top b jet momentum, py
   FPARAM_BJetLepPz_TOPLEP,    ///< Leptonic top b jet momentum, pz
   FPARAM_BJetLepE_TOPLEP,     ///< Leptonic top b jet energy
   FPARAM_LeptonPx_TOPLEP,     ///< Lepton momentum, px
   FPARAM_LeptonPy_TOPLEP,     ///< Lepton momentum, py
   FPARAM_LeptonPz_TOPLEP,     ///< Lepton momentum, pz
   FPARAM_LeptonE_TOPLEP,      ///< Lepton energy
   FPARAM_PhotonPx_TOPLEP,     ///< Photon momentum, px
   FPARAM_PhotonPy_TOPLEP,     ///< Photon momentum, py
   FPARAM_PhotonPz_TOPLEP,     ///< Photon momentum, pz
   FPARAM_PhotonE_TOPLEP,      ///< Photon energy
   FPARAM_N_TOPLEP             ///< Total number of parameters
};

/// NLL terms
enum 
{
   NLL_W1_TOPTOPLEPLEP,
   NLL_W2_TOPTOPLEPLEP,
   NLL_Top1_TOPTOPLEPLEP,
   NLL_Top2_TOPTOPLEPLEP,
   NLL_EtMissX_TOPTOPLEPLEP,
   NLL_EtMissY_TOPTOPLEPLEP,
   NLL_BJet1Px_TOPTOPLEPLEP,
   NLL_BJet1Py_TOPTOPLEPLEP,
   NLL_BJet1Pz_TOPTOPLEPLEP,
   NLL_BJet2Px_TOPTOPLEPLEP,
   NLL_BJet2Py_TOPTOPLEPLEP,
   NLL_BJet2Pz_TOPTOPLEPLEP,
   NLL_Lepton1Px_TOPTOPLEPLEP,
   NLL_Lepton1Py_TOPTOPLEPLEP,
   NLL_Lepton1Pz_TOPTOPLEPLEP,
   NLL_Lepton2Px_TOPTOPLEPLEP,
   NLL_Lepton2Py_TOPTOPLEPLEP,
   NLL_Lepton2Pz_TOPTOPLEPLEP,
   NLL_PhotonPx_TOPTOPLEPLEP,
   NLL_PhotonPy_TOPTOPLEPLEP,
   NLL_PhotonPz_TOPTOPLEPLEP,
   NLL_RacPT_TOPTOPLEPLEP,
   NLL_mWPT_TOPTOPLEPLEP,
   NLL_mTopPT_TOPTOPLEPLEP,
   NLL_TopTopMass_TOPTOPLEPLEP,
   NLL_N_TOPTOPLEPLEP
};

enum 
{
   NLL_WLep_TOPTOPLEPHAD,
   NLL_WHad_TOPTOPLEPHAD,
   NLL_TopLep_TOPTOPLEPHAD,
   NLL_TopHad_TOPTOPLEPHAD,
   NLL_EtMissX_TOPTOPLEPHAD,
   NLL_EtMissY_TOPTOPLEPHAD,
   NLL_BJetLepPx_TOPTOPLEPHAD,
   NLL_BJetLepPy_TOPTOPLEPHAD,
   NLL_BJetLepPz_TOPTOPLEPHAD,
   NLL_BJetHadPx_TOPTOPLEPHAD,
   NLL_BJetHadPy_TOPTOPLEPHAD,
   NLL_BJetHadPz_TOPTOPLEPHAD,
   NLL_NonBJet1Px_TOPTOPLEPHAD,
   NLL_NonBJet1Py_TOPTOPLEPHAD,
   NLL_NonBJet1Pz_TOPTOPLEPHAD,
   NLL_NonBJet2Px_TOPTOPLEPHAD,
   NLL_NonBJet2Py_TOPTOPLEPHAD,
   NLL_NonBJet2Pz_TOPTOPLEPHAD,
   NLL_LeptonPx_TOPTOPLEPHAD,
   NLL_LeptonPy_TOPTOPLEPHAD,
   NLL_LeptonPz_TOPTOPLEPHAD,
   NLL_PhotonPx_TOPTOPLEPHAD,
   NLL_PhotonPy_TOPTOPLEPHAD,
   NLL_PhotonPz_TOPTOPLEPHAD,
   NLL_RacPT_TOPTOPLEPHAD,
   NLL_mWPT_TOPTOPLEPHAD,
   NLL_mTopPT_TOPTOPLEPHAD,
   NLL_TopTopMass_TOPTOPLEPHAD,
   NLL_N_TOPTOPLEPHAD
};

enum 
{
   NLL_WLep_TOPLEP,
   NLL_TopLep_TOPLEP,
   NLL_EtMissX_TOPLEP,
   NLL_EtMissY_TOPLEP,
   NLL_BJetLepPx_TOPLEP,
   NLL_BJetLepPy_TOPLEP,
   NLL_BJetLepPz_TOPLEP,
   NLL_LeptonPx_TOPLEP,
   NLL_LeptonPy_TOPLEP,
   NLL_LeptonPz_TOPLEP,
   NLL_PhotonPx_TOPLEP,
   NLL_PhotonPy_TOPLEP,
   NLL_PhotonPz_TOPLEP,
   NLL_RacPT_TOPLEP,
   NLL_mWPT_TOPLEP,
   NLL_mTopPT_TOPLEP,
   NLL_N_TOPLEP
};

/// Process definitions
typedef enum 
{
   TOPTOPLEPLEP,
   TOPTOPLEPHAD,
   TOPLEP,
   HYPO_N
} HYPO;

/// Object mapping for the indices in the input collections
typedef enum 
{
   ELECTRON1_TOPTOPLEPLEP,    ///< Electron from the first top quark's W boson decay
   MUON1_TOPTOPLEPLEP,        ///< Muon from the first top quark's W boson decay
   ELECTRON2_TOPTOPLEPLEP,    ///< Electron from the second top quark's W boson decay
   MUON2_TOPTOPLEPLEP,        ///< Muon from the second top quark's W boson decay
   BJET1_TOPTOPLEPLEP,        ///< b jet from the first top quark decay
   BJET2_TOPTOPLEPLEP,        ///< b jet from the second top quark decay
   PHOTON_TOPTOPLEPLEP,       ///< photon radiated in the top quark decay

   ELECTRON_TOPTOPLEPHAD,     ///< Electron from the leptonic top quark's W boson decay
   MUON_TOPTOPLEPHAD,         ///< Muon from the leptonic top quark's W boson decay
   BJETLEP_TOPTOPLEPHAD,      ///< b jet from the leptonic top quark decay
   BJETHAD_TOPTOPLEPHAD,      ///< b jet from the hadronic top quark decay
   NONBJET1_TOPTOPLEPHAD,     ///< first non b jet from the hadronic top quark decay
   NONBJET2_TOPTOPLEPHAD,     ///< second non b jet from the hadronic top quark decay
   PHOTON_TOPTOPLEPHAD,       ///< photon radiated in the top quark decay

   ELECTRON_TOPLEP,     ///< Electron from the leptonic top quark's W boson decay
   MUON_TOPLEP,         ///< Muon from the leptonic top quark's W boson decay
   BJETLEP_TOPLEP,      ///< b jet from the leptonic top quark decay
   PHOTON_TOPLEP        ///< photon radiated in the top quark decay
} OBJ;

/// Transfer functions
typedef enum
{
   PDF_TopWMass,     ///< Reconstructed W boson mass from leptonic decay
   PDF_TopMass,      ///< Reconstructed top quark mass from t -> b W -> b l nu
   PDF_TopWHadMass,  ///< Reconstructed W boson mass from hadronic decay
   PDF_TopHadMass,   ///< Reconstructed top quark mass from t -> b W -> b q q'bar
//   PDF_TopTopMass,   ///< Reconstructed top quark mass from t -> b W -> b l nu
   PDF_MetPx,        ///< (Px(gen)-Px(rec))/Px(gen) for total missing transverse energy
   PDF_MetPy,        ///< (Py(gen)-Py(rec))/Py(gen) for total missing transverse energy
   PDF_BJetPx,       ///< (Px(gen)-Px(rec))/Px(gen) for b jets
   PDF_BJetPy,       ///< (Py(gen)-Py(rec))/Py(gen) for b jets
   PDF_BJetPz,       ///< (Pz(gen)-Pz(rec))/Pz(gen) for b jets
   PDF_NonBJetPx,    ///< (Px(gen)-Px(rec))/Px(gen) for non b jets
   PDF_NonBJetPy,    ///< (Py(gen)-Py(rec))/Py(gen) for non b jets
   PDF_NonBJetPz,    ///< (Pz(gen)-Pz(rec))/Pz(gen) for non b jets
   PDF_ElecPx,       ///< (Px(gen)-Px(rec))/Px(gen) for electrons
   PDF_ElecPy,       ///< (Py(gen)-Py(rec))/Py(gen) for electrons
   PDF_ElecPz,       ///< (Pz(gen)-Pz(rec))/Pz(gen) for electrons
   PDF_MuonPx,       ///< (Px(gen)-Px(rec))/Px(gen) for muons
   PDF_MuonPy,       ///< (Py(gen)-Py(rec))/Py(gen) for muons
   PDF_MuonPz,       ///< (Pz(gen)-Pz(rec))/Pz(gen) for muons
   PDF_PhotonPx,     ///< (Px(gen)-Px(rec))/Px(gen) for photons
   PDF_PhotonPy,     ///< (Py(gen)-Py(rec))/Py(gen) for photons
   PDF_PhotonPz,     ///< (Pz(gen)-Pz(rec))/Pz(gen) for photons
   PDF_N             ///< Total number of transfer function definitions
} PDF;

/// Origin of radiated photons
typedef enum
{
   PHOTON_FROM_TOP1_TOPTOPLEPLEP,
   PHOTON_FROM_LEPTON1_TOPTOPLEPLEP,
   PHOTON_FROM_W1_TOPTOPLEPLEP,
   PHOTON_FROM_TOP2_TOPTOPLEPLEP,
   PHOTON_FROM_LEPTON2_TOPTOPLEPLEP,
   PHOTON_FROM_W2_TOPTOPLEPLEP,
   PHOTON_FROM_BJET1_TOPTOPLEPLEP,
   PHOTON_FROM_BJET2_TOPTOPLEPLEP,
   PHOTON_FROM_ISR_TOPTOPLEPLEP,
   PHOTON_FROM_TOP1_COMB_TOPTOPLEPLEP, // top or b quark
   PHOTON_FROM_W1_COMB_TOPTOPLEPLEP, // w, lepton, or light quarks
   PHOTON_FROM_TOP2_COMB_TOPTOPLEPLEP, // top or b quark
   PHOTON_FROM_W2_COMB_TOPTOPLEPLEP, // w, lepton, or light quarks
   PHOTON_ORIGIN_N_TOPTOPLEPLEP
} PHOTON_ORIGIN_TOPTOPLEPLEP;

typedef enum
{
   PHOTON_FROM_TOPLEP_TOPTOPLEPHAD,
   PHOTON_FROM_LEPTON_TOPTOPLEPHAD,
   PHOTON_FROM_WLEP_TOPTOPLEPHAD,
   PHOTON_FROM_TOPHAD_TOPTOPLEPHAD,
   PHOTON_FROM_WHAD_TOPTOPLEPHAD,
   PHOTON_FROM_BJETLEP_TOPTOPLEPHAD,
   PHOTON_FROM_BJETHAD_TOPTOPLEPHAD,
   PHOTON_FROM_NONBJET1_TOPTOPLEPHAD,
   PHOTON_FROM_NONBJET2_TOPTOPLEPHAD,
   PHOTON_FROM_ISR_TOPTOPLEPHAD,
   PHOTON_FROM_TOPLEP_COMB_TOPTOPLEPHAD, // top or b quark
   PHOTON_FROM_WLEP_COMB_TOPTOPLEPHAD, // w, lepton, or light quarks
   PHOTON_FROM_TOPHAD_COMB_TOPTOPLEPHAD, // top or b quark
   PHOTON_FROM_WHAD_COMB_TOPTOPLEPHAD, // w, lepton, or light quarks
   PHOTON_ORIGIN_N_TOPTOPLEPHAD
} PHOTON_ORIGIN_TOPTOPLEPHAD;

typedef enum
{
   PHOTON_FROM_TOPLEP_TOPLEP,
   PHOTON_FROM_LEPTON_TOPLEP,
   PHOTON_FROM_WLEP_TOPLEP,
   PHOTON_FROM_BJETLEP_TOPLEP,
   PHOTON_FROM_ISR_TOPLEP,
   PHOTON_FROM_TOPLEP_COMB_TOPLEP, // top or b quark
   PHOTON_FROM_WLEP_COMB_TOPLEP, // w or lepton
   PHOTON_ORIGIN_N_TOPLEP
} PHOTON_ORIGIN_TOPLEP;

namespace KINFIT
{
   /// Main class
   class kfit
     {
	
      public:
	
	/// Default constructor
	kfit();
	
	/// Default destructor
	virtual ~kfit();

	/******************/
	/*  User methods  */
	/******************/
	
	/// Initialization of the method for specific process
	void Init(HYPO hypoMode=TOPTOPLEPLEP);

	/// Run the method
	void Run();

	/// Read input collection with reconstructed b jets
	void SetBJet(std::vector<float> pt,
		     std::vector<float> eta,
		     std::vector<float> phi,
		     std::vector<float> E);

	/// Read input collection with reconstructed non b jets
	void SetNonBJet(std::vector<float> pt,
			std::vector<float> eta,
			std::vector<float> phi,
			std::vector<float> E);

	/// Read input collection with reconstructed electrons
	void SetElectron(std::vector<float> pt,
			 std::vector<float> eta,
			 std::vector<float> phi,
			 std::vector<float> E,
			 std::vector<int> charge);

	/// Read input collection with reconstructed muons
	void SetMuon(std::vector<float> pt,
		     std::vector<float> eta,
		     std::vector<float> phi,
		     std::vector<float> E,
		     std::vector<int> charge);

	/// Read input collection with reconstructed photons
	void SetPhoton(std::vector<float> pt,
		       std::vector<float> eta,
		       std::vector<float> phi,
		       std::vector<float> E);
	
	/// Read reconstructed missing energy
	void SetMet(float px, float py);

	/// Manually set the absolute difference in reconstructed neutrino momenta
	void SetEtxEty(float etx, float ety);

	/// Manually set the W boson masses
	void SetWMass(float wm1, float wm2);

	/// Manually set the top quark masses
	void SetTopMass(float topm1, float topm2);
	
	/// Read transfer function from the input file
	void SetPDF(std::string obj, std::string fileName, std::string hName);

	/// Set the number of grid points for the first-layer minimization
	void SetNToy(int nMC) { NToy_ = nMC; };
	
	/// Set the maximal number of fits to perform in the second-layer minimization
	void SetNFitMax(int nFit) { NFitMax_ = nFit; };
	
	/// Set the maximal number of grid points to consider as input for the second-layer minimization
	void SetNGrid(int nGrid) { NGrid_ = nGrid; };
	
	/// Set the upper cut-off value on the NLL results obtained on grid points
	void SetLHMaxGeneric(float cut) { LHMaxGeneric_ = cut; };
	
	/// Set the upper threshold on the NLL minimal value in the second-layer minimization
	void SetLHMaxMinuit(float cut) { LHMaxMinuit_ = cut; };

	/// Set the maximal number of RMS in the variations on the transfer functions for missing transverse energy
	void SetNMetRMS(int nRMS) { NMetRMS_ = nRMS; };
	
	/// Set the maximal number of RMS in the variations on the transfer functions for b jet px
	void SetNBJetPxRMS(int nRMS) { NBJetPxRMS_ = nRMS; };
	
	/// Set the maximal number of RMS in the variations on the transfer functions for b jet py
	void SetNBJetPyRMS(int nRMS) { NBJetPyRMS_ = nRMS; };
	
	/// Set the maximal number of RMS in the variations on the transfer functions for b jet pz
	void SetNBJetPzRMS(int nRMS) { NBJetPzRMS_ = nRMS; };

	/// Set the maximal number of RMS in the variations on the transfer functions for non-b jet px
	void SetNNonBJetPxRMS(int nRMS) { NNonBJetPxRMS_ = nRMS; };
	
	/// Set the maximal number of RMS in the variations on the transfer functions for non-b jet py
	void SetNNonBJetPyRMS(int nRMS) { NNonBJetPyRMS_ = nRMS; };
	
	/// Set the maximal number of RMS in the variations on the transfer functions for non-b jet pz
	void SetNNonBJetPzRMS(int nRMS) { NNonBJetPzRMS_ = nRMS; };
	
	/// Set the maximal number of RMS in the variations on the transfer functions for electron px
	void SetNElecPxRMS(int nRMS) { NElecPxRMS_ = nRMS; };
	
	/// Set the maximal number of RMS in the variations on the transfer functions for electron py
	void SetNElecPyRMS(int nRMS) { NElecPyRMS_ = nRMS; };
	
	/// Set the maximal number of RMS in the variations on the transfer functions for electron pz
	void SetNElecPzRMS(int nRMS) { NElecPzRMS_ = nRMS; };

	/// Set the maximal number of RMS in the variations on the transfer functions for muon px
	void SetNMuonPxRMS(int nRMS) { NMuonPxRMS_ = nRMS; };
	
	/// Set the maximal number of RMS in the variations on the transfer functions for muon py
	void SetNMuonPyRMS(int nRMS) { NMuonPyRMS_ = nRMS; };
	
	/// Set the maximal number of RMS in the variations on the transfer functions for muon pz
	void SetNMuonPzRMS(int nRMS) { NMuonPzRMS_ = nRMS; };

	/// Set the maximal number of RMS in the variations on the transfer functions for photon px
	void SetNPhotonPxRMS(int nRMS) { NPhotonPxRMS_ = nRMS; };
	
	/// Set the maximal number of RMS in the variations on the transfer functions for photon py
	void SetNPhotonPyRMS(int nRMS) { NPhotonPyRMS_ = nRMS; };
	
	/// Set the maximal number of RMS in the variations on the transfer functions for photon pz
	void SetNPhotonPzRMS(int nRMS) { NPhotonPzRMS_ = nRMS; };
	
	/// Set the fit parameter range in terms of transfer function's RMS for lepton and jet kinematics
	void SetLimNRMS(int nRMS) { LimNRMS_ = nRMS; };

	/// Get the photon origin for a given permutation
	float GetPhotonOrigin(int idxPerm) { return PhotonOrigin_[idxMin_[idxPerm]]; };
	
	/// Get the reconstructed dR(top, top) for a given permutation
	float GetDrTopTop(int idxPerm) { return drTopTop_[idxMin_[idxPerm]]; };
	
	/// Get the reconstructed dEta(top, top) for a given permutation
	float GetDetaTopTop(int idxPerm) { return detaTopTop_[idxMin_[idxPerm]]; };
	
	/// Get the reconstructed dPhi(top, top) for a given permutation
	float GetDphiTopTop(int idxPerm) { return dphiTopTop_[idxMin_[idxPerm]]; };
	
	/// Get the reconstructed mass of the di-top system for a given permutation
	float GetMTopTop(int idxPerm) { return mTopTop_[idxMin_[idxPerm]]; };
	
	/// Get the reconstructed pt of the di-top system for a given permutation
	float GetPtTopTop(int idxPerm) { return ptTopTop_[idxMin_[idxPerm]]; };
	
	/// Get the reconstructed momentum of the di-top system for a given permutation
	float GetPTopTop(int idxPerm) { return pTopTop_[idxMin_[idxPerm]]; };
	
	/// Get the reconstructed eta of the di-top system for a given permutation
	float GetEtaTopTop(int idxPerm) { return etaTopTop_[idxMin_[idxPerm]]; };
	
	/// Get the reconstructed phi of the di-top system for a given permutation
	float GetPhiTopTop(int idxPerm) { return phiTopTop_[idxMin_[idxPerm]]; };

	/// Set minimal value for the EtX parameter in the generic minimization
	void SetEtxMin(float min) { EtxMin_ = min; };

	/// Set maximal value for the EtX parameter in the generic minimization
	void SetEtxMax(float max) { EtxMax_ = max; };
	
	/// Get wall time spent on the first-layer minimization for a given permutation
	double GetWallTimeFirst(int idxPerm) { return timerWallGeneric_[idxMin_[idxPerm]]; };

	/// Get CPU time spent on the first-layer minimization for a given permutation
	double GetCPUTimeFirst(int idxPerm) { return timerCPUGeneric_[idxMin_[idxPerm]]; };
	
	/// Get wall time spent on the second-layer minimization for a given permutation
	double GetWallTimeSecond(int idxPerm) { return timerWallMinuit_[idxMin_[idxPerm]]; };

	/// Get CPU time spent on the second-layer minimization for a given permutation
	double GetCPUTimeSecond(int idxPerm) { return timerCPUMinuit_[idxMin_[idxPerm]]; };

	/// Include individual photon permutations
	void CheckAllPhotonOrigins() { CheckAllPhotonOrigins_ = true; };
	
	/// Get the minimal NLL value for a given permutation
	float GetDisc(int idx)
	  {
	     checkIndex(idx, "permutation", NPerm_);
	     return chi_[idxMin_[idx]];
	  };

	/// Get a specific term contribution to the minimal NLL value
	float GetDiscTerm(int idxPerm, int idxTerm)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
	     checkIndex(idxTerm, "term", NTerm_);
	     return chiTerm_[idxMin_[idxPerm]][idxTerm];
	  };

	/// Get the NLL term name by its index
	std::string GetDiscTermName(int idxTerm)
	  {
	     checkIndex(idxTerm, "term", NTerm_);
	     return chiTermName_[idxTerm];
	  };

	/// Get the number of grid points in generic minimization
	int GetNGeneric(int idx)
	  {
	     checkIndex(idx, "permutation", NPerm_);
	     return NGeneric_[idxMin_[idx]];
	  };
	
	/// Get the total number of available parameters in the model
	int GetNPar() { return FPARAM_N; }

	/// Get the model parameter name by its index
	std::string GetParName(int idxPar)
	  {
	     checkIndex(idxPar, "parameter", FPARAM_N);
	     return FPARAM_NAME[idxPar];
	  }
	
	/// Get the fitted parameter value
	float GetPar(int idxPerm, int idxPar)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
	     checkIndex(idxPar, "parameter", FPARAM_N);
	     return par_[idxMin_[idxPerm]][idxPar];
	  };

	/// Get the lower thershold of the allowed range of the fit parameter
	float GetParMin(int idxPar)
	  {
	     checkIndex(idxPar, "parameter", FPARAM_N);
	     return (*ParMin)[idxPar];
	  };

	/// Get the upper thershold of the allowed range of the fit parameter
	float GetParMax(int idxPar)
	  {
	     checkIndex(idxPar, "parameter", FPARAM_N);
	     return (*ParMax)[idxPar];
	  };

	/// Check if parameter is fixed in the fit
	bool IsFixed(int idxPar)
	  {
	     checkIndex(idxPar, "parameter", FPARAM_N);
	     return (*IsParFixed)[idxPar];
	  };

	/// Fix the given model parameter in the fit
	void FixParameter(int idxPar)
	  {
	     checkIndex(idxPar, "parameter", FPARAM_N);
	     (*IsParFixed)[idxPar] = true;
	  };

	/// Free the given model parameter in the fit
	void FreeParameter(int idxPar)
	  {
	     checkIndex(idxPar, "parameter", FPARAM_N);
	     (*IsParFixed)[idxPar] = false;
	  };
	
	/// Get the reconstructed neutrino px
	float GetNuPx(int idxPerm, int idxNu)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_); 
	     checkIndex(idxNu, "neutrino", NNu_);
	     return nuPx_[idxMin_[idxPerm]][idxNu];
	  };

	/// Get the reconstructed neutrino py
	float GetNuPy(int idxPerm, int idxNu)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
	     checkIndex(idxNu, "neutrino", NNu_);
	     return nuPy_[idxMin_[idxPerm]][idxNu];
	  };

	/// Get the reconstructed neutrino pz
	float GetNuPz(int idxPerm, int idxNu)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
	     checkIndex(idxNu, "neutrino", NNu_);
	     return nuPz_[idxMin_[idxPerm]][idxNu];
	  };
	
	/// Get the reconstructed missing transverse energy px after minimization
	float GetMetX(int idxPerm)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
	     return MetPx_[idxMin_[idxPerm]];
	  };
	
	/// Get the reconstructed missing transverse energy py after minimization
	float GetMetY(int idxPerm)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
	     return MetPy_[idxMin_[idxPerm]];
	  };
	
	/// Get the reconstructed W boson mass after minimization
	float GetWMass(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
	     return WMass_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed W boson momentum
	float GetWP(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return WP_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed W boson pt
	float GetWPt(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return WPt_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed W boson px
	float GetWPx(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return WPx_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed W boson py
	float GetWPy(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return WPy_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed W boson pz
	float GetWPz(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return WPz_[idxMin_[idxPerm]][index];
	  };
	
	/// Get the reconstructed W boson eta
	float GetWEta(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return WEta_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed W boson phi
	float GetWPhi(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return WPhi_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed W boson energy
	float GetWE(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return WE_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed top quark mass
	float GetTopMass(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return TopMass_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed top quark pt
	float GetTopPt(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return TopPt_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed top quark px
	float GetTopPx(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return TopPx_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed top quark py
	float GetTopPy(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return TopPy_[idxMin_[idxPerm]][index];
	  };
	
	/// Get the reconstructed top quark pz
	float GetTopPz(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return TopPz_[idxMin_[idxPerm]][index];
	  };
	
	/// Get the reconstructed top quark momentum
	float GetTopP(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return TopP_[idxMin_[idxPerm]][index];
	  };
	
	/// Get the reconstructed top quark eta
	float GetTopEta(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return TopEta_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed top quark phi
	float GetTopPhi(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return TopPhi_[idxMin_[idxPerm]][index];
	  };

	/// Get the reconstructed top quark energy
	float GetTopE(int idxPerm, int index)
	  {
	     checkIndex(idxPerm, "permutation", NPerm_);
//	     checkIndex(index, "particle", NNu_);
	     return TopE_[idxMin_[idxPerm]][index];
	  };
	
	/// Get reconstructed object index in the input collection according to kinematic fit
	int GetIndex(OBJ objType=BJET1_TOPTOPLEPLEP, int idx = -1);
	
	/// Get the total number of permutations per event
	int GetNPerm() { return NPerm_; };
	
	/// Get the total number of NLL terms
	int GetNTerm() { return NTerm_; };
	
	/// Include additional term in the NLL with the reconstructed TopTopMass
//	void AddTopTopMassToNLL() { AddTopTopMassToNLL_ = true; };

	/// Sort best fits according to the TopTopMass value (in ascending order)
	void SortByTopTopMass(bool flag) { SortByTopTopMass_ = flag; };

	/// Sort best fits according to the TopTopPz value (in ascending order)
	void SortByTopTopPz(bool flag) { SortByTopTopPz_ = flag; };

	/// Sort best fits according to the TopTopPzAbs value (in ascending order)
	void SortByTopTopPzAbs(bool flag) { SortByTopTopPzAbs_ = flag; };
	
	/// Sort best fits according to the |PzNu1|+|PzNu2| (in ascending order)
	void SortByNeutrinoPz(bool flag) { SortByNeutrinoPz_ = flag; };

	/// Sort best fits according to the |PtNu1|+|PtNu2| (in ascending order)
	void SortByNeutrinoPt(bool flag) { SortByNeutrinoPt_ = flag; };
	
	/// Sort best fits according to the |PNu1|+|PNu2| (in ascending order)
	void SortByNeutrinoP(bool flag) { SortByNeutrinoP_ = flag; };
	
	/// Perform the runtime measurements
	void DoTiming() { DoTiming_ = true; };

	/// Do MINUIT fits
	void DoFit(bool flag = true) { DoFit_ = flag; };
	
	/**********************/
	/*  Internal methods  */
	/**********************/
	
      private:

	/// Clean all variables on each event
	void Reset();

	/// Clean fit-related derivatives on each event
	void ResetFitInfo();
	
	/// Run kinematic fit for the TOPTOPLEPLEP process
	void TopTopLepLep();

	/// Run kinematic fit for the TOPTOPLEPHAD process
	void TopTopLepHad();

	/// Run kinematic fit for the TOPLEP process
	void TopLep();
	
	/// Put electrons and muons into one collection
	void FillLeptons();
	
      protected:
	
	/// Calculate normalized probability for specific variable on a given transfer function
	static float getProb(TF1 *hPDF, float var, float max, double xmin, double xmax);
	
	/// Perform a randomized variation on a given transfer function
	static float getProbGaus(TF1 *hPDF, float max, float mean, float sigma, TRandom3 *rnd, float nSigma);
	
	/// Associate object index in kinematic fit to the input collection
	int getPermIndex(int idx, int *arr);
	
	/// Sort permutations according to the minimized NLL value
	void sortPermIndex();
	
	/// Sort fits per permutation according to the minimized NLL value
	void sortPermVector(std::vector<double> vRef, std::vector<double> &vSort);

	/// Calculate delta R variable
	float getDeltaR(float eta1, float phi1, float eta2, float phi2);
	
	/// Calculate delta Phi variable
	float getDeltaPhi(float phi1, float phi2);
	
	/// Calculate eta variable
	float getEta(float pt, float pz);
	
	/// Check if the transfer function is defined
	void checkPDF(TF1* hf, std::string tfname);

      private:
	
	void checkIndex(int idx, std::string name, int max)
	  {
	     if( idx < 0 )
	       {
		  std::cout << "Requested " << name << " index is negative" << std::endl;
		  exit(1);
	       }	
	     else if( idx >= max )
	       {
		  std::cout << "Requested index exceeds the total number of permutations (" << max-1 << ")" << std::endl;
		  exit(1);
	       }	  
	  };
	
	/*********************/
	/*  Class variables  */
	/*********************/

      protected:

	bool ISINITIALIZED_; ///< Flag for the successful tool initialization
	
	TRandom3 *rnd; ///< Random number generator
	
	int hypoMode; ///< Process hypothesis index

	///
	/// Input collections
	///
	
	float MetPx;  ///< Reconstructed missing transverse energy (px)
	float MetPy;  ///< Reconstructed missing transverse energy (py)
	
	int nBJet;                   ///< Number of reconstructed b jets
	std::vector<float> BJetPt;   ///< b jet pt
	std::vector<float> BJetEta;  ///< b jet eta
	std::vector<float> BJetPhi;  ///< b jet phi
	std::vector<float> BJetE;    ///< b jet energy
	std::vector<float> BJetPx;   ///< b jet px
	std::vector<float> BJetPy;   ///< b jet py
	std::vector<float> BJetPz;   ///< b jet pz

	int nNonBJet;                   ///< Number of reconstructed non-b jets
	std::vector<float> NonBJetPt;   ///< Non b jet pt
	std::vector<float> NonBJetEta;  ///< Non b jet eta
	std::vector<float> NonBJetPhi;  ///< Non b jet phi
	std::vector<float> NonBJetE;    ///< Non b jet energy
	std::vector<float> NonBJetPx;   ///< Non b jet px
	std::vector<float> NonBJetPy;   ///< Non b jet py
	std::vector<float> NonBJetPz;   ///< Non b jet pz

	int nElectron;                     ///< Number of reconstructed electrons
	std::vector<float> ElectronPt;     ///< Electron pt
	std::vector<float> ElectronEta;    ///< Electron eta
	std::vector<float> ElectronPhi;    ///< Electron phi
	std::vector<float> ElectronE;      ///< Electron energy
	std::vector<float> ElectronPx;     ///< Electron px
	std::vector<float> ElectronPy;     ///< Electron py
	std::vector<float> ElectronPz;     ///< Electron pz
	std::vector<int> ElectronCharge;   ///< Electron charge

	int nMuon;                       ///< Number of reconstructed muons
	std::vector<float> MuonPt;       ///< Muon pt
	std::vector<float> MuonEta;      ///< Muon eta
	std::vector<float> MuonPhi;      ///< Muon phi
	std::vector<float> MuonE;        ///< Muon E
	std::vector<float> MuonPx;       ///< Muon px
	std::vector<float> MuonPy;       ///< Muon py
	std::vector<float> MuonPz;       ///< Muon pz
	std::vector<int> MuonCharge;     ///< Muon charge

	int nLepton;                     ///< Number of reconstructed leptons
	std::vector<float> LeptonPt;     ///< Lepton pt
	std::vector<float> LeptonEta;    ///< Lepton eta
	std::vector<float> LeptonPhi;    ///< Lepton phi
	std::vector<float> LeptonE;      ///< Lepton E
	std::vector<float> LeptonPx;     ///< Lepton px
	std::vector<float> LeptonPy;     ///< Lepton py
	std::vector<float> LeptonPz;     ///< Lepton pz
	std::vector<int> LeptonCharge;   ///< Lepton charge
	std::vector<int> LeptonLabel;    ///< Lepton type
	std::vector<int> LeptonIdx;      ///< Lepton index in input collections

	int nPhoton;                     ///< Number of reconstructed photons
	std::vector<float> PhotonPt;     ///< Photon pt
	std::vector<float> PhotonEta;    ///< Photon eta
	std::vector<float> PhotonPhi;    ///< Photon phi
	std::vector<float> PhotonE;      ///< Photon E
	std::vector<float> PhotonPx;     ///< Photon px
	std::vector<float> PhotonPy;     ///< Photon py
	std::vector<float> PhotonPz;     ///< Photon pz
	
	///
	/// Results
	///
	
	float *MetPx_;  ///< Reconstructed missing transverse energy (px)
	float *MetPy_;  ///< Reconstructed missing transverse energy (py)
	
	int *PhotonOrigin_;  ///< Photon origin
	
	float **WMass_;  ///< Reconstructed mass of the W boson
	float **WP_;     ///< Reconstructed momentum of the W boson
	float **WPt_;    ///< Reconstructed pt of the W boson
	float **WPx_;    ///< Reconstructed px of the W boson
	float **WPy_;    ///< Reconstructed py of the W boson
	float **WPz_;    ///< Reconstructed pz of the W boson
	float **WEta_;   ///< Reconstructed eta of the W boson
	float **WPhi_;   ///< Reconstructed phi of the W boson
	float **WE_;     ///< Reconstructed energy of the W boson
	
	float **TopMass_;  ///< Reconstructed mass of the top quark
	float **TopP_;     ///< Reconstructed momentum of the top quark
	float **TopPt_;    ///< Reconstructed pt of the top quark
	float **TopPx_;    ///< Reconstructed px of the top quark
	float **TopPy_;    ///< Reconstructed py of the top quark
	float **TopPz_;    ///< Reconstructed pz of the top quark
	float **TopEta_;   ///< Reconstructed eta of the top quark
	float **TopPhi_;   ///< Reconstructed phi of the top quark
	float **TopE_;     ///< Reconstructed energy of the top quark

	float *drTopTop_;    ///< Reconstructed dR(top, top)
	float *detaTopTop_;  ///< Reconstructed dEta(top, top)
	float *dphiTopTop_;  ///< Reconstructed dPhi(top, top)
	float *mTopTop_;     ///< Reconstructed invariant mass of the two top quark system
	float *ptTopTop_;    ///< Reconstructed pt of the two top quark system
	float *pTopTop_;     ///< Reconstructed momentum of the two top quark system
	float *etaTopTop_;   ///< Reconstructed eta of the two top quark system
	float *phiTopTop_;   ///< Reconstructed phi of the two top quark system
	
	float **nuPx_;  ///< Reconstructed neutrino px
	float **nuPy_;  ///< Reconstructed neutrino py
	float **nuPz_;  ///< Reconstructed neutrino pz

	int *TopTopLepLep_Electron1Idx;  ///< Index of electron from the first top quark decay in the input collection
	int *TopTopLepLep_Muon1Idx;      ///< Index of muon from the first top quark decay in the input collection
	int *TopTopLepLep_Electron2Idx;  ///< Index of electron from the second top quark decay in the input collection
	int *TopTopLepLep_Muon2Idx;      ///< Index of muon from the second top quark decay in the input collection
	int *TopTopLepLep_BJet1Idx;      ///< Index of b jet from the first top quark decay in the input collection
	int *TopTopLepLep_BJet2Idx;      ///< Index of b jet from the second top quark decay in the input collection
	int *TopTopLepLep_PhotonIdx;     ///< Index of the radiated photon in the top quark decay in the input collection

	int *TopTopLepHad_ElectronIdx;   ///< Index of electron from the leptonic top quark decay in the input collection
	int *TopTopLepHad_MuonIdx;       ///< Index of muon from the leptonic top quark decay in the input collection
	int *TopTopLepHad_BJetLepIdx;    ///< Index of b jet from the leptonic top quark decay in the input collection
	int *TopTopLepHad_BJetHadIdx;    ///< Index of b jet from the hadronic top quark decay in the input collection
	int *TopTopLepHad_NonBJet1Idx;   ///< Index of the first non b jet from the hadronic top quark decay in the input collection
	int *TopTopLepHad_NonBJet2Idx;   ///< Index of the second non b jet from the hadronic top quark decay in the input collection
	int *TopTopLepHad_PhotonIdx;     ///< Index of the radiated photon in the top quark decay in the input collection

	int *TopLep_ElectronIdx;   ///< Index of electron from the leptonic top quark decay in the input collection
	int *TopLep_MuonIdx;       ///< Index of muon from the leptonic top quark decay in the input collection
	int *TopLep_BJetLepIdx;    ///< Index of b jet from the leptonic top quark decay in the input collection
	int *TopLep_PhotonIdx;     ///< Index of the radiated photon in the top quark decay in the input collection
	
	int NPERMEVENT;       ///< Number of permutations in the current event
	int NPERMEVENT_PREV;  ///< Number of permutations in the previous event
	
	static bool IncludePhotons_; ///< Include photons in permutations
	static bool CheckAllPhotonOrigins_; ///< Consider individual photon categories
	
	///
	/// Fit settings
	///
	
	int NToy_;
	float LHMaxGeneric_;
	float LHMaxMinuit_;
	int NFitMax_;
	int NGrid_;
	int NPerm_;
	int NTerm_;
//	static bool AddTopTopMassToNLL_;
	static bool SortByTopTopMass_;
	static bool SortByTopTopPz_;
	static bool SortByTopTopPzAbs_;
	static bool SortByNeutrinoPz_;
	static bool SortByNeutrinoPt_;
	static bool SortByNeutrinoP_;
	bool DoTiming_;
	bool DoFit_;
	int NNu_;
	
	static std::shared_ptr<std::vector<float> > ParMin;
	static std::shared_ptr<std::vector<float> > ParMax;
	static std::shared_ptr<std::vector<bool> > IsParFixed;

	///
	/// Fit parameters
	///

	static std::string FPARAM_NAME[NPARMAX];

	float *chi_;
	std::string *chiTermName_;
	float **chiTerm_;
	int *NGeneric_;
	float **par_;
	int *idxMin_;
	
	///
	/// Timers
	///

	double *timerWallGeneric_;
	double *timerCPUGeneric_;
	double *timerWallMinuit_;
	double *timerCPUMinuit_;
	
	///
	/// Transfer functions
	///
	
	static std::string PDF_NAME[PDF_N];

	int NMetRMS_;

	int NBJetPxRMS_; 
	int NBJetPyRMS_;
	int NBJetPzRMS_;
	int NBJetERMS_;

	int NNonBJetPxRMS_;
	int NNonBJetPyRMS_;
	int NNonBJetPzRMS_;
	int NNonBJetERMS_;
	
	int NElecPxRMS_;
	int NElecPyRMS_;
	int NElecPzRMS_;
	int NElecERMS_;

	int NMuonPxRMS_;
	int NMuonPyRMS_;
	int NMuonPzRMS_;
	int NMuonERMS_;

	int NPhotonPxRMS_;
	int NPhotonPyRMS_;
	int NPhotonPzRMS_;
	int NPhotonERMS_;
	
	int LimNRMS_;
	
	float EtxMin_;
	float EtxMax_;

	static float maxPDFTopWMass;
	static float meanPDFTopWMass;
	static float sigmaPDFTopWMass;
	static double xminPDFTopWMass;
	static double xmaxPDFTopWMass;
	
	static float maxPDFTopMass;
	static float meanPDFTopMass;
	static float sigmaPDFTopMass;
	static double xminPDFTopMass;
	static double xmaxPDFTopMass;

	static float maxPDFTopWHadMass;
	static float meanPDFTopWHadMass;
	static float sigmaPDFTopWHadMass;
	static double xminPDFTopWHadMass;
	static double xmaxPDFTopWHadMass;
	
	static float maxPDFTopHadMass;
	static float meanPDFTopHadMass;
	static float sigmaPDFTopHadMass;
	static double xminPDFTopHadMass;
	static double xmaxPDFTopHadMass;
	
//	static float maxPDFTopTopMass;
//	static float meanPDFTopTopMass;
//	static float sigmaPDFTopTopMass;
//	static double xminPDFTopTopMass;
//	static double xmaxPDFTopTopMass;
	
	static float maxPDFBJetPx;
	static float meanPDFBJetPx;
	static float sigmaPDFBJetPx;
	static double xminPDFBJetPx;
	static double xmaxPDFBJetPx;

	static float maxPDFBJetPy;
	static float meanPDFBJetPy;
	static float sigmaPDFBJetPy;
	static double xminPDFBJetPy;
	static double xmaxPDFBJetPy;

	static float maxPDFBJetPz;
	static float meanPDFBJetPz;
	static float sigmaPDFBJetPz;
	static double xminPDFBJetPz;
	static double xmaxPDFBJetPz;

	static float maxPDFMetPx;
	static float meanPDFMetPx;
	static float sigmaPDFMetPx;
	static double xminPDFMetPx;
	static double xmaxPDFMetPx;
	
	static float maxPDFMetPy;
	static float meanPDFMetPy;
	static float sigmaPDFMetPy;
	static double xminPDFMetPy;
	static double xmaxPDFMetPy;
	
	static float maxPDFElecPx;
	static float meanPDFElecPx;
	static float sigmaPDFElecPx;
	static double xminPDFElecPx;
	static double xmaxPDFElecPx;
	
	static float maxPDFElecPy;
	static float meanPDFElecPy;
	static float sigmaPDFElecPy;
	static double xminPDFElecPy;
	static double xmaxPDFElecPy;
	
	static float maxPDFElecPz;
	static float meanPDFElecPz;
	static float sigmaPDFElecPz;
	static double xminPDFElecPz;
	static double xmaxPDFElecPz;
	
	static float maxPDFMuonPx;
	static float meanPDFMuonPx;
	static float sigmaPDFMuonPx;
	static double xminPDFMuonPx;
	static double xmaxPDFMuonPx;
	
	static float maxPDFMuonPy;
	static float meanPDFMuonPy;
	static float sigmaPDFMuonPy;
	static double xminPDFMuonPy;
	static double xmaxPDFMuonPy;
	
	static float maxPDFMuonPz;
	static float meanPDFMuonPz;
	static float sigmaPDFMuonPz;
	static double xminPDFMuonPz;
	static double xmaxPDFMuonPz;

	static float maxPDFNonBJetPx;
	static float meanPDFNonBJetPx;
	static float sigmaPDFNonBJetPx;
	static double xminPDFNonBJetPx;
	static double xmaxPDFNonBJetPx;

	static float maxPDFNonBJetPy;
	static float meanPDFNonBJetPy;
	static float sigmaPDFNonBJetPy;
	static double xminPDFNonBJetPy;
	static double xmaxPDFNonBJetPy;

	static float maxPDFNonBJetPz;
	static float meanPDFNonBJetPz;
	static float sigmaPDFNonBJetPz;
	static double xminPDFNonBJetPz;
	static double xmaxPDFNonBJetPz;

	static float maxPDFPhotonPx;
	static float meanPDFPhotonPx;
	static float sigmaPDFPhotonPx;
	static double xminPDFPhotonPx;
	static double xmaxPDFPhotonPx;

	static float maxPDFPhotonPy;
	static float meanPDFPhotonPy;
	static float sigmaPDFPhotonPy;
	static double xminPDFPhotonPy;
	static double xmaxPDFPhotonPy;

	static float maxPDFPhotonPz;
	static float meanPDFPhotonPz;
	static float sigmaPDFPhotonPz;
	static double xminPDFPhotonPz;
	static double xmaxPDFPhotonPz;
	
	///
	/// Additional variables
	///
	
	float Etx;     ///< Fix Etx parameter to a certain value
	float Ety;     ///< Fix Ety parameter to a certain value
	float WMass1;  ///< Fix the first W mass to a certain value
	float WMass2;  ///< Fix the second W mass to a certain value

	///
	/// Pointers to fit variables created on every event and defined transfer functions
	///
	
	static std::shared_ptr<double> CHISQ;
	
	static std::shared_ptr<double> EtMissX;
	static std::shared_ptr<double> EtMissY;
	
	static std::shared_ptr<double> PxLepton1;
	static std::shared_ptr<double> PyLepton1;
	static std::shared_ptr<double> PzLepton1;
	static std::shared_ptr<double> ELepton1;
	static std::shared_ptr<double> MassLepton1;
	static std::shared_ptr<int> LabelLepton1;
	
	static std::shared_ptr<double> PxLepton2;
	static std::shared_ptr<double> PyLepton2;
	static std::shared_ptr<double> PzLepton2;
	static std::shared_ptr<double> ELepton2;
	static std::shared_ptr<double> MassLepton2;
	static std::shared_ptr<int> LabelLepton2;
	
	static std::shared_ptr<double> PxBJet1;
	static std::shared_ptr<double> PyBJet1;
	static std::shared_ptr<double> PzBJet1;
	static std::shared_ptr<double> EBJet1;
	static std::shared_ptr<double> MassBJet1;
	
	static std::shared_ptr<double> PxBJet2;
	static std::shared_ptr<double> PyBJet2;
	static std::shared_ptr<double> PzBJet2;
	static std::shared_ptr<double> EBJet2;
	static std::shared_ptr<double> MassBJet2;

	static std::shared_ptr<double> PxNonBJet1;
	static std::shared_ptr<double> PyNonBJet1;
	static std::shared_ptr<double> PzNonBJet1;
	static std::shared_ptr<double> ENonBJet1;
	static std::shared_ptr<double> MassNonBJet1;
	
	static std::shared_ptr<double> PxNonBJet2;
	static std::shared_ptr<double> PyNonBJet2;
	static std::shared_ptr<double> PzNonBJet2;
	static std::shared_ptr<double> ENonBJet2;
	static std::shared_ptr<double> MassNonBJet2;

	static std::shared_ptr<double> PxPhoton;
	static std::shared_ptr<double> PyPhoton;
	static std::shared_ptr<double> PzPhoton;
	static std::shared_ptr<double> EPhoton;
	static std::shared_ptr<int> PhotonOrigin;
	
	static std::shared_ptr<double> PxNu1;
	static std::shared_ptr<double> PyNu1;
	static std::shared_ptr<double> PzNu1;
	
	static std::shared_ptr<double> PxNu2;
	static std::shared_ptr<double> PyNu2;
	static std::shared_ptr<double> PzNu2;
	
	static std::shared_ptr<double> TopTopMass;
	static std::shared_ptr<double> TopTopPz;
	static std::shared_ptr<double> TopTopPzAbs;
	
	static std::shared_ptr<std::vector<double> > FitParam;
	static std::shared_ptr<std::vector<double> > ChiTerm;
	static std::shared_ptr<std::vector<std::string> > ChiTermName;
	
	static std::shared_ptr<TF1> hPDFTopWMass;
	static std::shared_ptr<TF1> hPDFTopMass;
	static std::shared_ptr<TF1> hPDFTopWHadMass;
	static std::shared_ptr<TF1> hPDFTopHadMass;
//	static std::shared_ptr<TF1> hPDFTopTopMass;
	
	static std::shared_ptr<TF1> hPDFMetPx;
	static std::shared_ptr<TF1> hPDFMetPy;
	
	static std::shared_ptr<TF1> hPDFBJetPx;
	static std::shared_ptr<TF1> hPDFBJetPy;
	static std::shared_ptr<TF1> hPDFBJetPz;
	
	static std::shared_ptr<TF1> hPDFElecPx;
	static std::shared_ptr<TF1> hPDFElecPy;
	static std::shared_ptr<TF1> hPDFElecPz;
	
	static std::shared_ptr<TF1> hPDFMuonPx;
	static std::shared_ptr<TF1> hPDFMuonPy;
	static std::shared_ptr<TF1> hPDFMuonPz;
	
	static std::shared_ptr<TF1> hPDFNonBJetPx;
	static std::shared_ptr<TF1> hPDFNonBJetPy;
	static std::shared_ptr<TF1> hPDFNonBJetPz;

	static std::shared_ptr<TF1> hPDFPhotonPx;
	static std::shared_ptr<TF1> hPDFPhotonPy;
	static std::shared_ptr<TF1> hPDFPhotonPz;
     };

}

#endif
