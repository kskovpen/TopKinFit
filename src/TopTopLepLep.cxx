// Top quark pair production in dilepton channel: t -> (W->lnu) b, tbar -> (W->lnu) b

#include "../include/TopTopLepLep.h"

ClassImp(KINFIT::TopTopLepLep)

KINFIT::TopTopLepLep::TopTopLepLep() {}

KINFIT::TopTopLepLep::~TopTopLepLep() {}

/// Main routine
void KINFIT::TopTopLepLep::TopTopLepLepRun()
{
   NPERMEVENT_PREV = NPERMEVENT;

   /// Calculate number of permutations in the current event
   CalcNPerm(); 

   /// Reinitialize all variables and release memory used in the previous event
   
   NPerm_ = 0;
   NTerm_ = 0;

   if( idxMin_ ) {delete[] idxMin_; idxMin_ = 0;}
   
   if( TopTopLepLep_Electron1Idx )  {delete[] TopTopLepLep_Electron1Idx; TopTopLepLep_Electron1Idx = 0;}
   if( TopTopLepLep_Muon1Idx )      {delete[] TopTopLepLep_Muon1Idx; TopTopLepLep_Muon1Idx = 0;}
   if( TopTopLepLep_Electron2Idx )  {delete[] TopTopLepLep_Electron2Idx; TopTopLepLep_Electron2Idx = 0;}
   if( TopTopLepLep_Muon2Idx )      {delete[] TopTopLepLep_Muon2Idx; TopTopLepLep_Muon2Idx = 0;}
   if( TopTopLepLep_BJet1Idx )      {delete[] TopTopLepLep_BJet1Idx; TopTopLepLep_BJet1Idx = 0;}
   if( TopTopLepLep_BJet2Idx )      {delete[] TopTopLepLep_BJet2Idx; TopTopLepLep_BJet2Idx = 0;}
   if( TopTopLepLep_PhotonIdx )     {delete[] TopTopLepLep_PhotonIdx; TopTopLepLep_PhotonIdx = 0;}
   
   if( chi_ ) {delete[] chi_; chi_ = 0;}
   if( NGeneric_ ) {delete[] NGeneric_; NGeneric_ = 0;}
   if( par_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] par_[i]; delete[] par_; par_ = 0;}   
   if( chiTerm_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] chiTerm_[i]; delete[] chiTerm_; chiTerm_ = 0;}
   if( chiTermName_ ) {delete[] chiTermName_; chiTermName_ = 0;}
   if( timerWallGeneric_ ) {delete[] timerWallGeneric_; timerWallGeneric_ = 0;}
   if( timerCPUGeneric_ ) {delete[] timerCPUGeneric_; timerCPUGeneric_ = 0;}
   if( timerWallMinuit_ ) {delete[] timerWallMinuit_; timerWallMinuit_ = 0;}
   if( timerCPUMinuit_ ) {delete[] timerCPUMinuit_; timerCPUMinuit_ = 0;}
   
   if( MetPx_ ) {delete[] MetPx_; MetPx_ = 0;}
   if( MetPy_ ) {delete[] MetPy_; MetPy_ = 0;}
   
   if( PhotonOrigin_ ) {delete[] PhotonOrigin_; PhotonOrigin_ = 0;}   

   if( drTopTop_ ) {delete[] drTopTop_; drTopTop_ = 0;}
   if( detaTopTop_ ) {delete[] detaTopTop_; detaTopTop_ = 0;}
   if( dphiTopTop_ ) {delete[] dphiTopTop_; dphiTopTop_ = 0;}
   if( mTopTop_ ) {delete[] mTopTop_; mTopTop_ = 0;}
   if( ptTopTop_ ) {delete[] ptTopTop_; ptTopTop_ = 0;}
   if( pTopTop_ ) {delete[] pTopTop_; pTopTop_ = 0;}
   if( etaTopTop_ ) {delete[] etaTopTop_; etaTopTop_ = 0;}
   if( phiTopTop_ ) {delete[] phiTopTop_; phiTopTop_ = 0;}
   
   if( WMass_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] WMass_[i]; delete[] WMass_; WMass_ = 0;}
   if( WPt_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] WPt_[i]; delete[] WPt_; WPt_ = 0;}
   if( WPx_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] WPx_[i]; delete[] WPx_; WPx_ = 0;}
   if( WPy_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] WPy_[i]; delete[] WPy_; WPy_ = 0;}
   if( WPz_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] WPz_[i]; delete[] WPz_; WPz_ = 0;}
   if( WP_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] WP_[i]; delete[] WP_; WP_ = 0;}
   if( WEta_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] WEta_[i]; delete[] WEta_; WEta_ = 0;}
   if( WPhi_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] WPhi_[i]; delete[] WPhi_; WPhi_ = 0;}
   if( WE_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] WE_[i]; delete[] WE_; WE_ = 0;}
   
   if( TopMass_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] TopMass_[i]; delete[] TopMass_; TopMass_ = 0;}
   if( TopPt_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] TopPt_[i]; delete[] TopPt_; TopPt_ = 0;}
   if( TopPx_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] TopPx_[i]; delete[] TopPx_; TopPx_ = 0;}
   if( TopPy_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] TopPy_[i]; delete[] TopPy_; TopPy_ = 0;}
   if( TopPz_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] TopPz_[i]; delete[] TopPz_; TopPz_ = 0;}
   if( TopP_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] TopP_[i]; delete[] TopP_; TopP_ = 0;}
   if( TopEta_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] TopEta_[i]; delete[] TopEta_; TopEta_ = 0;}
   if( TopPhi_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] TopPhi_[i]; delete[] TopPhi_; TopPhi_ = 0;}
   if( TopE_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] TopE_[i]; delete[] TopE_; TopE_ = 0;}
   
   if( nuPx_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] nuPx_[i]; delete[] nuPx_; nuPx_ = 0;}
   if( nuPy_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] nuPy_[i]; delete[] nuPy_; nuPy_ = 0;}
   if( nuPz_ ) {for(int i=0;i<NPERMEVENT_PREV;i++) delete[] nuPz_[i]; delete[] nuPz_; nuPz_ = 0;}
   
   /// Allocate memory for the current event
   
   idxMin_ = new int[NPERMEVENT];
   
   TopTopLepLep_Electron1Idx = new int[NPERMEVENT];
   TopTopLepLep_Muon1Idx = new int[NPERMEVENT];
   TopTopLepLep_Electron2Idx = new int[NPERMEVENT];
   TopTopLepLep_Muon2Idx = new int[NPERMEVENT];
   TopTopLepLep_BJet1Idx = new int[NPERMEVENT];
   TopTopLepLep_BJet2Idx = new int[NPERMEVENT];
   TopTopLepLep_PhotonIdx = new int[NPERMEVENT];
   
   chi_ = new float[NPERMEVENT];
   NGeneric_ = new int[NPERMEVENT];
   par_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) par_[i] = new float[FPARAM_N];
   chiTermName_ = new std::string[NTERMMAX];
   chiTerm_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) chiTerm_[i] = new float[NTERMMAX];
   timerWallGeneric_ = new double[NPERMEVENT];
   timerCPUGeneric_ = new double[NPERMEVENT];
   timerWallMinuit_ = new double[NPERMEVENT];
   timerCPUMinuit_ = new double[NPERMEVENT];
   
   MetPx_ = new float[NPERMEVENT];
   MetPy_ = new float[NPERMEVENT];
   
   PhotonOrigin_ = new int[NPERMEVENT];

   drTopTop_ = new float[NPERMEVENT];
   detaTopTop_ = new float[NPERMEVENT];
   dphiTopTop_ = new float[NPERMEVENT];
   mTopTop_ = new float[NPERMEVENT];
   ptTopTop_ = new float[NPERMEVENT];
   pTopTop_ = new float[NPERMEVENT];
   etaTopTop_ = new float[NPERMEVENT];
   phiTopTop_ = new float[NPERMEVENT];
   
   WMass_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WMass_[i] = new float[NNUMAX];
   WPt_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WPt_[i] = new float[NNUMAX];
   WPx_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WPx_[i] = new float[NNUMAX];
   WPy_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WPy_[i] = new float[NNUMAX];
   WPz_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WPz_[i] = new float[NNUMAX];
   WP_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WP_[i] = new float[NNUMAX];
   WEta_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WEta_[i] = new float[NNUMAX];
   WPhi_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WPhi_[i] = new float[NNUMAX];
   WE_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WE_[i] = new float[NNUMAX];
   
   TopMass_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopMass_[i] = new float[NNUMAX];
   TopPt_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopPt_[i] = new float[NNUMAX];
   TopPx_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopPx_[i] = new float[NNUMAX];
   TopPy_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopPy_[i] = new float[NNUMAX];
   TopPz_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopPz_[i] = new float[NNUMAX];
   TopP_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopP_[i] = new float[NNUMAX];
   TopEta_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopEta_[i] = new float[NNUMAX];
   TopPhi_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopPhi_[i] = new float[NNUMAX];
   TopE_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopE_[i] = new float[NNUMAX];
   
   nuPx_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) nuPx_[i] = new float[NNUMAX];
   nuPy_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) nuPy_[i] = new float[NNUMAX];
   nuPz_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) nuPz_[i] = new float[NNUMAX];

   /// Initialize output variables

   const unsigned long int INIT = 10E+10;
   
   for( int i=0;i<NPERMEVENT;i++ )
     {
	drTopTop_[i] = INIT;
	detaTopTop_[i] = INIT;
	dphiTopTop_[i] = INIT;
	mTopTop_[i] = INIT;
	ptTopTop_[i] = INIT;
	pTopTop_[i] = INIT;
	etaTopTop_[i] = INIT;
	phiTopTop_[i] = INIT;

	TopTopLepLep_Electron1Idx[i] = -1;
	TopTopLepLep_Muon1Idx[i] = -1;
	TopTopLepLep_Electron2Idx[i] = -1;
	TopTopLepLep_Muon2Idx[i] = -1;
	TopTopLepLep_BJet1Idx[i] = -1;
	TopTopLepLep_BJet2Idx[i] = -1;
	TopTopLepLep_PhotonIdx[i] = -1;
	
	chi_[i] = INIT;
	NGeneric_[i] = -1;
	
	timerWallGeneric_[i] = INIT;
	timerCPUGeneric_[i] = INIT;
	timerWallMinuit_[i] = INIT;
	timerCPUMinuit_[i] = INIT;

	MetPx_[i] = INIT;
	MetPy_[i] = INIT;
	
	PhotonOrigin_[i] = -1;
	
	for( int ipar=0;ipar<FPARAM_N;ipar++ ) par_[i][ipar] = INIT;

	for( int it=0;it<NTERMMAX;it++ ) chiTerm_[i][it] = INIT;
	
	for( int inu=0;inu<NNUMAX;inu++ )
	  {	     
	     WMass_[i][inu] = INIT;
	     WPt_[i][inu] = INIT;
	     WPx_[i][inu] = INIT;
	     WPy_[i][inu] = INIT;
	     WPz_[i][inu] = INIT;
	     WP_[i][inu] = INIT;
	     WEta_[i][inu] = INIT;
	     WPhi_[i][inu] = INIT;
	     WE_[i][inu] = INIT;
	     
	     TopMass_[i][inu] = INIT;
	     TopPt_[i][inu] = INIT;
	     TopPx_[i][inu] = INIT;
	     TopPy_[i][inu] = INIT;
	     TopPz_[i][inu] = INIT;
	     TopP_[i][inu] = INIT;
	     TopEta_[i][inu] = INIT;
	     TopPhi_[i][inu] = INIT;
	     TopE_[i][inu] = INIT;
	     
	     nuPx_[i][inu] = INIT;
	     nuPy_[i][inu] = INIT;
	     nuPz_[i][inu] = INIT;
	  }	
     }

   *EtMissX = MetPx;
   *EtMissY = MetPy;

   /// Set parameter ranges
   
   (*ParMin)[FPARAM_Etx_TOPTOPLEPLEP] = -1000.;
   (*ParMax)[FPARAM_Etx_TOPTOPLEPLEP] = 1000.;
   (*ParMin)[FPARAM_Ety_TOPTOPLEPLEP] = -1000.;
   (*ParMax)[FPARAM_Ety_TOPTOPLEPLEP] = 1000.;
   
   (*ParMin)[FPARAM_Sign1_TOPTOPLEPLEP] = -1.;
   (*ParMax)[FPARAM_Sign1_TOPTOPLEPLEP] = 1.;
   (*ParMin)[FPARAM_Sign2_TOPTOPLEPLEP] = -1.;
   (*ParMax)[FPARAM_Sign2_TOPTOPLEPLEP] = 1.;

   (*ParMin)[FPARAM_mW1_TOPTOPLEPLEP] = 50.;
   (*ParMax)[FPARAM_mW1_TOPTOPLEPLEP] = 110.;
   (*ParMin)[FPARAM_mW2_TOPTOPLEPLEP] = 50.;
   (*ParMax)[FPARAM_mW2_TOPTOPLEPLEP] = 110.;

   /// Loop through permutations

   for( int il=0;il<nLepton;il++ )
     {
	if( LeptonCharge[il] != 1 ) continue;
	
	for( int il2=0;il2<nLepton;il2++ )
	  {
	     if( il == il2 ) continue;
	     
	     if( LeptonCharge[il2] != -1 ) continue;	     
	     
	     (*ParMin)[FPARAM_EtRealX_TOPTOPLEPLEP] = (*EtMissX) - LimNRMS_*fabs(*EtMissX)*sigmaPDFMetPx;
	     (*ParMax)[FPARAM_EtRealX_TOPTOPLEPLEP] = (*EtMissX) + LimNRMS_*fabs(*EtMissX)*sigmaPDFMetPx;
	     (*ParMin)[FPARAM_EtRealY_TOPTOPLEPLEP] = (*EtMissY) - LimNRMS_*fabs(*EtMissY)*sigmaPDFMetPy;
	     (*ParMax)[FPARAM_EtRealY_TOPTOPLEPLEP] = (*EtMissY) + LimNRMS_*fabs(*EtMissY)*sigmaPDFMetPy;
	     
	     *PxLepton1 = LeptonPx[il];
	     *PyLepton1 = LeptonPy[il];
	     *PzLepton1 = LeptonPz[il];
	     *ELepton1 = LeptonE[il];
	     *MassLepton1 = (*ELepton1)*(*ELepton1) - (*PxLepton1)*(*PxLepton1) - (*PyLepton1)*(*PyLepton1) - (*PzLepton1)*(*PzLepton1);
	     *MassLepton1 = (*MassLepton1 > 0.) ? sqrt(*MassLepton1) : 0.;
	     *LabelLepton1 = LeptonLabel[il];

	     float sigmaPDFLeptonPx = sigmaPDFElecPx;
	     float sigmaPDFLeptonPy = sigmaPDFElecPy;
	     float sigmaPDFLeptonPz = sigmaPDFElecPz;
	     
	     if( *LabelLepton1 == 1 )
	       {
		  sigmaPDFLeptonPx = sigmaPDFMuonPx;
		  sigmaPDFLeptonPy = sigmaPDFMuonPy;
		  sigmaPDFLeptonPz = sigmaPDFMuonPz;
	       }	     

	     (*ParMin)[FPARAM_Lepton1Px_TOPTOPLEPLEP] = (*PxLepton1) - LimNRMS_*fabs(*PxLepton1)*sigmaPDFLeptonPx;
	     (*ParMax)[FPARAM_Lepton1Px_TOPTOPLEPLEP] = (*PxLepton1) + LimNRMS_*fabs(*PxLepton1)*sigmaPDFLeptonPx;
	     (*ParMin)[FPARAM_Lepton1Py_TOPTOPLEPLEP] = (*PyLepton1) - LimNRMS_*fabs(*PyLepton1)*sigmaPDFLeptonPy;
	     (*ParMax)[FPARAM_Lepton1Py_TOPTOPLEPLEP] = (*PyLepton1) + LimNRMS_*fabs(*PyLepton1)*sigmaPDFLeptonPy;
	     (*ParMin)[FPARAM_Lepton1Pz_TOPTOPLEPLEP] = (*PzLepton1) - LimNRMS_*fabs(*PzLepton1)*sigmaPDFLeptonPz;
	     (*ParMax)[FPARAM_Lepton1Pz_TOPTOPLEPLEP] = (*PzLepton1) + LimNRMS_*fabs(*PzLepton1)*sigmaPDFLeptonPz;
	     
	     int label1_l = LeptonLabel[il];
	     int idx1_l = LeptonIdx[il];

	     *PxLepton2 = LeptonPx[il2];
	     *PyLepton2 = LeptonPy[il2];
	     *PzLepton2 = LeptonPz[il2];
	     *ELepton2 = LeptonE[il2];
	     *MassLepton2 = (*ELepton2)*(*ELepton2) - (*PxLepton2)*(*PxLepton2) - (*PyLepton2)*(*PyLepton2) - (*PzLepton2)*(*PzLepton2);
	     *MassLepton2 = (*MassLepton2 > 0.) ? sqrt(*MassLepton2) : 0.;
	     *LabelLepton2 = LeptonLabel[il2];

	     sigmaPDFLeptonPx = sigmaPDFElecPx;
	     sigmaPDFLeptonPy = sigmaPDFElecPy;
	     sigmaPDFLeptonPz = sigmaPDFElecPz;

	     if( *LabelLepton2 == 1 )
	       {
		  sigmaPDFLeptonPx = sigmaPDFMuonPx;
		  sigmaPDFLeptonPy = sigmaPDFMuonPy;
		  sigmaPDFLeptonPz = sigmaPDFMuonPz;
	       }	     

	     (*ParMin)[FPARAM_Lepton2Px_TOPTOPLEPLEP] = (*PxLepton2) - LimNRMS_*fabs(*PxLepton2)*sigmaPDFLeptonPx;
	     (*ParMax)[FPARAM_Lepton2Px_TOPTOPLEPLEP] = (*PxLepton2) + LimNRMS_*fabs(*PxLepton2)*sigmaPDFLeptonPx;
	     (*ParMin)[FPARAM_Lepton2Py_TOPTOPLEPLEP] = (*PyLepton2) - LimNRMS_*fabs(*PyLepton2)*sigmaPDFLeptonPy;
	     (*ParMax)[FPARAM_Lepton2Py_TOPTOPLEPLEP] = (*PyLepton2) + LimNRMS_*fabs(*PyLepton2)*sigmaPDFLeptonPy;
	     (*ParMin)[FPARAM_Lepton2Pz_TOPTOPLEPLEP] = (*PzLepton2) - LimNRMS_*fabs(*PzLepton2)*sigmaPDFLeptonPz;
	     (*ParMax)[FPARAM_Lepton2Pz_TOPTOPLEPLEP] = (*PzLepton2) + LimNRMS_*fabs(*PzLepton2)*sigmaPDFLeptonPz;

	     int label2_l = LeptonLabel[il2];
	     int idx2_l = LeptonIdx[il2];
	     
	     for( int ib=0;ib<nBJet;ib++ )
	       {
		  *EBJet1 = BJetE[ib];
		  *PxBJet1 = BJetPx[ib];
		  *PyBJet1 = BJetPy[ib];
		  *PzBJet1 = BJetPz[ib];
		  *MassBJet1 = (*EBJet1)*(*EBJet1) - (*PxBJet1)*(*PxBJet1) - (*PyBJet1)*(*PyBJet1) - (*PzBJet1)*(*PzBJet1);
		  *MassBJet1 = (*MassBJet1 > 0.) ? sqrt(*MassBJet1) : 0.;

		  (*ParMin)[FPARAM_BJet1Px_TOPTOPLEPLEP] = (*PxBJet1) - LimNRMS_*fabs(*PxBJet1)*sigmaPDFBJetPx;
		  (*ParMax)[FPARAM_BJet1Px_TOPTOPLEPLEP] = (*PxBJet1) + LimNRMS_*fabs(*PxBJet1)*sigmaPDFBJetPx;
		  (*ParMin)[FPARAM_BJet1Py_TOPTOPLEPLEP] = (*PyBJet1) - LimNRMS_*fabs(*PyBJet1)*sigmaPDFBJetPy;
		  (*ParMax)[FPARAM_BJet1Py_TOPTOPLEPLEP] = (*PyBJet1) + LimNRMS_*fabs(*PyBJet1)*sigmaPDFBJetPy;
		  (*ParMin)[FPARAM_BJet1Pz_TOPTOPLEPLEP] = (*PzBJet1) - LimNRMS_*fabs(*PzBJet1)*sigmaPDFBJetPz;
		  (*ParMax)[FPARAM_BJet1Pz_TOPTOPLEPLEP] = (*PzBJet1) + LimNRMS_*fabs(*PzBJet1)*sigmaPDFBJetPz;

		  for( int ib2=0;ib2<nBJet;ib2++ )
		    {
		       if( ib == ib2 ) continue;
		       
		       chi_[NPerm_] = INIT;
		       
		       *EBJet2 = BJetE[ib2];
		       *PxBJet2 = BJetPx[ib2];
		       *PyBJet2 = BJetPy[ib2];
		       *PzBJet2 = BJetPz[ib2];
		       *MassBJet2 = (*EBJet2)*(*EBJet2) - (*PxBJet2)*(*PxBJet2) - (*PyBJet2)*(*PyBJet2) - (*PzBJet2)*(*PzBJet2);
		       *MassBJet2 = (*MassBJet2 > 0.) ? sqrt(*MassBJet2) : 0.;

		       (*ParMin)[FPARAM_BJet2Px_TOPTOPLEPLEP] = (*PxBJet2) - LimNRMS_*fabs(*PxBJet2)*sigmaPDFBJetPx;
		       (*ParMax)[FPARAM_BJet2Px_TOPTOPLEPLEP] = (*PxBJet2) + LimNRMS_*fabs(*PxBJet2)*sigmaPDFBJetPx;
		       (*ParMin)[FPARAM_BJet2Py_TOPTOPLEPLEP] = (*PyBJet2) - LimNRMS_*fabs(*PyBJet2)*sigmaPDFBJetPy;
		       (*ParMax)[FPARAM_BJet2Py_TOPTOPLEPLEP] = (*PyBJet2) + LimNRMS_*fabs(*PyBJet2)*sigmaPDFBJetPy;
		       (*ParMin)[FPARAM_BJet2Pz_TOPTOPLEPLEP] = (*PzBJet2) - LimNRMS_*fabs(*PzBJet2)*sigmaPDFBJetPz;
		       (*ParMax)[FPARAM_BJet2Pz_TOPTOPLEPLEP] = (*PzBJet2) + LimNRMS_*fabs(*PzBJet2)*sigmaPDFBJetPz;
		       
		       for( int ipho=nPhoton;ipho>=0;ipho-- )
			 {
			    if( ipho < nPhoton ) // permute according to the number of photons + permutation with no radiation (ISR)
			      {
				 *EPhoton  = PhotonE[ipho];
				 *PxPhoton = PhotonPx[ipho];
				 *PyPhoton = PhotonPy[ipho];
				 *PzPhoton = PhotonPz[ipho];
				 
				 (*ParMin)[FPARAM_PhotonPx_TOPTOPLEPLEP] = (*PxPhoton) - LimNRMS_*fabs(*PxPhoton)*sigmaPDFPhotonPx;
				 (*ParMax)[FPARAM_PhotonPx_TOPTOPLEPLEP] = (*PxPhoton) + LimNRMS_*fabs(*PxPhoton)*sigmaPDFPhotonPx;
				 (*ParMin)[FPARAM_PhotonPy_TOPTOPLEPLEP] = (*PyPhoton) - LimNRMS_*fabs(*PyPhoton)*sigmaPDFPhotonPy;
				 (*ParMax)[FPARAM_PhotonPy_TOPTOPLEPLEP] = (*PyPhoton) + LimNRMS_*fabs(*PyPhoton)*sigmaPDFPhotonPy;
				 (*ParMin)[FPARAM_PhotonPz_TOPTOPLEPLEP] = (*PzPhoton) - LimNRMS_*fabs(*PzPhoton)*sigmaPDFPhotonPz;
				 (*ParMax)[FPARAM_PhotonPz_TOPTOPLEPLEP] = (*PzPhoton) + LimNRMS_*fabs(*PzPhoton)*sigmaPDFPhotonPz;
			      }
			    
			    for( int ipo=0;ipo<PHOTON_ORIGIN_N_TOPTOPLEPLEP;ipo++ )
			      {
				 if( ipo != PHOTON_FROM_TOP1_COMB_TOPTOPLEPLEP &&
				     ipo != PHOTON_FROM_W1_COMB_TOPTOPLEPLEP &&
				     ipo != PHOTON_FROM_TOP2_COMB_TOPTOPLEPLEP &&
				     ipo != PHOTON_FROM_W2_COMB_TOPTOPLEPLEP &&
				     ipo != PHOTON_FROM_ISR_TOPTOPLEPLEP &&
				     !CheckAllPhotonOrigins_ ) continue;
				 
				 if( !IncludePhotons_ && ipo != PHOTON_FROM_ISR_TOPTOPLEPLEP ) continue;
				 else if( IncludePhotons_ && ipho == nPhoton ) continue;
				 
				 *PhotonOrigin = ipo;

				 std::vector<FRESULT> vp;

				 struct timespec walltime_calcNuGrid_begin, walltime_calcNuGrid_end;
				 struct timespec cputime_calcNuGrid_begin, cputime_calcNuGrid_end;
				 
				 if( DoTiming_ )
				   {			    
				      clock_gettime(CLOCK_REALTIME, &walltime_calcNuGrid_begin);
				      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cputime_calcNuGrid_begin);
				   }

				 /// Produce a grid with valid solutions (first-layer minimization)
				 calcNuGrid(vp);
		       
				 NGeneric_[NPerm_] = vp.size();
		       
				 /// Sort according to NLL
				 std::sort(vp.begin(), vp.end(), [](const FRESULT &lhs, const FRESULT &rhs) { return lhs.lh < rhs.lh; });

				 if( DoTiming_ )
				   {			    
				      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cputime_calcNuGrid_end);
				      clock_gettime(CLOCK_REALTIME, &walltime_calcNuGrid_end);
				      
				      long walltime_calcNuGrid_seconds = walltime_calcNuGrid_end.tv_sec - walltime_calcNuGrid_begin.tv_sec;
				      long walltime_calcNuGrid_nanoseconds = walltime_calcNuGrid_end.tv_nsec - walltime_calcNuGrid_begin.tv_nsec;
				      double walltime_calcNuGrid_elapsed = walltime_calcNuGrid_seconds + walltime_calcNuGrid_nanoseconds*1E-9;
				      
				      long cputime_calcNuGrid_seconds = cputime_calcNuGrid_end.tv_sec - cputime_calcNuGrid_begin.tv_sec;
				      long cputime_calcNuGrid_nanoseconds = cputime_calcNuGrid_end.tv_nsec - cputime_calcNuGrid_begin.tv_nsec;
				      double cputime_calcNuGrid_elapsed = cputime_calcNuGrid_seconds + cputime_calcNuGrid_nanoseconds*1E-9;
				      
				      timerWallGeneric_[NPerm_] = walltime_calcNuGrid_elapsed;
				      timerCPUGeneric_[NPerm_] = cputime_calcNuGrid_elapsed;
				   }
			    
				 float disc_ = INIT;

				 struct timespec walltime_fit_begin, walltime_fit_end;
				 struct timespec cputime_fit_begin, cputime_fit_end;
				 
				 if( DoTiming_ )
				   {			    
				      clock_gettime(CLOCK_REALTIME, &walltime_fit_begin);
				      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cputime_fit_begin);
				   }		       
				 
				 std::vector<FRESULT> lhres = vp;
				 
				 if( DoFit_ )
				   {
				      std::vector<FRESULT> vpp;
		       
				      /// Second-layer minimization
				      for( int ig=0;ig<vp.size();ig++ )
					{
					   if( ig >= NFitMax_ ) break;
			    
					   double par[FPARAM_N];
					   for( int ip=0;ip<FPARAM_N;ip++ ) par[ip] = vp[ig].par[ip];

					   fit(par, vpp);
					}

				      /// Sort according to NLL
				      std::sort(vpp.begin(), vpp.end(), [](const FRESULT &lhs, const FRESULT &rhs) { return lhs.lh < rhs.lh; });
				      
				      lhres = vpp;
				   }
				 		       
				 /// Apply additional sorting according to various kinematic variables (for good fits) split in different neutrino solutions, and NLL (for bad ones)
				 std::vector<FRESULT> lhgood_s1Plus_s2Plus, lhgood_s1Plus_s2Minus, lhgood_s1Minus_s2Plus, lhgood_s1Minus_s2Minus;
				 std::vector<FRESULT> lhbad;
		       
				 for( int ig=0;ig<lhres.size();ig++ )
				   {
				      FRESULT r = lhres[ig];
			    
				      if( r.lh < LHMaxMinuit_ ) 
					{
					   if( r.par[FPARAM_Sign1_TOPTOPLEPLEP] == 1 && r.par[FPARAM_Sign2_TOPTOPLEPLEP] == 1 ) lhgood_s1Plus_s2Plus.push_back(r);
					   else if( r.par[FPARAM_Sign1_TOPTOPLEPLEP] == 1 && r.par[FPARAM_Sign2_TOPTOPLEPLEP] == -1 ) lhgood_s1Plus_s2Minus.push_back(r);
					   else if( r.par[FPARAM_Sign1_TOPTOPLEPLEP] == -1 && r.par[FPARAM_Sign2_TOPTOPLEPLEP] == 1 ) lhgood_s1Minus_s2Plus.push_back(r);
					   else lhgood_s1Minus_s2Minus.push_back(r);
					}
				      else lhbad.push_back(r);
				   }
		       
				 std::vector<FRESULT> lhgood;
		       
				 /// Select best NLL for each neutrino solution and sort it according to NLL
				 if( lhgood_s1Plus_s2Plus.size() > 0) lhgood.push_back(lhgood_s1Plus_s2Plus[0]);
				 if( lhgood_s1Plus_s2Minus.size() > 0) lhgood.push_back(lhgood_s1Plus_s2Minus[0]);
				 if( lhgood_s1Minus_s2Plus.size() > 0) lhgood.push_back(lhgood_s1Minus_s2Plus[0]);
				 if( lhgood_s1Minus_s2Minus.size() > 0) lhgood.push_back(lhgood_s1Minus_s2Minus[0]);
		       
				 std::sort(lhgood.begin(), lhgood.end(), [](const FRESULT &lhs, const FRESULT &rhs) { return lhs.lh < rhs.lh; });
		       
				 /// Sort best solutions
				 bool sorting[6] = {SortByTopTopMass_, SortByTopTopPz_, SortByTopTopPzAbs_, SortByNeutrinoPz_, SortByNeutrinoPt_, SortByNeutrinoP_};
				 int nsort = 0;
				 for( int s=0;s<6;s++ )
				   if( sorting[s] ) nsort++;
				 if( nsort > 1 )
				   {
				      std::cout << "Multiple sorting options are chosen, please stick to only one option" << std::endl;
				      exit(1);
				   }
			 
				 if( SortByTopTopMass_ ) std::sort(lhgood.begin(), lhgood.end(), [](const FRESULT &lhs, const FRESULT &rhs) { return lhs.TopTopMass < rhs.TopTopMass; });
				 else if( SortByTopTopPz_ ) std::sort(lhgood.begin(), lhgood.end(), [](const FRESULT &lhs, const FRESULT &rhs) { return lhs.TopTopPz < rhs.TopTopPz; });
				 else if( SortByTopTopPzAbs_ ) std::sort(lhgood.begin(), lhgood.end(), [](const FRESULT &lhs, const FRESULT &rhs) { return lhs.TopTopPzAbs < rhs.TopTopPzAbs; });
				 else if( SortByNeutrinoPz_ ) std::sort(lhgood.begin(), lhgood.end(), [](const FRESULT &lhs, const FRESULT &rhs) { return lhs.PzNuSum < rhs.PzNuSum; });
				 else if( SortByNeutrinoPt_ ) std::sort(lhgood.begin(), lhgood.end(), [](const FRESULT &lhs, const FRESULT &rhs) { return lhs.PtNuSum < rhs.PtNuSum; });
				 else if( SortByNeutrinoP_ ) std::sort(lhgood.begin(), lhgood.end(), [](const FRESULT &lhs, const FRESULT &rhs) { return lhs.PNuSum < rhs.PNuSum; });
		       
				 /// Append the best NLL solution from bad fits
				 if( lhbad.size() > 0 ) lhgood.push_back(lhbad[0]);
				 
				 if( lhgood.size() > 0 )
				   {			    
				      FRESULT fbest = lhgood[0];
			    
				      disc_ = fbest.lh;
				      
				      TopTopLepLep_PhotonIdx[NPerm_] = (IncludePhotons_) ? ipho : -1;
				      if( IncludePhotons_ && ipho == nPhoton ) TopTopLepLep_PhotonIdx[NPerm_] = -1;
				      
				      TopTopLepLep_BJet1Idx[NPerm_] = ib;
				      TopTopLepLep_BJet2Idx[NPerm_] = ib2;
				 
				      if( label1_l == 0 ) 
					{
					   TopTopLepLep_Electron1Idx[NPerm_] = idx1_l;
					   TopTopLepLep_Muon1Idx[NPerm_] = -1;
					}
				      else if( label1_l == 1 ) 
					{
					   TopTopLepLep_Muon1Idx[NPerm_] = idx1_l;
					   TopTopLepLep_Electron1Idx[NPerm_] = -1;
					}
				      
				      if( label2_l == 0 )
					{
					   TopTopLepLep_Electron2Idx[NPerm_] = idx2_l;
					   TopTopLepLep_Muon2Idx[NPerm_] = -1;
					}
				      else if( label2_l == 1 ) 
					{
					   TopTopLepLep_Muon2Idx[NPerm_] = idx2_l;
					   TopTopLepLep_Electron2Idx[NPerm_] = -1;
					}
				      
				      chi_[NPerm_] = disc_;
			    
				      NTerm_ = fbest.NTERM;
				      for( int ip=0;ip<FPARAM_N;ip++ )
					{
					   par_[NPerm_][ip] = fbest.par[ip];
					}
				      for( int ip=0;ip<NTerm_;ip++ )
					{
					   chiTerm_[NPerm_][ip] = fbest.chiTerm[ip];
					   chiTermName_[ip] = fbest.chiTermName[ip];
					}
			    
				      nuPx_[NPerm_][0] = fbest.par[FPARAM_EtRealX_TOPTOPLEPLEP]/2.+fbest.par[FPARAM_Etx_TOPTOPLEPLEP];
				      nuPy_[NPerm_][0] = fbest.par[FPARAM_EtRealY_TOPTOPLEPLEP]/2.+fbest.par[FPARAM_Ety_TOPTOPLEPLEP];
				      nuPz_[NPerm_][0] = fbest.PzNu1;
				      
				      nuPx_[NPerm_][1] = fbest.par[FPARAM_EtRealX_TOPTOPLEPLEP]/2.-fbest.par[FPARAM_Etx_TOPTOPLEPLEP];
				      nuPy_[NPerm_][1] = fbest.par[FPARAM_EtRealY_TOPTOPLEPLEP]/2.-fbest.par[FPARAM_Ety_TOPTOPLEPLEP];
				      nuPz_[NPerm_][1] = fbest.PzNu2;
				      
				      MetPx_[NPerm_] = fbest.par[FPARAM_EtRealX_TOPTOPLEPLEP];
				      MetPy_[NPerm_] = fbest.par[FPARAM_EtRealY_TOPTOPLEPLEP];
				      
				      WMass_[NPerm_][0] = fbest.par[FPARAM_mW1_TOPTOPLEPLEP];
				      WMass_[NPerm_][1] = fbest.par[FPARAM_mW2_TOPTOPLEPLEP];
				      
				      PhotonOrigin_[NPerm_] = *PhotonOrigin;
				   }

				 if( DoTiming_ )
				   {			    
				      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cputime_fit_end);
				      clock_gettime(CLOCK_REALTIME, &walltime_fit_end);
				      
				      long walltime_fit_seconds = walltime_fit_end.tv_sec - walltime_fit_begin.tv_sec;
				      long walltime_fit_nanoseconds = walltime_fit_end.tv_nsec - walltime_fit_begin.tv_nsec;
				      double walltime_fit_elapsed = walltime_fit_seconds + walltime_fit_nanoseconds*1E-9;
				      
				      long cputime_fit_seconds = cputime_fit_end.tv_sec - cputime_fit_begin.tv_sec;
				      long cputime_fit_nanoseconds = cputime_fit_end.tv_nsec - cputime_fit_begin.tv_nsec;
				      double cputime_fit_elapsed = cputime_fit_seconds + cputime_fit_nanoseconds*1E-9;
				      
				      timerWallMinuit_[NPerm_] = walltime_fit_elapsed;
				      timerCPUMinuit_[NPerm_] = cputime_fit_elapsed;
				   }
				 
				 NPerm_++;
			      }
			 }
		    }		  
	       }
	  }
     } // end permutations

   // choose the best permutation
   for(int ip=0;ip<NPerm_;ip++)
     {
	calcVar(ip);
     }
}

/// Calculate the number of permutations per event
void KINFIT::TopTopLepLep::CalcNPerm()
{
   NPERMEVENT = 0;

   for(int ip=nPhoton;ip>=0;ip--)
     {	
	for(int il=0;il<nLepton;il++)
	  {
	     if( LeptonCharge[il] != 1 ) continue;
	     
	     for(int il2=0;il2<nLepton;il2++)
	       {
		  if( il == il2 ) continue;
		  
		  if( LeptonCharge[il2] != -1 ) continue;
		  
		  for(int ib=0;ib<nBJet;ib++)
		    {
		       for(int ib2=0;ib2<nBJet;ib2++)
			 {
			    if( ib == ib2 ) continue;
			    
			    for(int ipo=0;ipo<PHOTON_ORIGIN_N_TOPTOPLEPLEP;ipo++)
			      {
				 if( ipo != PHOTON_FROM_TOP1_COMB_TOPTOPLEPLEP &&
				     ipo != PHOTON_FROM_W1_COMB_TOPTOPLEPLEP &&
				     ipo != PHOTON_FROM_TOP2_COMB_TOPTOPLEPLEP &&
				     ipo != PHOTON_FROM_W2_COMB_TOPTOPLEPLEP &&
				     ipo != PHOTON_FROM_ISR_TOPTOPLEPLEP &&
				     !CheckAllPhotonOrigins_ ) continue;
				 
				 if( !IncludePhotons_ && ip > 0 ) continue;
			    
				 NPERMEVENT++;
			      }			    
			 }
		    }
	       }
	  }
     }
}

/// Derive solutions on a given set of inputs and calculate the resultant NLL
double KINFIT::TopTopLepLep::func(float PxLepton1, float PyLepton1, float PzLepton1, float ELepton1, int LabelLepton1,
				  float PxLepton2, float PyLepton2, float PzLepton2, float ELepton2, int LabelLepton2,
				  float PxBJet1, float PyBJet1, float PzBJet1, float EBJet1,
				  float PxBJet2, float PyBJet2, float PzBJet2, float EBJet2,
				  float PxPhoton, float PyPhoton, float PzPhoton, float EPhoton,
				  int photonOrigin,
				  std::vector<double> &chi2t, double *par)
{
   chi2t[NLL_RacPT_TOPTOPLEPLEP] = 0.;
   chi2t[NLL_mWPT_TOPTOPLEPLEP] = 0.;
   chi2t[NLL_mTopPT_TOPTOPLEPLEP] = 0.;
         
   float PxNu1 = par[FPARAM_EtRealX_TOPTOPLEPLEP]/2.+par[FPARAM_Etx_TOPTOPLEPLEP];
   float PyNu1 = par[FPARAM_EtRealY_TOPTOPLEPLEP]/2.+par[FPARAM_Ety_TOPTOPLEPLEP];

   *KINFIT::kfit::PxNu1 = PxNu1;
   *KINFIT::kfit::PyNu1 = PyNu1;
   
   float PxNu2 = par[FPARAM_EtRealX_TOPTOPLEPLEP]/2.-par[FPARAM_Etx_TOPTOPLEPLEP];
   float PyNu2 = par[FPARAM_EtRealY_TOPTOPLEPLEP]/2.-par[FPARAM_Ety_TOPTOPLEPLEP];

   *KINFIT::kfit::PxNu2 = PxNu2;
   *KINFIT::kfit::PyNu2 = PyNu2;
   
   float Lepton1Px = par[FPARAM_Lepton1Px_TOPTOPLEPLEP];
   float Lepton1Py = par[FPARAM_Lepton1Py_TOPTOPLEPLEP];
   float Lepton1Pz = par[FPARAM_Lepton1Pz_TOPTOPLEPLEP];
   float Lepton1E = par[FPARAM_Lepton1E_TOPTOPLEPLEP];

   float Lepton2Px = par[FPARAM_Lepton2Px_TOPTOPLEPLEP];
   float Lepton2Py = par[FPARAM_Lepton2Py_TOPTOPLEPLEP];
   float Lepton2Pz = par[FPARAM_Lepton2Pz_TOPTOPLEPLEP];
   float Lepton2E = par[FPARAM_Lepton2E_TOPTOPLEPLEP];
   
   float PhotonPx = par[FPARAM_PhotonPx_TOPTOPLEPLEP];
   float PhotonPy = par[FPARAM_PhotonPy_TOPTOPLEPLEP];
   float PhotonPz = par[FPARAM_PhotonPz_TOPTOPLEPLEP];
   float PhotonE = par[FPARAM_PhotonE_TOPTOPLEPLEP];
   
   if( (photonOrigin == PHOTON_FROM_W1_TOPTOPLEPLEP ||
	photonOrigin == PHOTON_FROM_W1_COMB_TOPTOPLEPLEP ||
	photonOrigin == PHOTON_FROM_LEPTON1_TOPTOPLEPLEP) &&
       IncludePhotons_ )
     {
	Lepton1Px += PhotonPx;
	Lepton1Py += PhotonPy;
	Lepton1Pz += PhotonPz;
	Lepton1E += PhotonE;
     }
   
   if( (photonOrigin == PHOTON_FROM_W2_TOPTOPLEPLEP ||
	photonOrigin == PHOTON_FROM_W2_COMB_TOPTOPLEPLEP ||
	photonOrigin == PHOTON_FROM_LEPTON2_TOPTOPLEPLEP) &&
       IncludePhotons_ )
     {
	Lepton2Px += PhotonPx;
	Lepton2Py += PhotonPy;
	Lepton2Pz += PhotonPz;
	Lepton2E += PhotonE;
     }
      
   float a1 = sqrt(Lepton1E*Lepton1E-Lepton1Pz*Lepton1Pz);
   float a2 = sqrt(Lepton2E*Lepton2E-Lepton2Pz*Lepton2Pz);
   float b1 = Lepton1Pz;
   float b2 = Lepton2Pz;
   float d1 = sqrt(PxNu1*PxNu1+PyNu1*PyNu1);
   float d2 = sqrt(PxNu2*PxNu2+PyNu2*PyNu2);
   float f1 = Lepton1E;
   float f2 = Lepton2E;

   float c1 = (par[FPARAM_mW1_TOPTOPLEPLEP]*par[FPARAM_mW1_TOPTOPLEPLEP]-(Lepton1E*Lepton1E-Lepton1Px*Lepton1Px-Lepton1Py*Lepton1Py-Lepton1Pz*Lepton1Pz))/2.+Lepton1Px*PxNu1+Lepton1Py*PyNu1;
   float c2 = (par[FPARAM_mW2_TOPTOPLEPLEP]*par[FPARAM_mW2_TOPTOPLEPLEP]-(Lepton2E*Lepton2E-Lepton2Px*Lepton2Px-Lepton2Py*Lepton2Py-Lepton2Pz*Lepton2Pz))/2.+Lepton2Px*PxNu2+Lepton2Py*PyNu2;
   
   float rac1 = c1*c1*b1*b1-a1*a1*(d1*d1*f1*f1-c1*c1);
   float rac2 = c2*c2*b2*b2-a2*a2*(d2*d2*f2*f2-c2*c2);
   
   float racAbs1 = fabs(rac1);
   float racAbs2 = fabs(rac2);
   
   float racPos1 = rac1;
   float racPos2 = rac2;
   
   if( racPos1 < 0. ) racPos1 = 0.;
   if( racPos2 < 0. ) racPos2 = 0.;

   float PzNu1 = (c1*b1+par[FPARAM_Sign1_TOPTOPLEPLEP]*sqrt(racPos1))/a1/a1;
   float PzNu2 = (c2*b2+par[FPARAM_Sign2_TOPTOPLEPLEP]*sqrt(racPos2))/a2/a2;
   
   *KINFIT::kfit::PzNu1 = PzNu1;
   *KINFIT::kfit::PzNu2 = PzNu2;
   
   float ENu1 = sqrt(PxNu1*PxNu1+PyNu1*PyNu1+PzNu1*PzNu1);
   float ENu2 = sqrt(PxNu2*PxNu2+PyNu2*PyNu2+PzNu2*PzNu2);
   
   float BJet1Px = par[FPARAM_BJet1Px_TOPTOPLEPLEP];
   float BJet1Py = par[FPARAM_BJet1Py_TOPTOPLEPLEP];
   float BJet1Pz = par[FPARAM_BJet1Pz_TOPTOPLEPLEP];
   float BJet1E = par[FPARAM_BJet1E_TOPTOPLEPLEP];
   
   if( photonOrigin == PHOTON_FROM_BJET1_TOPTOPLEPLEP && IncludePhotons_ )
     {
	BJet1Px += PhotonPx;
	BJet1Py += PhotonPy;
	BJet1Pz += PhotonPz;
	BJet1E += PhotonE;
     }   

   float BJet2Px = par[FPARAM_BJet2Px_TOPTOPLEPLEP];
   float BJet2Py = par[FPARAM_BJet2Py_TOPTOPLEPLEP];
   float BJet2Pz = par[FPARAM_BJet2Pz_TOPTOPLEPLEP];
   float BJet2E = par[FPARAM_BJet2E_TOPTOPLEPLEP];
   
   if( photonOrigin == PHOTON_FROM_BJET2_TOPTOPLEPLEP && IncludePhotons_ )
     {
	BJet2Px += PhotonPx;
	BJet2Py += PhotonPy;
	BJet2Pz += PhotonPz;
	BJet2E += PhotonE;
     }   
   
   float totPx1 = BJet1Px+PxNu1+Lepton1Px;
   float totPy1 = BJet1Py+PyNu1+Lepton1Py;
   float totPz1 = BJet1Pz+PzNu1+Lepton1Pz;
   float totE1 = BJet1E+ENu1+Lepton1E;
   
   if( (photonOrigin == PHOTON_FROM_TOP1_TOPTOPLEPLEP ||
	photonOrigin == PHOTON_FROM_TOP1_COMB_TOPTOPLEPLEP) && IncludePhotons_ )
     {
	totPx1 += PhotonPx;
	totPy1 += PhotonPy;
	totPz1 += PhotonPz;
	totE1 += PhotonE;
     }

   float totPx2 = BJet2Px+PxNu2+Lepton2Px;
   float totPy2 = BJet2Py+PyNu2+Lepton2Py;
   float totPz2 = BJet2Pz+PzNu2+Lepton2Pz;
   float totE2 = BJet2E+ENu2+Lepton2E;

   if( (photonOrigin == PHOTON_FROM_TOP2_TOPTOPLEPLEP ||
	photonOrigin == PHOTON_FROM_TOP2_COMB_TOPTOPLEPLEP) && IncludePhotons_ )
     {
	totPx2 += PhotonPx;
	totPy2 += PhotonPy;
	totPz2 += PhotonPz;
	totE2 += PhotonE;
     }
   
   float ENuL1 = ENu1+Lepton1E;
   float PxNuL1 = PxNu1+Lepton1Px;
   float PyNuL1 = PyNu1+Lepton1Py;
   float PzNuL1 = PzNu1+Lepton1Pz;
   float mW1 = ENuL1*ENuL1-PxNuL1*PxNuL1-PyNuL1*PyNuL1-PzNuL1*PzNuL1;
   float mtop1 = totE1*totE1-totPx1*totPx1-totPy1*totPy1-totPz1*totPz1;

   float ENuL2 = ENu2+Lepton2E;
   float PxNuL2 = PxNu2+Lepton2Px;
   float PyNuL2 = PyNu2+Lepton2Py;
   float PzNuL2 = PzNu2+Lepton2Pz;
   float mW2 = ENuL2*ENuL2-PxNuL2*PxNuL2-PyNuL2*PyNuL2-PzNuL2*PzNuL2;
   float mtop2 = totE2*totE2-totPx2*totPx2-totPy2*totPy2-totPz2*totPz2;
   
   float totE1E2 = totE1+totE2;
   float totPx1Px2 = totPx1+totPx2;
   float totPy1Py2 = totPy1+totPy2;
   float totPz1Pz2 = totPz1+totPz2;
  
   float mtoptop = totE1E2*totE1E2-totPx1Px2*totPx1Px2-totPy1Py2*totPy1Py2-totPz1Pz2*totPz1Pz2;
   float mtoptopAbs = sqrt(fabs(mtoptop));
   
   *KINFIT::kfit::TopTopMass = mtoptopAbs;
   *KINFIT::kfit::TopTopPz = fabs(totPz1Pz2);
   *KINFIT::kfit::TopTopPzAbs = fabs(totPz1) + fabs(totPz2);

   double val = 0.;

   /// Penalty term for negative roots
   if( rac1 < 0 ) 
     {
	float rac1PT = log(1.+racAbs1);
	val += rac1PT;
	chi2t[NLL_RacPT_TOPTOPLEPLEP] += rac1PT;
     }   
   if( rac2 < 0 )
     {
	float rac2PT = log(1.+racAbs2);
	val += rac2PT;
	chi2t[NLL_RacPT_TOPTOPLEPLEP] += rac2PT;
     }   

   float mW1Abs = fabs(mW1);
   float mW2Abs = fabs(mW2);
   float mtop1Abs = fabs(mtop1);
   float mtop2Abs = fabs(mtop2);
   
   /// Penalty terms for negative masses
   if( mW1 < 0 ) 
     {
	float mW1PT = log(1.+mW1Abs);
	val += mW1PT;
	chi2t[NLL_mWPT_TOPTOPLEPLEP] += mW1PT;
     }   
   if( mW2 < 0 )
     {
	float mW2PT = log(1.+mW2Abs);
	val += mW2PT;
	chi2t[NLL_mWPT_TOPTOPLEPLEP] += mW2PT;
     }   
   if( mtop1 < 0 )
     {
	float mTop1PT = log(1.+mtop1Abs);
	val += mTop1PT;
	chi2t[NLL_mTopPT_TOPTOPLEPLEP] += mTop1PT;
     }   
   if( mtop2 < 0 ) 
     {
	float mTop2PT = log(1.+mtop2Abs);
	val += mTop2PT;
	chi2t[NLL_mTopPT_TOPTOPLEPLEP] += mTop2PT;
     }   
   
   mtop1Abs = sqrt(mtop1Abs);
   mW1Abs = sqrt(mW1Abs);

   mtop2Abs = sqrt(mtop2Abs);
   mW2Abs = sqrt(mW2Abs);

   float mW1Prob = getProb(hPDFTopWMass.get(), mW1Abs, maxPDFTopWMass, xminPDFTopWMass, xmaxPDFTopWMass);
   float mW2Prob = getProb(hPDFTopWMass.get(), mW2Abs, maxPDFTopWMass, xminPDFTopWMass, xmaxPDFTopWMass);

   float mTop1Prob = getProb(hPDFTopMass.get(), mtop1Abs, maxPDFTopMass, xminPDFTopMass, xmaxPDFTopMass);
   float mTop2Prob = getProb(hPDFTopMass.get(), mtop2Abs, maxPDFTopMass, xminPDFTopMass, xmaxPDFTopMass);
   
//   float mTopTopProb = (AddTopTopMassToNLL_) ? getProb(hPDFTopTopMass.get(), mtoptopAbs, maxPDFTopTopMass, xminPDFTopTopMass, xmaxPDFTopTopMass) : 1.0;

   float MetPxProb = (! (*IsParFixed)[FPARAM_EtRealX_TOPTOPLEPLEP]) ? getProb(hPDFMetPx.get(), (par[FPARAM_EtRealX_TOPTOPLEPLEP]-*EtMissX)/par[FPARAM_EtRealX_TOPTOPLEPLEP], maxPDFMetPx, xminPDFMetPx, xmaxPDFMetPx) : 1.0;
   float MetPyProb = (! (*IsParFixed)[FPARAM_EtRealY_TOPTOPLEPLEP]) ? getProb(hPDFMetPy.get(), (par[FPARAM_EtRealY_TOPTOPLEPLEP]-*EtMissY)/par[FPARAM_EtRealY_TOPTOPLEPLEP], maxPDFMetPy, xminPDFMetPy, xmaxPDFMetPy) : 1.0;
   
   float BJet1PxProb = (! (*IsParFixed)[FPARAM_BJet1Px_TOPTOPLEPLEP]) ? getProb(hPDFBJetPx.get(), (par[FPARAM_BJet1Px_TOPTOPLEPLEP]-PxBJet1)/par[FPARAM_BJet1Px_TOPTOPLEPLEP], maxPDFBJetPx, xminPDFBJetPx, xmaxPDFBJetPx) : 1.0;
   float BJet1PyProb = (! (*IsParFixed)[FPARAM_BJet1Py_TOPTOPLEPLEP]) ? getProb(hPDFBJetPy.get(), (par[FPARAM_BJet1Py_TOPTOPLEPLEP]-PyBJet1)/par[FPARAM_BJet1Py_TOPTOPLEPLEP], maxPDFBJetPy, xminPDFBJetPy, xmaxPDFBJetPy) : 1.0;
   float BJet1PzProb = (! (*IsParFixed)[FPARAM_BJet1Pz_TOPTOPLEPLEP]) ? getProb(hPDFBJetPz.get(), (par[FPARAM_BJet1Pz_TOPTOPLEPLEP]-PzBJet1)/par[FPARAM_BJet1Pz_TOPTOPLEPLEP], maxPDFBJetPz, xminPDFBJetPz, xmaxPDFBJetPz) : 1.0;

   float BJet2PxProb = (! (*IsParFixed)[FPARAM_BJet2Px_TOPTOPLEPLEP]) ? getProb(hPDFBJetPx.get(), (par[FPARAM_BJet2Px_TOPTOPLEPLEP]-PxBJet2)/par[FPARAM_BJet2Px_TOPTOPLEPLEP], maxPDFBJetPx, xminPDFBJetPx, xmaxPDFBJetPx) : 1.0;
   float BJet2PyProb = (! (*IsParFixed)[FPARAM_BJet2Py_TOPTOPLEPLEP]) ? getProb(hPDFBJetPy.get(), (par[FPARAM_BJet2Py_TOPTOPLEPLEP]-PyBJet2)/par[FPARAM_BJet2Py_TOPTOPLEPLEP], maxPDFBJetPy, xminPDFBJetPy, xmaxPDFBJetPy) : 1.0;
   float BJet2PzProb = (! (*IsParFixed)[FPARAM_BJet2Pz_TOPTOPLEPLEP]) ? getProb(hPDFBJetPz.get(), (par[FPARAM_BJet2Pz_TOPTOPLEPLEP]-PzBJet2)/par[FPARAM_BJet2Pz_TOPTOPLEPLEP], maxPDFBJetPz, xminPDFBJetPz, xmaxPDFBJetPz) : 1.0;
   
   float Lepton1PxProb = 1.0;
   float Lepton1PyProb = 1.0;
   float Lepton1PzProb = 1.0;
   
   float Lepton2PxProb = 1.0;
   float Lepton2PyProb = 1.0;
   float Lepton2PzProb = 1.0;
   
   if (! (*IsParFixed)[FPARAM_Lepton1Px_TOPTOPLEPLEP]) Lepton1PxProb = (LabelLepton1 == 0) ? getProb(hPDFElecPx.get(), (par[FPARAM_Lepton1Px_TOPTOPLEPLEP]-PxLepton1)/par[FPARAM_Lepton1Px_TOPTOPLEPLEP], maxPDFElecPx, xminPDFElecPx, xmaxPDFElecPx) : getProb(hPDFMuonPx.get(), (par[FPARAM_Lepton1Px_TOPTOPLEPLEP]-PxLepton1)/par[FPARAM_Lepton1Px_TOPTOPLEPLEP], maxPDFMuonPx, xminPDFMuonPx, xmaxPDFMuonPx);
   if (! (*IsParFixed)[FPARAM_Lepton1Py_TOPTOPLEPLEP]) Lepton1PyProb = (LabelLepton1 == 0) ? getProb(hPDFElecPy.get(), (par[FPARAM_Lepton1Py_TOPTOPLEPLEP]-PyLepton1)/par[FPARAM_Lepton1Py_TOPTOPLEPLEP], maxPDFElecPy, xminPDFElecPy, xmaxPDFElecPy) : getProb(hPDFMuonPy.get(), (par[FPARAM_Lepton1Py_TOPTOPLEPLEP]-PyLepton1)/par[FPARAM_Lepton1Py_TOPTOPLEPLEP], maxPDFMuonPy, xminPDFMuonPy, xmaxPDFMuonPy);
   if (! (*IsParFixed)[FPARAM_Lepton1Pz_TOPTOPLEPLEP]) Lepton1PzProb = (LabelLepton1 == 0) ? getProb(hPDFElecPz.get(), (par[FPARAM_Lepton1Pz_TOPTOPLEPLEP]-PzLepton1)/par[FPARAM_Lepton1Pz_TOPTOPLEPLEP], maxPDFElecPz, xminPDFElecPz, xmaxPDFElecPz) : getProb(hPDFMuonPz.get(), (par[FPARAM_Lepton1Pz_TOPTOPLEPLEP]-PzLepton1)/par[FPARAM_Lepton1Pz_TOPTOPLEPLEP], maxPDFMuonPz, xminPDFMuonPz, xmaxPDFMuonPz);
   
   if (! (*IsParFixed)[FPARAM_Lepton2Px_TOPTOPLEPLEP]) Lepton2PxProb = (LabelLepton2 == 0) ? getProb(hPDFElecPx.get(), (par[FPARAM_Lepton2Px_TOPTOPLEPLEP]-PxLepton2)/par[FPARAM_Lepton2Px_TOPTOPLEPLEP], maxPDFElecPx, xminPDFElecPx, xmaxPDFElecPx) : getProb(hPDFMuonPx.get(), (par[FPARAM_Lepton2Px_TOPTOPLEPLEP]-PxLepton2)/par[FPARAM_Lepton2Px_TOPTOPLEPLEP], maxPDFMuonPx, xminPDFMuonPx, xmaxPDFMuonPx);
   if (! (*IsParFixed)[FPARAM_Lepton2Py_TOPTOPLEPLEP]) Lepton2PyProb = (LabelLepton2 == 0) ? getProb(hPDFElecPy.get(), (par[FPARAM_Lepton2Py_TOPTOPLEPLEP]-PyLepton2)/par[FPARAM_Lepton2Py_TOPTOPLEPLEP], maxPDFElecPy, xminPDFElecPy, xmaxPDFElecPy) : getProb(hPDFMuonPy.get(), (par[FPARAM_Lepton2Py_TOPTOPLEPLEP]-PyLepton2)/par[FPARAM_Lepton2Py_TOPTOPLEPLEP], maxPDFMuonPy, xminPDFMuonPy, xmaxPDFMuonPy);
   if (! (*IsParFixed)[FPARAM_Lepton2Pz_TOPTOPLEPLEP]) Lepton2PzProb = (LabelLepton2 == 0) ? getProb(hPDFElecPz.get(), (par[FPARAM_Lepton2Pz_TOPTOPLEPLEP]-PzLepton2)/par[FPARAM_Lepton2Pz_TOPTOPLEPLEP], maxPDFElecPz, xminPDFElecPz, xmaxPDFElecPz) : getProb(hPDFMuonPz.get(), (par[FPARAM_Lepton2Pz_TOPTOPLEPLEP]-PzLepton2)/par[FPARAM_Lepton2Pz_TOPTOPLEPLEP], maxPDFMuonPz, xminPDFMuonPz, xmaxPDFMuonPz);

   double minProb = 1E-20; // NLL = 92.1034
   
   chi2t[NLL_W1_TOPTOPLEPLEP] = (mW1Prob > minProb) ? mW1Prob : minProb;
   chi2t[NLL_W2_TOPTOPLEPLEP] = (mW2Prob > minProb) ? mW2Prob : minProb;
   chi2t[NLL_Top1_TOPTOPLEPLEP] = (mTop1Prob > minProb) ? mTop1Prob : minProb;
   chi2t[NLL_Top2_TOPTOPLEPLEP] = (mTop2Prob > minProb) ? mTop2Prob : minProb;
   chi2t[NLL_EtMissX_TOPTOPLEPLEP] = (MetPxProb > minProb) ? MetPxProb : minProb;
   chi2t[NLL_EtMissY_TOPTOPLEPLEP] = (MetPyProb > minProb) ? MetPyProb : minProb;
   chi2t[NLL_BJet1Px_TOPTOPLEPLEP] = (BJet1PxProb > minProb) ? BJet1PxProb : minProb;
   chi2t[NLL_BJet1Py_TOPTOPLEPLEP] = (BJet1PyProb > minProb) ? BJet1PyProb : minProb; 
   chi2t[NLL_BJet1Pz_TOPTOPLEPLEP] = (BJet1PzProb > minProb) ? BJet1PzProb : minProb;
   chi2t[NLL_BJet2Px_TOPTOPLEPLEP] = (BJet2PxProb > minProb) ? BJet2PxProb : minProb;
   chi2t[NLL_BJet2Py_TOPTOPLEPLEP] = (BJet2PyProb > minProb) ? BJet2PyProb : minProb;
   chi2t[NLL_BJet2Pz_TOPTOPLEPLEP] = (BJet2PzProb > minProb) ? BJet2PzProb : minProb;
   chi2t[NLL_Lepton1Px_TOPTOPLEPLEP] = (Lepton1PxProb > minProb) ? Lepton1PxProb : minProb;
   chi2t[NLL_Lepton1Py_TOPTOPLEPLEP] = (Lepton1PyProb > minProb) ? Lepton1PyProb : minProb;
   chi2t[NLL_Lepton1Pz_TOPTOPLEPLEP] = (Lepton1PzProb > minProb) ? Lepton1PzProb : minProb;
   chi2t[NLL_Lepton2Px_TOPTOPLEPLEP] = (Lepton2PxProb > minProb) ? Lepton2PxProb : minProb;
   chi2t[NLL_Lepton2Py_TOPTOPLEPLEP] = (Lepton2PyProb > minProb) ? Lepton2PyProb : minProb;
   chi2t[NLL_Lepton2Pz_TOPTOPLEPLEP] = (Lepton2PzProb > minProb) ? Lepton2PzProb : minProb;
   
//   chi2t[NLL_TopTopMass_TOPTOPLEPLEP] = (mTopTopProb > minProb) ? mTopTopProb : minProb;
   
   double lh = 1.0;
   
   lh *= chi2t[NLL_W1_TOPTOPLEPLEP];
   lh *= chi2t[NLL_W2_TOPTOPLEPLEP];
   lh *= chi2t[NLL_Top1_TOPTOPLEPLEP];
   lh *= chi2t[NLL_Top2_TOPTOPLEPLEP];

   /// Include additional terms only if is corresponding parameter is free
   if(! (*IsParFixed)[FPARAM_EtRealX_TOPTOPLEPLEP]) lh *= chi2t[NLL_EtMissX_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_EtRealY_TOPTOPLEPLEP]) lh *= chi2t[NLL_EtMissY_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_BJet1Px_TOPTOPLEPLEP]) lh *= chi2t[NLL_BJet1Px_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_BJet1Py_TOPTOPLEPLEP]) lh *= chi2t[NLL_BJet1Py_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_BJet1Pz_TOPTOPLEPLEP]) lh *= chi2t[NLL_BJet1Pz_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_BJet2Px_TOPTOPLEPLEP]) lh *= chi2t[NLL_BJet2Px_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_BJet2Py_TOPTOPLEPLEP]) lh *= chi2t[NLL_BJet2Py_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_BJet2Pz_TOPTOPLEPLEP]) lh *= chi2t[NLL_BJet2Pz_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_Lepton1Px_TOPTOPLEPLEP]) lh *= chi2t[NLL_Lepton1Px_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_Lepton1Py_TOPTOPLEPLEP]) lh *= chi2t[NLL_Lepton1Py_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_Lepton1Pz_TOPTOPLEPLEP]) lh *= chi2t[NLL_Lepton1Pz_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_Lepton2Px_TOPTOPLEPLEP]) lh *= chi2t[NLL_Lepton2Px_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_Lepton2Py_TOPTOPLEPLEP]) lh *= chi2t[NLL_Lepton2Py_TOPTOPLEPLEP];
   if(! (*IsParFixed)[FPARAM_Lepton2Pz_TOPTOPLEPLEP]) lh *= chi2t[NLL_Lepton2Pz_TOPTOPLEPLEP];
   
//   if( AddTopTopMassToNLL_ ) lh *= chi2t[NLL_TopTopMass_TOPTOPLEPLEP];

   val += -2.*log(lh);

   return val;
}

/// Calculate output variables for a given permutation
void KINFIT::TopTopLepLep::calcVar(int iPerm)
{
   if( chi_[iPerm] > 10E+9 ) return;

   int idxElec1 = TopTopLepLep_Electron1Idx[iPerm];
   int idxMuon1 = TopTopLepLep_Muon1Idx[iPerm];

   int idxElec2 = TopTopLepLep_Electron2Idx[iPerm];
   int idxMuon2 = TopTopLepLep_Muon2Idx[iPerm];
   
   int idxBJet1 = TopTopLepLep_BJet1Idx[iPerm];
   int idxBJet2 = TopTopLepLep_BJet2Idx[iPerm];
   
   int idxPhoton = TopTopLepLep_PhotonIdx[iPerm];
   int photonOrigin = PhotonOrigin_[iPerm];
   
   float Lepton1Px = 0;
   float Lepton1Py = 0;
   float Lepton1Pz = 0;
   float Lepton1E = 0;

   float Lepton2Px = 0;
   float Lepton2Py = 0;
   float Lepton2Pz = 0;
   float Lepton2E = 0;

   float BJet1Px = 0;
   float BJet1Py = 0;
   float BJet1Pz = 0;
   float BJet1E = 0;

   float BJet2Px = 0;
   float BJet2Py = 0;
   float BJet2Pz = 0;
   float BJet2E = 0;
   
   float RadPhotonPx = 0;
   float RadPhotonPy = 0;
   float RadPhotonPz = 0;
   float RadPhotonE = 0;
   
   if( idxPhoton >= 0 )
     {
	RadPhotonPx = PhotonPx[idxPhoton];
	RadPhotonPy = PhotonPy[idxPhoton];
	RadPhotonPz = PhotonPz[idxPhoton];
	RadPhotonE = PhotonE[idxPhoton];
     }
   
   if( idxElec1 >= 0 )
     {
	Lepton1Px = ElectronPx[idxElec1];
	Lepton1Py = ElectronPy[idxElec1];
	Lepton1Pz = ElectronPz[idxElec1];
	Lepton1E = ElectronE[idxElec1];
     }   
   else if( idxMuon1 >= 0 )
     {
	Lepton1Px = MuonPx[idxMuon1];
	Lepton1Py = MuonPy[idxMuon1];
	Lepton1Pz = MuonPz[idxMuon1];
	Lepton1E = MuonE[idxMuon1];
     }

   if( photonOrigin == PHOTON_FROM_LEPTON1_TOPTOPLEPLEP && (idxElec1 >= 0 || idxMuon1 >= 0) )
     {
	Lepton1Px += RadPhotonPx;
	Lepton1Py += RadPhotonPy;
	Lepton1Pz += RadPhotonPz;
	Lepton1E += RadPhotonE;
     }

   if( idxElec2 >= 0 )
     {
	Lepton2Px = ElectronPx[idxElec2];
	Lepton2Py = ElectronPy[idxElec2];
	Lepton2Pz = ElectronPz[idxElec2];
	Lepton2E = ElectronE[idxElec2];
     }   
   else if( idxMuon2 >= 0 )
     {
	Lepton2Px = MuonPx[idxMuon2];
	Lepton2Py = MuonPy[idxMuon2];
	Lepton2Pz = MuonPz[idxMuon2];
	Lepton2E = MuonE[idxMuon2];
     }   

   if( photonOrigin == PHOTON_FROM_LEPTON2_TOPTOPLEPLEP && (idxElec2 >= 0 || idxMuon2 >= 0) )
     {
	Lepton2Px += RadPhotonPx;
	Lepton2Py += RadPhotonPy;
	Lepton2Pz += RadPhotonPz;
	Lepton2E += RadPhotonE;
     }
   
   if( idxBJet1 >= 0 )
     {
	BJet1Px = BJetPx[idxBJet1];
	BJet1Py = BJetPy[idxBJet1];
	BJet1Pz = BJetPz[idxBJet1];
	BJet1E = BJetE[idxBJet1];
	
	if( photonOrigin == PHOTON_FROM_BJET1_TOPTOPLEPLEP )
	  {
	     BJet1Px += RadPhotonPx;
	     BJet1Py += RadPhotonPy;
	     BJet1Pz += RadPhotonPz;
	     BJet1E += RadPhotonE;
	  }
     }
   if( idxBJet2 >= 0 )
     {
	BJet2Px = BJetPx[idxBJet2];
	BJet2Py = BJetPy[idxBJet2];
	BJet2Pz = BJetPz[idxBJet2];
	BJet2E = BJetE[idxBJet2];
	
	if( photonOrigin == PHOTON_FROM_BJET2_TOPTOPLEPLEP )
	  {
	     BJet2Px += RadPhotonPx;
	     BJet2Py += RadPhotonPy;
	     BJet2Pz += RadPhotonPz;
	     BJet2E += RadPhotonE;
	  }
     }

   float nu1Px = nuPx_[iPerm][0];
   float nu1Py = nuPy_[iPerm][0];
   float nu1Pz = nuPz_[iPerm][0];

   float nu2Px = nuPx_[iPerm][1];
   float nu2Py = nuPy_[iPerm][1];
   float nu2Pz = nuPz_[iPerm][1];
   
   float nu1E = sqrt(nu1Px*nu1Px+nu1Py*nu1Py+nu1Pz*nu1Pz);
   float nu2E = sqrt(nu2Px*nu2Px+nu2Py*nu2Py+nu2Pz*nu2Pz);

   /// First W boson
   float W1E = nu1E+Lepton1E;
   float W1Px = nu1Px+Lepton1Px;
   float W1Py = nu1Py+Lepton1Py;
   float W1Pz = nu1Pz+Lepton1Pz;
   
   if( photonOrigin == PHOTON_FROM_W1_TOPTOPLEPLEP ||
       photonOrigin == PHOTON_FROM_W1_COMB_TOPTOPLEPLEP )
     {
	W1Px += RadPhotonPx;
	W1Py += RadPhotonPy;
	W1Pz += RadPhotonPz;
	W1E += RadPhotonE;
     }

   float W1Pt = sqrt(W1Px*W1Px+W1Py*W1Py);
   float W1Eta = getEta(W1Pt, W1Pz);
   float W1Phi = atan2(W1Py, W1Px);
   
   /// First top quark
   float top1E = W1E+BJet1E;
   float top1Px = W1Px+BJet1Px;
   float top1Py = W1Py+BJet1Py;
   float top1Pz = W1Pz+BJet1Pz;
   
   if( photonOrigin == PHOTON_FROM_TOP1_TOPTOPLEPLEP ||
       photonOrigin == PHOTON_FROM_TOP1_COMB_TOPTOPLEPLEP )
     {
	top1Px += RadPhotonPx;
	top1Py += RadPhotonPy;
	top1Pz += RadPhotonPz;
	top1E += RadPhotonE;
     }
   
   float top1Pt = sqrt(top1Px*top1Px+top1Py*top1Py);
   float top1Eta = getEta(top1Pt, top1Pz);
   float top1Phi = atan2(top1Py, top1Px);

   /// Second W boson
   float W2E = nu2E+Lepton2E;
   float W2Px = nu2Px+Lepton2Px;
   float W2Py = nu2Py+Lepton2Py;
   float W2Pz = nu2Pz+Lepton2Pz;

   if( photonOrigin == PHOTON_FROM_W2_TOPTOPLEPLEP ||
       photonOrigin == PHOTON_FROM_W2_COMB_TOPTOPLEPLEP )
     {
	W2Px += RadPhotonPx;
	W2Py += RadPhotonPy;
	W2Pz += RadPhotonPz;
	W2E += RadPhotonE;
     }
   
   float W2Pt = sqrt(W2Px*W2Px+W2Py*W2Py);
   float W2Eta = getEta(W2Pt, W2Pz);
   float W2Phi = atan2(W2Py, W2Px);

   /// Second top quark
   float top2E = W2E+BJet2E;
   float top2Px = W2Px+BJet2Px;
   float top2Py = W2Py+BJet2Py;
   float top2Pz = W2Pz+BJet2Pz;

   if( photonOrigin == PHOTON_FROM_TOP2_TOPTOPLEPLEP ||
       photonOrigin == PHOTON_FROM_TOP2_COMB_TOPTOPLEPLEP )
     {
	top2Px += RadPhotonPx;
	top2Py += RadPhotonPy;
	top2Pz += RadPhotonPz;
	top2E += RadPhotonE;
     }
   
   float top2Pt = sqrt(top2Px*top2Px+top2Py*top2Py);
   float top2Eta = getEta(top2Pt, top2Pz);
   float top2Phi = atan2(top2Py, top2Px);
   
   /// Di-top system
   float toptopE = top1E+top2E;
   float toptopPx = top1Px+top2Px;
   float toptopPy = top1Py+top2Py;
   float toptopPz = top1Pz+top2Pz;

   /// Fill the outputs
   
   drTopTop_[iPerm] = getDeltaR(top1Eta, top1Phi, top2Eta, top2Phi);
   detaTopTop_[iPerm] = top1Eta-top2Eta;
   dphiTopTop_[iPerm] = getDeltaPhi(top1Phi, top2Phi);

   mTopTop_[iPerm] = sqrt(toptopE*toptopE-toptopPx*toptopPx-toptopPy*toptopPy-toptopPz*toptopPz);
   ptTopTop_[iPerm] = sqrt(toptopPx*toptopPx+toptopPy*toptopPy);
   pTopTop_[iPerm] = sqrt(toptopPx*toptopPx+toptopPy*toptopPy+toptopPz*toptopPz);
   etaTopTop_[iPerm] = getEta(ptTopTop_[iPerm], toptopPz);
   phiTopTop_[iPerm] = atan2(toptopPy, toptopPx);
   
   TopMass_[iPerm][0] = sqrt(top1E*top1E-top1Px*top1Px-top1Py*top1Py-top1Pz*top1Pz);
   TopMass_[iPerm][1] = sqrt(top2E*top2E-top2Px*top2Px-top2Py*top2Py-top2Pz*top2Pz);
   TopPt_[iPerm][0] = sqrt(top1Px*top1Px+top1Py*top1Py);
   TopPt_[iPerm][1] = sqrt(top2Px*top2Px+top2Py*top2Py);
   TopP_[iPerm][0] = sqrt(top1Px*top1Px+top1Py*top1Py+top1Pz*top1Pz);
   TopP_[iPerm][1] = sqrt(top2Px*top2Px+top2Py*top2Py+top2Pz*top2Pz);
   TopEta_[iPerm][0] = getEta(TopPt_[iPerm][0],top1Pz);
   TopEta_[iPerm][1] = getEta(TopPt_[iPerm][1],top2Pz);
   TopPhi_[iPerm][0] = top1Phi;
   TopPhi_[iPerm][1] = top2Phi;
   TopE_[iPerm][0] = top1E;
   TopE_[iPerm][1] = top2E;
   TopPx_[iPerm][0] = top1Px;
   TopPx_[iPerm][1] = top2Px;
   TopPy_[iPerm][0] = top1Py;
   TopPy_[iPerm][1] = top2Py;
   TopPz_[iPerm][0] = top1Pz;
   TopPz_[iPerm][1] = top2Pz;

   WPt_[iPerm][0] = sqrt(W1Px*W1Px+W1Py*W1Py);
   WPt_[iPerm][1] = sqrt(W2Px*W2Px+W2Py*W2Py);
   WP_[iPerm][0] = sqrt(W1Px*W1Px+W1Py*W1Py+W1Pz*W1Pz);
   WP_[iPerm][1] = sqrt(W2Px*W2Px+W2Py*W2Py+W2Pz*W2Pz);
   WEta_[iPerm][0] = getEta(WPt_[iPerm][0],W1Pz);
   WEta_[iPerm][1] = getEta(WPt_[iPerm][1],W2Pz);
   WPhi_[iPerm][0] = atan2(W1Py, W1Px);
   WPhi_[iPerm][1] = atan2(W2Py, W2Px);
   WE_[iPerm][0] = W1E;
   WE_[iPerm][1] = W2E;
   WPx_[iPerm][0] = W1Px;
   WPx_[iPerm][1] = W2Px;
   WPy_[iPerm][0] = W1Py;
   WPy_[iPerm][1] = W2Py;
   WPz_[iPerm][0] = W1Pz;
   WPz_[iPerm][1] = W2Pz;
}

/// FCN definition
void KINFIT::TopTopLepLep::fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
//   std::string chi2tNames[NLL_N_TOPTOPLEPLEP] = {"W1", "W2", "Top1", "Top2", "EtMissX", "EtMissY", "BJet1Px", "BJet1Py", "BJet1Pz", "BJet2Px", "BJet2Py", "BJet2Pz", "Lepton1Px", "Lepton1Py", "Lepton1Pz", "Lepton2Px", "Lepton2Py", "Lepton2Pz", "RacPT", "mWPT", "mTopPT", "TopTopMass"};
   std::string chi2tNames[NLL_N_TOPTOPLEPLEP] = {"W1", "W2", "Top1", "Top2", "EtMissX", "EtMissY", "BJet1Px", "BJet1Py", "BJet1Pz", "BJet2Px", "BJet2Py", "BJet2Pz", "Lepton1Px", "Lepton1Py", "Lepton1Pz", "Lepton2Px", "Lepton2Py", "Lepton2Pz", "RacPT", "mWPT", "mTopPT"};
   std::vector<double> chi2t(NLL_N_TOPTOPLEPLEP);
   
   double lh = func(*PxLepton1, *PyLepton1, *PzLepton1, *ELepton1, *LabelLepton1,
		    *PxLepton2, *PyLepton2, *PzLepton2, *ELepton2, *LabelLepton2,
		    *PxBJet1, *PyBJet1, *PzBJet1, *EBJet1,
		    *PxBJet2, *PyBJet2, *PzBJet2, *EBJet2,
		    *PxPhoton, *PyPhoton, *PzPhoton, *EPhoton,
		    *PhotonOrigin,
		    chi2t, par);

   /// converged
   if( iflag == 3 )
     {
	*CHISQ = lh;

	for( int ip=0;ip<FPARAM_N;ip++ ) (*FitParam).push_back(par[ip]);
	for( int it=0;it<NLL_N_TOPTOPLEPLEP;it++ ) {(*ChiTerm).push_back(chi2t[it]); (*ChiTermName).push_back(chi2tNames[it]);}
     }
   
   f = lh;
}

/// Run the fit
void KINFIT::TopTopLepLep::fit(double *par, std::vector<FRESULT> &vpp)
{
   (*FitParam).clear();
   (*ChiTerm).clear();
   (*ChiTermName).clear();
   
   *CHISQ = 10E+10;
   
   double perr[FPARAM_N];
   
   TMinuit *gMinuit = new TMinuit(FPARAM_N);
   gMinuit->SetFCN(fcn);
   gMinuit->SetPrintLevel(-1);
   
   Double_t arglist[10];
   Int_t ierflg = 0;
   
   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
   
   arglist[0] = 5000; /// number of steps
   arglist[1] = 0.1; /// tolerance

   /// Set dummy ranges for fixed parameters
   (*ParMin)[FPARAM_BJet1E_TOPTOPLEPLEP] = par[FPARAM_BJet1E_TOPTOPLEPLEP]-1000.; (*ParMax)[FPARAM_BJet1E_TOPTOPLEPLEP] = par[FPARAM_BJet1E_TOPTOPLEPLEP]+1000.;
   (*ParMin)[FPARAM_BJet2E_TOPTOPLEPLEP] = par[FPARAM_BJet2E_TOPTOPLEPLEP]-1000.; (*ParMax)[FPARAM_BJet2E_TOPTOPLEPLEP] = par[FPARAM_BJet2E_TOPTOPLEPLEP]+1000.;
   (*ParMin)[FPARAM_Lepton1E_TOPTOPLEPLEP] = par[FPARAM_Lepton1E_TOPTOPLEPLEP]-1000.; (*ParMax)[FPARAM_Lepton1E_TOPTOPLEPLEP] = par[FPARAM_Lepton1E_TOPTOPLEPLEP]+1000.;
   (*ParMin)[FPARAM_Lepton2E_TOPTOPLEPLEP] = par[FPARAM_Lepton2E_TOPTOPLEPLEP]-1000.; (*ParMax)[FPARAM_Lepton2E_TOPTOPLEPLEP] = par[FPARAM_Lepton2E_TOPTOPLEPLEP]+1000.;
   (*ParMin)[FPARAM_PhotonE_TOPTOPLEPLEP] = par[FPARAM_PhotonE_TOPTOPLEPLEP]-1000.; (*ParMax)[FPARAM_PhotonE_TOPTOPLEPLEP] = par[FPARAM_PhotonE_TOPTOPLEPLEP]+1000.;
      
   gMinuit->mnparm(FPARAM_Etx_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_Etx_TOPTOPLEPLEP], par[FPARAM_Etx_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_Etx_TOPTOPLEPLEP], (*ParMax)[FPARAM_Etx_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_Ety_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_Ety_TOPTOPLEPLEP], par[FPARAM_Ety_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_Ety_TOPTOPLEPLEP], (*ParMax)[FPARAM_Ety_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_Sign1_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_Sign1_TOPTOPLEPLEP], par[FPARAM_Sign1_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_Sign1_TOPTOPLEPLEP], (*ParMax)[FPARAM_Sign1_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_Sign2_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_Sign2_TOPTOPLEPLEP], par[FPARAM_Sign2_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_Sign2_TOPTOPLEPLEP], (*ParMax)[FPARAM_Sign2_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_EtRealX_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_EtRealX_TOPTOPLEPLEP], par[FPARAM_EtRealX_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_EtRealX_TOPTOPLEPLEP], (*ParMax)[FPARAM_EtRealX_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_EtRealY_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_EtRealY_TOPTOPLEPLEP], par[FPARAM_EtRealY_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_EtRealY_TOPTOPLEPLEP], (*ParMax)[FPARAM_EtRealY_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_mW1_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_mW1_TOPTOPLEPLEP], par[FPARAM_mW1_TOPTOPLEPLEP], 1E-2, (*ParMin)[FPARAM_mW1_TOPTOPLEPLEP], (*ParMax)[FPARAM_mW1_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_mW2_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_mW2_TOPTOPLEPLEP], par[FPARAM_mW2_TOPTOPLEPLEP], 1E-2, (*ParMin)[FPARAM_mW2_TOPTOPLEPLEP], (*ParMax)[FPARAM_mW2_TOPTOPLEPLEP], ierflg);

   gMinuit->mnparm(FPARAM_BJet1Px_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_BJet1Px_TOPTOPLEPLEP], par[FPARAM_BJet1Px_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_BJet1Px_TOPTOPLEPLEP], (*ParMax)[FPARAM_BJet1Px_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_BJet1Py_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_BJet1Py_TOPTOPLEPLEP], par[FPARAM_BJet1Py_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_BJet1Py_TOPTOPLEPLEP], (*ParMax)[FPARAM_BJet1Py_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_BJet1Pz_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_BJet1Pz_TOPTOPLEPLEP], par[FPARAM_BJet1Pz_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_BJet1Pz_TOPTOPLEPLEP], (*ParMax)[FPARAM_BJet1Pz_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_BJet1E_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_BJet1E_TOPTOPLEPLEP], par[FPARAM_BJet1E_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_BJet1E_TOPTOPLEPLEP], (*ParMax)[FPARAM_BJet1E_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_BJet2Px_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_BJet2Px_TOPTOPLEPLEP], par[FPARAM_BJet2Px_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_BJet2Px_TOPTOPLEPLEP], (*ParMax)[FPARAM_BJet2Px_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_BJet2Py_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_BJet2Py_TOPTOPLEPLEP], par[FPARAM_BJet2Py_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_BJet2Py_TOPTOPLEPLEP], (*ParMax)[FPARAM_BJet2Py_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_BJet2Pz_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_BJet2Pz_TOPTOPLEPLEP], par[FPARAM_BJet2Pz_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_BJet2Pz_TOPTOPLEPLEP], (*ParMax)[FPARAM_BJet2Pz_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_BJet2E_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_BJet2E_TOPTOPLEPLEP], par[FPARAM_BJet2E_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_BJet2E_TOPTOPLEPLEP], (*ParMax)[FPARAM_BJet2E_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_Lepton1Px_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_Lepton1Px_TOPTOPLEPLEP], par[FPARAM_Lepton1Px_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_Lepton1Px_TOPTOPLEPLEP], (*ParMax)[FPARAM_Lepton1Px_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_Lepton1Py_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_Lepton1Py_TOPTOPLEPLEP], par[FPARAM_Lepton1Py_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_Lepton1Py_TOPTOPLEPLEP], (*ParMax)[FPARAM_Lepton1Py_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_Lepton1Pz_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_Lepton1Pz_TOPTOPLEPLEP], par[FPARAM_Lepton1Pz_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_Lepton1Pz_TOPTOPLEPLEP], (*ParMax)[FPARAM_Lepton1Pz_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_Lepton1E_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_Lepton1E_TOPTOPLEPLEP], par[FPARAM_Lepton1E_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_Lepton1E_TOPTOPLEPLEP], (*ParMax)[FPARAM_Lepton1E_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_Lepton2Px_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_Lepton2Px_TOPTOPLEPLEP], par[FPARAM_Lepton2Px_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_Lepton2Px_TOPTOPLEPLEP], (*ParMax)[FPARAM_Lepton2Px_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_Lepton2Py_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_Lepton2Py_TOPTOPLEPLEP], par[FPARAM_Lepton2Py_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_Lepton2Py_TOPTOPLEPLEP], (*ParMax)[FPARAM_Lepton2Py_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_Lepton2Pz_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_Lepton2Pz_TOPTOPLEPLEP], par[FPARAM_Lepton2Pz_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_Lepton2Pz_TOPTOPLEPLEP], (*ParMax)[FPARAM_Lepton2Pz_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_Lepton2E_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_Lepton2E_TOPTOPLEPLEP], par[FPARAM_Lepton2E_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_Lepton2E_TOPTOPLEPLEP], (*ParMax)[FPARAM_Lepton2E_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_PhotonPx_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_PhotonPx_TOPTOPLEPLEP], par[FPARAM_PhotonPx_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_PhotonPx_TOPTOPLEPLEP], (*ParMax)[FPARAM_PhotonPx_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_PhotonPy_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_PhotonPy_TOPTOPLEPLEP], par[FPARAM_PhotonPy_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_PhotonPy_TOPTOPLEPLEP], (*ParMax)[FPARAM_PhotonPy_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_PhotonPz_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_PhotonPz_TOPTOPLEPLEP], par[FPARAM_PhotonPz_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_PhotonPz_TOPTOPLEPLEP], (*ParMax)[FPARAM_PhotonPz_TOPTOPLEPLEP], ierflg);
   gMinuit->mnparm(FPARAM_PhotonE_TOPTOPLEPLEP, FPARAM_NAME[FPARAM_PhotonE_TOPTOPLEPLEP], par[FPARAM_PhotonE_TOPTOPLEPLEP], 1E-1, (*ParMin)[FPARAM_PhotonE_TOPTOPLEPLEP], (*ParMax)[FPARAM_PhotonE_TOPTOPLEPLEP], ierflg);
   
   for( int i=0;i<FPARAM_N;i++ ) if( (*IsParFixed)[i] ) gMinuit->FixParameter(i);

   /// Alternative option
//   gMinuit->mnexcm("MINI", arglist ,2, ierflg);
//   gMinuit->mnexcm("MINI", arglist ,2, ierflg);
//   gMinuit->mnexcm("MINI", arglist ,2, ierflg);
   
   gMinuit->mnexcm("MINI", arglist ,2, ierflg);
   gMinuit->mnexcm("MINOS", arglist ,2, ierflg);
   
   arglist[0] = 3;
   gMinuit->mnexcm("CALL FCN", arglist, 1, ierflg);

   FRESULT fres = {};
   
   /// Get the resultant fit parameters
   for( int ip=0;ip<FPARAM_N;ip++ )
     {	
	gMinuit->GetParameter(ip, par[ip], perr[ip]);
	fres.par[ip] = par[ip];
     }   

   for( int ic=0;ic<(*ChiTerm).size();ic++ )
     {	
	fres.chiTerm[ic] = (*ChiTerm)[ic];
	fres.chiTermName[ic] = (*ChiTermName)[ic];
     }   

   fres.NTERM = (*ChiTerm).size();
   fres.PzNu1 = (*PzNu1);
   fres.PzNu2 = (*PzNu2);
   fres.PzNuSum = fabs(*PzNu1) + fabs(*PzNu2);
   
   float PtNu1 = sqrt((*PxNu1)*(*PxNu1) + (*PyNu1)*(*PyNu1));
   float PNu1 = sqrt((PtNu1)*(PtNu1) + (*PzNu1)*(*PzNu1));
   float PtNu2 = sqrt((*PxNu2)*(*PxNu2) + (*PyNu2)*(*PyNu2));
   float PNu2 = sqrt((PtNu2)*(PtNu2) + (*PzNu2)*(*PzNu2));
   
   fres.PtNuSum = PtNu1 + PtNu2;
   fres.PNuSum = PNu1 + PNu2;   
   
   fres.TopTopMass = (*TopTopMass);
   fres.TopTopPz = (*TopTopPz);
   fres.TopTopPzAbs = (*TopTopPzAbs);
   fres.lh = *CHISQ;
   
   vpp.push_back(fres);
   
   delete gMinuit;
}

/// First-layer minimization (generic)
void KINFIT::TopTopLepLep::calcNuGrid(std::vector<FRESULT> &vp)
{
   double par[FPARAM_N];
   for(int p=0;p<FPARAM_N;p++) par[p] = 0.;
   
   std::vector<double> chi2t(NLL_N_TOPTOPLEPLEP);

   bool doToys = bool(NToy_ > 1);
   
   for(int it=0;it<NToy_;it++)
     {
	if( vp.size() >= NGrid_ ) break;
	
	double thres = -1E-6;

	do
	  {	     
	     double gv = (doToys) ? getProbGaus(hPDFBJetPx.get(), maxPDFBJetPx, meanPDFBJetPx, sigmaPDFBJetPx, rnd, NBJetPxRMS_) : 0.;
	     par[FPARAM_BJet1Px_TOPTOPLEPLEP] = *PxBJet1/(1.-gv);
	     gv = (doToys) ? getProbGaus(hPDFBJetPy.get(), maxPDFBJetPy, meanPDFBJetPy, sigmaPDFBJetPy, rnd, NBJetPyRMS_) : 0.;
	     par[FPARAM_BJet1Py_TOPTOPLEPLEP] = *PyBJet1/(1.-gv);
	     gv = (doToys) ? getProbGaus(hPDFBJetPz.get(), maxPDFBJetPz, meanPDFBJetPz, sigmaPDFBJetPz, rnd, NBJetPzRMS_) : 0.;
	     par[FPARAM_BJet1Pz_TOPTOPLEPLEP] = *PzBJet1/(1.-gv);
	     par[FPARAM_BJet1E_TOPTOPLEPLEP] = sqrt((*MassBJet1)*(*MassBJet1) + par[FPARAM_BJet1Px_TOPTOPLEPLEP]*par[FPARAM_BJet1Px_TOPTOPLEPLEP] + par[FPARAM_BJet1Py_TOPTOPLEPLEP]*par[FPARAM_BJet1Py_TOPTOPLEPLEP] + par[FPARAM_BJet1Pz_TOPTOPLEPLEP]*par[FPARAM_BJet1Pz_TOPTOPLEPLEP]);
	  } while( (par[FPARAM_BJet1E_TOPTOPLEPLEP]*par[FPARAM_BJet1E_TOPTOPLEPLEP]-par[FPARAM_BJet1Px_TOPTOPLEPLEP]*par[FPARAM_BJet1Px_TOPTOPLEPLEP]-par[FPARAM_BJet1Py_TOPTOPLEPLEP]*par[FPARAM_BJet1Py_TOPTOPLEPLEP]-par[FPARAM_BJet1Pz_TOPTOPLEPLEP]*par[FPARAM_BJet1Pz_TOPTOPLEPLEP]) < thres );

	do
	  {	     
	     double gv = (doToys) ? getProbGaus(hPDFBJetPx.get(), maxPDFBJetPx, meanPDFBJetPx, sigmaPDFBJetPx, rnd, NBJetPxRMS_) : 0.;
	     par[FPARAM_BJet2Px_TOPTOPLEPLEP] = *PxBJet2/(1.-gv);
	     gv = (doToys) ? getProbGaus(hPDFBJetPy.get(), maxPDFBJetPy, meanPDFBJetPy, sigmaPDFBJetPy, rnd, NBJetPyRMS_) : 0.;
	     par[FPARAM_BJet2Py_TOPTOPLEPLEP] = *PyBJet2/(1.-gv);
	     gv = (doToys) ? getProbGaus(hPDFBJetPz.get(), maxPDFBJetPz, meanPDFBJetPz, sigmaPDFBJetPz, rnd, NBJetPzRMS_) : 0.;
	     par[FPARAM_BJet2Pz_TOPTOPLEPLEP] = *PzBJet2/(1.-gv);
	     par[FPARAM_BJet2E_TOPTOPLEPLEP] = sqrt((*MassBJet2)*(*MassBJet2) + par[FPARAM_BJet2Px_TOPTOPLEPLEP]*par[FPARAM_BJet2Px_TOPTOPLEPLEP] + par[FPARAM_BJet2Py_TOPTOPLEPLEP]*par[FPARAM_BJet2Py_TOPTOPLEPLEP] + par[FPARAM_BJet2Pz_TOPTOPLEPLEP]*par[FPARAM_BJet2Pz_TOPTOPLEPLEP]);
	  } while( (par[FPARAM_BJet2E_TOPTOPLEPLEP]*par[FPARAM_BJet2E_TOPTOPLEPLEP]-par[FPARAM_BJet2Px_TOPTOPLEPLEP]*par[FPARAM_BJet2Px_TOPTOPLEPLEP]-par[FPARAM_BJet2Py_TOPTOPLEPLEP]*par[FPARAM_BJet2Py_TOPTOPLEPLEP]-par[FPARAM_BJet2Pz_TOPTOPLEPLEP]*par[FPARAM_BJet2Pz_TOPTOPLEPLEP]) < thres );
	
	if( IncludePhotons_ )
	  {
	     do
	       {
		  double gv = (doToys) ? getProbGaus(hPDFPhotonPx.get(), maxPDFPhotonPx, meanPDFPhotonPx, sigmaPDFPhotonPx, rnd, NPhotonPxRMS_) : 0.;
		  par[FPARAM_PhotonPx_TOPTOPLEPLEP] = *PxPhoton/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFPhotonPy.get(), maxPDFPhotonPy, meanPDFPhotonPy, sigmaPDFPhotonPy, rnd, NPhotonPyRMS_) : 0.;
		  par[FPARAM_PhotonPy_TOPTOPLEPLEP] = *PyPhoton/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFPhotonPz.get(), maxPDFPhotonPz, meanPDFPhotonPz, sigmaPDFPhotonPz, rnd, NPhotonPzRMS_) : 0.;
		  par[FPARAM_PhotonPz_TOPTOPLEPLEP] = *PzPhoton/(1.-gv);
		  par[FPARAM_PhotonE_TOPTOPLEPLEP] = sqrt(par[FPARAM_PhotonPx_TOPTOPLEPLEP]*par[FPARAM_PhotonPx_TOPTOPLEPLEP] + par[FPARAM_PhotonPy_TOPTOPLEPLEP]*par[FPARAM_PhotonPy_TOPTOPLEPLEP] + par[FPARAM_PhotonPz_TOPTOPLEPLEP]*par[FPARAM_PhotonPz_TOPTOPLEPLEP]);
	       } while( (par[FPARAM_PhotonE_TOPTOPLEPLEP]*par[FPARAM_PhotonE_TOPTOPLEPLEP]-par[FPARAM_PhotonPx_TOPTOPLEPLEP]*par[FPARAM_PhotonPx_TOPTOPLEPLEP]-par[FPARAM_PhotonPy_TOPTOPLEPLEP]*par[FPARAM_PhotonPy_TOPTOPLEPLEP]-par[FPARAM_PhotonPz_TOPTOPLEPLEP]*par[FPARAM_PhotonPz_TOPTOPLEPLEP]) < thres );
	  }	

	  {
	     double gv = (doToys) ? getProbGaus(hPDFMetPx.get(), maxPDFMetPx, meanPDFMetPx, sigmaPDFMetPx, rnd, NMetRMS_) : 0.;
	     par[FPARAM_EtRealX_TOPTOPLEPLEP] = *EtMissX/(1.-gv);
	     gv = (doToys) ? getProbGaus(hPDFMetPy.get(), maxPDFMetPy, meanPDFMetPy, sigmaPDFMetPy, rnd, NMetRMS_) : 0.;
	     par[FPARAM_EtRealY_TOPTOPLEPLEP] = *EtMissY/(1.-gv);
	  }
	
	par[FPARAM_mW1_TOPTOPLEPLEP] = (doToys) ? getProbGaus(hPDFTopWMass.get(), maxPDFTopWMass, meanPDFTopWMass, sigmaPDFTopWMass, rnd, 2) : WMass1;
	par[FPARAM_mW2_TOPTOPLEPLEP] = (doToys) ? getProbGaus(hPDFTopWMass.get(), maxPDFTopWMass, meanPDFTopWMass, sigmaPDFTopWMass, rnd, 2) : WMass2;

	if( *LabelLepton1 == 0 )
	  {
	     do
	       {
		  double gv = (doToys) ? getProbGaus(hPDFElecPx.get(), maxPDFElecPx, meanPDFElecPx, sigmaPDFElecPx, rnd, NElecPxRMS_) : 0.;
		  par[FPARAM_Lepton1Px_TOPTOPLEPLEP] = *PxLepton1/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFElecPy.get(), maxPDFElecPy, meanPDFElecPy, sigmaPDFElecPy, rnd, NElecPyRMS_) : 0.;
		  par[FPARAM_Lepton1Py_TOPTOPLEPLEP] = *PyLepton1/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFElecPz.get(), maxPDFElecPz, meanPDFElecPz, sigmaPDFElecPz, rnd, NElecPzRMS_) : 0.;
		  par[FPARAM_Lepton1Pz_TOPTOPLEPLEP] = *PzLepton1/(1.-gv);
		  par[FPARAM_Lepton1E_TOPTOPLEPLEP] = sqrt((*MassLepton1)*(*MassLepton1) + par[FPARAM_Lepton1Px_TOPTOPLEPLEP]*par[FPARAM_Lepton1Px_TOPTOPLEPLEP] + par[FPARAM_Lepton1Py_TOPTOPLEPLEP]*par[FPARAM_Lepton1Py_TOPTOPLEPLEP] + par[FPARAM_Lepton1Pz_TOPTOPLEPLEP]*par[FPARAM_Lepton1Pz_TOPTOPLEPLEP]);
	       } while( (par[FPARAM_Lepton1E_TOPTOPLEPLEP]*par[FPARAM_Lepton1E_TOPTOPLEPLEP]-par[FPARAM_Lepton1Px_TOPTOPLEPLEP]*par[FPARAM_Lepton1Px_TOPTOPLEPLEP]-par[FPARAM_Lepton1Py_TOPTOPLEPLEP]*par[FPARAM_Lepton1Py_TOPTOPLEPLEP]-par[FPARAM_Lepton1Pz_TOPTOPLEPLEP]*par[FPARAM_Lepton1Pz_TOPTOPLEPLEP]) < thres );
	  }
	else
	  {
	     do
	       {
		  double gv = (doToys) ? getProbGaus(hPDFMuonPx.get(), maxPDFMuonPx, meanPDFMuonPx, sigmaPDFMuonPx, rnd, NMuonPxRMS_) : 0.;
		  par[FPARAM_Lepton1Px_TOPTOPLEPLEP] = *PxLepton1/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFMuonPy.get(), maxPDFMuonPy, meanPDFMuonPy, sigmaPDFMuonPy, rnd, NMuonPyRMS_) : 0.;
		  par[FPARAM_Lepton1Py_TOPTOPLEPLEP] = *PyLepton1/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFMuonPz.get(), maxPDFMuonPz, meanPDFMuonPz, sigmaPDFMuonPz, rnd, NMuonPzRMS_) : 0.;
		  par[FPARAM_Lepton1Pz_TOPTOPLEPLEP] = *PzLepton1/(1.-gv);
		  par[FPARAM_Lepton1E_TOPTOPLEPLEP] = sqrt((*MassLepton1)*(*MassLepton1) + par[FPARAM_Lepton1Px_TOPTOPLEPLEP]*par[FPARAM_Lepton1Px_TOPTOPLEPLEP] + par[FPARAM_Lepton1Py_TOPTOPLEPLEP]*par[FPARAM_Lepton1Py_TOPTOPLEPLEP] + par[FPARAM_Lepton1Pz_TOPTOPLEPLEP]*par[FPARAM_Lepton1Pz_TOPTOPLEPLEP]);
	       } while( (par[FPARAM_Lepton1E_TOPTOPLEPLEP]*par[FPARAM_Lepton1E_TOPTOPLEPLEP]-par[FPARAM_Lepton1Px_TOPTOPLEPLEP]*par[FPARAM_Lepton1Px_TOPTOPLEPLEP]-par[FPARAM_Lepton1Py_TOPTOPLEPLEP]*par[FPARAM_Lepton1Py_TOPTOPLEPLEP]-par[FPARAM_Lepton1Pz_TOPTOPLEPLEP]*par[FPARAM_Lepton1Pz_TOPTOPLEPLEP]) < thres );
	  }

	if( *LabelLepton2 == 0 )
	  {
	     do
	       {
		  double gv = (doToys) ? getProbGaus(hPDFElecPx.get(), maxPDFElecPx, meanPDFElecPx, sigmaPDFElecPx, rnd, NElecPxRMS_) : 0.;
		  par[FPARAM_Lepton2Px_TOPTOPLEPLEP] = *PxLepton2/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFElecPy.get(), maxPDFElecPy, meanPDFElecPy, sigmaPDFElecPy, rnd, NElecPyRMS_) : 0.;
		  par[FPARAM_Lepton2Py_TOPTOPLEPLEP] = *PyLepton2/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFElecPz.get(), maxPDFElecPz, meanPDFElecPz, sigmaPDFElecPz, rnd, NElecPzRMS_) : 0.;
		  par[FPARAM_Lepton2Pz_TOPTOPLEPLEP] = *PzLepton2/(1.-gv);
		  par[FPARAM_Lepton2E_TOPTOPLEPLEP] = sqrt((*MassLepton2)*(*MassLepton2) + par[FPARAM_Lepton2Px_TOPTOPLEPLEP]*par[FPARAM_Lepton2Px_TOPTOPLEPLEP] + par[FPARAM_Lepton2Py_TOPTOPLEPLEP]*par[FPARAM_Lepton2Py_TOPTOPLEPLEP] + par[FPARAM_Lepton2Pz_TOPTOPLEPLEP]*par[FPARAM_Lepton2Pz_TOPTOPLEPLEP]);
	       } while( (par[FPARAM_Lepton2E_TOPTOPLEPLEP]*par[FPARAM_Lepton2E_TOPTOPLEPLEP]-par[FPARAM_Lepton2Px_TOPTOPLEPLEP]*par[FPARAM_Lepton2Px_TOPTOPLEPLEP]-par[FPARAM_Lepton2Py_TOPTOPLEPLEP]*par[FPARAM_Lepton2Py_TOPTOPLEPLEP]-par[FPARAM_Lepton2Pz_TOPTOPLEPLEP]*par[FPARAM_Lepton2Pz_TOPTOPLEPLEP]) < thres );
	  }
	else
	  {
	     do
	       {		  
		  double gv = (doToys) ? getProbGaus(hPDFMuonPx.get(), maxPDFMuonPx, meanPDFMuonPx, sigmaPDFMuonPx, rnd, NMuonPxRMS_) : 0.;
		  par[FPARAM_Lepton2Px_TOPTOPLEPLEP] = *PxLepton2/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFMuonPy.get(), maxPDFMuonPy, meanPDFMuonPy, sigmaPDFMuonPy, rnd, NMuonPyRMS_) : 0.;
		  par[FPARAM_Lepton2Py_TOPTOPLEPLEP] = *PyLepton2/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFMuonPz.get(), maxPDFMuonPz, meanPDFMuonPz, sigmaPDFMuonPz, rnd, NMuonPzRMS_) : 0.;
		  par[FPARAM_Lepton2Pz_TOPTOPLEPLEP] = *PzLepton2/(1.-gv);
		  par[FPARAM_Lepton2E_TOPTOPLEPLEP] = sqrt((*MassLepton2)*(*MassLepton2) + par[FPARAM_Lepton2Px_TOPTOPLEPLEP]*par[FPARAM_Lepton2Px_TOPTOPLEPLEP] + par[FPARAM_Lepton2Py_TOPTOPLEPLEP]*par[FPARAM_Lepton2Py_TOPTOPLEPLEP] + par[FPARAM_Lepton2Pz_TOPTOPLEPLEP]*par[FPARAM_Lepton2Pz_TOPTOPLEPLEP]);
	       } while( (par[FPARAM_Lepton2E_TOPTOPLEPLEP]*par[FPARAM_Lepton2E_TOPTOPLEPLEP]-par[FPARAM_Lepton2Px_TOPTOPLEPLEP]*par[FPARAM_Lepton2Px_TOPTOPLEPLEP]-par[FPARAM_Lepton2Py_TOPTOPLEPLEP]*par[FPARAM_Lepton2Py_TOPTOPLEPLEP]-par[FPARAM_Lepton2Pz_TOPTOPLEPLEP]*par[FPARAM_Lepton2Pz_TOPTOPLEPLEP]) < thres );
	  }

	if( doToys )
	  {	     
	     par[FPARAM_Etx_TOPTOPLEPLEP] = EtxMin_+EtxMax_*2.*(rnd->Rndm());
	     par[FPARAM_Ety_TOPTOPLEPLEP] = EtxMin_+EtxMax_*2.*(rnd->Rndm());
	  }
	else
	  {
	     par[FPARAM_Etx_TOPTOPLEPLEP] = Etx;
	     par[FPARAM_Ety_TOPTOPLEPLEP] = Ety;
	     
	     if( par[FPARAM_Etx_TOPTOPLEPLEP] > 1E+9 || par[FPARAM_Ety_TOPTOPLEPLEP] > 1E+9 )
	       {
		  std::cout << "Please call SetEtxEty without toy generation" << std::endl;
		  exit(1);
	       }	     
	  }	

	float mtoptopMin = 1E+10;
	std::vector<float> res;
	
	for( int is1=-1;is1<=1;is1++ )
	  {	
	     if( is1 == 0 ) continue;
	     
	     par[FPARAM_Sign1_TOPTOPLEPLEP] = is1;
	     
	     for( int is2=-1;is2<=1;is2++ )
	       {
		  if( is2 == 0 ) continue;
		  
		  par[FPARAM_Sign2_TOPTOPLEPLEP] = is2;

		  chi2t.clear();

		  double lh = func(*PxLepton1, *PyLepton1, *PzLepton1, *ELepton1, *LabelLepton1,
				   *PxLepton2, *PyLepton2, *PzLepton2, *ELepton2, *LabelLepton2,
				   *PxBJet1, *PyBJet1, *PzBJet1, *EBJet1,
				   *PxBJet2, *PyBJet2, *PzBJet2, *EBJet2,
				   *PxPhoton, *PyPhoton, *PzPhoton, *EPhoton,
				   *PhotonOrigin,
				   chi2t, par);

		  if( lh < LHMaxGeneric_ )
		    {
		       FRESULT fres = {};
		       for( int ip=0;ip<FPARAM_N;ip++ ) fres.par[ip] = par[ip];
		       fres.lh = lh;
		       
		       vp.push_back(fres);

		       if( vp.size() >= NGrid_ ) return;
		    }
	       }
	  }
     }   
}
