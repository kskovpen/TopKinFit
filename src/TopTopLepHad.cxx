// Top quark pair production in semileptonic channel: t -> (W->lnu) b, tbar -> (W->qq') b

#include "../include/TopTopLepHad.h"

ClassImp(KINFIT::TopTopLepHad)

KINFIT::TopTopLepHad::TopTopLepHad() {}

KINFIT::TopTopLepHad::~TopTopLepHad() {}

/// Main routine
void KINFIT::TopTopLepHad::TopTopLepHadRun()
{
   NPERMEVENT_PREV = NPERMEVENT;

   /// Calculate number of permutations in the current event
   CalcNPerm();

   /// Reinitialize all variables and release memory used in the previous event
   
   NPerm_ = 0;
   NTerm_ = 0;

   if( idxMin_ ) {delete[] idxMin_; idxMin_ = 0;}
   
   if( TopTopLepHad_ElectronIdx )     {delete[] TopTopLepHad_ElectronIdx; TopTopLepHad_ElectronIdx = 0;}
   if( TopTopLepHad_MuonIdx )         {delete[] TopTopLepHad_MuonIdx; TopTopLepHad_MuonIdx = 0;}
   if( TopTopLepHad_NonBJet1Idx )     {delete[] TopTopLepHad_NonBJet1Idx; TopTopLepHad_NonBJet1Idx = 0;}
   if( TopTopLepHad_NonBJet2Idx )     {delete[] TopTopLepHad_NonBJet2Idx; TopTopLepHad_NonBJet2Idx = 0;}
   if( TopTopLepHad_BJetLepIdx )      {delete[] TopTopLepHad_BJetLepIdx; TopTopLepHad_BJetLepIdx = 0;}
   if( TopTopLepHad_BJetHadIdx )      {delete[] TopTopLepHad_BJetHadIdx; TopTopLepHad_BJetHadIdx = 0;}
   if( TopTopLepHad_PhotonIdx )       {delete[] TopTopLepHad_PhotonIdx; TopTopLepHad_PhotonIdx = 0;}

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
   
   TopTopLepHad_ElectronIdx = new int[NPERMEVENT];
   TopTopLepHad_MuonIdx = new int[NPERMEVENT];
   TopTopLepHad_NonBJet1Idx = new int[NPERMEVENT];
   TopTopLepHad_NonBJet2Idx = new int[NPERMEVENT];
   TopTopLepHad_BJetLepIdx = new int[NPERMEVENT];
   TopTopLepHad_BJetHadIdx = new int[NPERMEVENT];
   TopTopLepHad_PhotonIdx = new int[NPERMEVENT];
   
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

	TopTopLepHad_ElectronIdx[i] = -1;
	TopTopLepHad_MuonIdx[i] = -1;
	TopTopLepHad_NonBJet1Idx[i] = -1;
	TopTopLepHad_NonBJet2Idx[i] = -1;
	TopTopLepHad_BJetLepIdx[i] = -1;
	TopTopLepHad_BJetHadIdx[i] = -1;
	TopTopLepHad_PhotonIdx[i] = -1;
	
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
   
   (*ParMin)[FPARAM_Sign_TOPTOPLEPHAD] = -1.;
   (*ParMax)[FPARAM_Sign_TOPTOPLEPHAD] = 1.;

   (*ParMin)[FPARAM_mWLep_TOPTOPLEPHAD] = 50.;
   (*ParMax)[FPARAM_mWLep_TOPTOPLEPHAD] = 110.;
   (*ParMin)[FPARAM_mWHad_TOPTOPLEPHAD] = 10.;
   (*ParMax)[FPARAM_mWHad_TOPTOPLEPHAD] = 250.;

   /// Loop through permutations

   for( int il=0;il<nLepton;il++ )
     {	
	(*ParMin)[FPARAM_EtRealX_TOPTOPLEPHAD] = (*EtMissX) - LimNRMS_*fabs(*EtMissX)*sigmaPDFMetPx;
	(*ParMax)[FPARAM_EtRealX_TOPTOPLEPHAD] = (*EtMissX) + LimNRMS_*fabs(*EtMissX)*sigmaPDFMetPx;
	(*ParMin)[FPARAM_EtRealY_TOPTOPLEPHAD] = (*EtMissY) - LimNRMS_*fabs(*EtMissY)*sigmaPDFMetPy;
	(*ParMax)[FPARAM_EtRealY_TOPTOPLEPHAD] = (*EtMissY) + LimNRMS_*fabs(*EtMissY)*sigmaPDFMetPy;	     
	
	*PxLepton1 = LeptonPx[il];
	*PyLepton1 = LeptonPy[il];
	*PzLepton1 = LeptonPz[il];
	*ELepton1 = LeptonE[il];
	*PtLepton1 = sqrt((*PxLepton1)*(*PxLepton1)+(*PyLepton1)*(*PyLepton1));
	*EtaLepton1 = getEta(*PtLepton1, *PzLepton1);
	*PhiLepton1 = atan2(*PyLepton1, *PxLepton1);
	*MassLepton1 = (*ELepton1)*(*ELepton1) - (*PxLepton1)*(*PxLepton1) - (*PyLepton1)*(*PyLepton1) - (*PzLepton1)*(*PzLepton1);
	*MassLepton1 = (*MassLepton1 > 0.) ? sqrt(*MassLepton1) : 0.;
	*LabelLepton1 = LeptonLabel[il];
	
	float sigmaPDFLeptonPx = sigmaPDFElecPx;
	float sigmaPDFLeptonPy = sigmaPDFElecPy;
	float sigmaPDFLeptonPz = sigmaPDFElecPz;
	float sigmaPDFLeptonPt = sigmaPDFElecPt;
	float sigmaPDFLeptonEta = sigmaPDFElecEta;
	float sigmaPDFLeptonPhi = sigmaPDFElecPhi;
	
	if( *LabelLepton1 == 1 )
	  {
	     sigmaPDFLeptonPx = sigmaPDFMuonPx;
	     sigmaPDFLeptonPy = sigmaPDFMuonPy;
	     sigmaPDFLeptonPz = sigmaPDFMuonPz;
	     sigmaPDFLeptonPt = sigmaPDFMuonPt;
	     sigmaPDFLeptonEta = sigmaPDFMuonEta;
	     sigmaPDFLeptonPhi = sigmaPDFMuonPhi;
	  }	     
	
	(*ParMin)[FPARAM_LeptonPx_TOPTOPLEPHAD] = (*PxLepton1) - LimNRMS_*fabs(*PxLepton1)*sigmaPDFLeptonPx;
	(*ParMax)[FPARAM_LeptonPx_TOPTOPLEPHAD] = (*PxLepton1) + LimNRMS_*fabs(*PxLepton1)*sigmaPDFLeptonPx;
	(*ParMin)[FPARAM_LeptonPy_TOPTOPLEPHAD] = (*PyLepton1) - LimNRMS_*fabs(*PyLepton1)*sigmaPDFLeptonPy;
	(*ParMax)[FPARAM_LeptonPy_TOPTOPLEPHAD] = (*PyLepton1) + LimNRMS_*fabs(*PyLepton1)*sigmaPDFLeptonPy;
	(*ParMin)[FPARAM_LeptonPz_TOPTOPLEPHAD] = (*PzLepton1) - LimNRMS_*fabs(*PzLepton1)*sigmaPDFLeptonPz;
	(*ParMax)[FPARAM_LeptonPz_TOPTOPLEPHAD] = (*PzLepton1) + LimNRMS_*fabs(*PzLepton1)*sigmaPDFLeptonPz;
	(*ParMin)[FPARAM_LeptonPt_TOPTOPLEPHAD] = (*PtLepton1) - LimNRMS_*fabs(*PtLepton1)*sigmaPDFLeptonPt;
	(*ParMax)[FPARAM_LeptonPt_TOPTOPLEPHAD] = (*PtLepton1) + LimNRMS_*fabs(*PtLepton1)*sigmaPDFLeptonPt;
	(*ParMin)[FPARAM_LeptonEta_TOPTOPLEPHAD] = (*EtaLepton1) - LimNRMS_*fabs(*EtaLepton1)*sigmaPDFLeptonEta;
	(*ParMax)[FPARAM_LeptonEta_TOPTOPLEPHAD] = (*EtaLepton1) + LimNRMS_*fabs(*EtaLepton1)*sigmaPDFLeptonEta;
	(*ParMin)[FPARAM_LeptonPhi_TOPTOPLEPHAD] = (*PhiLepton1) - LimNRMS_*fabs(*PhiLepton1)*sigmaPDFLeptonPhi;
	(*ParMax)[FPARAM_LeptonPhi_TOPTOPLEPHAD] = (*PhiLepton1) + LimNRMS_*fabs(*PhiLepton1)*sigmaPDFLeptonPhi;
	
	int label1_l = LeptonLabel[il];
	int idx1_l = LeptonIdx[il];
	
	for( int ib=0;ib<nBJet;ib++ )
	  {
	     *EBJet1 = BJetE[ib];
	     *PxBJet1 = BJetPx[ib];
	     *PyBJet1 = BJetPy[ib];
	     *PzBJet1 = BJetPz[ib];
	     *PtBJet1 = sqrt((*PxBJet1)*(*PxBJet1)+(*PyBJet1)*(*PyBJet1));
	     *EtaBJet1 = getEta(*PtBJet1, *PzBJet1);
	     *PhiBJet1 = atan2(*PyBJet1, *PxBJet1);
	     *MassBJet1 = (*EBJet1)*(*EBJet1) - (*PxBJet1)*(*PxBJet1) - (*PyBJet1)*(*PyBJet1) - (*PzBJet1)*(*PzBJet1);
	     *MassBJet1 = (*MassBJet1 > 0.) ? sqrt(*MassBJet1) : 0.;
	     
	     (*ParMin)[FPARAM_BJetLepPx_TOPTOPLEPHAD] = (*PxBJet1) - LimNRMS_*fabs(*PxBJet1)*sigmaPDFBJetPx;
	     (*ParMax)[FPARAM_BJetLepPx_TOPTOPLEPHAD] = (*PxBJet1) + LimNRMS_*fabs(*PxBJet1)*sigmaPDFBJetPx;
	     (*ParMin)[FPARAM_BJetLepPy_TOPTOPLEPHAD] = (*PyBJet1) - LimNRMS_*fabs(*PyBJet1)*sigmaPDFBJetPy;
	     (*ParMax)[FPARAM_BJetLepPy_TOPTOPLEPHAD] = (*PyBJet1) + LimNRMS_*fabs(*PyBJet1)*sigmaPDFBJetPy;
	     (*ParMin)[FPARAM_BJetLepPz_TOPTOPLEPHAD] = (*PzBJet1) - LimNRMS_*fabs(*PzBJet1)*sigmaPDFBJetPz;
	     (*ParMax)[FPARAM_BJetLepPz_TOPTOPLEPHAD] = (*PzBJet1) + LimNRMS_*fabs(*PzBJet1)*sigmaPDFBJetPz;
	     (*ParMin)[FPARAM_BJetLepPt_TOPTOPLEPHAD] = (*PtBJet1) - LimNRMS_*fabs(*PtBJet1)*sigmaPDFBJetPt;
	     (*ParMax)[FPARAM_BJetLepPt_TOPTOPLEPHAD] = (*PtBJet1) + LimNRMS_*fabs(*PtBJet1)*sigmaPDFBJetPt;
	     (*ParMin)[FPARAM_BJetLepEta_TOPTOPLEPHAD] = (*EtaBJet1) - LimNRMS_*fabs(*EtaBJet1)*sigmaPDFBJetEta;
	     (*ParMax)[FPARAM_BJetLepEta_TOPTOPLEPHAD] = (*EtaBJet1) + LimNRMS_*fabs(*EtaBJet1)*sigmaPDFBJetEta;
	     (*ParMin)[FPARAM_BJetLepPhi_TOPTOPLEPHAD] = (*PhiBJet1) - LimNRMS_*fabs(*PhiBJet1)*sigmaPDFBJetPhi;
	     (*ParMax)[FPARAM_BJetLepPhi_TOPTOPLEPHAD] = (*PhiBJet1) + LimNRMS_*fabs(*PhiBJet1)*sigmaPDFBJetPhi;
	     
	     for( int ib2=0;ib2<nBJet;ib2++ )
	       {
		  if( ib == ib2 ) continue;
		  
		  *EBJet2 = BJetE[ib2];
		  *PxBJet2 = BJetPx[ib2];
		  *PyBJet2 = BJetPy[ib2];
		  *PzBJet2 = BJetPz[ib2];
		  *PtBJet2 = sqrt((*PxBJet2)*(*PxBJet2)+(*PyBJet2)*(*PyBJet2));
		  *EtaBJet2 = getEta(*PtBJet2, *PzBJet2);
		  *PhiBJet2 = atan2(*PyBJet2, *PxBJet2);
		  *MassBJet2 = (*EBJet2)*(*EBJet2) - (*PxBJet2)*(*PxBJet2) - (*PyBJet2)*(*PyBJet2) - (*PzBJet2)*(*PzBJet2);
		  *MassBJet2 = (*MassBJet2 > 0.) ? sqrt(*MassBJet2) : 0.;
		       
		  (*ParMin)[FPARAM_BJetHadPx_TOPTOPLEPHAD] = (*PxBJet2) - LimNRMS_*fabs(*PxBJet2)*sigmaPDFBJetPx;
		  (*ParMax)[FPARAM_BJetHadPx_TOPTOPLEPHAD] = (*PxBJet2) + LimNRMS_*fabs(*PxBJet2)*sigmaPDFBJetPx;
		  (*ParMin)[FPARAM_BJetHadPy_TOPTOPLEPHAD] = (*PyBJet2) - LimNRMS_*fabs(*PyBJet2)*sigmaPDFBJetPy;
		  (*ParMax)[FPARAM_BJetHadPy_TOPTOPLEPHAD] = (*PyBJet2) + LimNRMS_*fabs(*PyBJet2)*sigmaPDFBJetPy;
		  (*ParMin)[FPARAM_BJetHadPz_TOPTOPLEPHAD] = (*PzBJet2) - LimNRMS_*fabs(*PzBJet2)*sigmaPDFBJetPz;
		  (*ParMax)[FPARAM_BJetHadPz_TOPTOPLEPHAD] = (*PzBJet2) + LimNRMS_*fabs(*PzBJet2)*sigmaPDFBJetPz;
		  (*ParMin)[FPARAM_BJetHadPt_TOPTOPLEPHAD] = (*PtBJet2) - LimNRMS_*fabs(*PtBJet2)*sigmaPDFBJetPt;
		  (*ParMax)[FPARAM_BJetHadPt_TOPTOPLEPHAD] = (*PtBJet2) + LimNRMS_*fabs(*PtBJet2)*sigmaPDFBJetPt;
		  (*ParMin)[FPARAM_BJetHadEta_TOPTOPLEPHAD] = (*EtaBJet2) - LimNRMS_*fabs(*EtaBJet2)*sigmaPDFBJetEta;
		  (*ParMax)[FPARAM_BJetHadEta_TOPTOPLEPHAD] = (*EtaBJet2) + LimNRMS_*fabs(*EtaBJet2)*sigmaPDFBJetEta;
		  (*ParMin)[FPARAM_BJetHadPhi_TOPTOPLEPHAD] = (*PhiBJet2) - LimNRMS_*fabs(*PhiBJet2)*sigmaPDFBJetPhi;
		  (*ParMax)[FPARAM_BJetHadPhi_TOPTOPLEPHAD] = (*PhiBJet2) + LimNRMS_*fabs(*PhiBJet2)*sigmaPDFBJetPhi;
		  
		  for( int ij=0;ij<nNonBJet;ij++ )
		    {
		       *ENonBJet1 = NonBJetE[ij];
		       *PxNonBJet1 = NonBJetPx[ij];
		       *PyNonBJet1 = NonBJetPy[ij];
		       *PzNonBJet1 = NonBJetPz[ij];
		       *PtNonBJet1 = sqrt((*PxNonBJet1)*(*PxNonBJet1)+(*PyNonBJet1)*(*PyNonBJet1));
		       *EtaNonBJet1 = getEta(*PtNonBJet1, *PzNonBJet1);
		       *PhiNonBJet1 = atan2(*PyNonBJet1, *PxNonBJet1);
		       *MassNonBJet1 = (*ENonBJet1)*(*ENonBJet1) - (*PxNonBJet1)*(*PxNonBJet1) - (*PyNonBJet1)*(*PyNonBJet1) - (*PzNonBJet1)*(*PzNonBJet1);
		       *MassNonBJet1 = (*MassNonBJet1 > 0.) ? sqrt(*MassNonBJet1) : 0.;
		       
		       (*ParMin)[FPARAM_NonBJet1Px_TOPTOPLEPHAD] = (*PxNonBJet1) - LimNRMS_*fabs(*PxNonBJet1)*sigmaPDFNonBJetPx;
		       (*ParMax)[FPARAM_NonBJet1Px_TOPTOPLEPHAD] = (*PxNonBJet1) + LimNRMS_*fabs(*PxNonBJet1)*sigmaPDFNonBJetPx;
		       (*ParMin)[FPARAM_NonBJet1Py_TOPTOPLEPHAD] = (*PyNonBJet1) - LimNRMS_*fabs(*PyNonBJet1)*sigmaPDFNonBJetPy;
		       (*ParMax)[FPARAM_NonBJet1Py_TOPTOPLEPHAD] = (*PyNonBJet1) + LimNRMS_*fabs(*PyNonBJet1)*sigmaPDFNonBJetPy;
		       (*ParMin)[FPARAM_NonBJet1Pz_TOPTOPLEPHAD] = (*PzNonBJet1) - LimNRMS_*fabs(*PzNonBJet1)*sigmaPDFNonBJetPz;
		       (*ParMax)[FPARAM_NonBJet1Pz_TOPTOPLEPHAD] = (*PzNonBJet1) + LimNRMS_*fabs(*PzNonBJet1)*sigmaPDFNonBJetPz;
		       (*ParMin)[FPARAM_NonBJet1Pt_TOPTOPLEPHAD] = (*PtNonBJet1) - LimNRMS_*fabs(*PtNonBJet1)*sigmaPDFNonBJetPt;
		       (*ParMax)[FPARAM_NonBJet1Pt_TOPTOPLEPHAD] = (*PtNonBJet1) + LimNRMS_*fabs(*PtNonBJet1)*sigmaPDFNonBJetPt;
		       (*ParMin)[FPARAM_NonBJet1Eta_TOPTOPLEPHAD] = (*EtaNonBJet1) - LimNRMS_*fabs(*EtaNonBJet1)*sigmaPDFNonBJetEta;
		       (*ParMax)[FPARAM_NonBJet1Eta_TOPTOPLEPHAD] = (*EtaNonBJet1) + LimNRMS_*fabs(*EtaNonBJet1)*sigmaPDFNonBJetEta;
		       (*ParMin)[FPARAM_NonBJet1Phi_TOPTOPLEPHAD] = (*PhiNonBJet1) - LimNRMS_*fabs(*PhiNonBJet1)*sigmaPDFNonBJetPhi;
		       (*ParMax)[FPARAM_NonBJet1Phi_TOPTOPLEPHAD] = (*PhiNonBJet1) + LimNRMS_*fabs(*PhiNonBJet1)*sigmaPDFNonBJetPhi;
		       
		       for( int ij2=ij+1;ij2<nNonBJet;ij2++ )
			 {
			    chi_[NPerm_] = INIT;
			    
			    *ENonBJet2 = NonBJetE[ij2];
			    *PxNonBJet2 = NonBJetPx[ij2];
			    *PyNonBJet2 = NonBJetPy[ij2];
			    *PzNonBJet2 = NonBJetPz[ij2];
			    *PtNonBJet2 = sqrt((*PxNonBJet2)*(*PxNonBJet2)+(*PyNonBJet2)*(*PyNonBJet2));
			    *EtaNonBJet2 = getEta(*PtNonBJet2, *PzNonBJet2);
			    *PhiNonBJet2 = atan2(*PyNonBJet2, *PxNonBJet2);
			    *MassNonBJet2 = (*ENonBJet2)*(*ENonBJet2) - (*PxNonBJet2)*(*PxNonBJet2) - (*PyNonBJet2)*(*PyNonBJet2) - (*PzNonBJet2)*(*PzNonBJet2);
			    *MassNonBJet2 = (*MassNonBJet2 > 0.) ? sqrt(*MassNonBJet2) : 0.;
			    
			    (*ParMin)[FPARAM_NonBJet2Px_TOPTOPLEPHAD] = (*PxNonBJet2) - LimNRMS_*fabs(*PxNonBJet2)*sigmaPDFNonBJetPx;
			    (*ParMax)[FPARAM_NonBJet2Px_TOPTOPLEPHAD] = (*PxNonBJet2) + LimNRMS_*fabs(*PxNonBJet2)*sigmaPDFNonBJetPx;
			    (*ParMin)[FPARAM_NonBJet2Py_TOPTOPLEPHAD] = (*PyNonBJet2) - LimNRMS_*fabs(*PyNonBJet2)*sigmaPDFNonBJetPy;
			    (*ParMax)[FPARAM_NonBJet2Py_TOPTOPLEPHAD] = (*PyNonBJet2) + LimNRMS_*fabs(*PyNonBJet2)*sigmaPDFNonBJetPy;
			    (*ParMin)[FPARAM_NonBJet2Pz_TOPTOPLEPHAD] = (*PzNonBJet2) - LimNRMS_*fabs(*PzNonBJet2)*sigmaPDFNonBJetPz;
			    (*ParMax)[FPARAM_NonBJet2Pz_TOPTOPLEPHAD] = (*PzNonBJet2) + LimNRMS_*fabs(*PzNonBJet2)*sigmaPDFNonBJetPz;
			    (*ParMin)[FPARAM_NonBJet2Pt_TOPTOPLEPHAD] = (*PtNonBJet2) - LimNRMS_*fabs(*PtNonBJet2)*sigmaPDFNonBJetPt;
			    (*ParMax)[FPARAM_NonBJet2Pt_TOPTOPLEPHAD] = (*PtNonBJet2) + LimNRMS_*fabs(*PtNonBJet2)*sigmaPDFNonBJetPt;
			    (*ParMin)[FPARAM_NonBJet2Eta_TOPTOPLEPHAD] = (*EtaNonBJet2) - LimNRMS_*fabs(*EtaNonBJet2)*sigmaPDFNonBJetEta;
			    (*ParMax)[FPARAM_NonBJet2Eta_TOPTOPLEPHAD] = (*EtaNonBJet2) + LimNRMS_*fabs(*EtaNonBJet2)*sigmaPDFNonBJetEta;
			    (*ParMin)[FPARAM_NonBJet2Phi_TOPTOPLEPHAD] = (*PhiNonBJet2) - LimNRMS_*fabs(*PhiNonBJet2)*sigmaPDFNonBJetPhi;
			    (*ParMax)[FPARAM_NonBJet2Phi_TOPTOPLEPHAD] = (*PhiNonBJet2) + LimNRMS_*fabs(*PhiNonBJet2)*sigmaPDFNonBJetPhi;
			    
			    for( int ipho=nPhoton;ipho>=0;ipho-- )
			      {
				 if( ipho < nPhoton ) // permute according to the number of photons + permutation with no radiation (ISR)
				   {	     
				      *EPhoton  = PhotonE[ipho];
				      *PxPhoton = PhotonPx[ipho];
				      *PyPhoton = PhotonPy[ipho];
				      *PzPhoton = PhotonPz[ipho];
				      *PtPhoton = sqrt((*PxPhoton)*(*PxPhoton)+(*PyPhoton)*(*PyPhoton));
				      *EtaPhoton = getEta(*PtPhoton, *PzPhoton);
				      *PhiPhoton = atan2(*PyPhoton, *PxPhoton);
				      
				      (*ParMin)[FPARAM_PhotonPx_TOPTOPLEPHAD] = (*PxPhoton) - LimNRMS_*fabs(*PxPhoton)*sigmaPDFPhotonPx;
				      (*ParMax)[FPARAM_PhotonPx_TOPTOPLEPHAD] = (*PxPhoton) + LimNRMS_*fabs(*PxPhoton)*sigmaPDFPhotonPx;
				      (*ParMin)[FPARAM_PhotonPy_TOPTOPLEPHAD] = (*PyPhoton) - LimNRMS_*fabs(*PyPhoton)*sigmaPDFPhotonPy;
				      (*ParMax)[FPARAM_PhotonPy_TOPTOPLEPHAD] = (*PyPhoton) + LimNRMS_*fabs(*PyPhoton)*sigmaPDFPhotonPy;
				      (*ParMin)[FPARAM_PhotonPz_TOPTOPLEPHAD] = (*PzPhoton) - LimNRMS_*fabs(*PzPhoton)*sigmaPDFPhotonPz;
				      (*ParMax)[FPARAM_PhotonPz_TOPTOPLEPHAD] = (*PzPhoton) + LimNRMS_*fabs(*PzPhoton)*sigmaPDFPhotonPz;
				      (*ParMin)[FPARAM_PhotonPt_TOPTOPLEPHAD] = (*PtPhoton) - LimNRMS_*fabs(*PtPhoton)*sigmaPDFPhotonPt;
				      (*ParMax)[FPARAM_PhotonPt_TOPTOPLEPHAD] = (*PtPhoton) + LimNRMS_*fabs(*PtPhoton)*sigmaPDFPhotonPt;
				      (*ParMin)[FPARAM_PhotonEta_TOPTOPLEPHAD] = (*EtaPhoton) - LimNRMS_*fabs(*EtaPhoton)*sigmaPDFPhotonEta;
				      (*ParMax)[FPARAM_PhotonEta_TOPTOPLEPHAD] = (*EtaPhoton) + LimNRMS_*fabs(*EtaPhoton)*sigmaPDFPhotonEta;
				      (*ParMin)[FPARAM_PhotonPhi_TOPTOPLEPHAD] = (*PhiPhoton) - LimNRMS_*fabs(*PhiPhoton)*sigmaPDFPhotonPhi;
				      (*ParMax)[FPARAM_PhotonPhi_TOPTOPLEPHAD] = (*PhiPhoton) + LimNRMS_*fabs(*PhiPhoton)*sigmaPDFPhotonPhi;
				   }
				 
				 for( int ipo=0;ipo<PHOTON_ORIGIN_N_TOPTOPLEPHAD;ipo++ )
				   {
				      if( ipo != PHOTON_FROM_TOPLEP_COMB_TOPTOPLEPHAD &&
					  ipo != PHOTON_FROM_WLEP_COMB_TOPTOPLEPHAD &&
					  ipo != PHOTON_FROM_TOPHAD_COMB_TOPTOPLEPHAD &&
					  ipo != PHOTON_FROM_WHAD_COMB_TOPTOPLEPHAD &&
					  ipo != PHOTON_FROM_ISR_TOPTOPLEPHAD &&
					  !CheckAllPhotonOrigins_ ) continue;
				      
				      if( !IncludePhotons_ && ipo != PHOTON_FROM_ISR_TOPTOPLEPHAD ) continue;
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
				      std::vector<FRESULT> lhgood_sPlus, lhgood_sMinus;
				      std::vector<FRESULT> lhbad;
				      
				      for( int ig=0;ig<lhres.size();ig++ )
					{
					   FRESULT r = lhres[ig];
						
					   if( r.lh < LHMaxMinuit_ ) 
					     {
						if( r.par[FPARAM_Sign_TOPTOPLEPHAD] == 1 ) lhgood_sPlus.push_back(r);
						else lhgood_sMinus.push_back(r);
					     }
					   else lhbad.push_back(r);
					}
					   
				      std::vector<FRESULT> lhgood;
					   
				      /// Select best NLL for each neutrino solution and sort it according to NLL
				      if( lhgood_sPlus.size() > 0) lhgood.push_back(lhgood_sPlus[0]);
				      if( lhgood_sMinus.size() > 0) lhgood.push_back(lhgood_sMinus[0]);
				      
				      std::sort(lhgood.begin(), lhgood.end(), [](const FRESULT &lhs, const FRESULT &rhs) { return lhs.lh < rhs.lh; });
					   
				      /// Sort best solutions
				      bool sorting[1] = {SortByNeutrinoPz_};
				      int nsort = 0;
				      for( int s=0;s<1;s++ )
					if( sorting[s] ) nsort++;
				      if( nsort > 1 )
					{
					   std::cout << "Multiple sorting options are chosen, please stick to only one option" << std::endl;
					   exit(1);
					}
				      
				      if( SortByNeutrinoPz_ ) std::sort(lhgood.begin(), lhgood.end(), [](const FRESULT &lhs, const FRESULT &rhs) { return lhs.PzNuSum < rhs.PzNuSum; });
				      
				      /// Append the best NLL solution from bad fits
				      if( lhbad.size() > 0 ) lhgood.push_back(lhbad[0]);

				      if( lhgood.size() > 0 )
					{
					   FRESULT fbest = lhgood[0];
					   
					   disc_ = fbest.lh;

					   TopTopLepHad_PhotonIdx[NPerm_] = (IncludePhotons_) ? ipho : -1;
					   if( IncludePhotons_ && ipho == nPhoton ) TopTopLepHad_PhotonIdx[NPerm_] = -1;
					   
					   TopTopLepHad_BJetLepIdx[NPerm_] = ib;
					   TopTopLepHad_BJetHadIdx[NPerm_] = ib2;
					   
					   TopTopLepHad_NonBJet1Idx[NPerm_] = ij;
					   TopTopLepHad_NonBJet2Idx[NPerm_] = ij2;
					   
					   if( label1_l == 0 )
					     {
						TopTopLepHad_ElectronIdx[NPerm_] = idx1_l;
						TopTopLepHad_MuonIdx[NPerm_] = -1;
					     }							  
					   else if( label1_l == 1 ) 
					     {
						TopTopLepHad_MuonIdx[NPerm_] = idx1_l;
						TopTopLepHad_ElectronIdx[NPerm_] = -1;
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
					   
					   nuPx_[NPerm_][0] = fbest.par[FPARAM_EtRealX_TOPTOPLEPHAD];
					   nuPy_[NPerm_][0] = fbest.par[FPARAM_EtRealY_TOPTOPLEPHAD];
					   nuPz_[NPerm_][0] = fbest.PzNu1;
				      
					   MetPx_[NPerm_] = fbest.par[FPARAM_EtRealX_TOPTOPLEPHAD];
					   MetPy_[NPerm_] = fbest.par[FPARAM_EtRealY_TOPTOPLEPHAD];
					   
					   WMass_[NPerm_][0] = fbest.par[FPARAM_mWLep_TOPTOPLEPHAD];
					   WMass_[NPerm_][1] = fbest.par[FPARAM_mWHad_TOPTOPLEPHAD];
					   
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
	  }	
     } // end permutations

   // choose the best permutation
   for(int ip=0;ip<NPerm_;ip++)
     {
	calcVar(ip);
     }
}

/// Calculate the number of permutations per event
void KINFIT::TopTopLepHad::CalcNPerm()
{
   NPERMEVENT = 0;

   for(int ip=nPhoton;ip>=0;ip--)
     {	
	for(int il=0;il<nLepton;il++)
	  {	
	     for(int ib=0;ib<nBJet;ib++)
	       {
		  for(int ib2=0;ib2<nBJet;ib2++)
		    {
		       if( ib == ib2 ) continue;
		       
		       for(int ij=0;ij<nNonBJet;ij++)
			 {
			    for(int ij2=0;ij2<nNonBJet;ij2++)
			      {
				 if( ij == ij2 ) continue;

				 for(int ipo=0;ipo<PHOTON_ORIGIN_N_TOPTOPLEPHAD;ipo++)
				   {
				      if( ipo != PHOTON_FROM_TOPLEP_COMB_TOPTOPLEPHAD &&
					  ipo != PHOTON_FROM_WLEP_COMB_TOPTOPLEPHAD &&
					  ipo != PHOTON_FROM_TOPHAD_COMB_TOPTOPLEPHAD &&
					  ipo != PHOTON_FROM_WHAD_COMB_TOPTOPLEPHAD &&
					  ipo != PHOTON_FROM_ISR_TOPTOPLEPHAD &&
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
}

/// Derive solutions on a given set of inputs and calculate the resultant NLL
double KINFIT::TopTopLepHad::func(float PxLepton, float PyLepton, float PzLepton, float PtLepton, float EtaLepton, float PhiLepton, float ELepton, int LabelLepton,
				  float PxBJet1, float PyBJet1, float PzBJet1, float PtBJet1, float EtaBJet1, float PhiBJet1, float EBJet1,
				  float PxBJet2, float PyBJet2, float PzBJet2, float PtBJet2, float EtaBJet2, float PhiBJet2, float EBJet2,
				  float PxNonBJet1, float PyNonBJet1, float PzNonBJet1, float PtNonBJet1, float EtaNonBJet1, float PhiNonBJet1, float ENonBJet1,
				  float PxNonBJet2, float PyNonBJet2, float PzNonBJet2, float PtNonBJet2, float EtaNonBJet2, float PhiNonBJet2, float ENonBJet2,
				  float PxPhoton, float PyPhoton, float PzPhoton, float PtPhoton, float EtaPhoton, float PhiPhoton, float EPhoton,
				  int photonOrigin,
				  std::vector<double> &chi2t, double *par)
{
   chi2t[NLL_RacPT_TOPTOPLEPHAD] = 0.;
   chi2t[NLL_mWPT_TOPTOPLEPHAD] = 0.;
   chi2t[NLL_mTopPT_TOPTOPLEPHAD] = 0.;
   
   float PxNu1 = par[FPARAM_EtRealX_TOPTOPLEPHAD];
   float PyNu1 = par[FPARAM_EtRealY_TOPTOPLEPHAD];

   *KINFIT::kfit::PxNu1 = PxNu1;
   *KINFIT::kfit::PyNu1 = PyNu1;
   
   float LeptonPx = par[FPARAM_LeptonPx_TOPTOPLEPHAD];
   float LeptonPy = par[FPARAM_LeptonPy_TOPTOPLEPHAD];
   float LeptonPz = par[FPARAM_LeptonPz_TOPTOPLEPHAD];
   float LeptonPt = par[FPARAM_LeptonPt_TOPTOPLEPHAD];
   float LeptonEta = par[FPARAM_LeptonEta_TOPTOPLEPHAD];
   float LeptonPhi = par[FPARAM_LeptonPhi_TOPTOPLEPHAD];
   float LeptonE = par[FPARAM_LeptonE_TOPTOPLEPHAD];

   float PhotonPx = par[FPARAM_PhotonPx_TOPTOPLEPHAD];
   float PhotonPy = par[FPARAM_PhotonPy_TOPTOPLEPHAD];
   float PhotonPz = par[FPARAM_PhotonPz_TOPTOPLEPHAD];
   float PhotonPt = par[FPARAM_PhotonPt_TOPTOPLEPHAD];
   float PhotonEta = par[FPARAM_PhotonEta_TOPTOPLEPHAD];
   float PhotonPhi = par[FPARAM_PhotonPhi_TOPTOPLEPHAD];
   float PhotonE = par[FPARAM_PhotonE_TOPTOPLEPHAD];   
   
   if( (photonOrigin == PHOTON_FROM_WLEP_TOPTOPLEPHAD ||
	photonOrigin == PHOTON_FROM_WLEP_COMB_TOPTOPLEPHAD ||
	photonOrigin == PHOTON_FROM_LEPTON_TOPTOPLEPHAD) &&
       IncludePhotons_ )
     {
	LeptonPx += PhotonPx;
	LeptonPy += PhotonPy;
	LeptonPz += PhotonPz;
	LeptonE += PhotonE;
     }
   
//   float a = sqrt(LeptonPx*LeptonPx+LeptonPy*LeptonPy);
   float a = sqrt(LeptonE*LeptonE-LeptonPz*LeptonPz);
   float b = LeptonPz;
   float d = sqrt(PxNu1*PxNu1+PyNu1*PyNu1);
   float f = LeptonE;

//   float c = par[FPARAM_mWLep_TOPTOPLEPHAD]*par[FPARAM_mWLep_TOPTOPLEPHAD]/2.+LeptonPx*PxNu1+LeptonPy*PyNu1;
   float c = (par[FPARAM_mWLep_TOPTOPLEPHAD]*par[FPARAM_mWLep_TOPTOPLEPHAD]-(LeptonE*LeptonE-LeptonPx*LeptonPx-LeptonPy*LeptonPy-LeptonPz*LeptonPz))/2.+LeptonPx*PxNu1+LeptonPy*PyNu1;
   
   float rac = c*c*b*b-a*a*(d*d*f*f-c*c);
   
   float racAbs = fabs(rac);
   
   float racPos = rac;
   
   if( racPos < 0. ) racPos = 0.;

   float PzNu1 = (c*b+par[FPARAM_Sign_TOPTOPLEPHAD]*sqrt(racPos))/a/a;
   
   *KINFIT::kfit::PzNu1 = PzNu1;
   
   float ENu1 = sqrt(PxNu1*PxNu1+PyNu1*PyNu1+PzNu1*PzNu1);

   float BJetLepPx = par[FPARAM_BJetLepPx_TOPTOPLEPHAD];
   float BJetLepPy = par[FPARAM_BJetLepPy_TOPTOPLEPHAD];
   float BJetLepPz = par[FPARAM_BJetLepPz_TOPTOPLEPHAD];
   float BJetLepPt = par[FPARAM_BJetLepPt_TOPTOPLEPHAD];
   float BJetLepEta = par[FPARAM_BJetLepEta_TOPTOPLEPHAD];
   float BJetLepPhi = par[FPARAM_BJetLepPhi_TOPTOPLEPHAD];
   float BJetLepE = par[FPARAM_BJetLepE_TOPTOPLEPHAD];

   if( photonOrigin == PHOTON_FROM_BJETLEP_TOPTOPLEPHAD && IncludePhotons_ )
     {
	BJetLepPx += PhotonPx;
	BJetLepPy += PhotonPy;
	BJetLepPz += PhotonPz;
	BJetLepE += PhotonE;
     }
   
   float BJetHadPx = par[FPARAM_BJetHadPx_TOPTOPLEPHAD];
   float BJetHadPy = par[FPARAM_BJetHadPy_TOPTOPLEPHAD];
   float BJetHadPz = par[FPARAM_BJetHadPz_TOPTOPLEPHAD];
   float BJetHadPt = par[FPARAM_BJetHadPt_TOPTOPLEPHAD];
   float BJetHadEta = par[FPARAM_BJetHadEta_TOPTOPLEPHAD];
   float BJetHadPhi = par[FPARAM_BJetHadPhi_TOPTOPLEPHAD];
   float BJetHadE = par[FPARAM_BJetHadE_TOPTOPLEPHAD];

   if( photonOrigin == PHOTON_FROM_BJETHAD_TOPTOPLEPHAD && IncludePhotons_ )
     {
	BJetHadPx += PhotonPx;
	BJetHadPy += PhotonPy;
	BJetHadPz += PhotonPz;
	BJetHadE += PhotonE;
     }

   float NonBJet1Px = par[FPARAM_NonBJet1Px_TOPTOPLEPHAD];
   float NonBJet1Py = par[FPARAM_NonBJet1Py_TOPTOPLEPHAD];
   float NonBJet1Pz = par[FPARAM_NonBJet1Pz_TOPTOPLEPHAD];
   float NonBJet1Pt = par[FPARAM_NonBJet1Pt_TOPTOPLEPHAD];
   float NonBJet1Eta = par[FPARAM_NonBJet1Eta_TOPTOPLEPHAD];
   float NonBJet1Phi = par[FPARAM_NonBJet1Phi_TOPTOPLEPHAD];
   float NonBJet1E = par[FPARAM_NonBJet1E_TOPTOPLEPHAD];

   if( photonOrigin == PHOTON_FROM_NONBJET1_TOPTOPLEPHAD && IncludePhotons_ )
     {
	NonBJet1Px += PhotonPx;
	NonBJet1Py += PhotonPy;
	NonBJet1Pz += PhotonPz;
	NonBJet1E += PhotonE;
     }

   float NonBJet2Px = par[FPARAM_NonBJet2Px_TOPTOPLEPHAD];
   float NonBJet2Py = par[FPARAM_NonBJet2Py_TOPTOPLEPHAD];
   float NonBJet2Pz = par[FPARAM_NonBJet2Pz_TOPTOPLEPHAD];
   float NonBJet2Pt = par[FPARAM_NonBJet2Pt_TOPTOPLEPHAD];
   float NonBJet2Eta = par[FPARAM_NonBJet2Eta_TOPTOPLEPHAD];
   float NonBJet2Phi = par[FPARAM_NonBJet2Phi_TOPTOPLEPHAD];
   float NonBJet2E = par[FPARAM_NonBJet2E_TOPTOPLEPHAD];
   
   if( photonOrigin == PHOTON_FROM_NONBJET2_TOPTOPLEPHAD && IncludePhotons_ )
     {
	NonBJet2Px += PhotonPx;
	NonBJet2Py += PhotonPy;
	NonBJet2Pz += PhotonPz;
	NonBJet2E += PhotonE;
     }
   
   float totPxLep = BJetLepPx+PxNu1+LeptonPx;
   float totPyLep = BJetLepPy+PyNu1+LeptonPy;
   float totPzLep = BJetLepPz+PzNu1+LeptonPz;
   float totELep = BJetLepE+ENu1+LeptonE;
   
   if( (photonOrigin == PHOTON_FROM_TOPLEP_TOPTOPLEPHAD ||
	photonOrigin == PHOTON_FROM_TOPLEP_COMB_TOPTOPLEPHAD) && IncludePhotons_ )
     {
	totPxLep += PhotonPx;
	totPyLep += PhotonPy;
	totPzLep += PhotonPz;
	totELep += PhotonE;
     }   
   
   float ENuL1 = ENu1+LeptonE;
   float PxNuL1 = PxNu1+LeptonPx;
   float PyNuL1 = PyNu1+LeptonPy;
   float PzNuL1 = PzNu1+LeptonPz;
   
   float mWLep = ENuL1*ENuL1-PxNuL1*PxNuL1-PyNuL1*PyNuL1-PzNuL1*PzNuL1;
   float mtopLep = totELep*totELep-totPxLep*totPxLep-totPyLep*totPyLep-totPzLep*totPzLep;
   
   float EJJ = NonBJet1E+NonBJet2E;
   float PxJJ = NonBJet1Px+NonBJet2Px;
   float PyJJ = NonBJet1Py+NonBJet2Py;
   float PzJJ = NonBJet1Pz+NonBJet2Pz;

   if( (photonOrigin == PHOTON_FROM_WHAD_TOPTOPLEPHAD ||
	photonOrigin == PHOTON_FROM_WHAD_COMB_TOPTOPLEPHAD) && IncludePhotons_ )
     {
	PxJJ += PhotonPx;
	PyJJ += PhotonPy;
	PzJJ += PhotonPz;
	EJJ += PhotonE;
     }
   
   float totPxHad = BJetHadPx+PxJJ;
   float totPyHad = BJetHadPy+PyJJ;
   float totPzHad = BJetHadPz+PzJJ;
   float totEHad = BJetHadE+EJJ;

   if( (photonOrigin == PHOTON_FROM_TOPHAD_TOPTOPLEPHAD ||
	photonOrigin == PHOTON_FROM_TOPHAD_COMB_TOPTOPLEPHAD) && IncludePhotons_ )
     {
	totPxHad += PhotonPx;
	totPyHad += PhotonPy;
	totPzHad += PhotonPz;
	totEHad += PhotonE;
     }   
   
   float mWHad = EJJ*EJJ-PxJJ*PxJJ-PyJJ*PyJJ-PzJJ*PzJJ;
   float mtopHad = totEHad*totEHad-totPxHad*totPxHad-totPyHad*totPyHad-totPzHad*totPzHad;
   
   float totELepEHad = totELep+totEHad;
   float totPxLepPxHad = totPxLep+totPxHad;
   float totPyLepPyHad = totPyLep+totPyHad;
   float totPzLepPzHad = totPzLep+totPzHad;
  
   float mtoptop = totELepEHad*totELepEHad-totPxLepPxHad*totPxLepPxHad-totPyLepPyHad*totPyLepPyHad-totPzLepPzHad*totPzLepPzHad;
   float mtoptopAbs = sqrt(fabs(mtoptop));
   
   *KINFIT::kfit::TopTopMass = mtoptopAbs;
   *KINFIT::kfit::TopTopPz = fabs(totPzLepPzHad);
   *KINFIT::kfit::TopTopPzAbs = fabs(totPzLep) + fabs(totPzHad);

   double val = 0.;

   /// Penalty term for negative roots
   if( rac < 0 ) 
     {
	float racPT = log(1.+racAbs);
	val += racPT;
	chi2t[NLL_RacPT_TOPTOPLEPHAD] += racPT;
     }   

   float mWLepAbs = fabs(mWLep);
   float mWHadAbs = fabs(mWHad);
   float mtopLepAbs = fabs(mtopLep);
   float mtopHadAbs = fabs(mtopHad);
   
   /// Penalty terms for negative masses
   if( mWLep < 0 )
     {
	float mWLepPT = log(1.+mWLepAbs);
	val += mWLepPT;
	chi2t[NLL_mWPT_TOPTOPLEPHAD] += mWLepPT;
     }   
   if( mWHad < 0 )
     {
	float mWHadPT = log(1.+mWHadAbs);
	val += mWHadPT;
	chi2t[NLL_mWPT_TOPTOPLEPHAD] += mWHadPT;
     }   
   if( mtopLep < 0 )
     {
	float mTopLepPT = log(1.+mtopLepAbs);
	val += mTopLepPT;
	chi2t[NLL_mTopPT_TOPTOPLEPHAD] += mTopLepPT;
     }   
   if( mtopHad < 0 )
     {
	float mTopHadPT = log(1.+mtopHadAbs);
	val += mTopHadPT;
	chi2t[NLL_mTopPT_TOPTOPLEPHAD] += mTopHadPT;
     }   
   
   mtopLepAbs = sqrt(mtopLepAbs);
   mWLepAbs = sqrt(mWLepAbs);

   mtopHadAbs = sqrt(mtopHadAbs);
   mWHadAbs = sqrt(mWHadAbs);

   float mWLepProb = (! isDeltaFuncPDFTopWMass) ? getProb(hPDFTopWMass.get(), mWLepAbs, maxPDFTopWMass, xminPDFTopWMass, xmaxPDFTopWMass) : 1.0;
   float mWHadProb = (! isDeltaFuncPDFTopWHadMass) ? getProb(hPDFTopWHadMass.get(), mWHadAbs, maxPDFTopWHadMass, xminPDFTopWHadMass, xmaxPDFTopWHadMass) : 1.0;

   float mTopLepProb = (! isDeltaFuncPDFTopMass) ? getProb(hPDFTopMass.get(), mtopLepAbs, maxPDFTopMass, xminPDFTopMass, xmaxPDFTopMass) : 1.0;
   float mTopHadProb = (! isDeltaFuncPDFTopHadMass) ? getProb(hPDFTopHadMass.get(), mtopHadAbs, maxPDFTopHadMass, xminPDFTopHadMass, xmaxPDFTopHadMass) : 1.0;

   float MetPxProb = (! (*IsParFixed)[FPARAM_EtRealX_TOPTOPLEPHAD] && ! isDeltaFuncPDFMetPx) ? getProb(hPDFMetPx.get(), (par[FPARAM_EtRealX_TOPTOPLEPHAD]-*EtMissX)/par[FPARAM_EtRealX_TOPTOPLEPHAD], maxPDFMetPx, xminPDFMetPx, xmaxPDFMetPx) : 1.0;
   float MetPyProb = (! (*IsParFixed)[FPARAM_EtRealY_TOPTOPLEPHAD] && ! isDeltaFuncPDFMetPy) ? getProb(hPDFMetPy.get(), (par[FPARAM_EtRealY_TOPTOPLEPHAD]-*EtMissY)/par[FPARAM_EtRealY_TOPTOPLEPHAD], maxPDFMetPy, xminPDFMetPy, xmaxPDFMetPy) : 1.0;
   
   float BJetLepPxProb = (! (*IsParFixed)[FPARAM_BJetLepPx_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPx) ? getProb(hPDFBJetPx.get(), (par[FPARAM_BJetLepPx_TOPTOPLEPHAD]-PxBJet1)/par[FPARAM_BJetLepPx_TOPTOPLEPHAD], maxPDFBJetPx, xminPDFBJetPx, xmaxPDFBJetPx) : 1.0;
   float BJetLepPyProb = (! (*IsParFixed)[FPARAM_BJetLepPy_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPy) ? getProb(hPDFBJetPy.get(), (par[FPARAM_BJetLepPy_TOPTOPLEPHAD]-PyBJet1)/par[FPARAM_BJetLepPy_TOPTOPLEPHAD], maxPDFBJetPy, xminPDFBJetPy, xmaxPDFBJetPy) : 1.0;
   float BJetLepPzProb = (! (*IsParFixed)[FPARAM_BJetLepPz_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPz) ? getProb(hPDFBJetPz.get(), (par[FPARAM_BJetLepPz_TOPTOPLEPHAD]-PzBJet1)/par[FPARAM_BJetLepPz_TOPTOPLEPHAD], maxPDFBJetPz, xminPDFBJetPz, xmaxPDFBJetPz) : 1.0;
   float BJetLepPtProb = (! (*IsParFixed)[FPARAM_BJetLepPt_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPt) ? getProb(hPDFBJetPt.get(), (par[FPARAM_BJetLepPt_TOPTOPLEPHAD]-PtBJet1)/par[FPARAM_BJetLepPt_TOPTOPLEPHAD], maxPDFBJetPt, xminPDFBJetPt, xmaxPDFBJetPt) : 1.0;
   float BJetLepEtaProb = (! (*IsParFixed)[FPARAM_BJetLepEta_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetEta) ? getProb(hPDFBJetEta.get(), (par[FPARAM_BJetLepEta_TOPTOPLEPHAD]-EtaBJet1)/par[FPARAM_BJetLepEta_TOPTOPLEPHAD], maxPDFBJetEta, xminPDFBJetEta, xmaxPDFBJetEta) : 1.0;
   float BJetLepPhiProb = (! (*IsParFixed)[FPARAM_BJetLepPhi_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPhi) ? getProb(hPDFBJetPhi.get(), (par[FPARAM_BJetLepPhi_TOPTOPLEPHAD]-PhiBJet1)/par[FPARAM_BJetLepPhi_TOPTOPLEPHAD], maxPDFBJetPhi, xminPDFBJetPhi, xmaxPDFBJetPhi) : 1.0;

   float BJetHadPxProb = (! (*IsParFixed)[FPARAM_BJetHadPx_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPx) ? getProb(hPDFBJetPx.get(), (par[FPARAM_BJetHadPx_TOPTOPLEPHAD]-PxBJet2)/par[FPARAM_BJetHadPx_TOPTOPLEPHAD], maxPDFBJetPx, xminPDFBJetPx, xmaxPDFBJetPx) : 1.0;
   float BJetHadPyProb = (! (*IsParFixed)[FPARAM_BJetHadPy_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPy) ? getProb(hPDFBJetPy.get(), (par[FPARAM_BJetHadPy_TOPTOPLEPHAD]-PyBJet2)/par[FPARAM_BJetHadPy_TOPTOPLEPHAD], maxPDFBJetPy, xminPDFBJetPy, xmaxPDFBJetPy) : 1.0;
   float BJetHadPzProb = (! (*IsParFixed)[FPARAM_BJetHadPz_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPz) ? getProb(hPDFBJetPz.get(), (par[FPARAM_BJetHadPz_TOPTOPLEPHAD]-PzBJet2)/par[FPARAM_BJetHadPz_TOPTOPLEPHAD], maxPDFBJetPz, xminPDFBJetPz, xmaxPDFBJetPz) : 1.0;
   float BJetHadPtProb = (! (*IsParFixed)[FPARAM_BJetHadPt_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPt) ? getProb(hPDFBJetPt.get(), (par[FPARAM_BJetHadPt_TOPTOPLEPHAD]-PtBJet2)/par[FPARAM_BJetHadPt_TOPTOPLEPHAD], maxPDFBJetPt, xminPDFBJetPt, xmaxPDFBJetPt) : 1.0;
   float BJetHadEtaProb = (! (*IsParFixed)[FPARAM_BJetHadEta_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetEta) ? getProb(hPDFBJetEta.get(), (par[FPARAM_BJetHadEta_TOPTOPLEPHAD]-EtaBJet2)/par[FPARAM_BJetHadEta_TOPTOPLEPHAD], maxPDFBJetEta, xminPDFBJetEta, xmaxPDFBJetEta) : 1.0;
   float BJetHadPhiProb = (! (*IsParFixed)[FPARAM_BJetHadPhi_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPhi) ? getProb(hPDFBJetPhi.get(), (par[FPARAM_BJetHadPhi_TOPTOPLEPHAD]-PhiBJet2)/par[FPARAM_BJetHadPhi_TOPTOPLEPHAD], maxPDFBJetPhi, xminPDFBJetPhi, xmaxPDFBJetPhi) : 1.0;

   float NonBJet1PxProb = (! (*IsParFixed)[FPARAM_NonBJet1Px_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPx) ? getProb(hPDFNonBJetPx.get(), (par[FPARAM_NonBJet1Px_TOPTOPLEPHAD]-PxNonBJet1)/par[FPARAM_NonBJet1Px_TOPTOPLEPHAD], maxPDFNonBJetPx, xminPDFNonBJetPx, xmaxPDFNonBJetPx) : 1.0;
   float NonBJet1PyProb = (! (*IsParFixed)[FPARAM_NonBJet1Py_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPy) ? getProb(hPDFNonBJetPy.get(), (par[FPARAM_NonBJet1Py_TOPTOPLEPHAD]-PyNonBJet1)/par[FPARAM_NonBJet1Py_TOPTOPLEPHAD], maxPDFNonBJetPy, xminPDFNonBJetPy, xmaxPDFNonBJetPy) : 1.0;
   float NonBJet1PzProb = (! (*IsParFixed)[FPARAM_NonBJet1Pz_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPz) ? getProb(hPDFNonBJetPz.get(), (par[FPARAM_NonBJet1Pz_TOPTOPLEPHAD]-PzNonBJet1)/par[FPARAM_NonBJet1Pz_TOPTOPLEPHAD], maxPDFNonBJetPz, xminPDFNonBJetPz, xmaxPDFNonBJetPz) : 1.0;
   float NonBJet1PtProb = (! (*IsParFixed)[FPARAM_NonBJet1Pt_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPt) ? getProb(hPDFNonBJetPt.get(), (par[FPARAM_NonBJet1Pt_TOPTOPLEPHAD]-PtNonBJet1)/par[FPARAM_NonBJet1Pt_TOPTOPLEPHAD], maxPDFNonBJetPt, xminPDFNonBJetPt, xmaxPDFNonBJetPt) : 1.0;
   float NonBJet1EtaProb = (! (*IsParFixed)[FPARAM_NonBJet1Eta_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetEta) ? getProb(hPDFNonBJetEta.get(), (par[FPARAM_NonBJet1Eta_TOPTOPLEPHAD]-EtaNonBJet1)/par[FPARAM_NonBJet1Eta_TOPTOPLEPHAD], maxPDFNonBJetEta, xminPDFNonBJetEta, xmaxPDFNonBJetEta) : 1.0;
   float NonBJet1PhiProb = (! (*IsParFixed)[FPARAM_NonBJet1Phi_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPhi) ? getProb(hPDFNonBJetPhi.get(), (par[FPARAM_NonBJet1Phi_TOPTOPLEPHAD]-PhiNonBJet1)/par[FPARAM_NonBJet1Phi_TOPTOPLEPHAD], maxPDFNonBJetPhi, xminPDFNonBJetPhi, xmaxPDFNonBJetPhi) : 1.0;

   float NonBJet2PxProb = (! (*IsParFixed)[FPARAM_NonBJet2Px_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPx) ? getProb(hPDFNonBJetPx.get(), (par[FPARAM_NonBJet2Px_TOPTOPLEPHAD]-PxNonBJet2)/par[FPARAM_NonBJet2Px_TOPTOPLEPHAD], maxPDFNonBJetPx, xminPDFNonBJetPx, xmaxPDFNonBJetPx) : 1.0;
   float NonBJet2PyProb = (! (*IsParFixed)[FPARAM_NonBJet2Py_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPy) ? getProb(hPDFNonBJetPy.get(), (par[FPARAM_NonBJet2Py_TOPTOPLEPHAD]-PyNonBJet2)/par[FPARAM_NonBJet2Py_TOPTOPLEPHAD], maxPDFNonBJetPy, xminPDFNonBJetPy, xmaxPDFNonBJetPy) : 1.0;
   float NonBJet2PzProb = (! (*IsParFixed)[FPARAM_NonBJet2Pz_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPz) ? getProb(hPDFNonBJetPz.get(), (par[FPARAM_NonBJet2Pz_TOPTOPLEPHAD]-PzNonBJet2)/par[FPARAM_NonBJet2Pz_TOPTOPLEPHAD], maxPDFNonBJetPz, xminPDFNonBJetPz, xmaxPDFNonBJetPz) : 1.0;
   float NonBJet2PtProb = (! (*IsParFixed)[FPARAM_NonBJet2Pt_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPt) ? getProb(hPDFNonBJetPt.get(), (par[FPARAM_NonBJet2Pt_TOPTOPLEPHAD]-PtNonBJet2)/par[FPARAM_NonBJet2Pt_TOPTOPLEPHAD], maxPDFNonBJetPt, xminPDFNonBJetPt, xmaxPDFNonBJetPt) : 1.0;
   float NonBJet2EtaProb = (! (*IsParFixed)[FPARAM_NonBJet2Eta_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetEta) ? getProb(hPDFNonBJetEta.get(), (par[FPARAM_NonBJet2Eta_TOPTOPLEPHAD]-EtaNonBJet2)/par[FPARAM_NonBJet2Eta_TOPTOPLEPHAD], maxPDFNonBJetEta, xminPDFNonBJetEta, xmaxPDFNonBJetEta) : 1.0;
   float NonBJet2PhiProb = (! (*IsParFixed)[FPARAM_NonBJet2Phi_TOPTOPLEPHAD] && ! isDeltaFuncPDFBJetPhi) ? getProb(hPDFNonBJetPhi.get(), (par[FPARAM_NonBJet2Phi_TOPTOPLEPHAD]-PhiNonBJet2)/par[FPARAM_NonBJet2Phi_TOPTOPLEPHAD], maxPDFNonBJetPhi, xminPDFNonBJetPhi, xmaxPDFNonBJetPhi) : 1.0;

   float PhotonPxProb = (! (*IsParFixed)[FPARAM_PhotonPx_TOPTOPLEPHAD] && IncludePhotons_ && ! isDeltaFuncPDFPhotonPx) ? getProb(hPDFPhotonPx.get(), (par[FPARAM_PhotonPx_TOPTOPLEPHAD]-PxPhoton)/par[FPARAM_PhotonPx_TOPTOPLEPHAD], maxPDFPhotonPx, xminPDFPhotonPx, xmaxPDFPhotonPx) : 1.0;
   float PhotonPyProb = (! (*IsParFixed)[FPARAM_PhotonPy_TOPTOPLEPHAD] && IncludePhotons_ && ! isDeltaFuncPDFPhotonPy) ? getProb(hPDFPhotonPy.get(), (par[FPARAM_PhotonPy_TOPTOPLEPHAD]-PyPhoton)/par[FPARAM_PhotonPy_TOPTOPLEPHAD], maxPDFPhotonPy, xminPDFPhotonPy, xmaxPDFPhotonPy) : 1.0;
   float PhotonPzProb = (! (*IsParFixed)[FPARAM_PhotonPz_TOPTOPLEPHAD] && IncludePhotons_ && ! isDeltaFuncPDFPhotonPz) ? getProb(hPDFPhotonPz.get(), (par[FPARAM_PhotonPz_TOPTOPLEPHAD]-PzPhoton)/par[FPARAM_PhotonPz_TOPTOPLEPHAD], maxPDFPhotonPz, xminPDFPhotonPz, xmaxPDFPhotonPz) : 1.0;
   float PhotonPtProb = (! (*IsParFixed)[FPARAM_PhotonPt_TOPTOPLEPHAD] && IncludePhotons_ && ! isDeltaFuncPDFPhotonPt) ? getProb(hPDFPhotonPt.get(), (par[FPARAM_PhotonPt_TOPTOPLEPHAD]-PtPhoton)/par[FPARAM_PhotonPt_TOPTOPLEPHAD], maxPDFPhotonPt, xminPDFPhotonPt, xmaxPDFPhotonPt) : 1.0;
   float PhotonEtaProb = (! (*IsParFixed)[FPARAM_PhotonEta_TOPTOPLEPHAD] && IncludePhotons_ && ! isDeltaFuncPDFPhotonEta) ? getProb(hPDFPhotonEta.get(), (par[FPARAM_PhotonEta_TOPTOPLEPHAD]-EtaPhoton)/par[FPARAM_PhotonEta_TOPTOPLEPHAD], maxPDFPhotonEta, xminPDFPhotonEta, xmaxPDFPhotonEta) : 1.0;
   float PhotonPhiProb = (! (*IsParFixed)[FPARAM_PhotonPhi_TOPTOPLEPHAD] && IncludePhotons_ && ! isDeltaFuncPDFPhotonPhi) ? getProb(hPDFPhotonPhi.get(), (par[FPARAM_PhotonPhi_TOPTOPLEPHAD]-PhiPhoton)/par[FPARAM_PhotonPhi_TOPTOPLEPHAD], maxPDFPhotonPhi, xminPDFPhotonPhi, xmaxPDFPhotonPhi) : 1.0;
   
   float LeptonPxProb = 1.0;
   float LeptonPyProb = 1.0;
   float LeptonPzProb = 1.0;
   float LeptonPtProb = 1.0;
   float LeptonEtaProb = 1.0;
   float LeptonPhiProb = 1.0;
   
   if (! (*IsParFixed)[FPARAM_LeptonPx_TOPTOPLEPHAD]) LeptonPxProb = (LabelLepton == 0) ? ((! isDeltaFuncPDFElecPx) ? getProb(hPDFElecPx.get(), (par[FPARAM_LeptonPx_TOPTOPLEPHAD]-PxLepton)/par[FPARAM_LeptonPx_TOPTOPLEPHAD], maxPDFElecPx, xminPDFElecPx, xmaxPDFElecPx) : 1.0) : ((! isDeltaFuncPDFMuonPx) ? getProb(hPDFMuonPx.get(), (par[FPARAM_LeptonPx_TOPTOPLEPHAD]-PxLepton)/par[FPARAM_LeptonPx_TOPTOPLEPHAD], maxPDFMuonPx, xminPDFMuonPx, xmaxPDFMuonPx) : 1.0);
   if (! (*IsParFixed)[FPARAM_LeptonPy_TOPTOPLEPHAD]) LeptonPyProb = (LabelLepton == 0) ? ((! isDeltaFuncPDFElecPy) ? getProb(hPDFElecPy.get(), (par[FPARAM_LeptonPy_TOPTOPLEPHAD]-PyLepton)/par[FPARAM_LeptonPy_TOPTOPLEPHAD], maxPDFElecPy, xminPDFElecPy, xmaxPDFElecPy) : 1.0) : ((! isDeltaFuncPDFMuonPy) ? getProb(hPDFMuonPy.get(), (par[FPARAM_LeptonPy_TOPTOPLEPHAD]-PyLepton)/par[FPARAM_LeptonPy_TOPTOPLEPHAD], maxPDFMuonPy, xminPDFMuonPy, xmaxPDFMuonPy) : 1.0);
   if (! (*IsParFixed)[FPARAM_LeptonPz_TOPTOPLEPHAD]) LeptonPzProb = (LabelLepton == 0) ? ((! isDeltaFuncPDFElecPz) ? getProb(hPDFElecPz.get(), (par[FPARAM_LeptonPz_TOPTOPLEPHAD]-PzLepton)/par[FPARAM_LeptonPz_TOPTOPLEPHAD], maxPDFElecPz, xminPDFElecPz, xmaxPDFElecPz) : 1.0) : ((! isDeltaFuncPDFMuonPz) ? getProb(hPDFMuonPz.get(), (par[FPARAM_LeptonPz_TOPTOPLEPHAD]-PzLepton)/par[FPARAM_LeptonPz_TOPTOPLEPHAD], maxPDFMuonPz, xminPDFMuonPz, xmaxPDFMuonPz) : 1.0);
   if (! (*IsParFixed)[FPARAM_LeptonPt_TOPTOPLEPHAD]) LeptonPtProb = (LabelLepton == 0) ? ((! isDeltaFuncPDFElecPt) ? getProb(hPDFElecPt.get(), (par[FPARAM_LeptonPt_TOPTOPLEPHAD]-PtLepton)/par[FPARAM_LeptonPt_TOPTOPLEPHAD], maxPDFElecPt, xminPDFElecPt, xmaxPDFElecPt) : 1.0) : ((! isDeltaFuncPDFMuonPt) ? getProb(hPDFMuonPt.get(), (par[FPARAM_LeptonPt_TOPTOPLEPHAD]-PtLepton)/par[FPARAM_LeptonPt_TOPTOPLEPHAD], maxPDFMuonPt, xminPDFMuonPt, xmaxPDFMuonPt) : 1.0);
   if (! (*IsParFixed)[FPARAM_LeptonEta_TOPTOPLEPHAD]) LeptonEtaProb = (LabelLepton == 0) ? ((! isDeltaFuncPDFElecEta) ? getProb(hPDFElecEta.get(), (par[FPARAM_LeptonEta_TOPTOPLEPHAD]-EtaLepton)/par[FPARAM_LeptonEta_TOPTOPLEPHAD], maxPDFElecEta, xminPDFElecEta, xmaxPDFElecEta) : 1.0) : ((! isDeltaFuncPDFMuonEta) ? getProb(hPDFMuonEta.get(), (par[FPARAM_LeptonEta_TOPTOPLEPHAD]-EtaLepton)/par[FPARAM_LeptonEta_TOPTOPLEPHAD], maxPDFMuonEta, xminPDFMuonEta, xmaxPDFMuonEta) : 1.0);
   if (! (*IsParFixed)[FPARAM_LeptonPhi_TOPTOPLEPHAD]) LeptonPhiProb = (LabelLepton == 0) ? ((! isDeltaFuncPDFElecPhi) ? getProb(hPDFElecPhi.get(), (par[FPARAM_LeptonPhi_TOPTOPLEPHAD]-PhiLepton)/par[FPARAM_LeptonPhi_TOPTOPLEPHAD], maxPDFElecPhi, xminPDFElecPhi, xmaxPDFElecPhi) : 1.0) : ((! isDeltaFuncPDFMuonPhi) ? getProb(hPDFMuonPhi.get(), (par[FPARAM_LeptonPhi_TOPTOPLEPHAD]-PhiLepton)/par[FPARAM_LeptonPhi_TOPTOPLEPHAD], maxPDFMuonPhi, xminPDFMuonPhi, xmaxPDFMuonPhi) : 1.0);

   double minProb = 1E-20; // NLL = 92.1034
   
   chi2t[NLL_WLep_TOPTOPLEPHAD] = (mWLepProb > minProb) ? mWLepProb : minProb;
   chi2t[NLL_WHad_TOPTOPLEPHAD] = (mWHadProb > minProb) ? mWHadProb : minProb;
   chi2t[NLL_TopLep_TOPTOPLEPHAD] = (mTopLepProb > minProb) ? mTopLepProb : minProb;
   chi2t[NLL_TopHad_TOPTOPLEPHAD] = (mTopHadProb > minProb) ? mTopHadProb : minProb;
   chi2t[NLL_EtMissX_TOPTOPLEPHAD] = (MetPxProb > minProb) ? MetPxProb : minProb;
   chi2t[NLL_EtMissY_TOPTOPLEPHAD] = (MetPyProb > minProb) ? MetPyProb : minProb;
   chi2t[NLL_BJetLepPx_TOPTOPLEPHAD] = (BJetLepPxProb > minProb) ? BJetLepPxProb : minProb;
   chi2t[NLL_BJetLepPy_TOPTOPLEPHAD] = (BJetLepPyProb > minProb) ? BJetLepPyProb : minProb; 
   chi2t[NLL_BJetLepPz_TOPTOPLEPHAD] = (BJetLepPzProb > minProb) ? BJetLepPzProb : minProb;
   chi2t[NLL_BJetLepPt_TOPTOPLEPHAD] = (BJetLepPtProb > minProb) ? BJetLepPtProb : minProb;
   chi2t[NLL_BJetLepEta_TOPTOPLEPHAD] = (BJetLepEtaProb > minProb) ? BJetLepEtaProb : minProb;
   chi2t[NLL_BJetLepPhi_TOPTOPLEPHAD] = (BJetLepPhiProb > minProb) ? BJetLepPhiProb : minProb;
   chi2t[NLL_BJetHadPx_TOPTOPLEPHAD] = (BJetHadPxProb > minProb) ? BJetHadPxProb : minProb;
   chi2t[NLL_BJetHadPy_TOPTOPLEPHAD] = (BJetHadPyProb > minProb) ? BJetHadPyProb : minProb;
   chi2t[NLL_BJetHadPz_TOPTOPLEPHAD] = (BJetHadPzProb > minProb) ? BJetHadPzProb : minProb;
   chi2t[NLL_BJetHadPt_TOPTOPLEPHAD] = (BJetHadPtProb > minProb) ? BJetHadPtProb : minProb;
   chi2t[NLL_BJetHadEta_TOPTOPLEPHAD] = (BJetHadEtaProb > minProb) ? BJetHadEtaProb : minProb;
   chi2t[NLL_BJetHadPhi_TOPTOPLEPHAD] = (BJetHadPhiProb > minProb) ? BJetHadPhiProb : minProb;
   chi2t[NLL_NonBJet1Px_TOPTOPLEPHAD] = (NonBJet1PxProb > minProb) ? NonBJet1PxProb : minProb;
   chi2t[NLL_NonBJet1Py_TOPTOPLEPHAD] = (NonBJet1PyProb > minProb) ? NonBJet1PyProb : minProb; 
   chi2t[NLL_NonBJet1Pz_TOPTOPLEPHAD] = (NonBJet1PzProb > minProb) ? NonBJet1PzProb : minProb;
   chi2t[NLL_NonBJet1Pt_TOPTOPLEPHAD] = (NonBJet1PtProb > minProb) ? NonBJet1PtProb : minProb;
   chi2t[NLL_NonBJet1Eta_TOPTOPLEPHAD] = (NonBJet1EtaProb > minProb) ? NonBJet1EtaProb : minProb;
   chi2t[NLL_NonBJet1Phi_TOPTOPLEPHAD] = (NonBJet1PhiProb > minProb) ? NonBJet1PhiProb : minProb;
   chi2t[NLL_NonBJet2Px_TOPTOPLEPHAD] = (NonBJet2PxProb > minProb) ? NonBJet2PxProb : minProb;
   chi2t[NLL_NonBJet2Py_TOPTOPLEPHAD] = (NonBJet2PyProb > minProb) ? NonBJet2PyProb : minProb;
   chi2t[NLL_NonBJet2Pz_TOPTOPLEPHAD] = (NonBJet2PzProb > minProb) ? NonBJet2PzProb : minProb;
   chi2t[NLL_NonBJet2Pt_TOPTOPLEPHAD] = (NonBJet2PtProb > minProb) ? NonBJet2PtProb : minProb;
   chi2t[NLL_NonBJet2Eta_TOPTOPLEPHAD] = (NonBJet2EtaProb > minProb) ? NonBJet2EtaProb : minProb;
   chi2t[NLL_NonBJet2Phi_TOPTOPLEPHAD] = (NonBJet2PhiProb > minProb) ? NonBJet2PhiProb : minProb;
   chi2t[NLL_LeptonPx_TOPTOPLEPHAD] = (LeptonPxProb > minProb) ? LeptonPxProb : minProb;
   chi2t[NLL_LeptonPy_TOPTOPLEPHAD] = (LeptonPyProb > minProb) ? LeptonPyProb : minProb;
   chi2t[NLL_LeptonPz_TOPTOPLEPHAD] = (LeptonPzProb > minProb) ? LeptonPzProb : minProb;
   chi2t[NLL_LeptonPt_TOPTOPLEPHAD] = (LeptonPtProb > minProb) ? LeptonPtProb : minProb;
   chi2t[NLL_LeptonEta_TOPTOPLEPHAD] = (LeptonEtaProb > minProb) ? LeptonEtaProb : minProb;
   chi2t[NLL_LeptonPhi_TOPTOPLEPHAD] = (LeptonPhiProb > minProb) ? LeptonPhiProb : minProb;
   chi2t[NLL_PhotonPx_TOPTOPLEPHAD] = (PhotonPxProb > minProb) ? PhotonPxProb : minProb;
   chi2t[NLL_PhotonPy_TOPTOPLEPHAD] = (PhotonPyProb > minProb) ? PhotonPyProb : minProb;
   chi2t[NLL_PhotonPz_TOPTOPLEPHAD] = (PhotonPzProb > minProb) ? PhotonPzProb : minProb;
   chi2t[NLL_PhotonPt_TOPTOPLEPHAD] = (PhotonPtProb > minProb) ? PhotonPtProb : minProb;
   chi2t[NLL_PhotonEta_TOPTOPLEPHAD] = (PhotonEtaProb > minProb) ? PhotonEtaProb : minProb;
   chi2t[NLL_PhotonPhi_TOPTOPLEPHAD] = (PhotonPhiProb > minProb) ? PhotonPhiProb : minProb;
   
   double lh = 1.0;
   
   lh *= chi2t[NLL_WLep_TOPTOPLEPHAD];
   lh *= chi2t[NLL_WHad_TOPTOPLEPHAD];
   lh *= chi2t[NLL_TopLep_TOPTOPLEPHAD];
   lh *= chi2t[NLL_TopHad_TOPTOPLEPHAD];

   /// Include additional terms only if is corresponding parameter is free
   if(! (*IsParFixed)[FPARAM_EtRealX_TOPTOPLEPHAD]) lh *= chi2t[NLL_EtMissX_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_EtRealY_TOPTOPLEPHAD]) lh *= chi2t[NLL_EtMissY_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_BJetLepPx_TOPTOPLEPHAD]) lh *= chi2t[NLL_BJetLepPx_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_BJetLepPy_TOPTOPLEPHAD]) lh *= chi2t[NLL_BJetLepPy_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_BJetLepPz_TOPTOPLEPHAD]) lh *= chi2t[NLL_BJetLepPz_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_BJetLepPt_TOPTOPLEPHAD]) lh *= chi2t[NLL_BJetLepPt_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_BJetLepEta_TOPTOPLEPHAD]) lh *= chi2t[NLL_BJetLepEta_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_BJetLepPhi_TOPTOPLEPHAD]) lh *= chi2t[NLL_BJetLepPhi_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_BJetHadPx_TOPTOPLEPHAD]) lh *= chi2t[NLL_BJetHadPx_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_BJetHadPy_TOPTOPLEPHAD]) lh *= chi2t[NLL_BJetHadPy_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_BJetHadPz_TOPTOPLEPHAD]) lh *= chi2t[NLL_BJetHadPz_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_BJetHadPt_TOPTOPLEPHAD]) lh *= chi2t[NLL_BJetHadPt_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_BJetHadEta_TOPTOPLEPHAD]) lh *= chi2t[NLL_BJetHadEta_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_BJetHadPhi_TOPTOPLEPHAD]) lh *= chi2t[NLL_BJetHadPhi_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_NonBJet1Px_TOPTOPLEPHAD]) lh *= chi2t[NLL_NonBJet1Px_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_NonBJet1Py_TOPTOPLEPHAD]) lh *= chi2t[NLL_NonBJet1Py_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_NonBJet1Pz_TOPTOPLEPHAD]) lh *= chi2t[NLL_NonBJet1Pz_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_NonBJet1Pt_TOPTOPLEPHAD]) lh *= chi2t[NLL_NonBJet1Pt_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_NonBJet1Eta_TOPTOPLEPHAD]) lh *= chi2t[NLL_NonBJet1Eta_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_NonBJet1Phi_TOPTOPLEPHAD]) lh *= chi2t[NLL_NonBJet1Phi_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_NonBJet2Px_TOPTOPLEPHAD]) lh *= chi2t[NLL_NonBJet2Px_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_NonBJet2Py_TOPTOPLEPHAD]) lh *= chi2t[NLL_NonBJet2Py_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_NonBJet2Pz_TOPTOPLEPHAD]) lh *= chi2t[NLL_NonBJet2Pz_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_NonBJet2Pt_TOPTOPLEPHAD]) lh *= chi2t[NLL_NonBJet2Pt_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_NonBJet2Eta_TOPTOPLEPHAD]) lh *= chi2t[NLL_NonBJet2Eta_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_NonBJet2Phi_TOPTOPLEPHAD]) lh *= chi2t[NLL_NonBJet2Phi_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_LeptonPx_TOPTOPLEPHAD]) lh *= chi2t[NLL_LeptonPx_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_LeptonPy_TOPTOPLEPHAD]) lh *= chi2t[NLL_LeptonPy_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_LeptonPz_TOPTOPLEPHAD]) lh *= chi2t[NLL_LeptonPz_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_LeptonPt_TOPTOPLEPHAD]) lh *= chi2t[NLL_LeptonPt_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_LeptonEta_TOPTOPLEPHAD]) lh *= chi2t[NLL_LeptonEta_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_LeptonPhi_TOPTOPLEPHAD]) lh *= chi2t[NLL_LeptonPhi_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_PhotonPx_TOPTOPLEPHAD]) lh *= chi2t[NLL_PhotonPx_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_PhotonPy_TOPTOPLEPHAD]) lh *= chi2t[NLL_PhotonPy_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_PhotonPz_TOPTOPLEPHAD]) lh *= chi2t[NLL_PhotonPz_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_PhotonPt_TOPTOPLEPHAD]) lh *= chi2t[NLL_PhotonPt_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_PhotonEta_TOPTOPLEPHAD]) lh *= chi2t[NLL_PhotonEta_TOPTOPLEPHAD];
   if(! (*IsParFixed)[FPARAM_PhotonPhi_TOPTOPLEPHAD]) lh *= chi2t[NLL_PhotonPhi_TOPTOPLEPHAD];

   val += -2.*log(lh);

   return val;
}

/// Calculate output variables for a given permutation
void KINFIT::TopTopLepHad::calcVar(int iPerm)
{
   if( chi_[iPerm] > 10E+9 ) return;

   int idxElec = TopTopLepHad_ElectronIdx[iPerm];
   int idxMuon = TopTopLepHad_MuonIdx[iPerm];
   
   int idxBJetLep = TopTopLepHad_BJetLepIdx[iPerm];
   int idxBJetHad = TopTopLepHad_BJetHadIdx[iPerm];   

   int idxNonBJet1 = TopTopLepHad_NonBJet1Idx[iPerm];
   int idxNonBJet2 = TopTopLepHad_NonBJet2Idx[iPerm];
   
   int idxPhoton = TopTopLepHad_PhotonIdx[iPerm];
   int photonOrigin = PhotonOrigin_[iPerm];
   
   float LeptonPx = 0;
   float LeptonPy = 0;
   float LeptonPz = 0;
   float LeptonE = 0;

   float BJetLepPx = 0;
   float BJetLepPy = 0;
   float BJetLepPz = 0;
   float BJetLepE = 0;

   float BJetHadPx = 0;
   float BJetHadPy = 0;
   float BJetHadPz = 0;
   float BJetHadE = 0;

   float NonBJet1Px = 0;
   float NonBJet1Py = 0;
   float NonBJet1Pz = 0;
   float NonBJet1E = 0;

   float NonBJet2Px = 0;
   float NonBJet2Py = 0;
   float NonBJet2Pz = 0;
   float NonBJet2E = 0;

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
   
   if( idxElec >= 0 )
     {
	LeptonPx = ElectronPx[idxElec];
	LeptonPy = ElectronPy[idxElec];
	LeptonPz = ElectronPz[idxElec];
	LeptonE = ElectronE[idxElec];
     }   
   else if( idxMuon >= 0 )
     {
	LeptonPx = MuonPx[idxMuon];
	LeptonPy = MuonPy[idxMuon];
	LeptonPz = MuonPz[idxMuon];
	LeptonE = MuonE[idxMuon];
     }
   
   if( photonOrigin == PHOTON_FROM_LEPTON_TOPTOPLEPHAD && (idxElec >= 0 || idxMuon >= 0) )
     {
	LeptonPx += RadPhotonPx;
	LeptonPy += RadPhotonPy;
	LeptonPz += RadPhotonPz;
	LeptonE += RadPhotonE;
     }

   if( idxBJetLep >= 0 )
     {
	BJetLepPx = BJetPx[idxBJetLep];
	BJetLepPy = BJetPy[idxBJetLep];
	BJetLepPz = BJetPz[idxBJetLep];
	BJetLepE = BJetE[idxBJetLep];
	
	if( photonOrigin == PHOTON_FROM_BJETLEP_TOPTOPLEPHAD )
	  {
	     BJetLepPx += RadPhotonPx;
	     BJetLepPy += RadPhotonPy;
	     BJetLepPz += RadPhotonPz;
	     BJetLepE += RadPhotonE;
	  }
     }   
   if( idxBJetHad >= 0 )
     {
	BJetHadPx = BJetPx[idxBJetHad];
	BJetHadPy = BJetPy[idxBJetHad];
	BJetHadPz = BJetPz[idxBJetHad];
	BJetHadE = BJetE[idxBJetHad];

	if( photonOrigin == PHOTON_FROM_BJETHAD_TOPTOPLEPHAD )
	  {
	     BJetHadPx += RadPhotonPx;
	     BJetHadPy += RadPhotonPy;
	     BJetHadPz += RadPhotonPz;
	     BJetHadE += RadPhotonE;
	  }
     }   

   if( idxNonBJet1 >= 0 )
     {
	NonBJet1Px = NonBJetPx[idxNonBJet1];
	NonBJet1Py = NonBJetPy[idxNonBJet1];
	NonBJet1Pz = NonBJetPz[idxNonBJet1];
	NonBJet1E = NonBJetE[idxNonBJet1];
	
	if( photonOrigin == PHOTON_FROM_NONBJET1_TOPTOPLEPHAD )
	  {
	     NonBJet1Px += RadPhotonPx;
	     NonBJet1Py += RadPhotonPy;
	     NonBJet1Pz += RadPhotonPz;
	     NonBJet1E += RadPhotonE;
	  }
     }   
   if( idxNonBJet2 >= 0 )
     {
	NonBJet2Px = NonBJetPx[idxNonBJet2];
	NonBJet2Py = NonBJetPy[idxNonBJet2];
	NonBJet2Pz = NonBJetPz[idxNonBJet2];
	NonBJet2E = NonBJetE[idxNonBJet2];

	if( photonOrigin == PHOTON_FROM_NONBJET2_TOPTOPLEPHAD )
	  {
	     NonBJet2Px += RadPhotonPx;
	     NonBJet2Py += RadPhotonPy;
	     NonBJet2Pz += RadPhotonPz;
	     NonBJet2E += RadPhotonE;
	  }
     }
   
   float nuPx = nuPx_[iPerm][0];
   float nuPy = nuPy_[iPerm][0];
   float nuPz = nuPz_[iPerm][0];
   
   float nuE = sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz);

   /// Leptonic W boson
   float WLepE = nuE+LeptonE;
   float WLepPx = nuPx+LeptonPx;
   float WLepPy = nuPy+LeptonPy;
   float WLepPz = nuPz+LeptonPz;
   
   if( photonOrigin == PHOTON_FROM_WLEP_TOPTOPLEPHAD ||
       photonOrigin == PHOTON_FROM_WLEP_COMB_TOPTOPLEPHAD )
     {
	WLepPx += RadPhotonPx;
	WLepPy += RadPhotonPy;
	WLepPz += RadPhotonPz;
	WLepE += RadPhotonE;
     }
   
   float WLepPt = sqrt(WLepPx*WLepPx+WLepPy*WLepPy);
   float WLepEta = getEta(WLepPt, WLepPz);
   float WLepPhi = atan2(WLepPy, WLepPx);
   
   /// Leptonic top quark
   float topLepE = WLepE+BJetLepE;
   float topLepPx = WLepPx+BJetLepPx;
   float topLepPy = WLepPy+BJetLepPy;
   float topLepPz = WLepPz+BJetLepPz;

   if( photonOrigin == PHOTON_FROM_TOPLEP_TOPTOPLEPHAD ||
       photonOrigin == PHOTON_FROM_TOPLEP_COMB_TOPTOPLEPHAD )
     {
	topLepPx += RadPhotonPx;
	topLepPy += RadPhotonPy;
	topLepPz += RadPhotonPz;
	topLepE += RadPhotonE;
     }
   
   float topLepPt = sqrt(topLepPx*topLepPx+topLepPy*topLepPy);
   float topLepEta = getEta(topLepPt, topLepPz);
   float topLepPhi = atan2(topLepPy, topLepPx);

   /// Hadronic W boson
   float WHadE = NonBJet1E+NonBJet2E;
   float WHadPx = NonBJet1Px+NonBJet2Px;
   float WHadPy = NonBJet1Py+NonBJet2Py;
   float WHadPz = NonBJet1Pz+NonBJet2Pz;

   if( photonOrigin == PHOTON_FROM_WHAD_TOPTOPLEPHAD ||
       photonOrigin == PHOTON_FROM_WHAD_COMB_TOPTOPLEPHAD )
     {
	WHadPx += RadPhotonPx;
	WHadPy += RadPhotonPy;
	WHadPz += RadPhotonPz;
	WHadE += RadPhotonE;
     }
   
   float WHadPt = sqrt(WHadPx*WHadPx+WHadPy*WHadPy);
   float WHadEta = getEta(WHadPt, WHadPz);
   float WHadPhi = atan2(WHadPy, WHadPx);

   /// Hadronic top quark
   float topHadE = WHadE+BJetHadE;
   float topHadPx = WHadPx+BJetHadPx;
   float topHadPy = WHadPy+BJetHadPy;
   float topHadPz = WHadPz+BJetHadPz;

   if( photonOrigin == PHOTON_FROM_TOPHAD_TOPTOPLEPHAD ||
       photonOrigin == PHOTON_FROM_TOPHAD_COMB_TOPTOPLEPHAD )
     {
	topHadPx += RadPhotonPx;
	topHadPy += RadPhotonPy;
	topHadPz += RadPhotonPz;
	topHadE += RadPhotonE;
     }
   
   float topHadPt = sqrt(topHadPx*topHadPx+topHadPy*topHadPy);
   float topHadEta = getEta(topHadPt, topHadPz);
   float topHadPhi = atan2(topHadPy, topHadPx);
   
   /// Di-top system
   float toptopE = topLepE+topHadE;
   float toptopPx = topLepPx+topHadPx;
   float toptopPy = topLepPy+topHadPy;
   float toptopPz = topLepPz+topHadPz;

   /// Fill the outputs
   
   drTopTop_[iPerm] = getDeltaR(topLepEta, topLepPhi, topHadEta, topHadPhi);
   detaTopTop_[iPerm] = topLepEta-topHadEta;
   dphiTopTop_[iPerm] = getDeltaPhi(topLepPhi, topHadPhi);

   mTopTop_[iPerm] = sqrt(toptopE*toptopE-toptopPx*toptopPx-toptopPy*toptopPy-toptopPz*toptopPz);
   ptTopTop_[iPerm] = sqrt(toptopPx*toptopPx+toptopPy*toptopPy);
   pTopTop_[iPerm] = sqrt(toptopPx*toptopPx+toptopPy*toptopPy+toptopPz*toptopPz);
   etaTopTop_[iPerm] = getEta(ptTopTop_[iPerm], toptopPz);
   phiTopTop_[iPerm] = atan2(toptopPy, toptopPx);

   TopMass_[iPerm][0] = sqrt(topLepE*topLepE-topLepPx*topLepPx-topLepPy*topLepPy-topLepPz*topLepPz);
   TopMass_[iPerm][1] = sqrt(topHadE*topHadE-topHadPx*topHadPx-topHadPy*topHadPy-topHadPz*topHadPz);   
   TopPt_[iPerm][0] = sqrt(topLepPx*topLepPx+topLepPy*topLepPy);
   TopPt_[iPerm][1] = sqrt(topHadPx*topHadPx+topHadPy*topHadPy);
   TopP_[iPerm][0] = sqrt(topLepPx*topLepPx+topLepPy*topLepPy+topLepPz*topLepPz);
   TopP_[iPerm][1] = sqrt(topHadPx*topHadPx+topHadPy*topHadPy+topHadPz*topHadPz);
   TopEta_[iPerm][0] = getEta(TopPt_[iPerm][0],topLepPz);
   TopEta_[iPerm][1] = getEta(TopPt_[iPerm][1],topHadPz);
   TopPhi_[iPerm][0] = topLepPhi;
   TopPhi_[iPerm][1] = topHadPhi;
   TopE_[iPerm][0] = topLepE;
   TopE_[iPerm][1] = topHadE;
   TopPx_[iPerm][0] = topLepPx;
   TopPx_[iPerm][1] = topHadPx;
   TopPy_[iPerm][0] = topLepPy;
   TopPy_[iPerm][1] = topHadPy;
   TopPz_[iPerm][0] = topLepPz;
   TopPz_[iPerm][1] = topHadPz;

   WMass_[iPerm][0] = sqrt(WLepE*WLepE-WLepPx*WLepPx-WLepPy*WLepPy-WLepPz*WLepPz);
   WMass_[iPerm][1] = sqrt(WHadE*WHadE-WHadPx*WHadPx-WHadPy*WHadPy-WHadPz*WHadPz);
   WPt_[iPerm][0] = sqrt(WLepPx*WLepPx+WLepPy*WLepPy);
   WPt_[iPerm][1] = sqrt(WHadPx*WHadPx+WHadPy*WHadPy);
   WP_[iPerm][0] = sqrt(WLepPx*WLepPx+WLepPy*WLepPy+WLepPz*WLepPz);
   WP_[iPerm][1] = sqrt(WHadPx*WHadPx+WHadPy*WHadPy+WHadPz*WHadPz);
   WEta_[iPerm][0] = getEta(WPt_[iPerm][0],WLepPz);
   WEta_[iPerm][1] = getEta(WPt_[iPerm][1],WHadPz);
   WPhi_[iPerm][0] = atan2(WLepPy, WLepPx);
   WPhi_[iPerm][1] = atan2(WHadPy, WHadPx);
   WE_[iPerm][0] = WLepE;
   WE_[iPerm][1] = WHadE;
   WPx_[iPerm][0] = WLepPx;
   WPx_[iPerm][1] = WHadPx;
   WPy_[iPerm][0] = WLepPy;
   WPy_[iPerm][1] = WHadPy;
   WPz_[iPerm][0] = WLepPz;
   WPz_[iPerm][1] = WHadPz;
}

/// FCN definition
void KINFIT::TopTopLepHad::fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
   std::string chi2tNames[NLL_N_TOPTOPLEPHAD] = {"WLep", "WHad", "TopLep", "TopHad", "EtMissX", "EtMissY", "BJetLepPx", "BJetLepPy", "BJetLepPz", "BJetLepPt", "BJetLepEta", "BJetLepPhi", "BJetHadPx", "BJetHadPy", "BJetHadPz", "BJetHadPt", "BJetHadEta", "BJetHadPhi", "NonBJet1Px", "NonBJet1Py", "NonBJet1Pz", "NonBJet1Pt", "NonBJet1Eta", "NonBJet1Phi", "NonBJet2Px", "NonBJet2Py", "NonBJet2Pz", "NonBJet2Pt", "NonBJet2Eta", "NonBJet2Phi", "LeptonPx", "LeptonPy", "LeptonPz", "LeptonPt", "LeptonEta", "LeptonPhi", "RacPT", "mWPT", "mTopPT"};
   std::vector<double> chi2t(NLL_N_TOPTOPLEPHAD);
   
   double lh = func(*PxLepton1, *PyLepton1, *PzLepton1, *PtLepton1, *EtaLepton1, *PhiLepton1, *ELepton1, *LabelLepton1,
		    *PxBJet1, *PyBJet1, *PzBJet1, *PtBJet1, *EtaBJet1, *PhiBJet1, *EBJet1,
		    *PxBJet2, *PyBJet2, *PzBJet2, *PtBJet2, *EtaBJet2, *PhiBJet2, *EBJet2,
		    *PxNonBJet1, *PyNonBJet1, *PzNonBJet1, *PtNonBJet1, *EtaNonBJet1, *PhiNonBJet1, *ENonBJet1,
		    *PxNonBJet2, *PyNonBJet2, *PzNonBJet2, *PtNonBJet2, *EtaNonBJet2, *PhiNonBJet2, *ENonBJet2,
		    *PxPhoton, *PyPhoton, *PzPhoton, *PtPhoton, *EtaPhoton, *PhiPhoton, *EPhoton,
		    *PhotonOrigin,
		    chi2t, par);

   /// converged
   if( iflag == 3 )
     {
	*CHISQ = lh;

	for( int ip=0;ip<FPARAM_N;ip++ ) (*FitParam).push_back(par[ip]);
	for( int it=0;it<NLL_N_TOPTOPLEPHAD;it++ ) {(*ChiTerm).push_back(chi2t[it]); (*ChiTermName).push_back(chi2tNames[it]);}
     }
   
   f = lh;
}

/// Run the fit
void KINFIT::TopTopLepHad::fit(double *par, std::vector<FRESULT> &vpp)
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
   (*ParMin)[FPARAM_BJetLepE_TOPTOPLEPHAD] = par[FPARAM_BJetLepE_TOPTOPLEPHAD]-1000.; (*ParMax)[FPARAM_BJetLepE_TOPTOPLEPHAD] = par[FPARAM_BJetLepE_TOPTOPLEPHAD]+1000.;
   (*ParMin)[FPARAM_BJetHadE_TOPTOPLEPHAD] = par[FPARAM_BJetHadE_TOPTOPLEPHAD]-1000.; (*ParMax)[FPARAM_BJetHadE_TOPTOPLEPHAD] = par[FPARAM_BJetHadE_TOPTOPLEPHAD]+1000.;
   (*ParMin)[FPARAM_NonBJet1E_TOPTOPLEPHAD] = par[FPARAM_NonBJet1E_TOPTOPLEPHAD]-1000.; (*ParMax)[FPARAM_NonBJet1E_TOPTOPLEPHAD] = par[FPARAM_NonBJet1E_TOPTOPLEPHAD]+1000.;
   (*ParMin)[FPARAM_NonBJet2E_TOPTOPLEPHAD] = par[FPARAM_NonBJet2E_TOPTOPLEPHAD]-1000.; (*ParMax)[FPARAM_NonBJet2E_TOPTOPLEPHAD] = par[FPARAM_NonBJet2E_TOPTOPLEPHAD]+1000.;
   (*ParMin)[FPARAM_LeptonE_TOPTOPLEPHAD] = par[FPARAM_LeptonE_TOPTOPLEPHAD]-1000.; (*ParMax)[FPARAM_LeptonE_TOPTOPLEPHAD] = par[FPARAM_LeptonE_TOPTOPLEPHAD]+1000.;
   (*ParMin)[FPARAM_PhotonE_TOPTOPLEPHAD] = par[FPARAM_PhotonE_TOPTOPLEPHAD]-1000.; (*ParMax)[FPARAM_PhotonE_TOPTOPLEPHAD] = par[FPARAM_PhotonE_TOPTOPLEPHAD]+1000.;
      
   gMinuit->mnparm(FPARAM_Sign_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_Sign_TOPTOPLEPHAD], par[FPARAM_Sign_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_Sign_TOPTOPLEPHAD], (*ParMax)[FPARAM_Sign_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_EtRealX_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_EtRealX_TOPTOPLEPHAD], par[FPARAM_EtRealX_TOPTOPLEPHAD], 1E-2, (*ParMin)[FPARAM_EtRealX_TOPTOPLEPHAD], (*ParMax)[FPARAM_EtRealX_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_EtRealY_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_EtRealY_TOPTOPLEPHAD], par[FPARAM_EtRealY_TOPTOPLEPHAD], 1E-2, (*ParMin)[FPARAM_EtRealY_TOPTOPLEPHAD], (*ParMax)[FPARAM_EtRealY_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_mWLep_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_mWLep_TOPTOPLEPHAD], par[FPARAM_mWLep_TOPTOPLEPHAD], 1E-4, (*ParMin)[FPARAM_mWLep_TOPTOPLEPHAD], (*ParMax)[FPARAM_mWLep_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_mWHad_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_mWHad_TOPTOPLEPHAD], par[FPARAM_mWHad_TOPTOPLEPHAD], 1E-4, (*ParMin)[FPARAM_mWHad_TOPTOPLEPHAD], (*ParMax)[FPARAM_mWHad_TOPTOPLEPHAD], ierflg);

   gMinuit->mnparm(FPARAM_BJetLepPx_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetLepPx_TOPTOPLEPHAD], par[FPARAM_BJetLepPx_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetLepPx_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetLepPx_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetLepPy_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetLepPy_TOPTOPLEPHAD], par[FPARAM_BJetLepPy_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetLepPy_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetLepPy_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetLepPz_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetLepPz_TOPTOPLEPHAD], par[FPARAM_BJetLepPz_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetLepPz_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetLepPz_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetLepPt_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetLepPt_TOPTOPLEPHAD], par[FPARAM_BJetLepPt_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetLepPt_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetLepPt_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetLepEta_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetLepEta_TOPTOPLEPHAD], par[FPARAM_BJetLepEta_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetLepEta_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetLepEta_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetLepPhi_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetLepPhi_TOPTOPLEPHAD], par[FPARAM_BJetLepPhi_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetLepPhi_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetLepPhi_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetLepE_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetLepE_TOPTOPLEPHAD], par[FPARAM_BJetLepE_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetLepE_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetLepE_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetHadPx_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetHadPx_TOPTOPLEPHAD], par[FPARAM_BJetHadPx_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetHadPx_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetHadPx_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetHadPy_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetHadPy_TOPTOPLEPHAD], par[FPARAM_BJetHadPy_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetHadPy_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetHadPy_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetHadPz_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetHadPz_TOPTOPLEPHAD], par[FPARAM_BJetHadPz_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetHadPz_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetHadPz_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetHadPt_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetHadPt_TOPTOPLEPHAD], par[FPARAM_BJetHadPt_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetHadPt_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetHadPt_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetHadEta_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetHadEta_TOPTOPLEPHAD], par[FPARAM_BJetHadEta_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetHadEta_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetHadEta_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetHadPhi_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetHadPhi_TOPTOPLEPHAD], par[FPARAM_BJetHadPhi_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetHadPhi_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetHadPhi_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_BJetHadE_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_BJetHadE_TOPTOPLEPHAD], par[FPARAM_BJetHadE_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_BJetHadE_TOPTOPLEPHAD], (*ParMax)[FPARAM_BJetHadE_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet1Px_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet1Px_TOPTOPLEPHAD], par[FPARAM_NonBJet1Px_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet1Px_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet1Px_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet1Py_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet1Py_TOPTOPLEPHAD], par[FPARAM_NonBJet1Py_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet1Py_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet1Py_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet1Pz_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet1Pz_TOPTOPLEPHAD], par[FPARAM_NonBJet1Pz_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet1Pz_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet1Pz_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet1Pt_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet1Pt_TOPTOPLEPHAD], par[FPARAM_NonBJet1Pt_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet1Pt_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet1Pt_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet1Eta_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet1Eta_TOPTOPLEPHAD], par[FPARAM_NonBJet1Eta_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet1Eta_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet1Eta_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet1Phi_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet1Phi_TOPTOPLEPHAD], par[FPARAM_NonBJet1Phi_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet1Phi_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet1Phi_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet1E_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet1E_TOPTOPLEPHAD], par[FPARAM_NonBJet1E_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet1E_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet1E_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet2Px_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet2Px_TOPTOPLEPHAD], par[FPARAM_NonBJet2Px_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet2Px_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet2Px_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet2Py_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet2Py_TOPTOPLEPHAD], par[FPARAM_NonBJet2Py_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet2Py_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet2Py_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet2Pz_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet2Pz_TOPTOPLEPHAD], par[FPARAM_NonBJet2Pz_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet2Pz_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet2Pz_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet2Pt_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet2Pt_TOPTOPLEPHAD], par[FPARAM_NonBJet2Pt_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet2Pt_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet2Pt_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet2Eta_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet2Eta_TOPTOPLEPHAD], par[FPARAM_NonBJet2Eta_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet2Eta_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet2Eta_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet2Phi_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet2Phi_TOPTOPLEPHAD], par[FPARAM_NonBJet2Phi_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet2Phi_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet2Phi_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_NonBJet2E_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_NonBJet2E_TOPTOPLEPHAD], par[FPARAM_NonBJet2E_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_NonBJet2E_TOPTOPLEPHAD], (*ParMax)[FPARAM_NonBJet2E_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_LeptonPx_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_LeptonPx_TOPTOPLEPHAD], par[FPARAM_LeptonPx_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_LeptonPx_TOPTOPLEPHAD], (*ParMax)[FPARAM_LeptonPx_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_LeptonPy_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_LeptonPy_TOPTOPLEPHAD], par[FPARAM_LeptonPy_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_LeptonPy_TOPTOPLEPHAD], (*ParMax)[FPARAM_LeptonPy_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_LeptonPz_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_LeptonPz_TOPTOPLEPHAD], par[FPARAM_LeptonPz_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_LeptonPz_TOPTOPLEPHAD], (*ParMax)[FPARAM_LeptonPz_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_LeptonPt_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_LeptonPt_TOPTOPLEPHAD], par[FPARAM_LeptonPt_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_LeptonPt_TOPTOPLEPHAD], (*ParMax)[FPARAM_LeptonPt_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_LeptonEta_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_LeptonEta_TOPTOPLEPHAD], par[FPARAM_LeptonEta_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_LeptonEta_TOPTOPLEPHAD], (*ParMax)[FPARAM_LeptonEta_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_LeptonPhi_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_LeptonPhi_TOPTOPLEPHAD], par[FPARAM_LeptonPhi_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_LeptonPhi_TOPTOPLEPHAD], (*ParMax)[FPARAM_LeptonPhi_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_LeptonE_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_LeptonE_TOPTOPLEPHAD], par[FPARAM_LeptonE_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_LeptonE_TOPTOPLEPHAD], (*ParMax)[FPARAM_LeptonE_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_PhotonPx_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_PhotonPx_TOPTOPLEPHAD], par[FPARAM_PhotonPx_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_PhotonPx_TOPTOPLEPHAD], (*ParMax)[FPARAM_PhotonPx_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_PhotonPy_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_PhotonPy_TOPTOPLEPHAD], par[FPARAM_PhotonPy_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_PhotonPy_TOPTOPLEPHAD], (*ParMax)[FPARAM_PhotonPy_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_PhotonPz_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_PhotonPz_TOPTOPLEPHAD], par[FPARAM_PhotonPz_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_PhotonPz_TOPTOPLEPHAD], (*ParMax)[FPARAM_PhotonPz_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_PhotonPt_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_PhotonPt_TOPTOPLEPHAD], par[FPARAM_PhotonPt_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_PhotonPt_TOPTOPLEPHAD], (*ParMax)[FPARAM_PhotonPt_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_PhotonEta_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_PhotonEta_TOPTOPLEPHAD], par[FPARAM_PhotonEta_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_PhotonEta_TOPTOPLEPHAD], (*ParMax)[FPARAM_PhotonEta_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_PhotonPhi_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_PhotonPhi_TOPTOPLEPHAD], par[FPARAM_PhotonPhi_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_PhotonPhi_TOPTOPLEPHAD], (*ParMax)[FPARAM_PhotonPhi_TOPTOPLEPHAD], ierflg);
   gMinuit->mnparm(FPARAM_PhotonE_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_PhotonE_TOPTOPLEPHAD], par[FPARAM_PhotonE_TOPTOPLEPHAD], 1E-1, (*ParMin)[FPARAM_PhotonE_TOPTOPLEPHAD], (*ParMax)[FPARAM_PhotonE_TOPTOPLEPHAD], ierflg);
   
   for( int i=0;i<FPARAM_N;i++ ) if( (*IsParFixed)[i] ) gMinuit->FixParameter(i);
   
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
   fres.PzNuSum = fabs(*PzNu1);
   
   float PtNu1 = sqrt((*PxNu1)*(*PxNu1) + (*PyNu1)*(*PyNu1));
   float PNu1 = sqrt((PtNu1)*(PtNu1) + (*PzNu1)*(*PzNu1));
   
   fres.PtNuSum = PtNu1;
   fres.PNuSum = PNu1;
   
   fres.TopTopMass = (*TopTopMass);
   fres.TopTopPz = (*TopTopPz);
   fres.TopTopPzAbs = (*TopTopPzAbs);
   fres.lh = *CHISQ;
   
   vpp.push_back(fres);
   
   delete gMinuit;
}

/// First-layer minimization (generic)
void KINFIT::TopTopLepHad::calcNuGrid(std::vector<FRESULT> &vp)
{
   double par[FPARAM_N];
   for(int p=0;p<FPARAM_N;p++) par[p] = 0.;
   
   std::vector<double> chi2t(NLL_N_TOPTOPLEPHAD);

   bool doToys = bool(NToy_ > 1);
   
   for(int it=0;it<NToy_;it++)
     {		     
	if( vp.size() >= NGrid_ ) break;
	
	double thres = -1E-6;

	if( *usePDFBJetPxPyPz )
	  {
	     do
	       {
		  if( ! isDeltaFuncPDFBJetPx ) {		       
		     double gv = (doToys) ? getProbGaus(hPDFBJetPx.get(), maxPDFBJetPx, meanPDFBJetPx, sigmaPDFBJetPx, rnd, NBJetPxRMS_) : 0.;
		     par[FPARAM_BJetLepPx_TOPTOPLEPHAD] = *PxBJet1/(1.-gv);
		  } else par[FPARAM_BJetLepPx_TOPTOPLEPHAD] = *PxBJet1;		  

		  if( ! isDeltaFuncPDFBJetPy ) {		       
		     double gv = (doToys) ? getProbGaus(hPDFBJetPy.get(), maxPDFBJetPy, meanPDFBJetPy, sigmaPDFBJetPy, rnd, NBJetPyRMS_) : 0.;
		     par[FPARAM_BJetLepPy_TOPTOPLEPHAD] = *PyBJet1/(1.-gv);
		  } else par[FPARAM_BJetLepPy_TOPTOPLEPHAD] = *PyBJet1;		  
		  
		  if( ! isDeltaFuncPDFBJetPz ) {		       
		     double gv = (doToys) ? getProbGaus(hPDFBJetPz.get(), maxPDFBJetPz, meanPDFBJetPz, sigmaPDFBJetPz, rnd, NBJetPzRMS_) : 0.;
		     par[FPARAM_BJetLepPz_TOPTOPLEPHAD] = *PzBJet1/(1.-gv);
		  } else par[FPARAM_BJetLepPz_TOPTOPLEPHAD] = *PzBJet1;		  
		  
		  par[FPARAM_BJetLepPt_TOPTOPLEPHAD] = sqrt(par[FPARAM_BJetLepPx_TOPTOPLEPHAD]*par[FPARAM_BJetLepPx_TOPTOPLEPHAD] + par[FPARAM_BJetLepPy_TOPTOPLEPHAD]*par[FPARAM_BJetLepPy_TOPTOPLEPHAD]);
		  par[FPARAM_BJetLepEta_TOPTOPLEPHAD] = getEta(par[FPARAM_BJetLepPt_TOPTOPLEPHAD], par[FPARAM_BJetLepPz_TOPTOPLEPHAD]);
		  par[FPARAM_BJetLepPhi_TOPTOPLEPHAD] = atan2(par[FPARAM_BJetLepPy_TOPTOPLEPHAD], par[FPARAM_BJetLepPx_TOPTOPLEPHAD]);
		  par[FPARAM_BJetLepE_TOPTOPLEPHAD] = sqrt((*MassBJet1)*(*MassBJet1) + par[FPARAM_BJetLepPx_TOPTOPLEPHAD]*par[FPARAM_BJetLepPx_TOPTOPLEPHAD] + par[FPARAM_BJetLepPy_TOPTOPLEPHAD]*par[FPARAM_BJetLepPy_TOPTOPLEPHAD] + par[FPARAM_BJetLepPz_TOPTOPLEPHAD]*par[FPARAM_BJetLepPz_TOPTOPLEPHAD]);
	       } while( (par[FPARAM_BJetLepE_TOPTOPLEPHAD]*par[FPARAM_BJetLepE_TOPTOPLEPHAD]-par[FPARAM_BJetLepPx_TOPTOPLEPHAD]*par[FPARAM_BJetLepPx_TOPTOPLEPHAD]-par[FPARAM_BJetLepPy_TOPTOPLEPHAD]*par[FPARAM_BJetLepPy_TOPTOPLEPHAD]-par[FPARAM_BJetLepPz_TOPTOPLEPHAD]*par[FPARAM_BJetLepPz_TOPTOPLEPHAD]) < thres );
	  }
	else {
	   do
	     {
		if( ! isDeltaFuncPDFBJetPt ) {
		   double gv = (doToys) ? getProbGaus(hPDFBJetPt.get(), maxPDFBJetPt, meanPDFBJetPt, sigmaPDFBJetPt, rnd, NBJetPtRMS_) : 0.;
		   par[FPARAM_BJetLepPt_TOPTOPLEPHAD] = *PtBJet1/(1.-gv);
		} else par[FPARAM_BJetLepPt_TOPTOPLEPHAD] = *PtBJet1;
		
		if( ! isDeltaFuncPDFBJetEta ) {		       
		   double gv = (doToys) ? getProbGaus(hPDFBJetEta.get(), maxPDFBJetEta, meanPDFBJetEta, sigmaPDFBJetEta, rnd, NBJetEtaRMS_) : 0.;
		   par[FPARAM_BJetLepEta_TOPTOPLEPHAD] = *EtaBJet1/(1.-gv);
		} else par[FPARAM_BJetLepEta_TOPTOPLEPHAD] = *EtaBJet1;		  
		  
		if( ! isDeltaFuncPDFBJetPhi ) {
		   double gv = (doToys) ? getProbGaus(hPDFBJetPhi.get(), maxPDFBJetPhi, meanPDFBJetPhi, sigmaPDFBJetPhi, rnd, NBJetPhiRMS_) : 0.;
		   par[FPARAM_BJetLepPhi_TOPTOPLEPHAD] = *PhiBJet1/(1.-gv);
		} else par[FPARAM_BJetLepPhi_TOPTOPLEPHAD] = *PhiBJet1;		  
		
		par[FPARAM_BJetLepPx_TOPTOPLEPHAD] = par[FPARAM_BJetLepPt_TOPTOPLEPHAD]*cos(par[FPARAM_BJetLepPhi_TOPTOPLEPHAD]);
		par[FPARAM_BJetLepPy_TOPTOPLEPHAD] = par[FPARAM_BJetLepPt_TOPTOPLEPHAD]*sin(par[FPARAM_BJetLepPhi_TOPTOPLEPHAD]);
		par[FPARAM_BJetLepPz_TOPTOPLEPHAD] = par[FPARAM_BJetLepPt_TOPTOPLEPHAD]*sinh(par[FPARAM_BJetLepEta_TOPTOPLEPHAD]);
		par[FPARAM_BJetLepE_TOPTOPLEPHAD] = sqrt((*MassBJet1)*(*MassBJet1) + par[FPARAM_BJetLepPx_TOPTOPLEPHAD]*par[FPARAM_BJetLepPx_TOPTOPLEPHAD] + par[FPARAM_BJetLepPy_TOPTOPLEPHAD]*par[FPARAM_BJetLepPy_TOPTOPLEPHAD] + par[FPARAM_BJetLepPz_TOPTOPLEPHAD]*par[FPARAM_BJetLepPz_TOPTOPLEPHAD]);
	     } while( (par[FPARAM_BJetLepE_TOPTOPLEPHAD]*par[FPARAM_BJetLepE_TOPTOPLEPHAD]-par[FPARAM_BJetLepPx_TOPTOPLEPHAD]*par[FPARAM_BJetLepPx_TOPTOPLEPHAD]-par[FPARAM_BJetLepPy_TOPTOPLEPHAD]*par[FPARAM_BJetLepPy_TOPTOPLEPHAD]-par[FPARAM_BJetLepPz_TOPTOPLEPHAD]*par[FPARAM_BJetLepPz_TOPTOPLEPHAD]) < thres );	   
	}

	if( *usePDFBJetPxPyPz )
	  {
	     do
	       {
		  if( ! isDeltaFuncPDFBJetPx ) {
		     double gv = (doToys) ? getProbGaus(hPDFBJetPx.get(), maxPDFBJetPx, meanPDFBJetPx, sigmaPDFBJetPx, rnd, NBJetPxRMS_) : 0.;
		     par[FPARAM_BJetHadPx_TOPTOPLEPHAD] = *PxBJet2/(1.-gv);
		  } else par[FPARAM_BJetHadPx_TOPTOPLEPHAD] = *PxBJet2;		  

		  if( ! isDeltaFuncPDFBJetPy ) {		       
		     double gv = (doToys) ? getProbGaus(hPDFBJetPy.get(), maxPDFBJetPy, meanPDFBJetPy, sigmaPDFBJetPy, rnd, NBJetPyRMS_) : 0.;
		     par[FPARAM_BJetHadPy_TOPTOPLEPHAD] = *PyBJet2/(1.-gv);
		  } else par[FPARAM_BJetHadPy_TOPTOPLEPHAD] = *PyBJet2;		  
		  
		  if( ! isDeltaFuncPDFBJetPz ) {		       
		     double gv = (doToys) ? getProbGaus(hPDFBJetPz.get(), maxPDFBJetPz, meanPDFBJetPz, sigmaPDFBJetPz, rnd, NBJetPzRMS_) : 0.;
		     par[FPARAM_BJetHadPz_TOPTOPLEPHAD] = *PzBJet2/(1.-gv);
		  } else par[FPARAM_BJetHadPz_TOPTOPLEPHAD] = *PzBJet2;		  
		  
		  par[FPARAM_BJetHadPt_TOPTOPLEPHAD] = sqrt(par[FPARAM_BJetHadPx_TOPTOPLEPHAD]*par[FPARAM_BJetHadPx_TOPTOPLEPHAD] + par[FPARAM_BJetHadPy_TOPTOPLEPHAD]*par[FPARAM_BJetHadPy_TOPTOPLEPHAD]);
		  par[FPARAM_BJetHadEta_TOPTOPLEPHAD] = getEta(par[FPARAM_BJetHadPt_TOPTOPLEPHAD], par[FPARAM_BJetHadPz_TOPTOPLEPHAD]);
		  par[FPARAM_BJetHadPhi_TOPTOPLEPHAD] = atan2(par[FPARAM_BJetHadPy_TOPTOPLEPHAD], par[FPARAM_BJetHadPx_TOPTOPLEPHAD]);
		  par[FPARAM_BJetHadE_TOPTOPLEPHAD] = sqrt((*MassBJet2)*(*MassBJet2) + par[FPARAM_BJetHadPx_TOPTOPLEPHAD]*par[FPARAM_BJetHadPx_TOPTOPLEPHAD] + par[FPARAM_BJetHadPy_TOPTOPLEPHAD]*par[FPARAM_BJetHadPy_TOPTOPLEPHAD] + par[FPARAM_BJetHadPz_TOPTOPLEPHAD]*par[FPARAM_BJetHadPz_TOPTOPLEPHAD]);
	       } while( (par[FPARAM_BJetHadE_TOPTOPLEPHAD]*par[FPARAM_BJetHadE_TOPTOPLEPHAD]-par[FPARAM_BJetHadPx_TOPTOPLEPHAD]*par[FPARAM_BJetHadPx_TOPTOPLEPHAD]-par[FPARAM_BJetHadPy_TOPTOPLEPHAD]*par[FPARAM_BJetHadPy_TOPTOPLEPHAD]-par[FPARAM_BJetHadPz_TOPTOPLEPHAD]*par[FPARAM_BJetHadPz_TOPTOPLEPHAD]) < thres );
	  }
	else {
	   do
	     {
		if( ! isDeltaFuncPDFBJetPt ) {
		   double gv = (doToys) ? getProbGaus(hPDFBJetPt.get(), maxPDFBJetPt, meanPDFBJetPt, sigmaPDFBJetPt, rnd, NBJetPtRMS_) : 0.;
		   par[FPARAM_BJetHadPt_TOPTOPLEPHAD] = *PtBJet2/(1.-gv);
		} else par[FPARAM_BJetHadPt_TOPTOPLEPHAD] = *PtBJet2;
		
		if( ! isDeltaFuncPDFBJetEta ) {		       
		   double gv = (doToys) ? getProbGaus(hPDFBJetEta.get(), maxPDFBJetEta, meanPDFBJetEta, sigmaPDFBJetEta, rnd, NBJetEtaRMS_) : 0.;
		   par[FPARAM_BJetHadEta_TOPTOPLEPHAD] = *EtaBJet2/(1.-gv);
		} else par[FPARAM_BJetHadEta_TOPTOPLEPHAD] = *EtaBJet2;		  
		  
		if( ! isDeltaFuncPDFBJetPhi ) {
		   double gv = (doToys) ? getProbGaus(hPDFBJetPhi.get(), maxPDFBJetPhi, meanPDFBJetPhi, sigmaPDFBJetPhi, rnd, NBJetPhiRMS_) : 0.;
		   par[FPARAM_BJetHadPhi_TOPTOPLEPHAD] = *PhiBJet2/(1.-gv);
		} else par[FPARAM_BJetHadPhi_TOPTOPLEPHAD] = *PhiBJet2;		  
		
		par[FPARAM_BJetHadPx_TOPTOPLEPHAD] = par[FPARAM_BJetHadPt_TOPTOPLEPHAD]*cos(par[FPARAM_BJetHadPhi_TOPTOPLEPHAD]);
		par[FPARAM_BJetHadPy_TOPTOPLEPHAD] = par[FPARAM_BJetHadPt_TOPTOPLEPHAD]*sin(par[FPARAM_BJetHadPhi_TOPTOPLEPHAD]);
		par[FPARAM_BJetHadPz_TOPTOPLEPHAD] = par[FPARAM_BJetHadPt_TOPTOPLEPHAD]*sinh(par[FPARAM_BJetHadEta_TOPTOPLEPHAD]);
		par[FPARAM_BJetHadE_TOPTOPLEPHAD] = sqrt((*MassBJet2)*(*MassBJet2) + par[FPARAM_BJetHadPx_TOPTOPLEPHAD]*par[FPARAM_BJetHadPx_TOPTOPLEPHAD] + par[FPARAM_BJetHadPy_TOPTOPLEPHAD]*par[FPARAM_BJetHadPy_TOPTOPLEPHAD] + par[FPARAM_BJetHadPz_TOPTOPLEPHAD]*par[FPARAM_BJetHadPz_TOPTOPLEPHAD]);
	     } while( (par[FPARAM_BJetHadE_TOPTOPLEPHAD]*par[FPARAM_BJetHadE_TOPTOPLEPHAD]-par[FPARAM_BJetHadPx_TOPTOPLEPHAD]*par[FPARAM_BJetHadPx_TOPTOPLEPHAD]-par[FPARAM_BJetHadPy_TOPTOPLEPHAD]*par[FPARAM_BJetHadPy_TOPTOPLEPHAD]-par[FPARAM_BJetHadPz_TOPTOPLEPHAD]*par[FPARAM_BJetHadPz_TOPTOPLEPHAD]) < thres );	   
	}

	if( *usePDFBJetPxPyPz )
	  {
	     do
	       {
		  if( ! isDeltaFuncPDFNonBJetPx ) {
		     double gv = (doToys) ? getProbGaus(hPDFNonBJetPx.get(), maxPDFNonBJetPx, meanPDFNonBJetPx, sigmaPDFNonBJetPx, rnd, NNonBJetPxRMS_) : 0.;
		     par[FPARAM_NonBJet1Px_TOPTOPLEPHAD] = *PxNonBJet1/(1.-gv);
		  } else par[FPARAM_NonBJet1Px_TOPTOPLEPHAD] = *PxNonBJet1;

		  if( ! isDeltaFuncPDFNonBJetPy ) {
		     double gv = (doToys) ? getProbGaus(hPDFNonBJetPy.get(), maxPDFNonBJetPy, meanPDFNonBJetPy, sigmaPDFNonBJetPy, rnd, NNonBJetPyRMS_) : 0.;
		     par[FPARAM_NonBJet1Py_TOPTOPLEPHAD] = *PyNonBJet1/(1.-gv);
		  } else par[FPARAM_NonBJet1Py_TOPTOPLEPHAD] = *PyNonBJet1;		  
		  
		  if( ! isDeltaFuncPDFNonBJetPz ) {
		     double gv = (doToys) ? getProbGaus(hPDFNonBJetPz.get(), maxPDFNonBJetPz, meanPDFNonBJetPz, sigmaPDFNonBJetPz, rnd, NNonBJetPzRMS_) : 0.;
		     par[FPARAM_NonBJet1Pz_TOPTOPLEPHAD] = *PzNonBJet1/(1.-gv);
		  } else par[FPARAM_NonBJet1Pz_TOPTOPLEPHAD] = *PzNonBJet1;
		  
		  par[FPARAM_NonBJet1Pt_TOPTOPLEPHAD] = sqrt(par[FPARAM_NonBJet1Px_TOPTOPLEPHAD]*par[FPARAM_NonBJet1Px_TOPTOPLEPHAD] + par[FPARAM_NonBJet1Py_TOPTOPLEPHAD]*par[FPARAM_NonBJet1Py_TOPTOPLEPHAD]);
		  par[FPARAM_NonBJet1Eta_TOPTOPLEPHAD] = getEta(par[FPARAM_NonBJet1Pt_TOPTOPLEPHAD], par[FPARAM_NonBJet1Pz_TOPTOPLEPHAD]);
		  par[FPARAM_NonBJet1Phi_TOPTOPLEPHAD] = atan2(par[FPARAM_NonBJet1Py_TOPTOPLEPHAD], par[FPARAM_NonBJet1Px_TOPTOPLEPHAD]);
		  par[FPARAM_NonBJet1E_TOPTOPLEPHAD] = sqrt((*MassNonBJet1)*(*MassNonBJet1) + par[FPARAM_NonBJet1Px_TOPTOPLEPHAD]*par[FPARAM_NonBJet1Px_TOPTOPLEPHAD] + par[FPARAM_NonBJet1Py_TOPTOPLEPHAD]*par[FPARAM_NonBJet1Py_TOPTOPLEPHAD] + par[FPARAM_NonBJet1Pz_TOPTOPLEPHAD]*par[FPARAM_NonBJet1Pz_TOPTOPLEPHAD]);
	       } while( (par[FPARAM_NonBJet1E_TOPTOPLEPHAD]*par[FPARAM_NonBJet1E_TOPTOPLEPHAD]-par[FPARAM_NonBJet1Px_TOPTOPLEPHAD]*par[FPARAM_NonBJet1Px_TOPTOPLEPHAD]-par[FPARAM_NonBJet1Py_TOPTOPLEPHAD]*par[FPARAM_NonBJet1Py_TOPTOPLEPHAD]-par[FPARAM_NonBJet1Pz_TOPTOPLEPHAD]*par[FPARAM_NonBJet1Pz_TOPTOPLEPHAD]) < thres );
	  }
	else {
	   do
	     {
		if( ! isDeltaFuncPDFNonBJetPt ) {
		   double gv = (doToys) ? getProbGaus(hPDFNonBJetPt.get(), maxPDFNonBJetPt, meanPDFNonBJetPt, sigmaPDFNonBJetPt, rnd, NNonBJetPtRMS_) : 0.;
		   par[FPARAM_NonBJet2Pt_TOPTOPLEPHAD] = *PtNonBJet2/(1.-gv);
		} else par[FPARAM_NonBJet2Pt_TOPTOPLEPHAD] = *PtNonBJet2;
		
		if( ! isDeltaFuncPDFNonBJetEta ) {		       
		   double gv = (doToys) ? getProbGaus(hPDFNonBJetEta.get(), maxPDFNonBJetEta, meanPDFNonBJetEta, sigmaPDFNonBJetEta, rnd, NNonBJetEtaRMS_) : 0.;
		   par[FPARAM_NonBJet2Eta_TOPTOPLEPHAD] = *EtaNonBJet2/(1.-gv);
		} else par[FPARAM_NonBJet2Eta_TOPTOPLEPHAD] = *EtaNonBJet2;
		  
		if( ! isDeltaFuncPDFNonBJetPhi ) {
		   double gv = (doToys) ? getProbGaus(hPDFNonBJetPhi.get(), maxPDFNonBJetPhi, meanPDFNonBJetPhi, sigmaPDFNonBJetPhi, rnd, NNonBJetPhiRMS_) : 0.;
		   par[FPARAM_NonBJet2Phi_TOPTOPLEPHAD] = *PhiNonBJet2/(1.-gv);
		} else par[FPARAM_NonBJet2Phi_TOPTOPLEPHAD] = *PhiNonBJet2;		  
		
		par[FPARAM_NonBJet2Px_TOPTOPLEPHAD] = par[FPARAM_NonBJet2Pt_TOPTOPLEPHAD]*cos(par[FPARAM_NonBJet2Phi_TOPTOPLEPHAD]);
		par[FPARAM_NonBJet2Py_TOPTOPLEPHAD] = par[FPARAM_NonBJet2Pt_TOPTOPLEPHAD]*sin(par[FPARAM_NonBJet2Phi_TOPTOPLEPHAD]);
		par[FPARAM_NonBJet2Pz_TOPTOPLEPHAD] = par[FPARAM_NonBJet2Pt_TOPTOPLEPHAD]*sinh(par[FPARAM_NonBJet2Eta_TOPTOPLEPHAD]);
		par[FPARAM_NonBJet2E_TOPTOPLEPHAD] = sqrt((*MassNonBJet2)*(*MassNonBJet2) + par[FPARAM_NonBJet2Px_TOPTOPLEPHAD]*par[FPARAM_NonBJet2Px_TOPTOPLEPHAD] + par[FPARAM_NonBJet2Py_TOPTOPLEPHAD]*par[FPARAM_NonBJet2Py_TOPTOPLEPHAD] + par[FPARAM_NonBJet2Pz_TOPTOPLEPHAD]*par[FPARAM_NonBJet2Pz_TOPTOPLEPHAD]);
	     } while( (par[FPARAM_NonBJet2E_TOPTOPLEPHAD]*par[FPARAM_NonBJet2E_TOPTOPLEPHAD]-par[FPARAM_NonBJet2Px_TOPTOPLEPHAD]*par[FPARAM_NonBJet2Px_TOPTOPLEPHAD]-par[FPARAM_NonBJet2Py_TOPTOPLEPHAD]*par[FPARAM_NonBJet2Py_TOPTOPLEPHAD]-par[FPARAM_NonBJet2Pz_TOPTOPLEPHAD]*par[FPARAM_NonBJet2Pz_TOPTOPLEPHAD]) < thres );
	}

	if( IncludePhotons_ )
	  {
	     if( *usePDFPhotonPxPyPz ) {		  
		do
		  {
		     if( ! isDeltaFuncPDFPhotonPx ) {
			double gv = (doToys) ? getProbGaus(hPDFPhotonPx.get(), maxPDFPhotonPx, meanPDFPhotonPx, sigmaPDFPhotonPx, rnd, NPhotonPxRMS_) : 0.;
			par[FPARAM_PhotonPx_TOPTOPLEPHAD] = *PxPhoton/(1.-gv);
		     } else par[FPARAM_PhotonPx_TOPTOPLEPHAD] = *PxPhoton;
		     
		     if( ! isDeltaFuncPDFPhotonPy ) {			  
			double gv = (doToys) ? getProbGaus(hPDFPhotonPy.get(), maxPDFPhotonPy, meanPDFPhotonPy, sigmaPDFPhotonPy, rnd, NPhotonPyRMS_) : 0.;
			par[FPARAM_PhotonPy_TOPTOPLEPHAD] = *PyPhoton/(1.-gv);
		     } else par[FPARAM_PhotonPy_TOPTOPLEPHAD] = *PyPhoton;
		     
		     if( ! isDeltaFuncPDFPhotonPz ) {			  
			double gv = (doToys) ? getProbGaus(hPDFPhotonPz.get(), maxPDFPhotonPz, meanPDFPhotonPz, sigmaPDFPhotonPz, rnd, NPhotonPzRMS_) : 0.;
			par[FPARAM_PhotonPz_TOPTOPLEPHAD] = *PzPhoton/(1.-gv);
		     } else par[FPARAM_PhotonPz_TOPTOPLEPHAD] = *PzPhoton;
		     
		     par[FPARAM_PhotonPt_TOPTOPLEPHAD] = sqrt(par[FPARAM_PhotonPx_TOPTOPLEPHAD]*par[FPARAM_PhotonPx_TOPTOPLEPHAD] + par[FPARAM_PhotonPy_TOPTOPLEPHAD]*par[FPARAM_PhotonPy_TOPTOPLEPHAD]);
		     par[FPARAM_PhotonEta_TOPTOPLEPHAD] = getEta(par[FPARAM_PhotonPt_TOPTOPLEPHAD], par[FPARAM_PhotonPz_TOPTOPLEPHAD]);
		     par[FPARAM_PhotonPhi_TOPTOPLEPHAD] = atan2(par[FPARAM_PhotonPy_TOPTOPLEPHAD], par[FPARAM_PhotonPx_TOPTOPLEPHAD]);
		     par[FPARAM_PhotonE_TOPTOPLEPHAD] = sqrt(par[FPARAM_PhotonPx_TOPTOPLEPHAD]*par[FPARAM_PhotonPx_TOPTOPLEPHAD] + par[FPARAM_PhotonPy_TOPTOPLEPHAD]*par[FPARAM_PhotonPy_TOPTOPLEPHAD] + par[FPARAM_PhotonPz_TOPTOPLEPHAD]*par[FPARAM_PhotonPz_TOPTOPLEPHAD]);
		  } while( (par[FPARAM_PhotonE_TOPTOPLEPHAD]*par[FPARAM_PhotonE_TOPTOPLEPHAD]-par[FPARAM_PhotonPx_TOPTOPLEPHAD]*par[FPARAM_PhotonPx_TOPTOPLEPHAD]-par[FPARAM_PhotonPy_TOPTOPLEPHAD]*par[FPARAM_PhotonPy_TOPTOPLEPHAD]-par[FPARAM_PhotonPz_TOPTOPLEPHAD]*par[FPARAM_PhotonPz_TOPTOPLEPHAD]) < thres );
	     }
	     else {
		do
		  {
		     if( ! isDeltaFuncPDFPhotonPt ) {
			double gv = (doToys) ? getProbGaus(hPDFPhotonPt.get(), maxPDFPhotonPt, meanPDFPhotonPt, sigmaPDFPhotonPt, rnd, NPhotonPtRMS_) : 0.;
			par[FPARAM_PhotonPt_TOPTOPLEPHAD] = *PtPhoton/(1.-gv);
		     } else par[FPARAM_PhotonPt_TOPTOPLEPHAD] = *PtPhoton;
		     
		     if( ! isDeltaFuncPDFPhotonEta ) {
			double gv = (doToys) ? getProbGaus(hPDFPhotonEta.get(), maxPDFPhotonEta, meanPDFPhotonEta, sigmaPDFPhotonEta, rnd, NPhotonEtaRMS_) : 0.;
			par[FPARAM_PhotonEta_TOPTOPLEPHAD] = *EtaPhoton/(1.-gv);
		     } else par[FPARAM_PhotonEta_TOPTOPLEPHAD] = *EtaPhoton;
		     
		     if( ! isDeltaFuncPDFPhotonPhi ) {
			double gv = (doToys) ? getProbGaus(hPDFPhotonPhi.get(), maxPDFPhotonPhi, meanPDFPhotonPhi, sigmaPDFPhotonPhi, rnd, NPhotonPhiRMS_) : 0.;
			par[FPARAM_PhotonPhi_TOPTOPLEPHAD] = *PhiPhoton/(1.-gv);
		     } else par[FPARAM_PhotonPhi_TOPTOPLEPHAD] = *PhiPhoton;
		     
		     par[FPARAM_PhotonPx_TOPTOPLEPHAD] = par[FPARAM_PhotonPt_TOPTOPLEPHAD]*cos(par[FPARAM_PhotonPhi_TOPTOPLEPHAD]);
		     par[FPARAM_PhotonPy_TOPTOPLEPHAD] = par[FPARAM_PhotonPt_TOPTOPLEPHAD]*sin(par[FPARAM_PhotonPhi_TOPTOPLEPHAD]);
		     par[FPARAM_PhotonPz_TOPTOPLEPHAD] = par[FPARAM_PhotonPt_TOPTOPLEPHAD]*sinh(par[FPARAM_PhotonEta_TOPTOPLEPHAD]);
		     par[FPARAM_PhotonE_TOPTOPLEPHAD] = sqrt(par[FPARAM_PhotonPx_TOPTOPLEPHAD]*par[FPARAM_PhotonPx_TOPTOPLEPHAD] + par[FPARAM_PhotonPy_TOPTOPLEPHAD]*par[FPARAM_PhotonPy_TOPTOPLEPHAD] + par[FPARAM_PhotonPz_TOPTOPLEPHAD]*par[FPARAM_PhotonPz_TOPTOPLEPHAD]);
		  } while( (par[FPARAM_PhotonE_TOPTOPLEPHAD]*par[FPARAM_PhotonE_TOPTOPLEPHAD]-par[FPARAM_PhotonPx_TOPTOPLEPHAD]*par[FPARAM_PhotonPx_TOPTOPLEPHAD]-par[FPARAM_PhotonPy_TOPTOPLEPHAD]*par[FPARAM_PhotonPy_TOPTOPLEPHAD]-par[FPARAM_PhotonPz_TOPTOPLEPHAD]*par[FPARAM_PhotonPz_TOPTOPLEPHAD]) < thres );		
	     }	     
	  }
	
	  {
	     if( ! isDeltaFuncPDFMetPx ) {		  
		double gv = (doToys) ? getProbGaus(hPDFMetPx.get(), maxPDFMetPx, meanPDFMetPx, sigmaPDFMetPx, rnd, NMetRMS_) : 0.;
		par[FPARAM_EtRealX_TOPTOPLEPHAD] = *EtMissX/(1.-gv);
	     } else par[FPARAM_EtRealX_TOPTOPLEPHAD] = *EtMissX;
	     
	     if( ! isDeltaFuncPDFMetPy ) {		  
		double gv = (doToys) ? getProbGaus(hPDFMetPy.get(), maxPDFMetPy, meanPDFMetPy, sigmaPDFMetPy, rnd, NMetRMS_) : 0.;
		par[FPARAM_EtRealY_TOPTOPLEPHAD] = *EtMissY/(1.-gv);
	     } else par[FPARAM_EtRealY_TOPTOPLEPHAD] = *EtMissY;	     
	  }	

	par[FPARAM_mWLep_TOPTOPLEPHAD] = (doToys) ? getProbGaus(hPDFTopWMass.get(), maxPDFTopWMass, meanPDFTopWMass, sigmaPDFTopWMass, rnd, 2) : WMass1;
	par[FPARAM_mWHad_TOPTOPLEPHAD] = (doToys) ? getProbGaus(hPDFTopWHadMass.get(), maxPDFTopWHadMass, meanPDFTopWHadMass, sigmaPDFTopWHadMass, rnd, 2) : WMass2;

	if( *LabelLepton1 == 0 )
	  {
	     if( *usePDFLeptonPxPyPz ) {
		do
		  {
		     if( ! isDeltaFuncPDFElecPx ) {
			double gv = (doToys) ? getProbGaus(hPDFElecPx.get(), maxPDFElecPx, meanPDFElecPx, sigmaPDFElecPx, rnd, NElecPxRMS_) : 0.;
			par[FPARAM_LeptonPx_TOPTOPLEPHAD] = *PxLepton1/(1.-gv);
		     } else par[FPARAM_LeptonPx_TOPTOPLEPHAD] = *PxLepton1;		     
		     
		     if( ! isDeltaFuncPDFElecPy ) {			  
			double gv = (doToys) ? getProbGaus(hPDFElecPy.get(), maxPDFElecPy, meanPDFElecPy, sigmaPDFElecPy, rnd, NElecPyRMS_) : 0.;
			par[FPARAM_LeptonPy_TOPTOPLEPHAD] = *PyLepton1/(1.-gv);
		     } else par[FPARAM_LeptonPy_TOPTOPLEPHAD] = *PyLepton1;		     
		     
		     if( ! isDeltaFuncPDFElecPz ) {			  
			double gv = (doToys) ? getProbGaus(hPDFElecPz.get(), maxPDFElecPz, meanPDFElecPz, sigmaPDFElecPz, rnd, NElecPzRMS_) : 0.;
			par[FPARAM_LeptonPz_TOPTOPLEPHAD] = *PzLepton1/(1.-gv);
		     } else par[FPARAM_LeptonPz_TOPTOPLEPHAD] = *PzLepton1;
		     
		     par[FPARAM_LeptonPt_TOPTOPLEPHAD] = sqrt(par[FPARAM_LeptonPx_TOPTOPLEPHAD]*par[FPARAM_LeptonPx_TOPTOPLEPHAD] + par[FPARAM_LeptonPy_TOPTOPLEPHAD]*par[FPARAM_LeptonPy_TOPTOPLEPHAD]);
		     par[FPARAM_LeptonEta_TOPTOPLEPHAD] = getEta(par[FPARAM_LeptonPt_TOPTOPLEPHAD], par[FPARAM_LeptonPz_TOPTOPLEPHAD]);
		     par[FPARAM_LeptonPhi_TOPTOPLEPHAD] = atan2(par[FPARAM_LeptonPy_TOPTOPLEPHAD], par[FPARAM_LeptonPx_TOPTOPLEPHAD]);
		     par[FPARAM_LeptonE_TOPTOPLEPHAD] = sqrt((*MassLepton1)*(*MassLepton1) + par[FPARAM_LeptonPx_TOPTOPLEPHAD]*par[FPARAM_LeptonPx_TOPTOPLEPHAD] + par[FPARAM_LeptonPy_TOPTOPLEPHAD]*par[FPARAM_LeptonPy_TOPTOPLEPHAD] + par[FPARAM_LeptonPz_TOPTOPLEPHAD]*par[FPARAM_LeptonPz_TOPTOPLEPHAD]);
		  } while( (par[FPARAM_LeptonE_TOPTOPLEPHAD]*par[FPARAM_LeptonE_TOPTOPLEPHAD]-par[FPARAM_LeptonPx_TOPTOPLEPHAD]*par[FPARAM_LeptonPx_TOPTOPLEPHAD]-par[FPARAM_LeptonPy_TOPTOPLEPHAD]*par[FPARAM_LeptonPy_TOPTOPLEPHAD]-par[FPARAM_LeptonPz_TOPTOPLEPHAD]*par[FPARAM_LeptonPz_TOPTOPLEPHAD]) < thres );
	     }
	     else {
		do
		  {
		     if( ! isDeltaFuncPDFElecPt ) {
			double gv = (doToys) ? getProbGaus(hPDFElecPt.get(), maxPDFElecPt, meanPDFElecPt, sigmaPDFElecPt, rnd, NElecPtRMS_) : 0.;
			par[FPARAM_LeptonPt_TOPTOPLEPHAD] = *PtLepton1/(1.-gv);
		     } else par[FPARAM_LeptonPt_TOPTOPLEPHAD] = *PtLepton1;		     
		     
		     if( ! isDeltaFuncPDFElecEta ) {
			double gv = (doToys) ? getProbGaus(hPDFElecEta.get(), maxPDFElecEta, meanPDFElecEta, sigmaPDFElecEta, rnd, NElecEtaRMS_) : 0.;
			par[FPARAM_LeptonEta_TOPTOPLEPHAD] = *EtaLepton1/(1.-gv);
		     } else par[FPARAM_LeptonEta_TOPTOPLEPHAD] = *EtaLepton1;		     
		     
		     if( ! isDeltaFuncPDFElecPhi ) {
			double gv = (doToys) ? getProbGaus(hPDFElecPhi.get(), maxPDFElecPhi, meanPDFElecPhi, sigmaPDFElecPhi, rnd, NElecPhiRMS_) : 0.;
			par[FPARAM_LeptonPhi_TOPTOPLEPHAD] = *PhiLepton1/(1.-gv);
		     } else par[FPARAM_LeptonPhi_TOPTOPLEPHAD] = *PhiLepton1;
		     
		     par[FPARAM_LeptonPx_TOPTOPLEPHAD] = par[FPARAM_LeptonPt_TOPTOPLEPHAD]*cos(par[FPARAM_LeptonPhi_TOPTOPLEPHAD]);
		     par[FPARAM_LeptonPy_TOPTOPLEPHAD] = par[FPARAM_LeptonPt_TOPTOPLEPHAD]*sin(par[FPARAM_LeptonPhi_TOPTOPLEPHAD]);
		     par[FPARAM_LeptonPz_TOPTOPLEPHAD] = par[FPARAM_LeptonPt_TOPTOPLEPHAD]*sinh(par[FPARAM_LeptonEta_TOPTOPLEPHAD]);
		     par[FPARAM_LeptonE_TOPTOPLEPHAD] = sqrt((*MassLepton1)*(*MassLepton1) + par[FPARAM_LeptonPx_TOPTOPLEPHAD]*par[FPARAM_LeptonPx_TOPTOPLEPHAD] + par[FPARAM_LeptonPy_TOPTOPLEPHAD]*par[FPARAM_LeptonPy_TOPTOPLEPHAD] + par[FPARAM_LeptonPz_TOPTOPLEPHAD]*par[FPARAM_LeptonPz_TOPTOPLEPHAD]);
		  } while( (par[FPARAM_LeptonE_TOPTOPLEPHAD]*par[FPARAM_LeptonE_TOPTOPLEPHAD]-par[FPARAM_LeptonPx_TOPTOPLEPHAD]*par[FPARAM_LeptonPx_TOPTOPLEPHAD]-par[FPARAM_LeptonPy_TOPTOPLEPHAD]*par[FPARAM_LeptonPy_TOPTOPLEPHAD]-par[FPARAM_LeptonPz_TOPTOPLEPHAD]*par[FPARAM_LeptonPz_TOPTOPLEPHAD]) < thres );		
	     }
	  }
	else
	  {
	     if( *usePDFLeptonPxPyPz ) {
		do
		  {
		     if( ! isDeltaFuncPDFMuonPx ) {
			double gv = (doToys) ? getProbGaus(hPDFMuonPx.get(), maxPDFMuonPx, meanPDFMuonPx, sigmaPDFMuonPx, rnd, NMuonPxRMS_) : 0.;
			par[FPARAM_LeptonPx_TOPTOPLEPHAD] = *PxLepton1/(1.-gv);
		     } else par[FPARAM_LeptonPx_TOPTOPLEPHAD] = *PxLepton1;		     
		     
		     if( ! isDeltaFuncPDFMuonPy ) {			  
			double gv = (doToys) ? getProbGaus(hPDFMuonPy.get(), maxPDFMuonPy, meanPDFMuonPy, sigmaPDFMuonPy, rnd, NMuonPyRMS_) : 0.;
			par[FPARAM_LeptonPy_TOPTOPLEPHAD] = *PyLepton1/(1.-gv);
		     } else par[FPARAM_LeptonPy_TOPTOPLEPHAD] = *PyLepton1;		     
		     
		     if( ! isDeltaFuncPDFMuonPz ) {			  
			double gv = (doToys) ? getProbGaus(hPDFMuonPz.get(), maxPDFMuonPz, meanPDFMuonPz, sigmaPDFMuonPz, rnd, NMuonPzRMS_) : 0.;
			par[FPARAM_LeptonPz_TOPTOPLEPHAD] = *PzLepton1/(1.-gv);
		     } else par[FPARAM_LeptonPz_TOPTOPLEPHAD] = *PzLepton1;
		     
		     par[FPARAM_LeptonPt_TOPTOPLEPHAD] = sqrt(par[FPARAM_LeptonPx_TOPTOPLEPHAD]*par[FPARAM_LeptonPx_TOPTOPLEPHAD] + par[FPARAM_LeptonPy_TOPTOPLEPHAD]*par[FPARAM_LeptonPy_TOPTOPLEPHAD]);
		     par[FPARAM_LeptonEta_TOPTOPLEPHAD] = getEta(par[FPARAM_LeptonPt_TOPTOPLEPHAD], par[FPARAM_LeptonPz_TOPTOPLEPHAD]);
		     par[FPARAM_LeptonPhi_TOPTOPLEPHAD] = atan2(par[FPARAM_LeptonPy_TOPTOPLEPHAD], par[FPARAM_LeptonPx_TOPTOPLEPHAD]);
		     par[FPARAM_LeptonE_TOPTOPLEPHAD] = sqrt((*MassLepton1)*(*MassLepton1) + par[FPARAM_LeptonPx_TOPTOPLEPHAD]*par[FPARAM_LeptonPx_TOPTOPLEPHAD] + par[FPARAM_LeptonPy_TOPTOPLEPHAD]*par[FPARAM_LeptonPy_TOPTOPLEPHAD] + par[FPARAM_LeptonPz_TOPTOPLEPHAD]*par[FPARAM_LeptonPz_TOPTOPLEPHAD]);
		  } while( (par[FPARAM_LeptonE_TOPTOPLEPHAD]*par[FPARAM_LeptonE_TOPTOPLEPHAD]-par[FPARAM_LeptonPx_TOPTOPLEPHAD]*par[FPARAM_LeptonPx_TOPTOPLEPHAD]-par[FPARAM_LeptonPy_TOPTOPLEPHAD]*par[FPARAM_LeptonPy_TOPTOPLEPHAD]-par[FPARAM_LeptonPz_TOPTOPLEPHAD]*par[FPARAM_LeptonPz_TOPTOPLEPHAD]) < thres );
	     }
	     else {
		do
		  {
		     if( ! isDeltaFuncPDFMuonPt ) {
			double gv = (doToys) ? getProbGaus(hPDFMuonPt.get(), maxPDFMuonPt, meanPDFMuonPt, sigmaPDFMuonPt, rnd, NMuonPtRMS_) : 0.;
			par[FPARAM_LeptonPt_TOPTOPLEPHAD] = *PtLepton1/(1.-gv);
		     } else par[FPARAM_LeptonPt_TOPTOPLEPHAD] = *PtLepton1;		     
		     
		     if( ! isDeltaFuncPDFMuonEta ) {
			double gv = (doToys) ? getProbGaus(hPDFMuonEta.get(), maxPDFMuonEta, meanPDFMuonEta, sigmaPDFMuonEta, rnd, NMuonEtaRMS_) : 0.;
			par[FPARAM_LeptonEta_TOPTOPLEPHAD] = *EtaLepton1/(1.-gv);
		     } else par[FPARAM_LeptonEta_TOPTOPLEPHAD] = *EtaLepton1;		     
		     
		     if( ! isDeltaFuncPDFMuonPhi ) {
			double gv = (doToys) ? getProbGaus(hPDFMuonPhi.get(), maxPDFMuonPhi, meanPDFMuonPhi, sigmaPDFMuonPhi, rnd, NMuonPhiRMS_) : 0.;
			par[FPARAM_LeptonPhi_TOPTOPLEPHAD] = *PhiLepton1/(1.-gv);
		     } else par[FPARAM_LeptonPhi_TOPTOPLEPHAD] = *PhiLepton1;
		     
		     par[FPARAM_LeptonPx_TOPTOPLEPHAD] = par[FPARAM_LeptonPt_TOPTOPLEPHAD]*cos(par[FPARAM_LeptonPhi_TOPTOPLEPHAD]);
		     par[FPARAM_LeptonPy_TOPTOPLEPHAD] = par[FPARAM_LeptonPt_TOPTOPLEPHAD]*sin(par[FPARAM_LeptonPhi_TOPTOPLEPHAD]);
		     par[FPARAM_LeptonPz_TOPTOPLEPHAD] = par[FPARAM_LeptonPt_TOPTOPLEPHAD]*sinh(par[FPARAM_LeptonEta_TOPTOPLEPHAD]);
		     par[FPARAM_LeptonE_TOPTOPLEPHAD] = sqrt((*MassLepton1)*(*MassLepton1) + par[FPARAM_LeptonPx_TOPTOPLEPHAD]*par[FPARAM_LeptonPx_TOPTOPLEPHAD] + par[FPARAM_LeptonPy_TOPTOPLEPHAD]*par[FPARAM_LeptonPy_TOPTOPLEPHAD] + par[FPARAM_LeptonPz_TOPTOPLEPHAD]*par[FPARAM_LeptonPz_TOPTOPLEPHAD]);
		  } while( (par[FPARAM_LeptonE_TOPTOPLEPHAD]*par[FPARAM_LeptonE_TOPTOPLEPHAD]-par[FPARAM_LeptonPx_TOPTOPLEPHAD]*par[FPARAM_LeptonPx_TOPTOPLEPHAD]-par[FPARAM_LeptonPy_TOPTOPLEPHAD]*par[FPARAM_LeptonPy_TOPTOPLEPHAD]-par[FPARAM_LeptonPz_TOPTOPLEPHAD]*par[FPARAM_LeptonPz_TOPTOPLEPHAD]) < thres );
	     }
	  }
			
	for( int is1=-1;is1<=1;is1++ )
	  {	
	     if( is1 == 0 ) continue;
	     
	     par[FPARAM_Sign_TOPTOPLEPHAD] = is1;

	     chi2t.clear();

	     double lh = func(*PxLepton1, *PyLepton1, *PzLepton1, *PtLepton1, *EtaLepton1, *PhiLepton1, *ELepton1, *LabelLepton1,
			      *PxBJet1, *PyBJet1, *PzBJet1, *PtBJet1, *EtaBJet1, *PhiBJet1, *EBJet1,
			      *PxBJet2, *PyBJet2, *PzBJet2, *PtBJet2, *EtaBJet2, *PhiBJet2, *EBJet2,
			      *PxNonBJet1, *PyNonBJet1, *PzNonBJet1, *PtNonBJet1, *EtaNonBJet1, *PhiNonBJet1, *ENonBJet1,
			      *PxNonBJet2, *PyNonBJet2, *PzNonBJet2, *PtNonBJet2, *EtaNonBJet2, *PhiNonBJet2, *ENonBJet2,
			      *PxPhoton, *PyPhoton, *PzPhoton, *PtPhoton, *EtaPhoton, *PhiPhoton, *EPhoton,
			      *PhotonOrigin,
			      chi2t, par);

	     if( lh < LHMaxGeneric_ )
	       {
		  FRESULT fres = {};
		  for( int ip=0;ip<FPARAM_N;ip++ ) fres.par[ip] = par[ip];
		  fres.lh = lh;
		  fres.PzNu1 = *KINFIT::kfit::PzNu1;
		  
		  vp.push_back(fres);
		  
		  if( vp.size() >= NGrid_ ) return;
	       }
	  }
     }   
}
