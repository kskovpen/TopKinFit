// Top quark pair production in single lepton channel: t -> (W->lnu) b

#include "../include/TopLep.h"

ClassImp(KINFIT::TopLep)

KINFIT::TopLep::TopLep() {}

KINFIT::TopLep::~TopLep() {}

/// Main routine
void KINFIT::TopLep::TopLepRun()
{
   NPERMEVENT_PREV = NPERMEVENT;

   /// Calculate number of permutations in the current event
   CalcNPerm();

   /// Reinitialize all variables and release memory used in the previous event
   
   NPerm_ = 0;
   NTerm_ = 0;

   if( idxMin_ ) {delete[] idxMin_; idxMin_ = 0;}
   
   if( TopLep_ElectronIdx )     {delete[] TopLep_ElectronIdx; TopLep_ElectronIdx = 0;}
   if( TopLep_MuonIdx )         {delete[] TopLep_MuonIdx; TopLep_MuonIdx = 0;}
   if( TopLep_BJetLepIdx )      {delete[] TopLep_BJetLepIdx; TopLep_BJetLepIdx = 0;}
   if( TopLep_PhotonIdx )       {delete[] TopLep_PhotonIdx; TopLep_PhotonIdx = 0;}

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
   
   TopLep_ElectronIdx = new int[NPERMEVENT];
   TopLep_MuonIdx = new int[NPERMEVENT];
   TopLep_BJetLepIdx = new int[NPERMEVENT];
   TopLep_PhotonIdx = new int[NPERMEVENT];
   
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
	TopLep_ElectronIdx[i] = -1;
	TopLep_MuonIdx[i] = -1;
	TopLep_BJetLepIdx[i] = -1;
	TopLep_PhotonIdx[i] = -1;
	
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
   
   (*ParMin)[FPARAM_Sign_TOPLEP] = -1.;
   (*ParMax)[FPARAM_Sign_TOPLEP] = 1.;

   (*ParMin)[FPARAM_mWLep_TOPLEP] = 50.;
   (*ParMax)[FPARAM_mWLep_TOPLEP] = 110.;

   /// Loop through permutations

   for( int il=0;il<nLepton;il++ )
     {	
	(*ParMin)[FPARAM_EtRealX_TOPLEP] = (*EtMissX) - LimNRMS_*fabs(*EtMissX)*sigmaPDFMetPx;
	(*ParMax)[FPARAM_EtRealX_TOPLEP] = (*EtMissX) + LimNRMS_*fabs(*EtMissX)*sigmaPDFMetPx;
	(*ParMin)[FPARAM_EtRealY_TOPLEP] = (*EtMissY) - LimNRMS_*fabs(*EtMissY)*sigmaPDFMetPy;
	(*ParMax)[FPARAM_EtRealY_TOPLEP] = (*EtMissY) + LimNRMS_*fabs(*EtMissY)*sigmaPDFMetPy;
	
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
	
	(*ParMin)[FPARAM_LeptonPx_TOPLEP] = (*PxLepton1) - LimNRMS_*fabs(*PxLepton1)*sigmaPDFLeptonPx;
	(*ParMax)[FPARAM_LeptonPx_TOPLEP] = (*PxLepton1) + LimNRMS_*fabs(*PxLepton1)*sigmaPDFLeptonPx;
	(*ParMin)[FPARAM_LeptonPy_TOPLEP] = (*PyLepton1) - LimNRMS_*fabs(*PyLepton1)*sigmaPDFLeptonPy;
	(*ParMax)[FPARAM_LeptonPy_TOPLEP] = (*PyLepton1) + LimNRMS_*fabs(*PyLepton1)*sigmaPDFLeptonPy;
	(*ParMin)[FPARAM_LeptonPz_TOPLEP] = (*PzLepton1) - LimNRMS_*fabs(*PzLepton1)*sigmaPDFLeptonPz;
	(*ParMax)[FPARAM_LeptonPz_TOPLEP] = (*PzLepton1) + LimNRMS_*fabs(*PzLepton1)*sigmaPDFLeptonPz;
	
	int label1_l = LeptonLabel[il];
	int idx1_l = LeptonIdx[il];
	
	for( int ib=0;ib<nBJet;ib++ )
	  {
	     *EBJet1 = BJetE[ib];
	     *PxBJet1 = BJetPx[ib];
	     *PyBJet1 = BJetPy[ib];
	     *PzBJet1 = BJetPz[ib];
	     *MassBJet1 = (*EBJet1)*(*EBJet1) - (*PxBJet1)*(*PxBJet1) - (*PyBJet1)*(*PyBJet1) - (*PzBJet1)*(*PzBJet1);
	     *MassBJet1 = (*MassBJet1 > 0.) ? sqrt(*MassBJet1) : 0.;
	     
	     (*ParMin)[FPARAM_BJetLepPx_TOPLEP] = (*PxBJet1) - LimNRMS_*fabs(*PxBJet1)*sigmaPDFBJetPx;
	     (*ParMax)[FPARAM_BJetLepPx_TOPLEP] = (*PxBJet1) + LimNRMS_*fabs(*PxBJet1)*sigmaPDFBJetPx;
	     (*ParMin)[FPARAM_BJetLepPy_TOPLEP] = (*PyBJet1) - LimNRMS_*fabs(*PyBJet1)*sigmaPDFBJetPy;
	     (*ParMax)[FPARAM_BJetLepPy_TOPLEP] = (*PyBJet1) + LimNRMS_*fabs(*PyBJet1)*sigmaPDFBJetPy;
	     (*ParMin)[FPARAM_BJetLepPz_TOPLEP] = (*PzBJet1) - LimNRMS_*fabs(*PzBJet1)*sigmaPDFBJetPz;
	     (*ParMax)[FPARAM_BJetLepPz_TOPLEP] = (*PzBJet1) + LimNRMS_*fabs(*PzBJet1)*sigmaPDFBJetPz;

	     chi_[NPerm_] = INIT;
			    
	     for( int ipho=nPhoton;ipho>=0;ipho-- )
	       {
		  if( ipho < nPhoton ) // permute according to the number of photons + permutation with no radiation (ISR)
		    {	     
		       *EPhoton  = PhotonE[ipho];
		       *PxPhoton = PhotonPx[ipho];
		       *PyPhoton = PhotonPy[ipho];
		       *PzPhoton = PhotonPz[ipho];
				      
		       (*ParMin)[FPARAM_PhotonPx_TOPLEP] = (*PxPhoton) - LimNRMS_*fabs(*PxPhoton)*sigmaPDFPhotonPx;
		       (*ParMax)[FPARAM_PhotonPx_TOPLEP] = (*PxPhoton) + LimNRMS_*fabs(*PxPhoton)*sigmaPDFPhotonPx;
		       (*ParMin)[FPARAM_PhotonPy_TOPLEP] = (*PyPhoton) - LimNRMS_*fabs(*PyPhoton)*sigmaPDFPhotonPy;
		       (*ParMax)[FPARAM_PhotonPy_TOPLEP] = (*PyPhoton) + LimNRMS_*fabs(*PyPhoton)*sigmaPDFPhotonPy;
		       (*ParMin)[FPARAM_PhotonPz_TOPLEP] = (*PzPhoton) - LimNRMS_*fabs(*PzPhoton)*sigmaPDFPhotonPz;
		       (*ParMax)[FPARAM_PhotonPz_TOPLEP] = (*PzPhoton) + LimNRMS_*fabs(*PzPhoton)*sigmaPDFPhotonPz;
		    }
				 
		  for( int ipo=0;ipo<PHOTON_ORIGIN_N_TOPLEP;ipo++ )
		    {
		       if( ipo != PHOTON_FROM_TOPLEP_COMB_TOPLEP &&
			   ipo != PHOTON_FROM_WLEP_COMB_TOPLEP &&
			   ipo != PHOTON_FROM_ISR_TOPLEP &&
			   !CheckAllPhotonOrigins_ ) continue;
		       
		       if( !IncludePhotons_ && ipo != PHOTON_FROM_ISR_TOPLEP ) continue;
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
				 if( r.par[FPARAM_Sign_TOPLEP] == 1 ) lhgood_sPlus.push_back(r);
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
			    
			    TopLep_PhotonIdx[NPerm_] = (IncludePhotons_) ? ipho : -1;
			    if( IncludePhotons_ && ipho == nPhoton ) TopLep_PhotonIdx[NPerm_] = -1;
					   
			    TopLep_BJetLepIdx[NPerm_] = ib;
					   
			    if( label1_l == 0 )
			      {
				 TopLep_ElectronIdx[NPerm_] = idx1_l;
				 TopLep_MuonIdx[NPerm_] = -1;
			      }
			    else if( label1_l == 1 ) 
			      {
				 TopLep_MuonIdx[NPerm_] = idx1_l;
				 TopLep_ElectronIdx[NPerm_] = -1;
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
					   
			    nuPx_[NPerm_][0] = fbest.par[FPARAM_EtRealX_TOPLEP];
			    nuPy_[NPerm_][0] = fbest.par[FPARAM_EtRealY_TOPLEP];
			    nuPz_[NPerm_][0] = fbest.PzNu1;
			    
			    MetPx_[NPerm_] = fbest.par[FPARAM_EtRealX_TOPLEP];
			    MetPy_[NPerm_] = fbest.par[FPARAM_EtRealY_TOPLEP];
			    
			    WMass_[NPerm_][0] = fbest.par[FPARAM_mWLep_TOPLEP];
					   
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
     } // end permutations

   // choose the best permutation
   for(int ip=0;ip<NPerm_;ip++)
     {
	calcVar(ip);
     }
}

/// Calculate the number of permutations per event
void KINFIT::TopLep::CalcNPerm()
{
   NPERMEVENT = 0;

   for(int ip=nPhoton;ip>=0;ip--)
     {	
	for(int il=0;il<nLepton;il++)
	  {	
	     for(int ib=0;ib<nBJet;ib++)
	       {
		  for(int ipo=0;ipo<PHOTON_ORIGIN_N_TOPLEP;ipo++)
		    {
		       if( ipo != PHOTON_FROM_TOPLEP_COMB_TOPLEP &&
			   ipo != PHOTON_FROM_WLEP_COMB_TOPLEP &&
			   ipo != PHOTON_FROM_ISR_TOPLEP &&
			   !CheckAllPhotonOrigins_ ) continue;
				      
		       if( !IncludePhotons_ && ip > 0 ) continue;
		       
		       NPERMEVENT++;
		    }
	       }
	  }
     }
}

/// Derive solutions on a given set of inputs and calculate the resultant NLL
double KINFIT::TopLep::func(float PxLepton, float PyLepton, float PzLepton, float ELepton, int LabelLepton,
			    float PxBJet1, float PyBJet1, float PzBJet1, float EBJet1,
			    float PxPhoton, float PyPhoton, float PzPhoton, float EPhoton,
			    int photonOrigin,
			    std::vector<double> &chi2t, double *par)
{
   chi2t[NLL_RacPT_TOPLEP] = 0.;
   chi2t[NLL_mWPT_TOPLEP] = 0.;
   chi2t[NLL_mTopPT_TOPLEP] = 0.;
   
   float PxNu1 = par[FPARAM_EtRealX_TOPLEP];
   float PyNu1 = par[FPARAM_EtRealY_TOPLEP];

   *KINFIT::kfit::PxNu1 = PxNu1;
   *KINFIT::kfit::PyNu1 = PyNu1;
   
   float LeptonPx = par[FPARAM_LeptonPx_TOPLEP];
   float LeptonPy = par[FPARAM_LeptonPy_TOPLEP];
   float LeptonPz = par[FPARAM_LeptonPz_TOPLEP];
   float LeptonE = par[FPARAM_LeptonE_TOPLEP];

   float PhotonPx = par[FPARAM_PhotonPx_TOPLEP];
   float PhotonPy = par[FPARAM_PhotonPy_TOPLEP];
   float PhotonPz = par[FPARAM_PhotonPz_TOPLEP];
   float PhotonE = par[FPARAM_PhotonE_TOPLEP];
   
   if( (photonOrigin == PHOTON_FROM_WLEP_TOPLEP ||
	photonOrigin == PHOTON_FROM_WLEP_COMB_TOPLEP ||
	photonOrigin == PHOTON_FROM_LEPTON_TOPLEP) &&
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

//   float c = par[FPARAM_mWLep_TOPLEP]*par[FPARAM_mWLep_TOPLEP]/2.+LeptonPx*PxNu1+LeptonPy*PyNu1;
   float c = (par[FPARAM_mWLep_TOPLEP]*par[FPARAM_mWLep_TOPLEP]-(LeptonE*LeptonE-LeptonPx*LeptonPx-LeptonPy*LeptonPy-LeptonPz*LeptonPz))/2.+LeptonPx*PxNu1+LeptonPy*PyNu1;
   
   float rac = c*c*b*b-a*a*(d*d*f*f-c*c);
   
   float racAbs = fabs(rac);
   
   float racPos = rac;
   
   if( racPos < 0. ) racPos = 0.;

   float PzNu1 = (c*b+par[FPARAM_Sign_TOPLEP]*sqrt(racPos))/a/a;
   
   *KINFIT::kfit::PzNu1 = PzNu1;
   
   float ENu1 = sqrt(PxNu1*PxNu1+PyNu1*PyNu1+PzNu1*PzNu1);

   float BJetLepPx = par[FPARAM_BJetLepPx_TOPLEP];
   float BJetLepPy = par[FPARAM_BJetLepPy_TOPLEP];
   float BJetLepPz = par[FPARAM_BJetLepPz_TOPLEP];
   float BJetLepE = par[FPARAM_BJetLepE_TOPLEP];

   if( photonOrigin == PHOTON_FROM_BJETLEP_TOPLEP && IncludePhotons_ )
     {
	BJetLepPx += PhotonPx;
	BJetLepPy += PhotonPy;
	BJetLepPz += PhotonPz;
	BJetLepE += PhotonE;
     }
   
   float totPxLep = BJetLepPx+PxNu1+LeptonPx;
   float totPyLep = BJetLepPy+PyNu1+LeptonPy;
   float totPzLep = BJetLepPz+PzNu1+LeptonPz;
   float totELep = BJetLepE+ENu1+LeptonE;
   
   if( (photonOrigin == PHOTON_FROM_TOPLEP_TOPLEP ||
	photonOrigin == PHOTON_FROM_TOPLEP_COMB_TOPLEP) && IncludePhotons_ )
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

   double val = 0.;

   /// Penalty term for negative roots
   if( rac < 0 ) 
     {
	float racPT = log(1.+racAbs);
	val += racPT;
	chi2t[NLL_RacPT_TOPLEP] += racPT;
     }   

   float mWLepAbs = fabs(mWLep);
   float mtopLepAbs = fabs(mtopLep);
   
   /// Penalty terms for negative masses
   if( mWLep < 0 )
     {
	float mWLepPT = log(1.+mWLepAbs);
	val += mWLepPT;
	chi2t[NLL_mWPT_TOPLEP] += mWLepPT;
     }   
   if( mtopLep < 0 )
     {
	float mTopLepPT = log(1.+mtopLepAbs);
	val += mTopLepPT;
	chi2t[NLL_mTopPT_TOPLEP] += mTopLepPT;
     }   
   
   mtopLepAbs = sqrt(mtopLepAbs);
   mWLepAbs = sqrt(mWLepAbs);

   float mWLepProb = getProb(hPDFTopWMass.get(), mWLepAbs, maxPDFTopWMass, xminPDFTopWMass, xmaxPDFTopWMass);

   float mTopLepProb = getProb(hPDFTopMass.get(), mtopLepAbs, maxPDFTopMass, xminPDFTopMass, xmaxPDFTopMass);

   float MetPxProb = (! (*IsParFixed)[FPARAM_EtRealX_TOPLEP]) ? getProb(hPDFMetPx.get(), (par[FPARAM_EtRealX_TOPLEP]-*EtMissX)/par[FPARAM_EtRealX_TOPLEP], maxPDFMetPx, xminPDFMetPx, xmaxPDFMetPx) : 1.0;
   float MetPyProb = (! (*IsParFixed)[FPARAM_EtRealY_TOPLEP]) ? getProb(hPDFMetPy.get(), (par[FPARAM_EtRealY_TOPLEP]-*EtMissY)/par[FPARAM_EtRealY_TOPLEP], maxPDFMetPy, xminPDFMetPy, xmaxPDFMetPy) : 1.0;
   
   float BJetLepPxProb = (! (*IsParFixed)[FPARAM_BJetLepPx_TOPLEP]) ? getProb(hPDFBJetPx.get(), (par[FPARAM_BJetLepPx_TOPLEP]-PxBJet1)/par[FPARAM_BJetLepPx_TOPLEP], maxPDFBJetPx, xminPDFBJetPx, xmaxPDFBJetPx) : 1.0;
   float BJetLepPyProb = (! (*IsParFixed)[FPARAM_BJetLepPy_TOPLEP]) ? getProb(hPDFBJetPy.get(), (par[FPARAM_BJetLepPy_TOPLEP]-PyBJet1)/par[FPARAM_BJetLepPy_TOPLEP], maxPDFBJetPy, xminPDFBJetPy, xmaxPDFBJetPy) : 1.0;
   float BJetLepPzProb = (! (*IsParFixed)[FPARAM_BJetLepPz_TOPLEP]) ? getProb(hPDFBJetPz.get(), (par[FPARAM_BJetLepPz_TOPLEP]-PzBJet1)/par[FPARAM_BJetLepPz_TOPLEP], maxPDFBJetPz, xminPDFBJetPz, xmaxPDFBJetPz) : 1.0;

   float PhotonPxProb = (! (*IsParFixed)[FPARAM_PhotonPx_TOPLEP] && IncludePhotons_) ? getProb(hPDFPhotonPx.get(), (par[FPARAM_PhotonPx_TOPLEP]-PxPhoton)/par[FPARAM_PhotonPx_TOPLEP], maxPDFPhotonPx, xminPDFPhotonPx, xmaxPDFPhotonPx) : 1.0;
   float PhotonPyProb = (! (*IsParFixed)[FPARAM_PhotonPy_TOPLEP] && IncludePhotons_) ? getProb(hPDFPhotonPy.get(), (par[FPARAM_PhotonPy_TOPLEP]-PyPhoton)/par[FPARAM_PhotonPy_TOPLEP], maxPDFPhotonPy, xminPDFPhotonPy, xmaxPDFPhotonPy) : 1.0;
   float PhotonPzProb = (! (*IsParFixed)[FPARAM_PhotonPz_TOPLEP] && IncludePhotons_) ? getProb(hPDFPhotonPz.get(), (par[FPARAM_PhotonPz_TOPLEP]-PzPhoton)/par[FPARAM_PhotonPz_TOPLEP], maxPDFPhotonPz, xminPDFPhotonPz, xmaxPDFPhotonPz) : 1.0;
   
   float LeptonPxProb = 1.0;
   float LeptonPyProb = 1.0;
   float LeptonPzProb = 1.0;
   
   if (! (*IsParFixed)[FPARAM_LeptonPx_TOPLEP]) LeptonPxProb = (LabelLepton == 0) ? getProb(hPDFElecPx.get(), (par[FPARAM_LeptonPx_TOPLEP]-PxLepton)/par[FPARAM_LeptonPx_TOPLEP], maxPDFElecPx, xminPDFElecPx, xmaxPDFElecPx) : getProb(hPDFMuonPx.get(), (par[FPARAM_LeptonPx_TOPLEP]-PxLepton)/par[FPARAM_LeptonPx_TOPLEP], maxPDFMuonPx, xminPDFMuonPx, xmaxPDFMuonPx);
   if (! (*IsParFixed)[FPARAM_LeptonPy_TOPLEP]) LeptonPyProb = (LabelLepton == 0) ? getProb(hPDFElecPy.get(), (par[FPARAM_LeptonPy_TOPLEP]-PyLepton)/par[FPARAM_LeptonPy_TOPLEP], maxPDFElecPy, xminPDFElecPy, xmaxPDFElecPy) : getProb(hPDFMuonPy.get(), (par[FPARAM_LeptonPy_TOPLEP]-PyLepton)/par[FPARAM_LeptonPy_TOPLEP], maxPDFMuonPy, xminPDFMuonPy, xmaxPDFMuonPy);
   if (! (*IsParFixed)[FPARAM_LeptonPz_TOPLEP]) LeptonPzProb = (LabelLepton == 0) ? getProb(hPDFElecPz.get(), (par[FPARAM_LeptonPz_TOPLEP]-PzLepton)/par[FPARAM_LeptonPz_TOPLEP], maxPDFElecPz, xminPDFElecPz, xmaxPDFElecPz) : getProb(hPDFMuonPz.get(), (par[FPARAM_LeptonPz_TOPLEP]-PzLepton)/par[FPARAM_LeptonPz_TOPLEP], maxPDFMuonPz, xminPDFMuonPz, xmaxPDFMuonPz);

   double minProb = 1E-20; // NLL = 92.1034
   
   chi2t[NLL_WLep_TOPLEP] = (mWLepProb > minProb) ? mWLepProb : minProb;
   chi2t[NLL_TopLep_TOPLEP] = (mTopLepProb > minProb) ? mTopLepProb : minProb;
   chi2t[NLL_EtMissX_TOPLEP] = (MetPxProb > minProb) ? MetPxProb : minProb;
   chi2t[NLL_EtMissY_TOPLEP] = (MetPyProb > minProb) ? MetPyProb : minProb;
   chi2t[NLL_BJetLepPx_TOPLEP] = (BJetLepPxProb > minProb) ? BJetLepPxProb : minProb;
   chi2t[NLL_BJetLepPy_TOPLEP] = (BJetLepPyProb > minProb) ? BJetLepPyProb : minProb; 
   chi2t[NLL_BJetLepPz_TOPLEP] = (BJetLepPzProb > minProb) ? BJetLepPzProb : minProb;
   chi2t[NLL_LeptonPx_TOPLEP] = (LeptonPxProb > minProb) ? LeptonPxProb : minProb;
   chi2t[NLL_LeptonPy_TOPLEP] = (LeptonPyProb > minProb) ? LeptonPyProb : minProb;
   chi2t[NLL_LeptonPz_TOPLEP] = (LeptonPzProb > minProb) ? LeptonPzProb : minProb;
   chi2t[NLL_PhotonPx_TOPLEP] = (PhotonPxProb > minProb) ? PhotonPxProb : minProb;
   chi2t[NLL_PhotonPy_TOPLEP] = (PhotonPyProb > minProb) ? PhotonPyProb : minProb;
   chi2t[NLL_PhotonPz_TOPLEP] = (PhotonPzProb > minProb) ? PhotonPzProb : minProb;
   
   double lh = 1.0;
   
   lh *= chi2t[NLL_WLep_TOPLEP];
   lh *= chi2t[NLL_TopLep_TOPLEP];

   /// Include additional terms only if is corresponding parameter is free
   if(! (*IsParFixed)[FPARAM_EtRealX_TOPLEP]) lh *= chi2t[NLL_EtMissX_TOPLEP];
   if(! (*IsParFixed)[FPARAM_EtRealY_TOPLEP]) lh *= chi2t[NLL_EtMissY_TOPLEP];
   if(! (*IsParFixed)[FPARAM_BJetLepPx_TOPLEP]) lh *= chi2t[NLL_BJetLepPx_TOPLEP];
   if(! (*IsParFixed)[FPARAM_BJetLepPy_TOPLEP]) lh *= chi2t[NLL_BJetLepPy_TOPLEP];
   if(! (*IsParFixed)[FPARAM_BJetLepPz_TOPLEP]) lh *= chi2t[NLL_BJetLepPz_TOPLEP];
   if(! (*IsParFixed)[FPARAM_LeptonPx_TOPLEP]) lh *= chi2t[NLL_LeptonPx_TOPLEP];
   if(! (*IsParFixed)[FPARAM_LeptonPy_TOPLEP]) lh *= chi2t[NLL_LeptonPy_TOPLEP];
   if(! (*IsParFixed)[FPARAM_LeptonPz_TOPLEP]) lh *= chi2t[NLL_LeptonPz_TOPLEP];
   if(! (*IsParFixed)[FPARAM_PhotonPx_TOPLEP]) lh *= chi2t[NLL_PhotonPx_TOPLEP];
   if(! (*IsParFixed)[FPARAM_PhotonPy_TOPLEP]) lh *= chi2t[NLL_PhotonPy_TOPLEP];
   if(! (*IsParFixed)[FPARAM_PhotonPz_TOPLEP]) lh *= chi2t[NLL_PhotonPz_TOPLEP];

   val += -2.*log(lh);

   return val;
}

/// Calculate output variables for a given permutation
void KINFIT::TopLep::calcVar(int iPerm)
{
   if( chi_[iPerm] > 10E+9 ) return;

   int idxElec = TopLep_ElectronIdx[iPerm];
   int idxMuon = TopLep_MuonIdx[iPerm];
   
   int idxBJetLep = TopLep_BJetLepIdx[iPerm];
   
   int idxPhoton = TopLep_PhotonIdx[iPerm];
   int photonOrigin = PhotonOrigin_[iPerm];
   
   float LeptonPx = 0;
   float LeptonPy = 0;
   float LeptonPz = 0;
   float LeptonE = 0;

   float BJetLepPx = 0;
   float BJetLepPy = 0;
   float BJetLepPz = 0;
   float BJetLepE = 0;

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
   
   if( photonOrigin == PHOTON_FROM_LEPTON_TOPLEP && (idxElec >= 0 || idxMuon >= 0) )
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
	
	if( photonOrigin == PHOTON_FROM_BJETLEP_TOPLEP )
	  {
	     BJetLepPx += RadPhotonPx;
	     BJetLepPy += RadPhotonPy;
	     BJetLepPz += RadPhotonPz;
	     BJetLepE += RadPhotonE;
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
   
   if( photonOrigin == PHOTON_FROM_WLEP_TOPLEP ||
       photonOrigin == PHOTON_FROM_WLEP_COMB_TOPLEP )
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

   if( photonOrigin == PHOTON_FROM_TOPLEP_TOPLEP ||
       photonOrigin == PHOTON_FROM_TOPLEP_COMB_TOPLEP )
     {
	topLepPx += RadPhotonPx;
	topLepPy += RadPhotonPy;
	topLepPz += RadPhotonPz;
	topLepE += RadPhotonE;
     }
   
   float topLepPt = sqrt(topLepPx*topLepPx+topLepPy*topLepPy);
   float topLepEta = getEta(topLepPt, topLepPz);
   float topLepPhi = atan2(topLepPy, topLepPx);

   TopMass_[iPerm][0] = sqrt(topLepE*topLepE-topLepPx*topLepPx-topLepPy*topLepPy-topLepPz*topLepPz);
   TopPt_[iPerm][0] = sqrt(topLepPx*topLepPx+topLepPy*topLepPy);
   TopP_[iPerm][0] = sqrt(topLepPx*topLepPx+topLepPy*topLepPy+topLepPz*topLepPz);
   TopEta_[iPerm][0] = getEta(TopPt_[iPerm][0],topLepPz);
   TopPhi_[iPerm][0] = topLepPhi;
   TopE_[iPerm][0] = topLepE;
   TopPx_[iPerm][0] = topLepPx;
   TopPy_[iPerm][0] = topLepPy;
   TopPz_[iPerm][0] = topLepPz;

   WMass_[iPerm][0] = sqrt(WLepE*WLepE-WLepPx*WLepPx-WLepPy*WLepPy-WLepPz*WLepPz);
   WPt_[iPerm][0] = sqrt(WLepPx*WLepPx+WLepPy*WLepPy);
   WP_[iPerm][0] = sqrt(WLepPx*WLepPx+WLepPy*WLepPy+WLepPz*WLepPz);
   WEta_[iPerm][0] = getEta(WPt_[iPerm][0],WLepPz);
   WPhi_[iPerm][0] = atan2(WLepPy, WLepPx);
   WE_[iPerm][0] = WLepE;
   WPx_[iPerm][0] = WLepPx;
   WPy_[iPerm][0] = WLepPy;
   WPz_[iPerm][0] = WLepPz;
}

/// FCN definition
void KINFIT::TopLep::fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
   std::string chi2tNames[NLL_N_TOPLEP] = {"WLep", "TopLep", "EtMissX", "EtMissY", "BJetLepPx", "BJetLepPy", "BJetLepPz", "LeptonPx", "LeptonPy", "LeptonPz", "RacPT", "mWPT", "mTopPT"};
   std::vector<double> chi2t(NLL_N_TOPLEP);
   
   double lh = func(*PxLepton1, *PyLepton1, *PzLepton1, *ELepton1, *LabelLepton1,
		    *PxBJet1, *PyBJet1, *PzBJet1, *EBJet1,
		    *PxPhoton, *PyPhoton, *PzPhoton, *EPhoton,
		    *PhotonOrigin,
		    chi2t, par);

   /// converged
   if( iflag == 3 )
     {
	*CHISQ = lh;

	for( int ip=0;ip<FPARAM_N;ip++ ) (*FitParam).push_back(par[ip]);
	for( int it=0;it<NLL_N_TOPLEP;it++ ) {(*ChiTerm).push_back(chi2t[it]); (*ChiTermName).push_back(chi2tNames[it]);}
     }
   
   f = lh;
}

/// Run the fit
void KINFIT::TopLep::fit(double *par, std::vector<FRESULT> &vpp)
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
   (*ParMin)[FPARAM_BJetLepE_TOPLEP] = par[FPARAM_BJetLepE_TOPLEP]-1000.; (*ParMax)[FPARAM_BJetLepE_TOPLEP] = par[FPARAM_BJetLepE_TOPLEP]+1000.;
   (*ParMin)[FPARAM_LeptonE_TOPLEP] = par[FPARAM_LeptonE_TOPLEP]-1000.; (*ParMax)[FPARAM_LeptonE_TOPLEP] = par[FPARAM_LeptonE_TOPLEP]+1000.;
   (*ParMin)[FPARAM_PhotonE_TOPLEP] = par[FPARAM_PhotonE_TOPLEP]-1000.; (*ParMax)[FPARAM_PhotonE_TOPLEP] = par[FPARAM_PhotonE_TOPLEP]+1000.;
      
   gMinuit->mnparm(FPARAM_Sign_TOPTOPLEPHAD, FPARAM_NAME[FPARAM_Sign_TOPLEP], par[FPARAM_Sign_TOPLEP], 1E-1, (*ParMin)[FPARAM_Sign_TOPLEP], (*ParMax)[FPARAM_Sign_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_EtRealX_TOPLEP, FPARAM_NAME[FPARAM_EtRealX_TOPLEP], par[FPARAM_EtRealX_TOPLEP], 1E-2, (*ParMin)[FPARAM_EtRealX_TOPLEP], (*ParMax)[FPARAM_EtRealX_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_EtRealY_TOPLEP, FPARAM_NAME[FPARAM_EtRealY_TOPLEP], par[FPARAM_EtRealY_TOPLEP], 1E-2, (*ParMin)[FPARAM_EtRealY_TOPLEP], (*ParMax)[FPARAM_EtRealY_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_mWLep_TOPLEP, FPARAM_NAME[FPARAM_mWLep_TOPLEP], par[FPARAM_mWLep_TOPLEP], 1E-4, (*ParMin)[FPARAM_mWLep_TOPLEP], (*ParMax)[FPARAM_mWLep_TOPLEP], ierflg);

   gMinuit->mnparm(FPARAM_BJetLepPx_TOPLEP, FPARAM_NAME[FPARAM_BJetLepPx_TOPLEP], par[FPARAM_BJetLepPx_TOPLEP], 1E-1, (*ParMin)[FPARAM_BJetLepPx_TOPLEP], (*ParMax)[FPARAM_BJetLepPx_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_BJetLepPy_TOPLEP, FPARAM_NAME[FPARAM_BJetLepPy_TOPLEP], par[FPARAM_BJetLepPy_TOPLEP], 1E-1, (*ParMin)[FPARAM_BJetLepPy_TOPLEP], (*ParMax)[FPARAM_BJetLepPy_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_BJetLepPz_TOPLEP, FPARAM_NAME[FPARAM_BJetLepPz_TOPLEP], par[FPARAM_BJetLepPz_TOPLEP], 1E-1, (*ParMin)[FPARAM_BJetLepPz_TOPLEP], (*ParMax)[FPARAM_BJetLepPz_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_BJetLepE_TOPLEP, FPARAM_NAME[FPARAM_BJetLepE_TOPLEP], par[FPARAM_BJetLepE_TOPLEP], 1E-1, (*ParMin)[FPARAM_BJetLepE_TOPLEP], (*ParMax)[FPARAM_BJetLepE_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_LeptonPx_TOPLEP, FPARAM_NAME[FPARAM_LeptonPx_TOPLEP], par[FPARAM_LeptonPx_TOPLEP], 1E-1, (*ParMin)[FPARAM_LeptonPx_TOPLEP], (*ParMax)[FPARAM_LeptonPx_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_LeptonPy_TOPLEP, FPARAM_NAME[FPARAM_LeptonPy_TOPLEP], par[FPARAM_LeptonPy_TOPLEP], 1E-1, (*ParMin)[FPARAM_LeptonPy_TOPLEP], (*ParMax)[FPARAM_LeptonPy_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_LeptonPz_TOPLEP, FPARAM_NAME[FPARAM_LeptonPz_TOPLEP], par[FPARAM_LeptonPz_TOPLEP], 1E-1, (*ParMin)[FPARAM_LeptonPz_TOPLEP], (*ParMax)[FPARAM_LeptonPz_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_LeptonE_TOPLEP, FPARAM_NAME[FPARAM_LeptonE_TOPLEP], par[FPARAM_LeptonE_TOPLEP], 1E-1, (*ParMin)[FPARAM_LeptonE_TOPLEP], (*ParMax)[FPARAM_LeptonE_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_PhotonPx_TOPLEP, FPARAM_NAME[FPARAM_PhotonPx_TOPLEP], par[FPARAM_PhotonPx_TOPLEP], 1E-1, (*ParMin)[FPARAM_PhotonPx_TOPLEP], (*ParMax)[FPARAM_PhotonPx_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_PhotonPy_TOPLEP, FPARAM_NAME[FPARAM_PhotonPy_TOPLEP], par[FPARAM_PhotonPy_TOPLEP], 1E-1, (*ParMin)[FPARAM_PhotonPy_TOPLEP], (*ParMax)[FPARAM_PhotonPy_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_PhotonPz_TOPLEP, FPARAM_NAME[FPARAM_PhotonPz_TOPLEP], par[FPARAM_PhotonPz_TOPLEP], 1E-1, (*ParMin)[FPARAM_PhotonPz_TOPLEP], (*ParMax)[FPARAM_PhotonPz_TOPLEP], ierflg);
   gMinuit->mnparm(FPARAM_PhotonE_TOPLEP, FPARAM_NAME[FPARAM_PhotonE_TOPLEP], par[FPARAM_PhotonE_TOPLEP], 1E-1, (*ParMin)[FPARAM_PhotonE_TOPLEP], (*ParMax)[FPARAM_PhotonE_TOPLEP], ierflg);
   
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

   fres.lh = *CHISQ;
   
   vpp.push_back(fres);
   
   delete gMinuit;
}

/// First-layer minimization (generic)
void KINFIT::TopLep::calcNuGrid(std::vector<FRESULT> &vp)
{
   double par[FPARAM_N];
   for(int p=0;p<FPARAM_N;p++) par[p] = 0.;
   
   std::vector<double> chi2t(NLL_N_TOPLEP);

   bool doToys = bool(NToy_ > 1);
   
   for(int it=0;it<NToy_;it++)
     {		     
	if( vp.size() >= NGrid_ ) break;
	
	double thres = -1E-6;

	do
	  {
	     double gv = (doToys) ? getProbGaus(hPDFBJetPx.get(), maxPDFBJetPx, meanPDFBJetPx, sigmaPDFBJetPx, rnd, NBJetPxRMS_) : 0.;
	     par[FPARAM_BJetLepPx_TOPLEP] = *PxBJet1/(1.-gv);
	     gv = (doToys) ? getProbGaus(hPDFBJetPy.get(), maxPDFBJetPy, meanPDFBJetPy, sigmaPDFBJetPy, rnd, NBJetPyRMS_) : 0.;
	     par[FPARAM_BJetLepPy_TOPLEP] = *PyBJet1/(1.-gv);
	     gv = (doToys) ? getProbGaus(hPDFBJetPz.get(), maxPDFBJetPz, meanPDFBJetPz, sigmaPDFBJetPz, rnd, NBJetPzRMS_) : 0.;
	     par[FPARAM_BJetLepPz_TOPLEP] = *PzBJet1/(1.-gv);
	     par[FPARAM_BJetLepE_TOPLEP] = sqrt((*MassBJet1)*(*MassBJet1) + par[FPARAM_BJetLepPx_TOPLEP]*par[FPARAM_BJetLepPx_TOPLEP] + par[FPARAM_BJetLepPy_TOPLEP]*par[FPARAM_BJetLepPy_TOPLEP] + par[FPARAM_BJetLepPz_TOPLEP]*par[FPARAM_BJetLepPz_TOPLEP]);
	  } while( (par[FPARAM_BJetLepE_TOPLEP]*par[FPARAM_BJetLepE_TOPLEP]-par[FPARAM_BJetLepPx_TOPLEP]*par[FPARAM_BJetLepPx_TOPLEP]-par[FPARAM_BJetLepPy_TOPLEP]*par[FPARAM_BJetLepPy_TOPLEP]-par[FPARAM_BJetLepPz_TOPLEP]*par[FPARAM_BJetLepPz_TOPLEP]) < thres );

	if( IncludePhotons_ )
	  {	     
	     do
	       {
		  double gv = (doToys) ? getProbGaus(hPDFPhotonPx.get(), maxPDFPhotonPx, meanPDFPhotonPx, sigmaPDFPhotonPx, rnd, NPhotonPxRMS_) : 0.;
		  par[FPARAM_PhotonPx_TOPLEP] = *PxPhoton/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFPhotonPy.get(), maxPDFPhotonPy, meanPDFPhotonPy, sigmaPDFPhotonPy, rnd, NPhotonPyRMS_) : 0.;
		  par[FPARAM_PhotonPy_TOPLEP] = *PyPhoton/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFPhotonPz.get(), maxPDFPhotonPz, meanPDFPhotonPz, sigmaPDFPhotonPz, rnd, NPhotonPzRMS_) : 0.;
		  par[FPARAM_PhotonPz_TOPLEP] = *PzPhoton/(1.-gv);
		  par[FPARAM_PhotonE_TOPLEP] = sqrt(par[FPARAM_PhotonPx_TOPLEP]*par[FPARAM_PhotonPx_TOPLEP] + par[FPARAM_PhotonPy_TOPLEP]*par[FPARAM_PhotonPy_TOPLEP] + par[FPARAM_PhotonPz_TOPLEP]*par[FPARAM_PhotonPz_TOPLEP]);
	       } while( (par[FPARAM_PhotonE_TOPLEP]*par[FPARAM_PhotonE_TOPLEP]-par[FPARAM_PhotonPx_TOPLEP]*par[FPARAM_PhotonPx_TOPLEP]-par[FPARAM_PhotonPy_TOPLEP]*par[FPARAM_PhotonPy_TOPLEP]-par[FPARAM_PhotonPz_TOPLEP]*par[FPARAM_PhotonPz_TOPLEP]) < thres );
	  }
	
	  {
	     double gv = (doToys) ? getProbGaus(hPDFMetPx.get(), maxPDFMetPx, meanPDFMetPx, sigmaPDFMetPx, rnd, NMetRMS_) : 0.;
	     par[FPARAM_EtRealX_TOPLEP] = *EtMissX/(1.-gv);
	     gv = (doToys) ? getProbGaus(hPDFMetPy.get(), maxPDFMetPy, meanPDFMetPy, sigmaPDFMetPy, rnd, NMetRMS_) : 0.;
	     par[FPARAM_EtRealY_TOPLEP] = *EtMissY/(1.-gv);
	  }	

	par[FPARAM_mWLep_TOPLEP] = (doToys) ? getProbGaus(hPDFTopWMass.get(), maxPDFTopWMass, meanPDFTopWMass, sigmaPDFTopWMass, rnd, 2) : WMass1;

	if( *LabelLepton1 == 0 )
	  {
	     do
	       {
		  double gv = (doToys) ? getProbGaus(hPDFElecPx.get(), maxPDFElecPx, meanPDFElecPx, sigmaPDFElecPx, rnd, NElecPxRMS_) : 0.;
		  par[FPARAM_LeptonPx_TOPLEP] = *PxLepton1/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFElecPy.get(), maxPDFElecPy, meanPDFElecPy, sigmaPDFElecPy, rnd, NElecPyRMS_) : 0.;
		  par[FPARAM_LeptonPy_TOPLEP] = *PyLepton1/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFElecPz.get(), maxPDFElecPz, meanPDFElecPz, sigmaPDFElecPz, rnd, NElecPzRMS_) : 0.;
		  par[FPARAM_LeptonPz_TOPLEP] = *PzLepton1/(1.-gv);
		  par[FPARAM_LeptonE_TOPLEP] = sqrt((*MassLepton1)*(*MassLepton1) + par[FPARAM_LeptonPx_TOPLEP]*par[FPARAM_LeptonPx_TOPLEP] + par[FPARAM_LeptonPy_TOPLEP]*par[FPARAM_LeptonPy_TOPLEP] + par[FPARAM_LeptonPz_TOPLEP]*par[FPARAM_LeptonPz_TOPLEP]);
	       } while( (par[FPARAM_LeptonE_TOPLEP]*par[FPARAM_LeptonE_TOPLEP]-par[FPARAM_LeptonPx_TOPLEP]*par[FPARAM_LeptonPx_TOPLEP]-par[FPARAM_LeptonPy_TOPLEP]*par[FPARAM_LeptonPy_TOPLEP]-par[FPARAM_LeptonPz_TOPLEP]*par[FPARAM_LeptonPz_TOPLEP]) < thres );
	  }
	else
	  {
	     do
	       {
		  double gv = (doToys) ? getProbGaus(hPDFMuonPx.get(), maxPDFMuonPx, meanPDFMuonPx, sigmaPDFMuonPx, rnd, NMuonPxRMS_) : 0.;
		  par[FPARAM_LeptonPx_TOPLEP] = *PxLepton1/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFMuonPy.get(), maxPDFMuonPy, meanPDFMuonPy, sigmaPDFMuonPy, rnd, NMuonPyRMS_) : 0.;
		  par[FPARAM_LeptonPy_TOPLEP] = *PyLepton1/(1.-gv);
		  gv = (doToys) ? getProbGaus(hPDFMuonPz.get(), maxPDFMuonPz, meanPDFMuonPz, sigmaPDFMuonPz, rnd, NMuonPzRMS_) : 0.;
		  par[FPARAM_LeptonPz_TOPLEP] = *PzLepton1/(1.-gv);
		  par[FPARAM_LeptonE_TOPLEP] = sqrt((*MassLepton1)*(*MassLepton1) + par[FPARAM_LeptonPx_TOPLEP]*par[FPARAM_LeptonPx_TOPLEP] + par[FPARAM_LeptonPy_TOPLEP]*par[FPARAM_LeptonPy_TOPLEP] + par[FPARAM_LeptonPz_TOPLEP]*par[FPARAM_LeptonPz_TOPLEP]);
	       } while( (par[FPARAM_LeptonE_TOPLEP]*par[FPARAM_LeptonE_TOPLEP]-par[FPARAM_LeptonPx_TOPLEP]*par[FPARAM_LeptonPx_TOPLEP]-par[FPARAM_LeptonPy_TOPLEP]*par[FPARAM_LeptonPy_TOPLEP]-par[FPARAM_LeptonPz_TOPLEP]*par[FPARAM_LeptonPz_TOPLEP]) < thres );
	  }
	
	for( int is1=-1;is1<=1;is1++ )
	  {	
	     if( is1 == 0 ) continue;
	     
	     par[FPARAM_Sign_TOPLEP] = is1;

	     chi2t.clear();

	     double lh = func(*PxLepton1, *PyLepton1, *PzLepton1, *ELepton1, *LabelLepton1,
			      *PxBJet1, *PyBJet1, *PzBJet1, *EBJet1,
			      *PxPhoton, *PyPhoton, *PzPhoton, *EPhoton,
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
