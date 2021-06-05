#ifndef TOPTOPLEPLEP_H
#define TOPTOPLEPLEP_H

#include "kinfit.h"

namespace KINFIT
{         
   class TopTopLepLep : public KINFIT::kfit
     {
	
      public:
	
	TopTopLepLep();
	virtual ~TopTopLepLep();
	
	static void fcn(int &npar, double *gin, double &f, double *par, int iflag);

	static double func(float PxLepton1, float PyLepton1, float PzLepton1, float PtLepton1, float EtaLepton1, float PhiLepton1, float ELepton1, int LabelLepton1,
			   float PxLepton2, float PyLepton2, float PzLepton2, float PtLepton2, float EtaLepton2, float PhiLepton2, float ELepton2, int LabelLepton2,
			   float PxBJet1, float PyBJet1, float PzBJet1, float PtBJet1, float EtaBJet1, float PhiBJet1, float EBJet1,
			   float PxBJet2, float PyBJet2, float PzBJet2, float PtBJet2, float EtaBJet2, float PhiBJet2, float EBJet2,
			   float PxPhoton, float PyPhoton, float PzPhoton, float PtPhoton, float EtaPhoton, float PhiPhoton, float EPhoton,
			   int pho,
			   std::vector<double> &chi2t, double *par);
      public:
	
	void TopTopLepLepRun();

	void fit(double *par, std::vector<FRESULT> &vpp);
	
	void calcVar(int iPerm);

	void calcNuGrid(std::vector<FRESULT> &vp);
	
      private:

	void CalcNPerm();

     };
}

#endif
