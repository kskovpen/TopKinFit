#ifndef TOPTOPLEPHAD_H
#define TOPTOPLEPHAD_H

#include "kinfit.h"

namespace KINFIT
{         
   class TopTopLepHad : public KINFIT::kfit
     {
	
      public:
	
	TopTopLepHad();
	virtual ~TopTopLepHad();
	
	static void fcn(int &npar, double *gin, double &f, double *par, int iflag);

	static double func(
			   float PxLepton1, float PyLepton1, float PzLepton1, float ELepton1, int LabelLepton1,
			   float PxBJet1, float PyBJet1, float PzBJet1, float EBJet1,
			   float PxBJet2, float PyBJet2, float PzBJet2, float EBJet2,
			   float PxNonBJet1, float PyNonBJet1, float PzNonBJet1, float ENonBJet1,
			   float PxNonBJet2, float PyNonBJet2, float PzNonBJet2, float ENonBJet2,
			   float PxPhoton, float PyPhoton, float PzPhoton, float EPhoton,
			   int pho,
			   std::vector<double> &chi2t, double *par);
      public:
	
	void TopTopLepHadRun();

	void fit(double *par, std::vector<FRESULT> &vpp);
	
	void calcVar(int iPerm);

	void calcNuGrid(std::vector<FRESULT> &vp);
	
      private:

	void CalcNPerm();

     };
}

#endif
