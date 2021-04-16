#ifndef TOPLEP_H
#define TOPLEP_H

#include "kinfit.h"

namespace KINFIT
{         
   class TopLep : public KINFIT::kfit
     {
	
      public:
	
	TopLep();
	virtual ~TopLep();
	
	static void fcn(int &npar, double *gin, double &f, double *par, int iflag);

	static double func(
			   float PxLepton1, float PyLepton1, float PzLepton1, float ELepton1, int LabelLepton1,
			   float PxBJet1, float PyBJet1, float PzBJet1, float EBJet1,
			   float PxPhoton, float PyPhoton, float PzPhoton, float EPhoton,
			   int pho,
			   std::vector<double> &chi2t, double *par);
      public:
	
	void TopLepRun();

	void fit(double *par, std::vector<FRESULT> &vpp);
	
	void calcVar(int iPerm);

	void calcNuGrid(std::vector<FRESULT> &vp);
	
      private:

	void CalcNPerm();

     };
}

#endif
