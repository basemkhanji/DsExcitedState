#include "TStyle.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TChain.h"
using namespace RooFit;

void D_mass(){
//First we open the input dataset.root file:
TChain* tree_data = new TChain("tree_data","tree_data");
 tree_data -> Add("/eos/lhcb/user/b/bkhanji/data/data2016.root/Ds1Tuple/Ds1Tuple");
	RooRealVar d("D_M","D_M",1860,2050);
     RooDataSet* data = new RooDataSet("data","data",tree_data,RooArgList(d));	//plot only mass
   RooPlot* Eframe = d.frame();
   data -> plotOn(Eframe);
    Eframe -> Draw();
 }
