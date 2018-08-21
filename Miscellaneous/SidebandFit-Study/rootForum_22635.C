#include "tdrstyle.C"
#include <iostream>
using namespace std;
using namespace RooFit ;

void rootForum_22635(bool extend = false, bool reduceData = false) \
{ 
   setTDRStyle();

   double minL=0;
   double MAXR=10; 
   double minR=6;
   double MAXL=4;
   double sigL=4;
   double sigR=6;


   // --- Define variable and its ranges
   //RooRealVar *mHNL = new RooRealVar("Lambda0_M","m_{HNL}",1500,6000,"MeV/c^{2}");
   RooRealVar *mHNL = new RooRealVar("Lambda0_M","m_{HNL}",minL,MAXR,"MeV/c^{2}");
   mHNL->setRange("R",minL,MAXR) ;
   mHNL->setRange("R1",minL,MAXL) ;
   mHNL->setRange("R2",minR,MAXR) ;
   mHNL->setRange("Sig",sigL,sigR) ;

   // --- Define fitting function
   RooRealVar *a0 = new RooRealVar("a0", "coefficient of x^1 term",-1.,1.);
   //RooChebychev *p0 = new RooChebychev("p0","p0",*mHNL,RooArgSet(*a0));
   //p0->selectNormalizationRange("R",true);
   RooPolynomial * p0 = new RooPolynomial("p0","p0",*mHNL,RooArgSet(*a0));
   // --- Extend
   RooRealVar *nbkg = new RooRealVar("nbkg","number of bkg in whole region",0,1e8);
   RooExtendPdf *ep0 = new RooExtendPdf("ep0","extended cheb PDF",*p0,*nbkg,"R");

   // --- Load file to dataset
   RooDataSet* mHNL_dataset_all = p0->generate(*mHNL,2000) ;

   mHNL_dataset_all->Print("V");

   RooAbsData * mHNL_dataset = mHNL_dataset_all;
   if (reduceData) { 
      TString dataCut = TString::Format(" Lambda0_M <= %f || Lambda0_M >= %f",MAXL,minR);
      mHNL_dataset = (RooDataSet*) mHNL_dataset_all->reduce(RooFit::Cut(dataCut) );
   }
   std::cout << "data " << mHNL_dataset->numEntries() << std::endl;

   if (mHNL_dataset->numEntries() == 0) return;
  
   nbkg->setVal(mHNL_dataset->numEntries());
   //nbkg->setConstant(true);

   RooAbsPdf * pdf = p0; 
   if (extend) pdf = ep0; 
  
   // --- Fit
   RooFitResult* fitep0 = pdf->fitTo(*mHNL_dataset,RooFit::Range("R1,R2"),RooFit::Save());
   cout<<"  +++  Fit results  +++  "<<endl;
   fitep0->Print();
   fitep0->correlationMatrix().Print();
   cout<<"Fit status= "<<fitep0->status()<<"; covQual= "<<fitep0->covQual()<<endl;

   // --- Plot
   TCanvas* c = new TCanvas( "c", "", 800, 700);
   gPad->SetRightMargin(0.05);
   gPad->SetTopMargin(0.08);
   RooPlot* frame = mHNL->frame(RooFit::Bins(25)) ;
   mHNL_dataset->plotOn(frame,XErrorSize(0.0), MarkerStyle(20), MarkerSize(1.5)) ;
   pdf->plotOn(frame);
   frame->SetMinimum(0.00001);
   frame->GetXaxis()->SetTickSize(0.02);
   frame->GetYaxis()->SetTickSize(0.02);
   frame->GetXaxis()->SetTitle("m_{HNL}");
   frame->GetXaxis()->SetTitleSize(0.05);
   frame->GetYaxis()->SetTitleSize(0.05);
   frame->GetXaxis()->SetLabelSize(0.04);
   frame->GetYaxis()->SetLabelSize(0.04);
   frame->GetXaxis()->SetTitleOffset(1.1);
   frame->GetYaxis()->SetRangeUser(0.1,150);
   c->cd();
   frame->Draw();
   const char* exttex1;
   if(extend) 
      exttex1 = "extended";
   else
      exttex1 = "nonextended";

   c->SaveAs(Form("rootForum_22635_%s.png",exttex1));
   c->SaveAs(Form("rootForum_22635_%s.pdf",exttex1));

   int nevt = mHNL_dataset->numEntries();

   ROOT::Fit::DataRange r;
   r.AddRange(0,4);
   r.AddRange(6,10);
   ROOT::Fit::DataOptions opt; 
   ROOT::Fit::UnBinData d(opt,r,nevt);
//std::cout << "size " << d.Size() << std::endl;

   mHNL_dataset->Print("V");
   mHNL_dataset->get(0)->Print("V");


   std::cout << mHNL_dataset->get(0)->getRealValue("Lambda0_M") << std::endl;

//d.Resize(nevt);
   for (int i =0; i < nevt; ++i) {
      double val = mHNL_dataset->get(i)->getRealValue("Lambda0_M");
      if (r.IsInside(val))
         d.Add(val);
   }

   std::cout << "number of entries " << d.Size() << std::endl;

   TF1 * f1 = new TF1("f1","[Const]*(1.+[A]*x)/(10.+50*[A])",0,10);

   f1->SetParameters(1,1);

//f1->SetNormalized(true);

   ROOT::Fit::Fitter fitter; 

   ROOT::Math::WrappedMultiTF1 wf(*f1,1);

   fitter.SetFunction(wf,false);
   fitter.Config().SetMinimizer("Minuit2");
   fitter.Config().MinimizerOptions().SetPrintLevel(1);
   if (extend) { 
      fitter.Config().ParamsSettings()[1].SetValue(2000.);
   }
   else {    
      fitter.Config().ParamsSettings()[1].Fix();
      fitter.Config().ParamsSettings()[0].SetValue(1.);
   }  
   fitter.LikelihoodFit(d,extend);

   fitter.Result().Print(std::cout);


   std::cout << "\nROOFIT result " << std::endl;
   fitep0->Print();

}
