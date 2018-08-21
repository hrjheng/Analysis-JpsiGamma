#include "tdrstyle.C"

using namespace RooFit;

void Excess_ToyExp(float Nsig, float Nbkg=10000.)
{
	setTDRStyle();

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    RooMsgService::instance().setSilentMode(true);

    TFile *f = new TFile("minitree_2016_HJpsiG_AllRuns.root","READ");
    TTree *t = (TTree*) f->Get("outTree");
    float m;
    t->SetBranchAddress("m_mmg", &m);
    TH1F *h_data = new TH1F("h_data"," ",25,100,150);
    for (Long64_t ev = 0; ev < t->GetEntriesFast(); ev++){
        t->GetEntry(ev);
        h_data->Fill(m);
    }

    RooRealVar *CMS_hmmg_mass = new RooRealVar("m_mmg","m_mmg", 100., 150., "GeV");

    RooDataSet dataset("dataset", " ", t, RooArgSet(*CMS_hmmg_mass), "m_mmg>100.&&m_mmg<150.", 0);

    RooRealVar *gen_p0 = new RooRealVar("gen_p0"," ",2.33867e+01); gen_p0->setConstant(true);
    RooRealVar *gen_p1 = new RooRealVar("gen_p1"," ",4.69699e+00); gen_p1->setConstant(true);
    RooRealVar *gen_p2 = new RooRealVar("gen_p2"," ",4.27736e+00); gen_p2->setConstant(true);
    RooBernstein *Pdf_data = new RooBernstein("Pdf_data","",*CMS_hmmg_mass,RooArgList(*gen_p0,*gen_p1,*gen_p2));
    RooGaussian *Pdf_Sig = new RooGaussian("Pdf_Sig","Pdf_Sig",*CMS_hmmg_mass,RooConst(125.),RooConst(2.0)) ;

    RooRealVar Nsig_gen("Nsig_gen"," ",Nsig);
    RooRealVar Ndata_gen("Ndata_gen"," ",Nbkg);
    RooAddPdf Pdf_totalgen("Pdf_totalgen", " ",RooArgList(*Pdf_data,*Pdf_Sig),RooArgList(Ndata_gen,Nsig_gen));
    Nsig_gen.setConstant(true);
    Ndata_gen.setConstant(true);

    RooDataSet *toy = Pdf_totalgen.generate(RooArgSet(*CMS_hmmg_mass),(Nbkg+Nsig)) ;

    RooRealVar *fit_p0 = new RooRealVar("fit_p0"," ",0,1); //fit_p0->setConstant(false);
    RooRealVar *fit_p1 = new RooRealVar("fit_p1"," ",0,1); //fit_p1->setConstant(false);
    RooRealVar *fit_p2 = new RooRealVar("fit_p2"," ",0,1); //fit_p2->setConstant(false);
    RooBernstein *BkgFit = new RooBernstein("BkgFit","",*CMS_hmmg_mass,RooArgList(*fit_p0,*fit_p1,*fit_p2));

    RooFitResult* r = BkgFit->fitTo(*toy, Save(kTRUE));

    TCanvas* c = new TCanvas( "c", "", 800, 700);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    RooPlot* xframe = CMS_hmmg_mass->frame() ;
    toy->plotOn(xframe,Binning(25),XErrorSize(0.0));
    Pdf_totalgen.plotOn(xframe,RooFit::Name(Pdf_totalgen.GetName()),LineColor(TColor::GetColor("#F6318C")), LineStyle(2), LineWidth(3));
    BkgFit->plotOn(xframe,RooFit::Name(BkgFit->GetName()),LineColor(TColor::GetColor("#0099FF")), LineWidth(3));
    toy->plotOn(xframe,Binning(25),RooFit::Name("Pseudo-data"), MarkerStyle(20), MarkerSize(1.5),XErrorSize(0.0));
    xframe->SetMinimum(0.00001);
    xframe->GetXaxis()->SetTickSize(0.02);
    xframe->GetYaxis()->SetTickSize(0.02);
    xframe->GetXaxis()->SetTitle("M_{#mu#mu#gamma} (GeV)");
    xframe->GetXaxis()->SetTitleSize(0.04);
    xframe->GetYaxis()->SetTitleSize(0.04);
    xframe->GetXaxis()->SetLabelSize(0.04);
    xframe->GetYaxis()->SetLabelSize(0.04);
    xframe->GetXaxis()->SetTitleOffset(1.1);
    xframe->GetYaxis()->SetTitleOffset(1.4);
    xframe->GetYaxis()->SetRangeUser(0.1,xframe->GetMaximum()*1.5);
    c->cd();
    xframe->Draw();
    TLatex* ltx1 = new TLatex();
    ltx1->SetNDC();
    ltx1->SetTextSize(0.04);
    ltx1->DrawLatex(0.22, 0.82, Form("N_{Bkg} = %.1f",Nbkg));
    ltx1->DrawLatex(0.22, 0.76, Form("N_{Signal} = %.1f", Nsig));
    TLegend* leg = new TLegend(0.5,0.65,0.85,0.87);
    leg->AddEntry(xframe->findObject("Pseudo-data"), "Pseudo-data", "lep");
    leg->AddEntry(xframe->findObject(Pdf_totalgen.GetName()), "Toy pdf", "l");
    leg->AddEntry(xframe->findObject(BkgFit->GetName()), "Background-only fit", "l");
    leg->Draw("same");
    c->SaveAs(Form("Fig/Excess_ToyExperiment/Excess_ToyExperiment_Nsig%.0f.png",Nsig));
    c->SaveAs(Form("Fig/Excess_ToyExperiment/Excess_ToyExperiment_Nsig%.0f.pdf",Nsig));
}