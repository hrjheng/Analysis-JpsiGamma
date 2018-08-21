#include "tdrstyle.C"

using namespace RooFit;

bool constparam = false;

void fitExclude_RooFitUnbinnedFit()
{
	setTDRStyle();

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
    CMS_hmmg_mass->setRange("low",100,120);
   	CMS_hmmg_mass->setRange("high",130,150);

    RooDataSet dataset("dataset", " ", t, RooArgSet(*CMS_hmmg_mass), "m_mmg>100.&&m_mmg<150.", 0);

    RooRealVar *fitexclude_p0 = new RooRealVar("fitexclude_p0"," ",2.33867e+01); fitexclude_p0->setConstant(constparam);
    RooRealVar *fitexclude_p1 = new RooRealVar("fitexclude_p1"," ",4.69699e+00); fitexclude_p1->setConstant(constparam);
    RooRealVar *fitexclude_p2 = new RooRealVar("fitexclude_p2"," ",4.27736e+00); fitexclude_p2->setConstant(constparam);

    RooRealVar *fitinclude_p0 = new RooRealVar("fitinclude_p0"," ",2.29339e+01); fitinclude_p0->setConstant(constparam);
    RooRealVar *fitinclude_p1 = new RooRealVar("fitinclude_p1"," ",6.81455e+00); fitinclude_p1->setConstant(constparam);
    RooRealVar *fitinclude_p2 = new RooRealVar("fitinclude_p2"," ",3.73678e+00); fitinclude_p2->setConstant(constparam);

    RooBernstein *Bern_exclude = new RooBernstein("Bern_exclude","",*CMS_hmmg_mass,RooArgList(*fitexclude_p0,*fitexclude_p1,*fitexclude_p2));
    RooBernstein *Bern_include = new RooBernstein("Bern_include","",*CMS_hmmg_mass,RooArgList(*fitinclude_p0,*fitinclude_p1,*fitinclude_p2));

    RooRealVar NBkg_exclude("NBkg_exclude","",t->GetEntriesFast(),t->GetEntriesFast()*0.25,t->GetEntriesFast()*1.75);
    RooRealVar NBkg_include("NBkg_include","",t->GetEntriesFast(),t->GetEntriesFast()*0.25,t->GetEntriesFast()*1.75);

    RooAddPdf model_exclude("model_exclude","",RooArgList(*Bern_exclude),RooArgList(NBkg_exclude));
    RooAddPdf model_include("model_include","",RooArgList(*Bern_include),RooArgList(NBkg_include));

    RooFitResult* r_exclude = model_exclude.fitTo(dataset,RooFit::Range("low,high"), Save(kTRUE));
    RooFitResult* r_include = model_include.fitTo(dataset, Save(kTRUE));

    TCanvas* c = new TCanvas( "c", "", 800, 700);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    RooPlot* xframe = CMS_hmmg_mass->frame() ;
    dataset.plotOn(xframe,Binning(25),XErrorSize(0.0));
    model_exclude.plotOn(xframe,RooFit::Name(model_exclude.GetName()),LineColor(TColor::GetColor("#F6318C")), LineWidth(3));
    model_include.plotOn(xframe,RooFit::Name(model_include.GetName()),LineColor(TColor::GetColor("#0099FF")), LineWidth(3));
    dataset.plotOn(xframe,Binning(25),RooFit::Name("data"), MarkerStyle(20), MarkerSize(1.5),XErrorSize(0.0));
    xframe->SetMinimum(0.00001);
    xframe->GetXaxis()->SetTickSize(0.02);
    xframe->GetYaxis()->SetTickSize(0.02);
    xframe->GetXaxis()->SetTitle("M_{#mu#mu#gamma} (GeV)");
    xframe->GetXaxis()->SetTitleSize(0.04);
    xframe->GetYaxis()->SetTitleSize(0.04);
    xframe->GetXaxis()->SetLabelSize(0.04);
    xframe->GetYaxis()->SetLabelSize(0.04);
    xframe->GetXaxis()->SetTitleOffset(1.1);
    xframe->GetYaxis()->SetRangeUser(0.1,(h_data->GetMaximum())*1.7);
    c->cd();
    xframe->Draw();
    TLine *L1 = new TLine(120,0,120,35);
    L1->SetLineWidth(3);
    L1->SetLineColor(TColor::GetColor("#547C66"));
    L1->Draw("same");
    TLine *L2 = new TLine(130,0,130,35);
    L2->SetLineWidth(3);
    L2->SetLineColor(TColor::GetColor("#547C66"));
    L2->Draw("same");
    TLegend* leg = new TLegend(0.45,0.70,0.85,0.9);
    leg->AddEntry(xframe->findObject("data"), "Data", "lep");
    leg->AddEntry(xframe->findObject(model_exclude.GetName()), "B-only fit (sideband only)", "l");
    leg->AddEntry(xframe->findObject(model_include.GetName()), "B-only fit (full range)", "l");
    leg->Draw("same");
    c->SaveAs("fitExclude_RooFitUnbinnedFit.png");
    c->SaveAs("fitExclude_RooFitUnbinnedFit.pdf");
}