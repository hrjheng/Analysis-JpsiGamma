// Illustrate how to fit excluding points in a given range
// Author: Rene Brun
#include "TH1.h"
#include "TF1.h"
#include "TList.h"
#include "tdrstyle.C"

Bool_t reject;
Double_t fline(Double_t *x, Double_t *par)
{
    if (reject && x[0] > (120.-100.)/(150.-100.) && x[0] < (130.-100.)/(150.-100.)) {
        TF1::RejectPoint();
        return 0;
    }
    //return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
    //return par[0]*par[0]*(1.-x[0])*(1.-x[0])+par[1]*par[1]*(2*x[0]*(1.-x[0]))*(2*x[0]*(1.-x[0]))+par[2]*par[2]*x[0]*x[0];
    return ((par[0]-2*par[1]+par[2])*x[0] + 2*(par[1]-par[0]))*x[0] + par[0];
}

void fitExclude_TF1BinnedFit()
{
    setTDRStyle();
    
    //Create a source function
    TFile *f = new TFile("minitree_2016_HJpsiG_AllRuns.root","READ");
    TTree *t = (TTree*) f->Get("outTree");
    float m;
    t->SetBranchAddress("m_mmg", &m);
    TH1F *h_data = new TH1F("h_data"," ",25,0,1);
    for (Long64_t ev = 0; ev < t->GetEntriesFast(); ev++){
        t->GetEntry(ev);
        h_data->Fill((m-100.)/(150.-100.));
    }
    
    TCanvas *c = new TCanvas("c","",800,700);
    c->cd();
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gStyle->SetErrorX(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    //TF1 *f1 = new TF1("f1","[0] +[1]*x +gaus(2)",0,5);
    //f1->SetParameters(6,-1,5,3,0.2);
    // create and fill histogram according to the source function
    //TH1F *h = new TH1F("h","background + signal",100,0,5);
    //h->FillRandom("f1",2000);
    TF1 *fl = new TF1("fl",fline,0,1,3);
    fl->SetLineColor(kRed);
    //fit only the linear background excluding the signal area
    reject = kTRUE;
    h_data->Fit(fl,"0L");
    reject = kFALSE;
    //store 2 separate functions for visualization
    TF1 *fleft = new TF1("fleft",fline,0,(120.-100.)/(150.-100.),3);
    fleft->SetParameters(fl->GetParameters());
    h_data->GetListOfFunctions()->Add(fleft);
    gROOT->GetListOfFunctions()->Remove(fleft);
    TF1 *fright = new TF1("fright",fline,(130.-100.)/(150.-100.),1,3);
    fright->SetParameters(fl->GetParameters());
    h_data->GetListOfFunctions()->Add(fright);
    gROOT->GetListOfFunctions()->Remove(fright);
    fleft->SetLineWidth(3);
    fright->SetLineWidth(3);
    TF1 *fl2 = new TF1("fl2",fline,0,1,3);
    fl2->SetLineColor(kBlue);
    fl2->SetLineWidth(3);
    h_data->Fit(fl2,"+L");
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(2);
    h_data->SetLineWidth(2);
    h_data->GetXaxis()->SetTitle("Transformed M_{#mu#mu#gamma} (GeV)");
    h_data->GetYaxis()->SetTitle(Form("Events / (%.2f GeV)",(float) 1./25.));
    h_data->GetYaxis()->SetRangeUser(0,(h_data->GetMaximum())*1.7);
    h_data->Draw("ep");
    TLegend* leg = new TLegend(0.35,0.6,0.85,0.9);
    leg->AddEntry(h_data, "Data", "lep");
    leg->AddEntry(fleft, Form("#splitline{Sideband-only fit}{(#chi^{2}/NDF = %.1f / %d)}", fl->GetChisquare(),fl->GetNDF()), "l");
    leg->AddEntry(fl2, Form("#splitline{Sideband+Signal region fit}{(#chi^{2}/NDF = %.1f / %d)}",fl2->GetChisquare(),fl2->GetNDF()), "l");
    leg->Draw("same");
    c->SaveAs("fitExclude_TF1BinnedFit.pdf");
    c->SaveAs("fitExclude_TF1BinnedFit.png");
}
