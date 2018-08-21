//Example taken from :  http://hep1.phys.ntu.edu.tw/~kfjack/lecture/hepstat/

#include "tdrstyle.C"
TH1D *hist = 0;
TNtupleD *nt = 0;
Bool_t reject;
Bool_t splusbfit;

double model(double x, double *par)
{
    double bprime = par[0];
    double fs     = par[1]/nt->GetEntries();
    double mu     = par[2];
    double sigma  = par[3];
    
    double norm  = 1./sqrt(2.*TMath::Pi())/sigma;
    double G     = norm*exp(-0.5 * pow((x-mu)/sigma,2));

    if(splusbfit)
        return fs * G + (1.-fs) * (1 + bprime*x)/(2. + 2.*bprime);
    else
    {
        if (reject && x > 0.8 && x < 1.2) 
        {
            TF1::RejectPoint();
            return 0;
        }
        return (1.-fs) * (1 + bprime*x)/(2. + 2.*bprime);
    }
    
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
    f = 0.;
    for (int i=0;i<nt->GetEntries();i++) {
        nt->GetEntry(i);
        double *mass = nt->GetArgs();
        
        double L = model(*mass,par);
        if (L>0.) f -= 2.*log(L);
        else { f = HUGE; return; }
    }
}

void example_06_modified()
{
    setTDRStyle();

    TFile *fin = new TFile("example_data.root");
    hist = (TH1D *)fin->Get("hist");
    nt = (TNtupleD *)fin->Get("nt");
    
    hist->SetStats(false);
    
    // Sideband-only fit
    TMinuit *gMinuit = new TMinuit(4);
    reject = kTRUE;
    splusbfit = kFALSE;
    gMinuit->SetFCN(fcn);
    gMinuit->DefineParameter(0, "bprime",-0.3, 1.,  -10.,   10.);
    gMinuit->DefineParameter(1, "area", 2000., 1.,    0.,20000.);
    gMinuit->DefineParameter(2, "mean",  1.00, 1.,   0.5,   1.5);
    gMinuit->DefineParameter(3, "width", 0.05, 1., 0.001,  0.15);
    gMinuit->Command("MIGRAD");
    gMinuit->Command("MIGRAD");
    gMinuit->Command("MINOS");

    double par[4],err[4];
    for(int i=0;i<4;i++) gMinuit->GetParameter(i,par[i],err[i]);
    
    TH1F* curve = new TH1F("curve","curve",hist->GetNbinsX()*100,
                           hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    
    for(int i=1;i<=curve->GetNbinsX();i++) {
        reject = kFALSE;
        double x = curve->GetBinCenter(i);
        double f = model(x,par);
        double BinWidth = hist->GetBinWidth(1);
        curve->SetBinContent(i,f*nt->GetEntries()*BinWidth);
    }
    TCanvas *c = new TCanvas("c","",800,700);
    c->cd();
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gStyle->SetErrorX(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    hist->SetLineColor(1);
    hist->SetMarkerStyle(21);
    hist->SetMarkerSize(0.5);
    hist->SetLineWidth(2);
    //hist->Draw("ep");
    curve->SetLineWidth(3);
    curve->SetLineColor(TColor::GetColor("#F96C48"));
    curve->SetFillColor(TColor::GetColor("#F96C48"));
    curve->GetXaxis()->SetTitle("Mass");
    curve->GetYaxis()->SetTitle("Entries");
    curve->GetYaxis()->SetRangeUser(0.,hist->GetMaximum()*1.5);
    curve->Draw();

    //Include signal component & signal+background fit
    TMinuit *gMinuit2 = new TMinuit(4);
    reject = kFALSE;
    splusbfit = kTRUE;
    gMinuit2->SetFCN(fcn);
    gMinuit2->DefineParameter(0, "bprime",-0.3, 1.,  -10.,   10.);
    gMinuit2->DefineParameter(1, "area", 2000., 1.,    0.,20000.);
    gMinuit2->DefineParameter(2, "mean",  1.00, 1.,   0.5,   1.5);
    gMinuit2->DefineParameter(3, "width", 0.05, 1., 0.001,  0.15);
    gMinuit2->Command("MIGRAD");
    gMinuit2->Command("MIGRAD");
    gMinuit2->Command("MINOS");

    //double par[4],err[4];
    for(int i=0;i<4;i++) gMinuit2->GetParameter(i,par[i],err[i]);
    
    TH1F* curve2 = new TH1F("curve","curve",hist->GetNbinsX()*5,
                           hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    
    for(int i=1;i<=curve2->GetNbinsX();i++) {
        reject = kFALSE;
        double x = curve2->GetBinCenter(i);
        double f = model(x,par);
        double BinWidth = hist->GetBinWidth(1);
        curve2->SetBinContent(i,f*nt->GetEntries()*BinWidth);
    }
    curve2->SetLineWidth(3);
    curve2->SetLineColor(TColor::GetColor("#466C95"));
    curve2->Draw("Csame");
    hist->Draw("epsame");

    TLegend* leg = new TLegend(0.5,0.72,0.9,0.85);
    leg->AddEntry(hist, "Pseudo-data", "lep");
    leg->AddEntry(curve, "Sideband-only fit", "fl");
    leg->AddEntry(curve2, "Signal+Background fit", "l");
    leg->Draw("same");
    c->SaveAs("example_06.png");
    c->SaveAs("example_06.pdf");
}