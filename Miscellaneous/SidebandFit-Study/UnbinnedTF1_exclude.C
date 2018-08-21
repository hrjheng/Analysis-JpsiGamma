#include "tdrstyle.C"
//TH1D *hist = 0;
TTree *t = 0;
Bool_t reject;
Bool_t ZAna;
float Mmmg_lower;
float Mmmg_upper;
float Signalrange_lower;
float Signalrange_upper;
float Norm_SB;

double model(double x, double *par) // Bernstein 2nd-order polynomial
{
    double a0 = par[0];
    double a1 = 2.*(par[1]-par[0]);
    double a2 = par[0]-2.*par[1]+par[2];
    double N = 1./(a0+0.5*a1+(1./3.)*a2);

    //cout << x << " " << a0 << " " << a1 << " " << a2 << " " << N << " " << (a0 + a1*x + a2*x*x)*N << endl;

    //cout << (Signalrange_lower-Mmmg_lower)/(Mmmg_upper-Mmmg_lower)  << " " << (Signalrange_upper-Mmmg_lower)/(Mmmg_upper-Mmmg_lower) << endl;
    if (reject && x > (Signalrange_lower-Mmmg_lower)/(Mmmg_upper-Mmmg_lower) && x < (Signalrange_upper-Mmmg_lower)/(Mmmg_upper-Mmmg_lower))
    {
        TF1::RejectPoint();
        return 0;
    }

    return (a0 + a1*x + a2*x*x)/(a0+0.5*a1+(1./3.)*a2);
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
    f = 0.;
    float m,trans_m;
    t->SetBranchAddress("m_mmg", &m);
    for (Long64_t ev = 0; ev < t->GetEntriesFast(); ev++) 
    {
        t->GetEntry(ev);

        trans_m = (m-Mmmg_lower)/(Mmmg_upper-Mmmg_lower);

        double L = model(trans_m,par);
        if (L>0.) f -= 2.*log(L);
        else { f = HUGE; return; }
    }
}

void UnbinnedTF1_exclude()
{
    setTDRStyle();

    //Create a source function
    TFile *f = new TFile("minitree_2016_HJpsiG_AllRuns.root","READ");
    t = (TTree*) f->Get("outTree");
    float m;
    t->SetBranchAddress("m_mmg", &m);

    ZAna=kFALSE;
    if(ZAna)
    {
        Mmmg_lower = 70.;
        Mmmg_upper = 120.;
        Signalrange_lower = 85.;
        Signalrange_upper = 95.;
    }
    else
    {
        Mmmg_lower = 100.;
        Mmmg_upper = 150.;
        Signalrange_lower = 120.;
        Signalrange_upper = 130.;
    }

    TH1F *h_data = new TH1F("h_data"," ",25,0,1);
    for (Long64_t ev = 0; ev < t->GetEntriesFast(); ev++){
        t->GetEntry(ev);
        h_data->Fill((m-Mmmg_lower)/(Mmmg_upper-Mmmg_lower));
        
        if(m<Signalrange_lower || m>Signalrange_upper)
            Norm_SB++;
    }
    Norm_SB = Norm_SB/(t->GetEntriesFast());
    
    // Sideband-only fit
    TMinuit *gMinuit = new TMinuit(3);
    reject = kTRUE;
    gMinuit->SetFCN(fcn);
    gMinuit->DefineParameter(0, "c0", 0.3, 0.1, 0, 1);
    gMinuit->DefineParameter(1, "c1", 0.2, 0.1, 0, 1);
    gMinuit->DefineParameter(2, "c2", 0.1, 0.1, 0, 1);
    gMinuit->Command("MIGRAD");
    //gMinuit->Command("MIGRAD");
    gMinuit->Command("MINOS");

    double par[3],err[3];
    for(int i=0;i<3;i++) gMinuit->GetParameter(i,par[i],err[i]);
    
    TH1F* curve = new TH1F("curve","curve",h_data->GetNbinsX()*100,
                           h_data->GetXaxis()->GetXmin(),h_data->GetXaxis()->GetXmax());
    
    for(int i=1;i<=curve->GetNbinsX();i++) 
    {
        reject = kFALSE;
        double x = curve->GetBinCenter(i);
        double f = model(x,par);
        double BinWidth = h_data->GetBinWidth(1);
        curve->SetBinContent(i,f*t->GetEntries()*BinWidth);
        //cout << f << " " << f*t->GetEntries()*BinWidth << endl;
    }

    //Include-the-signal-region fit
    TMinuit *gMinuit2 = new TMinuit(3);
    reject = kFALSE;
    gMinuit2->SetFCN(fcn);
    gMinuit2->DefineParameter(0, "c0_2", 0.3, 0.1, 0, 1);
    gMinuit2->DefineParameter(1, "c1_2", 0.2, 0.1, 0, 1);
    gMinuit2->DefineParameter(2, "c2_2", 0.1, 0.1, 0, 1);
    gMinuit2->Command("MIGRAD");
    //gMinuit2->Command("MIGRAD");
    //gMinuit2->Command("MINOS");

    for(int i=0;i<3;i++) gMinuit2->GetParameter(i,par[i],err[i]);
    
    TH1F* curve2 = new TH1F("curve","curve",h_data->GetNbinsX()*100,
                           h_data->GetXaxis()->GetXmin(),h_data->GetXaxis()->GetXmax());
    
    for(int i=1;i<=curve2->GetNbinsX();i++) 
    {
        reject = kFALSE;
        double x = curve2->GetBinCenter(i);
        double f = model(x,par);
        double BinWidth = h_data->GetBinWidth(1);
        curve2->SetBinContent(i,f*t->GetEntries()*BinWidth);
    }

    //Start making plot
    TCanvas *c = new TCanvas("c","",800,700);
    c->cd();
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gStyle->SetErrorX(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    h_data->SetLineColor(1);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(2);
    h_data->SetLineWidth(2);
    //hist->Draw("ep");
    curve->SetLineWidth(3);
    curve->SetLineColor(TColor::GetColor("#F96C48"));
    curve->SetFillColor(TColor::GetColor("#F96C48"));
    curve->GetXaxis()->SetTitle("Transformed M_{#mu#mu#gamma} (GeV)");
    curve->GetYaxis()->SetTitle(Form("Events / (%.2f GeV)",(float) 1./25.));
    curve->GetYaxis()->SetRangeUser(0.,h_data->GetMaximum()*1.7);
    curve->Draw();
    curve2->SetLineWidth(3);
    curve2->SetLineColor(TColor::GetColor("#466C95"));
    curve2->Draw("Csame");
    h_data->Draw("epsame");
    TLegend* leg = new TLegend(0.62,0.65,0.9,0.85);
    leg->AddEntry(h_data, "Data (H#rightarrowJ/#psi #gamma channel)", "lep");
    leg->AddEntry(curve, "Sideband-only fit", "fl");
    leg->AddEntry(curve2, "Signal+Background fit", "l");
    leg->Draw("same");
    c->SaveAs("UnbinnedTF1_exclude.pdf");
    c->SaveAs("UnbinnedTF1_exclude.png");
}
