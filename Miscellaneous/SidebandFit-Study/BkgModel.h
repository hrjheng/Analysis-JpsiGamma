#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"

RooAbsPdf* GetBernstein(RooRealVar *var, const char* prefix, int order)
{
    RooArgList *coeffList = new RooArgList();
    map<const char*,RooRealVar*> params;
    map<const char*,RooFormulaVar*> prods;
    map<const char*,RooAbsPdf*> utilities;
    
    for (int i=0; i<order; i++)
    {
        const char* name = Form("%s_p%d",prefix,i);

        RooRealVar *param = new RooRealVar(name,name,0.1*(i+1),-5.,5.);
        RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name),Form("%s_sq",name),"@0*@0",RooArgList(*param));
        params.insert(pair<const char*,RooRealVar*>(name,param));
        prods.insert(pair<const char*,RooFormulaVar*>(name,form));
        coeffList->add(*prods[name]);
    }
    
    if (order==1) {
        RooBernsteinFast<1> *bern = new RooBernsteinFast<1>(prefix,prefix,*var,*coeffList);
        return bern;
    } else if (order==2) {
        RooBernsteinFast<2> *bern = new RooBernsteinFast<2>(prefix,prefix,*var,*coeffList);
        return bern;
    } else if (order==3) {
        RooBernsteinFast<3> *bern = new RooBernsteinFast<3>(prefix,prefix,*var,*coeffList);
        return bern;
    } else if (order==4) {
        RooBernsteinFast<4> *bern = new RooBernsteinFast<4>(prefix,prefix,*var,*coeffList);
        return bern;
    } else if (order==5) {
        RooBernsteinFast<5> *bern = new RooBernsteinFast<5>(prefix,prefix,*var,*coeffList);
        return bern;
    } else if (order==6) {
        RooBernsteinFast<6> *bern = new RooBernsteinFast<6>(prefix,prefix,*var,*coeffList);
        return bern;
    } else {
        return NULL;
    }
}

RooAbsPdf* GetExponential(RooRealVar *var, const char* prefix, int order)
{
    map<const char*,RooRealVar*> params;
    map<const char*,RooFormulaVar*> prods;
    map<const char*,RooAbsPdf*> utilities;
    
    if (order%2==0)
    {
        cout << "ERROR -- addExponential -- only odd number of params allowed" << endl;
        return NULL;
    }
    else
    {
        int nfracs=(order-1)/2;
        int nexps=order-nfracs;
        assert(nfracs==nexps-1);
        RooArgList *fracs = new RooArgList();
        RooArgList *exps = new RooArgList();
        for (int i=1; i<=nfracs; i++){
            const char* name =  Form("%s_f%d",prefix,i);
            params.insert(pair<const char*,RooRealVar*>(name, new RooRealVar(name,name,0.9-float(i-1)*1./nfracs,0.0001,0.9999)));
            fracs->add(*params[name]);
        }
        for (int i=1; i<=nexps; i++){
            const char* name =  Form("%s_p%d",prefix,i);
            const char* ename =  Form("%s_e%d",prefix,i);
            params.insert(pair<const char*,RooRealVar*>(name, new RooRealVar(name,name,TMath::Max(-1.,-0.04*(i+1)),-1.,0.)));
            utilities.insert(pair<const char*,RooAbsPdf*>(ename, new RooExponential(ename,ename,*var,*params[name])));
            exps->add(*utilities[ename]);
        }

        RooAbsPdf *exp = new RooAddPdf(prefix,prefix,*exps,*fracs,true);

        return exp;
    }
}

RooAbsPdf* GetPowerLaw(RooRealVar *var, const char* prefix, int order)
{
    map<const char*,RooRealVar*> params;
    map<const char*,RooFormulaVar*> prods;
    map<const char*,RooAbsPdf*> utilities;
    
    if (order%2==0){
        cerr << "ERROR -- addPowerLaw -- only odd number of params allowed" << endl;
        return NULL;
    }
    else {
        int nfracs=(order-1)/2;
        int npows=order-nfracs;
        assert(nfracs==npows-1);
        RooArgList *fracs = new RooArgList();
        RooArgList *pows = new RooArgList();
        for (int i=1; i<=nfracs; i++){
            const char* name =  Form("%s_f%d",prefix,i);
            params.insert(pair<const char*,RooRealVar*>(name, new RooRealVar(name,name,0.9-float(i-1)*1./nfracs,0.,1.)));

            fracs->add(*params[name]);
        }
        for (int i=1; i<=npows; i++){
            const char* name =  Form("%s_p%d",prefix,i);
            const char* ename =  Form("%s_e%d",prefix,i);
            params.insert(pair<const char*,RooRealVar*>(name, new RooRealVar(name,name,TMath::Max(-9.,-1.*(i+1)),-9.,1.)));

            utilities.insert(pair<const char*,RooAbsPdf*>(ename, new RooPower(ename,ename,*var,*params[name])));
            pows->add(*utilities[ename]);
        }

        RooAbsPdf *pow = new RooAddPdf(prefix,prefix,*pows,*fracs,true);

        return pow;
    }
}

RooAbsPdf* GetLaurent(RooRealVar *var, const char* prefix, int order)
{
    map<const char*,RooRealVar*> params;
    map<const char*,RooFormulaVar*> prods;
    map<const char*,RooAbsPdf*> utilities;
    
    int nlower=int(ceil(order/2.));
    int nhigher=order-nlower;
    // first do 0th order
    RooArgList *pows = new RooArgList();
    RooArgList *plist = new RooArgList();
    const char* pname =  Form("%s_pow0",prefix);
    utilities.insert(pair<const char*,RooAbsPdf*>(pname, new RooPower(pname,pname,*var,RooConst(-4.))));
    pows->add(*utilities[pname]);
    
    // even terms
    for (int i=1; i<=nlower; i++){
        const char* name = Form("%s_l%d",prefix,i);
        params.insert(pair<const char*,RooRealVar*>(name, new RooRealVar(name,name,0.25/order,0,1)));
        plist->add(*params[name]);
        const char* pname =  Form("%s_powl%d",prefix,i);
        utilities.insert(pair<const char*,RooAbsPdf*>(pname, new RooPower(pname,pname,*var,RooConst(-4.-i))));
        pows->add(*utilities[pname]);
    }
    // odd terms
    for (int i=1; i<=nhigher; i++){
        const char* name = Form("%s_h%d",prefix,i);
        params.insert(pair<const char*,RooRealVar*>(name, new RooRealVar(name,name,0.25/order,0,1)));
        plist->add(*params[name]);
        const char* pname =  Form("%s_powh%d",prefix,i);
        utilities.insert(pair<const char*,RooAbsPdf*>(pname, new RooPower(pname,pname,*var,RooConst(-4.+i))));
        pows->add(*utilities[pname]);
    }
    RooAddPdf *pdf = new RooAddPdf(prefix,prefix,*pows,*plist,true);
    return pdf;
    //bkgPdfs.insert(pair<const char*,RooAbsPdf*>(pdf->GetName(),pdf));
}
