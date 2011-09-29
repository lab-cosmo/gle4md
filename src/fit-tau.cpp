#include "color.hpp"
#include "clparser.hpp"
#include "ioparser.hpp"
#include "matrix-io.hpp"
#include "minsearch.hpp"
#include <fstream>
using namespace std;
using namespace toolbox;
 
void banner() 
{
    std::cerr
            << " USAGE: fit-tau (-n nthermo | -i initial ) -t target  [-h] [-nmin minstep]      \n"
            << " compute best-fit parameters for the autocorrelation time of a set of           \n"
            << " thermostats.                                                                   \n"
            << "                                                                                \n";
}


ofstream fa("A.mat");
ofstream fb("B.mat");
ofstream fc("C.mat");
ofstream fthermo("thermo");

double chi2(const std::valarray<ThermoOptions>& thermo, const std::valarray<double>& lx, const std::valarray<double>& ly, const std::valarray<double>& lw)
{
    unsigned long nthermo=thermo.size(), n=nthermo+2;
    
    FMatrix<double> A(n,n),B(n,n),C(n,n),X00(n,n),X11(n,n);
    build_A(thermo,A);
    build_B(thermo,B);
    A(1,0)=-1.;
    for (int it=0; it<nthermo; ++it)
        std::cerr<<it<<":"<<thermo[it];
    double ichi=0;
    for (int i=0; i<lx.size(); ++i)
    {
        double w=lx[i];
        A(0,1)=w*w;
        static_covariance(A,B,C);
        get_Xab(A,C,0,0,X00);
        get_Xab(A,C,1,1,X11);
        double pppp=int_2acf(X00,X00,0,0,0,0);
        double ppqq=int_2acf(X00,X11,0,0,1,1)*w*w;
        double qqpp=int_2acf(X11,X00,1,1,0,0)*w*w;
        double qqqq=int_2acf(X11,X11,1,1,1,1)*w*w*w*w;
        double pppp0=(C(0,0)*C(0,0)+C(0,0)*C(0,0));
        double ppqq0=(C(0,1)*C(0,1)+C(0,1)*C(0,1))*w*w;
        double qqpp0=(C(1,0)*C(1,0)+C(1,0)*C(1,0))*w*w;
        double qqqq0=(C(1,1)*C(1,1)+C(1,1)*C(1,1))*w*w*w*w;
        double tau=(pppp+qqqq+ppqq+qqpp)/(pppp0+qqqq0+ppqq0+qqpp0);
        //std::cerr<<lx[i]<<","<<ly[i]<<","<<1/tau<<std::endl;
        ichi+=lw[i]*((1./tau)-(ly[i]))*((1./tau)-(ly[i]));
        
    }
    /*for (int it=0; it<nthermo; ++it)
        ichi+=thermo[it].taut*thermo[it].taut;
    */
    return ichi;
}

#define MAX(a,b) ((a)<=(b)?(b):(a))
class FTError
{
    public:
        std::valarray<ThermoOptions> thermo;
        double beta;
        std::valarray<double> lx, ly, lw;
        inline unsigned long size() { return thermo.size()*3; }

        inline void set_vars(const std::valarray<double>& rv) 
        {
            if (thermo.size()*3!=rv.size()) thermo.resize((rv.size())/3);
            
            thermo[0].beta=beta;
            for (unsigned long i=0; i<thermo.size();++i)
            {
                //thermo[i].tauw=MAX(rv[4*i+3], 1e-10);
                thermo[i].taut=exp(rv[3*i]);
                thermo[i].tauf=exp(rv[3*i+1]);
                thermo[i].tauw=exp(rv[3*i+2]);
            }
        } 
        
        inline void get_vars(std::valarray<double>& rv) const 
        { 
            if (thermo.size()*3!=rv.size()) rv.resize(thermo.size()*3);
            for (unsigned long i=0; i<thermo.size();++i)
            {
                rv[3*i]=log(thermo[i].taut);
                rv[3*i+1]=log(thermo[i].tauf);
                rv[3*i+2]=log(thermo[i].tauw);
            }
        }
        
        inline void get_value(double& rv) { rv=chi2(thermo, lx, ly, lw); }
        inline void get_gradient(std::valarray<double>& rv) const 
        { ERROR("Minimization function class does not define a gradient function"); }
        inline void get_hessian(std::valarray<std::valarray<double> >& rv) const 
        { ERROR("Minimization function class does not define a hessian function"); } 
        bool bailout() { return false; }
};


int main(int argc, char **argv)
{
    CLParser clp(argc, argv);
    bool fhelp, frest;
    std::string iname="", tname="";
    unsigned long nthermo, nsteps;
    FTError chi2min;
    bool fok=
            clp.getoption(nsteps,"nmin",(unsigned long) 500u) &&
            (clp.getoption(nthermo,"n") | (frest=clp.getoption(iname,"i"))) &&
            clp.getoption(tname,"t") &&
            clp.getoption(fhelp,"h",false);
    
    if ( fhelp || ! fok) { banner(); return 0; }
    
    std::valarray<double> ll, tw, tt, tp;
    std::ifstream itarget(tname.c_str());
    IField<std::valarray<double > > iftarget(ll,"");
    itarget>> iftarget;
    tw.resize(ll.size()/3);
    tt.resize(ll.size()/3);
    tp.resize(ll.size()/3);
    for (int it=0; it<ll.size()/3; it++)
    { tw[it]=ll[3*it]; tt[it]=ll[3*it+1]; tp[it]=ll[3*it+2];}
    itarget.close();
    
    std::valarray<ThermoOptions> thermo;
    if (frest) 
    {
        std::ifstream ithermo(iname.c_str());
        IField<std::valarray<ThermoOptions > > ifthermo(thermo,"");
        ithermo >> ifthermo;
        nthermo=thermo.size();
        ithermo.close();
    }
    else
    {
        //hardcoded init parameters
        thermo.resize(nthermo);
        for (unsigned long i=0; i<nthermo;++i)
        {
            thermo[i].tauw=1./(i+1.);
            thermo[i].tauf=1.;
            thermo[i].taut=1.;
            thermo[i].beta=1.;
        }
    }
    
    std::valarray<double> pars;
    std::cerr<<"init min object "<<"\n";
    chi2min.lx.resize(tw.size()); chi2min.ly.resize(tw.size()); chi2min.lw.resize(tw.size());
    chi2min.lx=tw; chi2min.ly=tt; chi2min.lw=tp;
    chi2min.thermo.resize(nthermo); chi2min.thermo=thermo;
    chi2min.get_vars(pars); chi2min.beta=1.;

    double ierr; 
    chi2min.get_value(ierr);
    std::cerr<<"error is "<<ierr<<"\n";
    double ferr;
    
    /*
    simplex
    IterOptions<double,2> iops=IterOptions<double,2>(
            nsteps,
            fixarray<double,2>(1e-5, 1e-5),
            fixarray<double,2>(0.,0.),
            fixarray<double,2>(ichk_change, ichk_default));
    min_simplex(chi2min,make_simplex(pars,0.002,1.1),pars,ferr, iops);
    */
    
    AnnealingOptions ao;
    ao.temp_init=ierr; ao.temp_final=ierr*1e-3;
    ao.steps=nsteps; ao.mc_step=1e-1;
    sim_annealing(chi2min,pars,pars,ferr,ao);
    
    
    std::cerr<<"final error is "<<ferr<<"\n";
    std::cerr<<"final thermo is \n";
    for (unsigned long i=0; i<nthermo;++i)
        std::cerr<<i<<":"<<chi2min.thermo[i];
    
    thermo=chi2min.thermo;
    IField<std::valarray<ThermoOptions > > othermo(thermo,"");
    fthermo<<"# thermostat parameters with beta set to unit\n";
    fthermo<<othermo;
    
    unsigned long n=nthermo+2;
    FMatrix<double> A(n,n), B(n,n), C(n,n), X00, X11, T;
    build_A(chi2min.thermo,A);
    build_B(chi2min.thermo,B);
    A(1,0)=-1.;
    fa<<A;
    fb<<B;
    for (double w=tw.min(); w<tw.max()*10; w+=(10*tw.max()-tw.min())/10000.)
    {
        A(0,1)=w*w;
        static_covariance(A,B,C);
        get_Xab(A,C,0,0,X00);
        get_Xab(A,C,1,1,X11);
        double pppp=int_2acf(X00,X00,0,0,0,0);
        double ppqq=int_2acf(X00,X11,0,0,1,1)*w*w;
        double qqpp=int_2acf(X11,X00,1,1,0,0)*w*w;
        double qqqq=int_2acf(X11,X11,1,1,1,1)*w*w*w*w;
        double pppp0=(C(0,0)*C(0,0)+C(0,0)*C(0,0));
        double ppqq0=(C(0,1)*C(0,1)+C(0,1)*C(0,1))*w*w;
        double qqpp0=(C(1,0)*C(1,0)+C(1,0)*C(1,0))*w*w;
        double qqqq0=(C(1,1)*C(1,1)+C(1,1)*C(1,1))*w*w*w*w;
        
        std::cout<<w<<"  "<<1./((pppp+qqqq+ppqq+qqpp)/(pppp0+qqqq0+ppqq0+qqpp0))
                <<std::endl;
    }
    
    return 0;
}
