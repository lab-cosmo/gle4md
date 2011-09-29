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
            << " USAGE: fit-debye (-n nthermo | -i initial ) -f max-freq  [-hbar hbar]          \n"
            << "                  [-beta beta] [-t tau0] [-nmin minstep] [-h]                   \n"
            << " compute best-fit parameters for a set of thermostats enforcing bose-like       \n"
            << " energy distribution for harmonic oscillators with frequencies up to <max-freq>.\n"
            << "                                                                                \n";
}


ofstream fa("A.mat");
ofstream fb("B.mat");
ofstream fc("C.mat");
ofstream fthermo("thermo");

double bose(double beta, double hbar, double w)
{
    return hbar*w*(0.5+1./(exp(beta*hbar*w)-1));
}

double chi2(const std::valarray<ThermoOptions> thermo, double beta, double hbar, double maxfreq)
{
    unsigned long nthermo=thermo.size(), n=nthermo+2;
    
    FMatrix<double> A(n,n),B(n,n),C(n,n);
    build_A(thermo,A);
    build_B(thermo,B);
    A(1,0)=-1.;
    for (int it=0; it<nthermo; ++it)
        std::cerr<<it<<":"<<thermo[it];
    double ichi=0;
    for (double t=0; t<=50; t+=1)
    {
        double w=maxfreq/pow(1.1,t);
        A(0,1)=w*w;
        static_covariance(A,B,C);
        double e=C(0,0)+C(1,1)*w*w,et=bose(beta,hbar,w);
        ichi+=(e-et)*(e-et);
        //std::cout<<w<<"  "<<bose(beta,hbar,w)<<"  "<< C(0,0)<<"  "<<C(1,1)*w*w<<"  "<<C(0,0)+C(1,1)*w*w<<std::endl;
    }
    
    return ichi;
}

#define MAX(a,b) ((a)<=(b)?(b):(a))
class FDError
{
    public:
        std::valarray<ThermoOptions> thermo;
        double beta,hbar;
        double maxfreq;
        inline unsigned long size() { return thermo.size()*3-1; }

        inline void set_vars(const std::valarray<double>& rv) 
        {
            if (thermo.size()*3-1!=rv.size()) thermo.resize((rv.size()+1)/3);
            
            thermo[0].beta=fabs(1./rv[0]);
            thermo[0].tauf=MAX(fabs(1./rv[1]),1e-10);
            for (unsigned long i=1; i<thermo.size();++i)
            {
                //thermo[i].tauw=MAX(rv[4*i+3], 1e-10);
                thermo[i].tauf=MAX(1./rv[3*i+1], 1e-10);
                thermo[i].taut=MAX(1./rv[3*i+0], 1e-10);
                thermo[i].beta=fabs(1./rv[3*i-1]);
            }
        } 
        
        inline void get_vars(std::valarray<double>& rv) const 
        { 
            if (thermo.size()*3-1!=rv.size()) rv.resize(thermo.size()*3-1);
            rv[0]=1./thermo[0].beta;
            rv[1]=1./thermo[0].tauf;
            for (unsigned long i=1; i<thermo.size();++i)
            {
                rv[3*i+1]=1./thermo[i].tauf;
                rv[3*i+0]=1./thermo[i].taut;
                rv[3*i-1]=1./thermo[i].beta;
            }
        }
        
        inline void get_value(double& rv) { rv=chi2(thermo, beta, hbar, maxfreq); }
        inline void get_gradient(std::valarray<double>& rv) const 
        { ERROR("Minimization function class does not define a gradient function"); }
        inline void get_hessian(std::valarray<std::valarray<double> >& rv) const 
        { ERROR("Minimization function class does not define a hessian function"); } 
        bool bailout() { return false; }
};

/*
class FDError
{
    public:
        std::valarray<ThermoOptions> thermo;
        double beta,hbar;
        double maxfreq;
        inline unsigned long size() { return thermo.size()*4; }

        inline void set_vars(const std::valarray<double>& rv) 
        {
            if (thermo.size()!=rv.size()/2) thermo.resize(rv.size()/2);
            
            for (unsigned long i=0; i<thermo.size();++i)
            {
                /*thermo[i].tauw=MAX(rv[4*i+3], 1e-10);
                thermo[i].tauf=MAX(rv[4*i+2], 1e-10);
                thermo[i].taut=MAX(rv[2*i+1], 1e-10);
                thermo[i].beta=MAX(rv[2*i], 1e-10);
            }
        } 
        
        inline void get_vars(std::valarray<double>& rv) const 
        { 
            if (thermo.size()!=rv.size()/2) rv.resize(thermo.size()*2);
            for (unsigned long i=0; i<thermo.size();++i)
            {
                /*rv[4*i+3]=thermo[i].tauw;
                rv[4*i+2]=thermo[i].tauf;
                rv[2*i+1]=thermo[i].taut;
                rv[2*i]=thermo[i].beta;
            }
        }
        
        inline void get_value(double& rv) { rv=chi2(thermo, beta, hbar, maxfreq); }
        inline void get_gradient(std::valarray<double>& rv) const 
        { ERROR("Minimization function class does not define a gradient function"); }
        inline void get_hessian(std::valarray<std::valarray<double> >& rv) const 
        { ERROR("Minimization function class does not define a hessian function"); } 
        bool bailout() { return false; }
};*/

int main(int argc, char **argv)
{
    CLParser clp(argc, argv);
    bool fhelp, frest;
    std::string iname="";
    unsigned long nthermo, nsteps;
    double taut0;
    FDError chi2min;
    bool fok=
            clp.getoption(chi2min.beta,"beta",1.) &&
            clp.getoption(chi2min.hbar,"hbar",1.) &&
            clp.getoption(taut0,"t",1.) &&
            clp.getoption(nsteps,"nmin",(unsigned long) 500u) &&
            (clp.getoption(nthermo,"n") | (frest=clp.getoption(iname,"i"))) &&
            clp.getoption(chi2min.maxfreq,"f") &&
            clp.getoption(fhelp,"h",false);
    
    if ( fhelp || ! fok) { banner(); return 0; }
    
    std::valarray<ThermoOptions> thermo;
    if (frest) 
    {
        std::ifstream ithermo(iname.c_str());
        IField<std::valarray<ThermoOptions > > ifthermo(thermo,"");
        ithermo >> ifthermo;
        nthermo=thermo.size();
    }
    else
    {
        //hardcoded init parameters
        thermo.resize(nthermo);
        for (unsigned long i=0; i<nthermo;++i)
        {
            thermo[i].tauw=1./(chi2min.maxfreq*1.5*i/(nthermo-1.)+1e-10);
            thermo[i].tauf=1./(chi2min.maxfreq*1.5/(nthermo-1.));
            thermo[i].taut=1;
            thermo[i].beta=2./(bose(chi2min.beta,chi2min.hbar,1./thermo[i].tauw));
        }
        thermo[0].taut=taut0; thermo[0].tauf/=3.;
    }
    unsigned long n=nthermo+2;
    std::valarray<double> pars;
    
    std::cerr<<"init min object "<<"\n";
    chi2min.thermo.resize(nthermo); chi2min.thermo=thermo;
    chi2min.get_vars(pars);
    
    double ierr; 
    chi2min.get_value(ierr);
    std::cerr<<"error is "<<ierr<<"\n";
    double ferr;
    
    IterOptions<double,2> iops=IterOptions<double,2>(
            nsteps,
            fixarray<double,2>(1e-5, 1e-5),
            fixarray<double,2>(0.,0.),
            fixarray<double,2>(ichk_change, ichk_default));
    min_simplex(chi2min,make_simplex(pars,0.001,1.1),pars,ferr, iops);
    std::cerr<<"final error is "<<ferr<<"\n";
    std::cerr<<"final thermo is \n";
    for (unsigned long i=0; i<nthermo;++i)
        std::cerr<<i<<":"<<chi2min.thermo[i];
     
    thermo=chi2min.thermo;
    for (unsigned long i=0; i<nthermo;++i)
    {
        thermo[i].beta/=chi2min.beta;
        thermo[i].taut/=(chi2min.beta*chi2min.hbar);
        thermo[i].tauf/=(chi2min.beta*chi2min.hbar);
        thermo[i].tauw/=(chi2min.beta*chi2min.hbar);
    }
    IField<std::valarray<ThermoOptions > > othermo(thermo,"");
    fthermo<<"# thermostat parameters with beta and hbar set to unit\n";
    fthermo<<othermo;
    
    FMatrix<double> A(n,n), B(n,n), C(n,n), X00, X11, T;
    build_A(chi2min.thermo,A);
    build_B(chi2min.thermo,B);
    A(1,0)=-1.;
    fa<<A;
    fb<<B;
    for (double w=1.; w<chi2min.maxfreq; w+=1.)
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
        
        std::cout<<w<<"  "<<bose(chi2min.beta,chi2min.hbar,w)<<"  "<< C(0,0)<<"  "<<C(1,1)*w*w<<"  "<<C(0,0)+C(1,1)*w*w
                <<" "<< pppp<<"  "<<ppqq<<"  "<<qqpp<<"  "<<qqqq
                <<"  "<<(pppp+qqqq+ppqq+qqpp)/(pppp0+qqqq0+ppqq0+qqpp0)
                <<std::endl;
    }
    
    return 0;
}
