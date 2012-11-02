#define TB_FFTAC 1
 
#include "color.hpp"
#include "clparser.hpp"
#include "conv-units.hpp"
#include "tools-autocorr.hpp"
#include "tools-histogram.hpp"
#include "rndgen.hpp"
#include "ioparser.hpp"
#include "matrix-io.hpp"
#include "minsearch.hpp"
#include <fstream>
/*******************************************************************
  INTERNAL UNITS ARE ATOMIC UNITS
  length=1 bohr=5.29177 10^-11 m
*******************************************************************/


using namespace std;
using namespace toolbox;
enum Potential { PHarmonic, PLennardJones, PQuasiHarmo, PFlexiWell, PDoubleWell, PSquareWell, PMorse, PFree};
enum Integrator { PBeads, PNormalModes };
class odops 
{
public:
    FMatrix<double> A, C;
    std::vector<FMatrix<double> > fA, fC;
    UnitConv units;
    std::string prefix;
    bool f_print_traj, f_histogram; 
    unsigned long h_nwin;
    long seed;
    double dt;
    unsigned long nstep, stride, nbeads;
    unsigned long maxlag;
    unsigned long nseries;
    unsigned long drop;
    Integrator prop;
    Potential pot;
    double temp, mass;
    std::valarray<double> pars;
    std::vector<double> ux; unsigned long nlin;
    unsigned long npderiv;
    
};

std::ostream& operator<< (std::ostream& ostr, const Integrator& pp)
{
    switch (pp) 
    {
    case PBeads:
        ostr<<"beads";
        break;
    case PNormalModes:
        ostr<<"nmodes";
        break;
    }
    return ostr;
}

std::istream& operator>> (std::istream& istr,  Integrator& pp)
{
    std::string spp;
    istr>>spp;
    if (spp=="beads") pp=PBeads;
    else if (spp=="nmodes") pp=PNormalModes;
    else ERROR(spp<<" is not a supported integrator [beads|nmodes].");
    std::cerr<<"INTEGRATOR READ "<<spp<<"\n";
    return istr;
}

std::ostream& operator<< (std::ostream& ostr, const odops& oo)
{
    std::vector<string> ppars;
    switch (oo.pot)
    {
    case PHarmonic:
        if (oo.pars.size()<1)
            ERROR("Wrong number of parameters for harmonic potential.");
        ppars.push_back("harmonic");
        ppars.push_back(float2str(oo.pars[0]));
        break;
    case PLennardJones:
        if (oo.pars.size()<2)
            ERROR("Wrong number of parameters for Lennard-Jones potential.");
        ppars.push_back("lennard-jones");
        ppars.push_back(float2str(oo.pars[0]));
        ppars.push_back(float2str(oo.pars[1]));
        break;
    case PMorse:
        if (oo.pars.size()<2)
            ERROR("Wrong number of parameters for Morse potential.");
        ppars.push_back("morse");
        ppars.push_back(float2str(oo.pars[0]));
        ppars.push_back(float2str(oo.pars[1]));
        break;
    case PQuasiHarmo:
        if (oo.pars.size()<2)
            ERROR("Wrong number of parameters for quasi-harmonic potential.");
        ppars.push_back("quasi-harmo");
        ppars.push_back(float2str(oo.pars[0]));
        ppars.push_back(float2str(oo.pars[1]));
        break;
    case PDoubleWell:
        if (oo.pars.size()<2)
            ERROR("Wrong number of parameters for double well potential.");
        ppars.push_back("double-well");
        ppars.push_back(float2str(oo.pars[0]));
        ppars.push_back(float2str(oo.pars[1]));
        break;
    case PFlexiWell:
            if (oo.pars.size()<3)
                ERROR("Wrong number of parameters for 'flexible' well potential.");
            ppars.push_back("flexi-well");
            ppars.push_back(float2str(oo.pars[0]));
            ppars.push_back(float2str(oo.pars[1]));
            ppars.push_back(float2str(oo.pars[2]));
            break;
    case PSquareWell:
        if (oo.pars.size()<1)
            ERROR("Wrong number of parameters for harmonic potential.");
        ppars.push_back("square-well");
        ppars.push_back(float2str(oo.pars[0]));
        break;
    case PFree:
        if (oo.pars.size()>0)
            ERROR("Wrong number of parameters for free particle.");
        ppars.push_back("free");
        break;
    }

    IOMap iom;
    iom.insert(oo.ux,"lin-x");
    iom.insert(oo.nlin,"lin-n");
    iom.insert(oo.A,"a_matrix"); 
    iom.insert(oo.C,"c_matrix");
    iom.insert(oo.units, "units");
    iom.insert(oo.dt,"timestep");
    iom.insert(oo.drop,"discard");
    iom.insert(ppars,"potential");
    iom.insert(oo.temp,"temperature");
    iom.insert(oo.mass,"mass");
    iom.insert(oo.nstep,"steps");
    iom.insert(oo.stride,"stride");
    iom.insert(oo.nbeads,"beads");
    iom.insert(oo.maxlag,"max_lag");
    iom.insert(oo.nseries,"max_data");
    iom.insert(oo.prefix,"prefix");
    iom.insert(oo.seed,"seed");
    iom.insert(oo.prop,"integrator");
    iom.insert(oo.f_print_traj,"print_traj");
    iom.insert(oo.f_histogram,"histogram");
    iom.insert(oo.h_nwin,"histo-bins");
    iom.insert(oo.npderiv,"p-derivs");    
    ostr<<iom;
    
    return ostr;
}

std::istream& operator>> (std::istream& istr, odops& oo)
{
    IOMap iom;
    std::vector<string> ppars;
    FMatrix<double> Z(0,0);
    std::string famat, fcmat;
    iom.insert(famat,"a_matrix_f",std::string("")); 
    iom.insert(fcmat,"c_matrix_f",std::string(""));
    iom.insert(oo.ux,"lin-x",std::vector<double> (0));
    iom.insert(oo.nlin,"lin-n",(unsigned long) 0);
    iom.insert(oo.A,"a_matrix");
    iom.insert(oo.C,"c_matrix");
    iom.insert(oo.fA,"f_a_matrix",std::vector<FMatrix<double> >(0));
    iom.insert(oo.fC,"f_c_matrix",std::vector<FMatrix<double> >(0));
    iom.insert(oo.dt,"timestep");
    iom.insert(oo.nstep,"steps");
    iom.insert(oo.stride,"stride",(unsigned long) 1);
    iom.insert(ppars,"potential");
    iom.insert(oo.temp,"temperature",0.);
    iom.insert(oo.mass,"mass",1.);
    iom.insert(oo.maxlag,"max_lag");
    iom.insert(oo.nseries,"max_data");
    iom.insert(oo.nbeads,"beads",(unsigned long) 1);
    iom.insert(oo.drop,"discard",(unsigned long) 0);
    iom.insert(oo.prefix,"prefix",std::string("oned"));
    iom.insert(oo.seed,"seed",(long)1234321);
    iom.insert(oo.prop,"integrator",PBeads);
    iom.insert(oo.f_print_traj,"print_traj",true);
    iom.insert(oo.f_histogram,"histogram",true);
    iom.insert(oo.h_nwin,"histo-bins",(unsigned long) 0);
    iom.insert(oo.units, "units",UnitConv());
    iom.insert(oo.npderiv,"p-derivs",(unsigned long) 0);    
    istr>>iom;
   
    //sets INTERNAL units
    oo.units.ITime=UnitConv::AtomicTime; oo.units.IEnergy=UnitConv::Hartree;
    oo.units.IFrequency=UnitConv::AtomicFrequency; oo.units.IMass=UnitConv::AtomicMass;
    
    if (oo.h_nwin==0) oo.h_nwin=sqrt((oo.nstep-oo.drop)/100.);
 
    if (ppars[0]=="harmonic")
    {
        oo.pot=PHarmonic;
        if (ppars.size()<2)
            ERROR("Wrong number of parameters for harmonic potential.");
        oo.pars.resize(1); oo.pars[0]=str2float(ppars[1]);
    }
    else if (ppars[0]=="lennard-jones")
    {
        oo.pot=PLennardJones;
        if (ppars.size()<3)
            ERROR("Wrong number of parameters for Lennard-Jones potential.");
        oo.pars.resize(2); oo.pars[0]=str2float(ppars[1]); oo.pars[1]=str2float(ppars[2]);
    }
    else if (ppars[0]=="morse")
    {
        oo.pot=PMorse;
        if (ppars.size()<3)
            ERROR("Wrong number of parameters for Morse potential.");
        oo.pars.resize(2); oo.pars[0]=str2float(ppars[1]); oo.pars[1]=str2float(ppars[2]);
    }
    else if (ppars[0]=="quasi-harmo")
    {
        oo.pot=PQuasiHarmo;
        if (ppars.size()<3)
            ERROR("Wrong number of parameters for quasi-harmonic potential.");
        oo.pars.resize(2); oo.pars[0]=str2float(ppars[1]); oo.pars[1]=str2float(ppars[2]);
    }
    else if (ppars[0]=="double-well")
    {
        oo.pot=PDoubleWell;
        if (ppars.size()<3)
            ERROR("Wrong number of parameters for double well potential.");
        oo.pars.resize(2); oo.pars[0]=str2float(ppars[1]); oo.pars[1]=str2float(ppars[2]);
    }
    else if (ppars[0]=="flexi-well")
    {
        oo.pot=PFlexiWell;
        if (ppars.size()<4)
            ERROR("Wrong number of parameters for 'flexible' well potential.");
        oo.pars.resize(3); oo.pars[0]=str2float(ppars[1]); oo.pars[1]=str2float(ppars[2]); oo.pars[2]=str2float(ppars[3]);
    }
    else if (ppars[0]=="square-well")
    {
        oo.pot=PSquareWell;
        if (ppars.size()<2)
            ERROR("Wrong number of parameters for square well potential.");
        oo.pars.resize(1); oo.pars[0]=str2float(ppars[1]);
    }
    else if (ppars[0]=="free")
    {
        oo.pot=PFree;
        if (ppars.size()>0)
            ERROR("Wrong number of parameters for free particle.");
        oo.pars.resize(0);
    }
    
    if (oo.A.rows()==0)
    {
        if (famat!="") { std::ifstream afile(famat.c_str()); afile>>oo.A; }
    }
    
    if (oo.C.rows()==0) if (fcmat!="") { std::ifstream cfile(fcmat.c_str()); cfile>>oo.C; }

    if (oo.C.rows()==0 && oo.fC.size() != oo.nbeads)
    {
        if (oo.temp<=0.) ERROR("C matrix or temperature should be given"); 
        //classical matrix for the target temperature
        oo.C.resize(oo.A.rows(),oo.A.cols());
        for (unsigned int i=0; i<oo.A.rows(); ++i) oo.C(i,i)=oo.temp*oo.nbeads;
    }
    
    return istr;
}

double (*pma) (double), (*pe)(double);
std::valarray<double> pl;

/* 
  HARMONIC POTENTIAL:
  E=1/2w^2 x^2
  pars[0]=w (frequency)
  pl[0]=w^2
*/
double aHarmonic(double mq)
{
    return -pl[0]*mq;
}

double eHarmonic(double mq)
{
    return pl[0]*mq*mq*0.5;
}

/* 
  QUASI-HARMONIC POTENTIAL:
  E=w^2 x *(1-exp(-k x)) / (2k)
  pars[0]=w (frequency)
  pars[1]=k (anharmonicity factor)
  pl[0]=w^2/2k
  pl[1]=k
*/
double aQuasiHarmo(double mq)
{
    double ek=exp(-pl[1]*mq);
    return -pl[0]*((pl[1]*mq-1.)*ek+1.);
}

double eQuasiHarmo(double mq)
{
    return pl[0]*mq*(1.-exp(-pl[1]*mq));
}


/* 
  DOUBLE WELL POTENTIAL:
  E=w[(x/d)^4-2(x/d)^2]
  pars[0]=d  (mass-scaled x coord. of each minimum)
  pars[1]=w  (well depth)
  pl[0]=d
  pl[1]=4w/d
*/
double aDW(double mq)
{
    double sq=mq/pl[0];
    return -pl[1]*(sq*sq-1)*sq;
}
double eDW(double mq)
{
    double sq=mq/pl[0]; sq*=sq;
    return pl[1]*pl[0]/4.*((sq-2.)*sq+1);
}

/* 
  FLEXIBLE WELL POTENTIAL:
  E=w[(x/d)^4-2(x/d)^2]^(1+k exp(-2(x/d)^2))
  pars[0]=d  (mass-scaled x coord. of each minimum)
  pars[1]=q  (well depth)
  pars[2]=k  (well "flatness")
  pl[0]=d
  pl[1]=4w/d
  pl[2]=k
*/

double aFW(double mq)
{
    double u=mq/pl[0], u2=u*u, u2m1=u2-1., u2m12=u2m1*u2m1, eu=exp(-2.*u2), k=pl[2];
    if (u2m12<1e-40) return 0.;
  //  std::cerr<<"force: "<< -pl[1]*u*u2m1*pow(u2m12,eu*(k-1.))*(1./eu-1.+k-(k-1.)*u2m1*log(u2m12));
    return -pl[1]*u*u2m1*pow(u2m12,eu*k)*(1.+k*eu*(1.-u2m1*log(u2m12)));
}
double eFW(double mq)
{
    double sq=mq/pl[0]; sq*=sq;
  //  std::cerr<<"energy: "<<pl[1]*pl[0]/4.*pow(((sq-2.)*sq+1),1.+(pl[2]-1.) *exp(-2. *sq));
    return pl[1]*pl[0]/4.*pow(((sq-2.)*sq+1),1.+pl[2] *exp(-2. *sq));
}


double aSW(double mq)
{    return 0.; }
double eSW(double mq)
{    return 1.; }


double aLJ(double mq)
{
    return 0;
}

double aMorse(double mq)
{
    double z=exp(-pl[1]*mq);
    return -2.*pl[0]*pl[1]*z*(1-z);
}

double eMorse(double mq)
{
    double z=(1.-exp(-pl[1]*mq));
    return pl[0]*z*z;
}

double aFree(double mq)
{
    return 0.;
}

double eFree(double mq)
{
    return 0.;
}

double piwn, lastpiv; 
void addpia(const std::valarray<double>& mq, std::valarray<double>& pa)
{
    double f;
    lastpiv=0.; //this global variable makes my eyes bleed! 
    f=(mq[mq.size()-1]-mq[0])*piwn; pa[0]+=f; pa[mq.size()-1]-=f; lastpiv+=f*0.5*(mq[mq.size()-1]-mq[0]);
    for (int b=1; b<mq.size(); ++b) { f=(mq[b-1]-mq[b])*piwn; pa[b]+=f; pa[b-1]-=f; lastpiv+=f*0.5*(mq[b-1]-mq[b]); }
} 

void gtrans_fw(double gamma, const std::valarray<double>& mq, std::valarray<double>& mq1)
{
    if(mq1.size()!=mq.size())    mq1.resize(mq.size());
    double n=mq.size(), g=(n-sqrt(n+n*(n-1)*gamma))/(n*(n-1)), dq=mq[0]-mq[mq.size()-1];
    for (int i=0; i<mq.size(); ++i)
        mq1[i]=mq[i]+g*(i-(n-1)/2)*dq;
}

void gtrans_bw(double gamma, const std::valarray<double>& mq, std::valarray<double>& mq1)
{
    if(mq1.size()!=mq.size())    mq1.resize(mq.size());
    double n=mq.size(), g=(n-sqrt(n+n*(n-1)*gamma))/(n*(n-1)), h=g/((n-1)*g-1), dq=mq[0]-mq[mq.size()-1];
    for (int i=0; i<mq.size(); ++i)
        mq1[i]=mq[i]+h*(i-(n-1)/2)*dq;
}

FMatrix<double> nmC; std::valarray<double> nmt;
FMatrix<double> nmev;
void mknmc(unsigned long nb, double dt)
{
    nmC.resize(nb,nb); nmt.resize(nb); unsigned long i,j,k;
    nmC.col(0)=1.0;
    for (k=1; k<nb/2; k++) for (j=0; j<nb; j++) nmC(j,k)=sqrt(2.)*cos(2.0/nb*toolbox::constant::pi*j*k);
    if (nb%2==0) for (j=0; j<nb; j++) nmC(j,nb/2)=(j%2==0?1:-1);
    for (k=nb/2+1; k<nb; k++) for (j=0; j<nb; j++) nmC(j,k)=sqrt(2.)*sin(2.0/nb*toolbox::constant::pi*j*k);
    nmC*=sqrt(1.0/nb);

    nmev.resize(2*nb,2); double wk;
    for (k=0; k<nb; k++) 
    {
        wk=2*sqrt(piwn)*sin(k*constant::pi/nb);
        std::cerr<<"freq "<<k<<" "<<wk<<"\n";
        nmev(2*k,0)=nmev(2*k+1,1)=cos(wk *dt);
        nmev(2*k,1)=nmev(2*k+1,0)=sin(wk *dt);
        nmev(2*k,1)*=1.0/wk;
        nmev(2*k+1,0)*=-wk;
    }
    nmev(0,1)=dt; //w_0 -> 0
    std::cerr<<"NM EVOLUTION MATRIX "<<nmev<<"\n";
}

void nm2bd(std::valarray<double>& mq, std::valarray<double>& mp)
{
    unsigned long i,j,nb=mq.size();
    nmt=0.0; for (i=0;  i<nb; i++)  for (j=0;  j<nb; j++) nmt[i]+=nmC(i,j)*mq[j]; mq=nmt; 
    nmt=0.0; for (i=0;  i<nb; i++)  for (j=0;  j<nb; j++) nmt[i]+=nmC(i,j)*mp[j]; mp=nmt; 
}

void bd2nm(std::valarray<double>& mq, std::valarray<double>& mp)
{
    unsigned long i,j,nb=mq.size();
    nmt=0.0; for (i=0;  i<nb; i++)  for (j=0;  j<nb; j++) nmt[i]+=nmC(j,i)*mq[j]; mq=nmt; 
    nmt=0.0; for (i=0;  i<nb; i++)  for (j=0;  j<nb; j++) nmt[i]+=nmC(j,i)*mp[j]; mp=nmt; 
}

void verlet(std::valarray<double>& mq, std::valarray<double>& mp, double& dt, Integrator prop=PBeads)
{
    static std::valarray<double> olda(0);

    switch(prop)
    {
    case PBeads:
        if (olda.size()==0) { olda.resize(mq.size()); for (int b=0; b<olda.size(); ++b) olda[b]=pma(mq[b]); addpia(mq, olda); }
        mp+=0.5*dt*olda;
        mq+=dt*mp;
        for (int b=0; b<olda.size(); ++b) olda[b]=pma(mq[b]); addpia(mq, olda); 
        mp+=0.5*dt*olda;

        break;
    case PNormalModes:
        nm2bd(mq,mp);
        if (olda.size()==0) { olda.resize(mq.size()); for (int b=0; b<olda.size(); ++b) olda[b]=pma(mq[b]); }
        mp+=0.5*dt*olda;
        bd2nm(mq,mp);

        double nq,np;
        for (unsigned long k=0; k<mq.size(); ++k)
        {
            nq=nmev(2*k,0)*mq[k]+nmev(2*k,1)*mp[k]; 
            np=nmev(2*k+1,0)*mq[k]+nmev(2*k+1,1)*mp[k]; 
            mq[k]=nq; mp[k]=np; 
        }
        
        //mq+=dt*mp;

        nm2bd(mq,mp);
        for (int b=0; b<olda.size(); ++b) olda[b]=pma(mq[b]);
        mp+=0.5*dt*olda;
        bd2nm(mq,mp);

        break;
    }
}

double binomial(int n, int k)
{
   if (k==0) return 1.0;
   if (n==0) return 0.0;
   return binomial(n-1,k-1) + binomial(n-1,k);
}

int main(int argc, char **argv)
{
    RndGaussian<double, MTRndUniform> rgauss(RGPars<double>(0.,1.));
    
    FMatrix<double>  S, T, rC;
    std::vector<FMatrix<double> > fS, fT;
    
    odops ops;
    std::cin>>ops;
    std::cerr<<"SEED "<<ops.seed<<"\n";
    rgauss.RUGenerator().seed(ops.seed);
    //we convert the options in internal units
    ops.C*=ops.units.u2i(UnitConv::UEnergy,1.); 
    ops.A*=1./ops.units.u2i(UnitConv::UTime,1.); 
    std::cerr<<"QH: "<<ops.temp<<","<<ops.mass<<std::endl;
    ops.temp=ops.units.u2i(UnitConv::UEnergy,ops.temp); 
    piwn=pow(ops.nbeads*ops.temp,2.);
    ops.mass=ops.units.u2i(UnitConv::UMass,ops.mass); 
    std::cerr<<"QH: "<<ops.temp<<","<<ops.mass<<std::endl;
    ops.dt=ops.units.u2i(UnitConv::UTime,ops.dt); 
    std::cerr<<"QH: "<<ops.dt<<std::endl;
    unsigned long n=ops.A.rows();
    /*if (ops.temp>0.) 
    {
        //temperature control!
        std::cerr<<"RESCALING MATRICES FOR TEMPERATURE CONTROL!"<<std::endl;
        ops.A*=ops.temp*ops.nbeads/ops.C(0,0);
        ops.C*=ops.temp*ops.nbeads/ops.C(0,0);
    }*/
    StabCholesky(ops.C,rC);
    
    std::valarray<double> q(0.0,ops.nbeads), p(0.0,ops.nbeads);
    pl.resize(ops.pars.size());
    std::cerr<<"integrator: "<<ops.prop<<"\n";
    switch(ops.pot)
    {
        case PHarmonic:
            //pars[0] is the oscillator FREQUENCY (nu)
            ops.pars[0]=ops.units.u2i(UnitConv::UFrequency,ops.pars[0]); 
            pl[0]=ops.pars[0]*ops.pars[0];
            pma=&aHarmonic; pe=&eHarmonic;
            for (unsigned long b=0; b<ops.nbeads; ++b) q[b]=rgauss()/ops.pars[0]*rC(0,0);
            break;
        case PQuasiHarmo:
            //pars[0] is the oscillator FREQUENCY (nu)
            ops.pars[0]=ops.units.u2i(UnitConv::UFrequency,ops.pars[0]); 
            //pars[1] is an inverse length, and we should also consider the mass-scaling
            ops.pars[1]*=1./(ops.units.u2i(UnitConv::ULength,1.)*sqrt(ops.mass)); 
            pl[0]=ops.pars[0]*ops.pars[0]/(2.*ops.pars[1]);
            pl[1]=ops.pars[1];
            pma=&aQuasiHarmo; pe=&eQuasiHarmo;
            q=fabs(rgauss()/ops.pars[0]*rC(0,0));
            break;
        case PDoubleWell:
            //pars[0] is the well spacing in units of length (and we should include mass-scaling)
            //pars[1] is the barrier depth in units of length (and we should include mass-scaling)
            pma=&aDW; pe=&eDW;
            
            ops.pars[0]=ops.units.u2i(UnitConv::ULength,ops.pars[0])*sqrt(ops.mass);
            ops.pars[1]=ops.units.u2i(UnitConv::UEnergy,ops.pars[1]);
            pl[0]=ops.pars[0];
            pl[1]=4.*ops.pars[1]/ops.pars[0];
            for (unsigned long b=0; b<ops.nbeads; ++b) q[b]=ops.pars[0]*(b%2==0?1:-1); //it will decorrelate afterwards, we start in the bottom of a well
            break;
        case PFlexiWell:
            pma=&aFW; pe=&eFW;
            pl[0]=ops.pars[0];
            pl[1]=4*ops.pars[1]/ops.pars[0];
            pl[2]=ops.pars[2];
            q=ops.pars[0]; //it will decorrelate afterwards, we start in the bottom of a well
            break;
            case PMorse:   //input parameters are equivalent harmonic frequency and well depth. here we convert to D and a
                pl[0]=ops.pars[0];
                pl[1]=ops.pars[1]/sqrt(2.*ops.pars[0]);
                pma=&aMorse; pe=&eMorse;
                q=rgauss()/ops.pars[1]*rC(0,0);
                break;
        case PSquareWell:
            pma=&aSW; pe=&eSW;
            ops.pars[0]=ops.units.u2i(UnitConv::ULength,ops.pars[0])*sqrt(ops.mass);
            std::cerr<<"UNIT LENGTH "<<ops.pars[0]<<"\n";
            pl[0]=ops.pars[0];
            q=0.;
            break;
        case PFree:
            pma=&aFree; pe=&eFree;
            q=0.;
            break;
        default:
            ERROR("Potential not yet implemented!\n");
    }

    std::ofstream acout((ops.prefix+std::string(".acf")).c_str());
    std::ofstream tout, lout, laout, npout;
    if (ops.f_print_traj)  tout.open((ops.prefix+std::string(".traj")).c_str());
    
    std::valarray<double> uxavg(ops.ux.size()); uxavg=0.; unsigned long nux=0;
    if (ops.ux.size()>0) laout.open((ops.prefix+std::string(".lin")).c_str());
    if (ops.f_print_traj && ops.ux.size()>0) { 
        lout.open((ops.prefix+std::string(".lin-traj")).c_str());
        lout <<" # displaced path estimator computed at x values\n # "; 
        for (int u=0; u<ops.ux.size(); ++u)
        {
            for (int u=0; u<ops.ux.size(); ++u) lout<<ops.ux[u]<<"  ";
        }
        lout <<std::endl;
    }
    for (int u=0; u<ops.ux.size(); ++u) ops.ux[u]=ops.units.u2i(UnitConv::ULength,ops.ux[u])*sqrt(ops.mass);

    if (ops.f_print_traj && ops.npderiv>0) npout.open((ops.prefix+std::string(".npd")).c_str());
    npout << "# Dnf/Dx^n | x=0\n";

    acout.precision(6); tout.precision(6); lout.precision(6);
    std::cerr<<ops.dt*0.5<<" TIMESTEP IN INT. UNITS\n";

    
    if (ops.prop==PNormalModes && ops.fA.size()>0 && ops.fC.size()>0 )
    {
        if (ops.fA.size()!=ops.nbeads) ERROR("Number of matrices "<<ops.fA.size()<<" does not match bead n.");
        
        fS.resize(ops.nbeads); fT.resize(ops.nbeads);
        for (unsigned long b=0; b<ops.nbeads; ++b)
        {
            ops.fC[b]*=ops.units.u2i(UnitConv::UEnergy,1.); 
            ops.fA[b]*=1./ops.units.u2i(UnitConv::UTime,1.); 

            get_TS(ops.fA[b], ops.fC[b], ops.dt*0.5, fT[b], fS[b]);
            std::cerr<<b<<" fT "<<fT[b]<<" fS "<<fS[b];
        }
    }       
    else
    {
    std::cerr<<"################# MATRICES ###################\n";
    std::cerr<<"****************** A *************************\n"<<ops.A<<"\n";
    std::cerr<<"****************** C *************************\n"<<ops.C<<"\n";

        get_TS(ops.A, ops.C, ops.dt*0.5, T, S);
    std::cerr<<"****************** T *************************\n"<<T<<"\n";
    std::cerr<<"****************** S *************************\n"<<S<<"\n";
    FMatrix<double> chk,Tt,St,tmp;
    transpose(T,Tt); transpose(S,St);
    mult(ops.C,Tt,chk); mult(T,chk,Tt); chk=Tt; chk-=ops.C;
    mult(S,St,tmp);
    chk+=tmp;
    std::cerr<<"******************chk*************************\n"<<chk<<"\n";

    }
    
    
    
    
    std::ofstream mout;
    mout.setf(std::ios::scientific); mout.width(14); mout.precision(6);
    mout.open((ops.prefix+std::string(".A")).c_str());
    mout <<ops.A;
    mout.close();
    mout.open((ops.prefix+std::string(".C")).c_str());
    mout <<ops.C;
    mout.close();
    
    //creates matrix for NM transformation
    double isqrt2=sqrt(0.5);
    mknmc(ops.nbeads, ops.dt);

    std::valarray<std::valarray<double> > pR(std::valarray<double>(n),ops.nbeads);
    std::valarray<double> npR(n), rnd(n);
    unsigned long i, j, k;
    //init 'noise' terms
      
    std::valarray<double> tp(ops.nbeads), tq(ops.nbeads); 
    if (ops.prop==PNormalModes && ops.fA.size()>0 && ops.fC.size()>0 )
    {
        for (unsigned long b=0; b<ops.nbeads; ++b)
        {
            StabCholesky(ops.fC[b],rC);
            pR[b]=0.;
            for (i=0; i<n; ++i) rnd[i]=rgauss();
            for (i=0; i<n; ++i) for (j=0; j<n; ++j)
                pR[b][i]+=rC(i,j)*rnd[j];
        }
        //converts initial moments to NM representation       
        tq=0.0;
        for (i=0; i<n; ++i)
        {
            for (unsigned long b=0; b<ops.nbeads; ++b) tp[b]=pR[b][i];
            nm2bd(tq,tp);
            for (unsigned long b=0; b<ops.nbeads; ++b) pR[b][i]=tp[b];
        }
        std::cerr<<"Initial moments"<< pR<<"\n";
    }
    else
    {
        for (unsigned long b=0; b<ops.nbeads; ++b)
        {
            pR[b]=0.;
            for (i=0; i<n; ++i) rnd[i]=rgauss();
            for (i=0; i<n; ++i) for (j=0; j<n; ++j)
                pR[b][i]+=rC(i,j)*rnd[j];
        }
    }
    for (unsigned long b=0; b<ops.nbeads; ++b) p[b]=pR[b][0];
        
    unsigned long is=0;
    AutoCorrelation<double> acfpot(ops.maxlag), acfkin(ops.maxlag),
                          acftot(ops.maxlag), acftz(ops.maxlag);
    ACOptions<AutoCorrelation<double> > aco;
    acfpot.get_options(aco); aco.timestep=ops.units.i2u(UnitConv::UTime,ops.dt); 
    acfpot.set_options(aco); acfkin.set_options(aco);
    acftot.set_options(aco); acftz.set_options(aco);
    
    HGOptions<Histogram<double> > hgo;
    hgo.window=HGWTriangle;
    
    hgo.boundaries.resize(ops.h_nwin+1);
    Histogram<double> hisp2(hgo), hisq2(hgo), hisp(hgo), hisq(hgo);
    double hpa=pR[0][0],hpb=pR[0][0],hqa=q[0],hqb=q[0],hp2=pR[0][0]*pR[0][0],hq2=pe(q[0])*2.;
    double dvdq;
    double ecns, eth=0;

    for (is=0; is< ops.nstep; ++is)
    {

        if (ops.prop==PBeads)
        {
            //primitive PI integrator
            //thermostat
            for (unsigned long b=0; b<ops.nbeads; ++b) //OK. this is the MOST STUPID PI propagator. 
            {
                npR=0.; 
                eth+=pR[b][0]*pR[b][0]/2.;
                
                for (i=0; i<n; ++i) rnd[i]=rgauss();
                for (i=0; i<n; ++i)
                    for (j=0; j<n; ++j)
                        npR[i]+=T(i,j)*pR[b][j]+S(i,j)*rnd[j];
                pR[b]=npR;
                eth-=pR[b][0]*pR[b][0]/2.;
            }
            
            
            //AGAIN. this is the MOST STUPID PI propagator. 
            for (unsigned long b=0; b<ops.nbeads; ++b) p[b]=pR[b][0];
            verlet(q,p,ops.dt);
            for (unsigned long b=0; b<ops.nbeads; ++b) pR[b][0]=p[b];
            
            for (unsigned long b=0; b<ops.nbeads; ++b) //Did I say it? this is the MOST STUPID PI propagator. Pardon me for my laziness
            {
                eth+=pR[b][0]*pR[b][0]/2.;
                npR=0;
                for (i=0; i<n; ++i) rnd[i]=rgauss();
                for (i=0; i<n; ++i)   
                    for (j=0; j<n; ++j)
                        npR[i]+=T(i,j)*pR[b][j]+S(i,j)*rnd[j];
                pR[b]=npR;
                eth-=pR[b][0]*pR[b][0]/2.;
            }
        }
        else if (ops.prop==PNormalModes)
        {
            bd2nm(q,p);
            for (unsigned long b=0; b<ops.nbeads; ++b) pR[b][0]=p[b];
            for (unsigned long b=0; b<ops.nbeads; ++b) 
            {
                npR=0.; 
                eth+=pR[b][0]*pR[b][0]/2.;
                
                for (i=0; i<n; ++i) rnd[i]=rgauss();

                if (fT.size()!=ops.nbeads)
                    for (i=0; i<n; ++i)
                        for (j=0; j<n; ++j)
                            npR[i]+=T(i,j)*pR[b][j]+S(i,j)*rnd[j];
                else
                    for (i=0; i<n; ++i)
                        for (j=0; j<n; ++j)
                            npR[i]+=fT[b](i,j)*pR[b][j]+fS[b](i,j)*rnd[j];
                    
                pR[b]=npR;
                eth-=pR[b][0]*pR[b][0]/2.;
            }

//            nm2bd(q,p);
            for (unsigned long b=0; b<ops.nbeads; ++b) p[b]=pR[b][0];
            verlet(q,p,ops.dt,PNormalModes);
            for (unsigned long b=0; b<ops.nbeads; ++b) pR[b][0]=p[b];
//            bd2nm(q,p);

            for (unsigned long b=0; b<ops.nbeads; ++b) 
            {
                npR=0.; 
                eth+=pR[b][0]*pR[b][0]/2.;
                
                for (i=0; i<n; ++i) rnd[i]=rgauss();

                if (fT.size()!=ops.nbeads)
                    for (i=0; i<n; ++i)
                        for (j=0; j<n; ++j)
                            npR[i]+=T(i,j)*pR[b][j]+S(i,j)*rnd[j];
                else
                    for (i=0; i<n; ++i)
                        for (j=0; j<n; ++j)
                            npR[i]+=fT[b](i,j)*pR[b][j]+fS[b](i,j)*rnd[j];
                    
                pR[b]=npR;
                eth-=pR[b][0]*pR[b][0]/2.;
            }
            for (unsigned long b=0; b<ops.nbeads; ++b) p[b]=pR[b][0];
            nm2bd(q,p);
        }
        
        
        double epot, ekin, etz, vest, kest, qc;
        vest=kest=ekin=epot=etz=0;
        for (unsigned long b=0; b<ops.nbeads; ++b) ekin+=ops.units.i2u(UnitConv::UEnergy,pR[b][0]*pR[b][0]/2.); 
        for (unsigned long b=0; b<ops.nbeads; ++b) epot+=ops.units.i2u(UnitConv::UEnergy,pe(q[b]));
        vest=epot/ops.nbeads; 
        
        //virial kinetic energy estimator
        qc=0.0; for (unsigned long b=0; b<ops.nbeads; ++b) qc+=q[b]; qc/=ops.nbeads;
        for (unsigned long b=0; b<ops.nbeads; ++b) kest+=(q[b]-qc)*pma(q[b]);
        kest=ops.units.i2u(UnitConv::UEnergy,(ops.temp-kest/ops.nbeads)*0.5);
        
        etz=ekin+epot;
        ecns=etz+ops.units.i2u(UnitConv::UEnergy,eth+lastpiv);
        for (unsigned long b=0; b<ops.nbeads; ++b) for (i=1; i<n; ++i) etz+=ops.units.i2u(UnitConv::UEnergy,pR[b][i]*pR[b][i]/2.);
        
        if (ops.f_histogram && (is<ops.drop || is<100)) //minimum number of data to accumulate before histogram
        {
            for (unsigned long b=0; b<ops.nbeads; ++b) {
                double nq=ops.units.i2u(UnitConv::ULength,q[b])/sqrt(ops.mass);
                double np=ops.units.i2u(UnitConv::ULength,pR[b][0])*sqrt(ops.mass)*ops.units.i2u(UnitConv::UMass,1.)/ops.units.i2u(UnitConv::UTime,1.);
                hpa=(hpa<np?hpa:np);
                hpb=(hpb>np?hpb:np);
                hqa=(hqa<nq?hqa:nq);
                hqb=(hqb>nq?hqb:nq);
            }
            hp2=(hp2>ekin?hp2:ekin);
            hq2=(hq2>epot?hq2:epot);
        }
        if (ops.f_histogram && is==(ops.drop<100?100:ops.drop))
        {
            
            for (int ij=0;ij<=ops.h_nwin;++ij)
            {
                hgo.boundaries[ij]=(hqa+hqb)/2.-(hqb-hqa)*2.*(1.-(2.*ij)/ops.h_nwin);
            }
            hgo.window=HGWTriangle;
            hgo.window_width=(hqb-hqa)/(ops.h_nwin)*0.5;
            hisq.set_options(hgo);
            for (int ij=0;ij<=ops.h_nwin;++ij)
                hgo.boundaries[ij]=(hpa+hpb)/2.-(hpb-hpa)*2.*(1.-(2.*ij)/ops.h_nwin);
            hgo.window_width=(hpb-hpa)/(ops.h_nwin)*0.5;
            hisp.set_options(hgo);
            for (int ij=0;ij<=ops.h_nwin;++ij)
                hgo.boundaries[ij]=(hp2)*(3.*ij)/ops.h_nwin;
            hgo.window_width=hp2/(ops.h_nwin)*0.5;
            hisp2.set_options(hgo);
            for (int ij=0;ij<=ops.h_nwin;++ij)
                hgo.boundaries[ij]=(hq2)*(3.*ij)/ops.h_nwin;
            hgo.window_width=hq2/(ops.h_nwin)*0.5;
            hisq2.set_options(hgo);
        }
        if (is>ops.drop)
        {
            if ((is-ops.drop) % ops.nseries==0)
            {
                std::cerr<<is<<" steps done.\n";
                //after the first series, we collapse on series zero
                if (ops.maxlag>0) 
                {
                    if (is-ops.drop>ops.nseries) 
                    {
                        acfpot.series_collapse();
                        acfkin.series_collapse();
                        acftot.series_collapse();
                        acftz.series_collapse();
                    }
                    ++acfkin;
                    ++acfpot;
                    ++acftot;
                    ++acftz;
                }
            }
                
            
            if (ops.maxlag>0)
            {
                acfpot<<epot;
                acfkin<<ekin;
                acftot<<(ekin+epot);
                acftz<<etz;
            }
            
            if (ops.f_histogram && is>100)
            {
                for (int b=0; b<ops.nbeads; ++b) 
                {
                    hisp<<ops.units.i2u(UnitConv::ULength,pR[b][0])*sqrt(ops.mass)*ops.units.i2u(UnitConv::UMass,1.)/ops.units.i2u(UnitConv::UTime,1.);
                    hisq<<ops.units.i2u(UnitConv::ULength,q[b])/sqrt(ops.mass);
                }
                hisp2<<ekin;
                hisq2<<epot;
            }
            
            if (is%ops.stride==0 && ops.f_print_traj)
            { 
                
                tout<<setw(12)<<std::scientific<<is*ops.units.i2u(UnitConv::UTime,ops.dt)<<" ";
                tout<< setw(12)<<vest<<" " <<setw(12)<< kest<<" ";
                tq=q; bd2nm(tq,tp);
//                std::cerr<<sqrt(pl[0])<<" w " <<sqrt(piwn)<<" wn\n";
                tout<<setw(12)<<ops.units.i2u(UnitConv::UEnergy,pl[0]*tq[0]*tq[0]/2*0.5)<<" "<<setw(12)<<ops.units.i2u(UnitConv::UEnergy,tq[1]*tq[1]/2*(pl[0]+4*piwn)*0.5)<<" ";;
//                  tout<<setw(12)<<ops.units.i2u(UnitConv::UEnergy,pl[0]*q[0]*q[0]*0.5)<<" ";
//                  tout<<setw(12)<<q[0]<<" "<<setw(12)<<q[1]<<" "<<setw(12)<<tq[0]<<" "<<setw(12)<<tq[1]<<" ";
                tout<< setw(12)<<epot<<" " <<setw(12)<< ekin<<" "<<setw(12)<< ecns<<" ";
                for (int b=0; b<ops.nbeads; ++b) tout<<setw(12)<<ops.units.i2u(UnitConv::ULength,q[b]/sqrt(ops.mass))<<" ";
                for (int b=0; b<ops.nbeads; ++b) tout<<setw(12)<<ops.units.i2u(UnitConv::ULength,pR[b][0]*sqrt(ops.mass))*ops.units.i2u(UnitConv::UMass,1.)/ops.units.i2u(UnitConv::UTime,1.)<<" ";
                
                tout<<std::endl;
            }

            if (is%ops.stride==0 && ops.npderiv>0)
            {
               double dx=1e-6;               
               std::valarray<double> y(ops.nbeads), q1(ops.nbeads), pre(ops.nbeads), fi(ops.npderiv+1), npd(ops.npderiv+1); int b;
               for (b=0; b<ops.nbeads; ++b) pre[b]=pe(q[b]);
               npd=0.0;

               //OK, this is a REALLY BAD way of computing derivatives. a lot of re-using could be done, but here I want to be quick&dirty               
               for (int ud=0; ud<=ops.npderiv; ++ud)
               {
                    for (int u=0; u<=ud; ++u)
                    {
                      double dispe=0; 
                      double du=dx*(ud*0.5-u);
                    
                      gtrans_bw(1.0+du,q,q1);
                      for (b=0; b<ops.nbeads; ++b) 
                      { y[b]=pe(q1[b])-pre[b]; dispe+=y[b]; }
                      dispe=exp(-dispe/(ops.temp*ops.nbeads));      
                      fi[u]=dispe;
                    }
                    // now compute the finite-differences derivative
                    for (int i=0; i<=ud; ++i) npd[ud]+=(i%2==0?1:-1)*binomial(ud,i)*fi[i];
                    npd[ud]*=pow(1.0/dx, ud);
               }
               
               //screw it, directly print out the mmoments estimators
               
               //for (int ud=0; ud<=ops.npderiv; ++ud) 
               
               npout<<setw(12)<<std::scientific<<ops.units.i2u(UnitConv::UEnergy,ops.temp*(0.5+ops.nbeads*npd[1]))    <<"  "
                       <<setw(12)<<std::scientific<<pow(ops.units.i2u(UnitConv::UEnergy,1.0),2)*ops.temp*ops.temp*(0.75+ops.nbeads*(2*ops.nbeads+1)*npd[1]+ops.nbeads*ops.nbeads*npd[2])    <<"  ";
               npout<<"\n";
               
            }
            
            if (is%ops.stride==0 && ops.ux.size()>0)
            { 
                double dispe;
                for (int u=0; u< ops.ux.size(); ++u)
                {
                    dispe=0; 
                    if (ops.nlin==0) 
                    {
                        for (int b=0; b<ops.nbeads; ++b) dispe+=(pe(q[b]+(0.5-b*1.0/ops.nbeads)*ops.ux[u])-pe(q[b]))*(b==0?0.5:1.0);
                        //for (int b=0; b<ops.nbeads; ++b) std::cerr<<(pe(q[b]+(0.5-b*1.0/ops.nbeads)*ops.ux[u])-pe(q[b]))<<"  ";
                        //std::cerr<<std::endl;
                        dispe+=(pe(q[0]-0.5*ops.ux[u])-pe(q[0]))*0.5;
                        dispe*=1.0/ops.nbeads; 
                    }
                    else
                    {
                        //for (int b=0; b<ops.nlin; ++b) dispe+=(pe(q[0]+(0.5-b*1.0/ops.nlin)*ops.ux[u])-pe(q[0]))*(b==0?0.5:1.0);
                        //dispe+=(pe(q[0]-0.5*ops.ux[u])-pe(q[0]))*0.5;
/*                        for (int bb=0; bb<ops.nlin; ++bb)
                        for (int b=0; b<ops.nlin; ++b) { 
                            dispe+=(pe(q[bb]+(0.5-((bb+b)%ops.nlin)*1.0/ops.nlin)*ops.ux[u])-pe(q[bb]))*(((bb+b)%ops.nlin)==0?0.5:1.0);
                            dispe+=(pe(q[ops.nlin-1-bb]-0.5*ops.ux[u])-pe(q[ops.nlin-1-bb]))*0.5;
                        }*/
                        std::valarray<double> y(ops.nbeads+1), x(ops.nbeads+1); int b;
                        for (b=0; b<ops.nbeads; ++b) { x[b]=0.5-b*1.0/ops.nbeads; y[b]=pe(q[b]+x[b]*ops.ux[u])-pe(q[b]); }
                        x[b]=-0.5; y[b]=pe(q[0]+x[b]*ops.ux[u])-pe(q[0]);
                        //for (b=0; b<=ops.nbeads; ++b)  std::cerr<<x[b]<<","<<y[b]<<" "; std::cerr<<"\n";
                        //!trapezoids
                        //for (b=0; b<ops.nbeads; ++b) dispe+=-(x[b+1]-x[b])*0.5*(y[b+1]+y[b]); 
                        //!parabolas
                    /*
                        for (b=1; b<ops.nbeads; ++b) dispe+=y[b-1]+22*y[b]+y[b+1];
                        dispe+=8*y[0]+5*y[1]-y[2];dispe+=8*y[ops.nbeads]+5*y[ops.nbeads-1]-y[ops.nbeads-2];
                        dispe*=1.0/(ops.nbeads*24.);
                      */
                         //! modified simpson (numerical recipes)  
                        for (b=3; b<=ops.nbeads-3; ++b) dispe+=y[b];
                        dispe+=3./8.*(y[0]+y[ops.nbeads])+7./6.*(y[1]+y[ops.nbeads-1])+23./24.*(y[2]+y[ops.nbeads-2]);
                        dispe*=1.0/ops.nbeads;
                    }
                    dispe=exp(-dispe/ops.temp);
                    
                    uxavg[u]+=dispe;
                    if (ops.f_print_traj) lout<<setw(12)<<std::scientific<<dispe<<"  ";
                }
                nux++;
                /*tout<<setw(12)<<std::scientific<<is*ops.units.i2u(UnitConv::UTime,ops.dt)<<" ";
                tout<< setw(12)<<epot<<" " <<setw(12)<< ekin<<" "<<setw(12)<< ecns<<" ";
                for (int b=0; b<ops.nbeads; ++b) tout<<setw(12)<<ops.units.i2u(UnitConv::ULength,q[b]/sqrt(ops.mass))<<" ";
                for (int b=0; b<ops.nbeads; ++b) tout<<setw(12)<<ops.units.i2u(UnitConv::ULength,pR[b][0]*sqrt(ops.mass))*ops.units.i2u(UnitConv::UMass,1.)/ops.units.i2u(UnitConv::UTime,1.)<<" ";*/
                
                if (ops.f_print_traj) lout<<std::endl;
            }
        }
    }
    
    if (ops.ux.size()>0)
    {
        laout<<"#   x      U(x)        n(x)\n";
        
        for (int u=0; u< ops.ux.size(); ++u)
            laout<<ops.units.i2u(UnitConv::ULength,ops.ux[u]/sqrt(ops.mass))<<"  "<<-log(uxavg[u]/nux)<<"  "
                <<(uxavg[u]/nux*exp(-ops.ux[u]*ops.ux[u]*ops.temp*0.5))<<std::endl;
        
    }
    if (ops.f_histogram)
    {
        std::ofstream hout;
        hout.open((ops.prefix+std::string(".hisq")).c_str());
        hout<<hisq;
        hout.close();
        hout.open((ops.prefix+std::string(".hisp")).c_str());
        hout<<hisp;
        hout.close();
        hout.open((ops.prefix+std::string(".hisq2")).c_str());
        hout<<hisq2;
        hout.close();
        hout.open((ops.prefix+std::string(".hisp2")).c_str());
        hout<<hisp2;
        hout.close();
    }
    
    if (ops.f_print_traj) tout.close();
    if (ops.maxlag>0) 
    {
        acout
                <<"#######################################################################\n"
                <<"# OUTPUT DATA:                                                        #\n"
                <<"# total trajectory time: "<<std::setw(14)<<ops.units.i2u(UnitConv::UTime,(ops.nstep-ops.drop)*ops.dt)<<"                               #\n"
                <<"# <U> :  "<<std::setw(14)<<acfpot.mean()
                <<" s(U) :  "<<std::setw(14)<<acfpot.sigma()
                <<" t(U) :  "<<std::setw(14)<<acfpot.actime()<<" #\n"
                <<"# <K> :  "<<std::setw(14)<<acfkin.mean()
                <<" s(K) :  "<<std::setw(14)<<acfkin.sigma()
                <<" t(K) :  "<<std::setw(14)<<acfkin.actime()<<" #\n"
                <<"# <H> :  "<<std::setw(14)<<acftot.mean()
                <<" s(H) :  "<<std::setw(14)<<acftot.sigma()
                <<" t(H) :  "<<std::setw(14)<<acftot.actime()<<" #\n"
                <<"# <X> :  "<<std::setw(14)<<acftz.mean()
                <<" s(X) :  "<<std::setw(14)<<acftz.sigma()
                <<" t(X) :  "<<std::setw(14)<<acftz.actime()<<" #\n"
                <<"#######################################################################\n";
        for (unsigned long h=0; h<ops.maxlag; ++h)
            acout<<std::setw(14)<<(h*aco.timestep)<<" "
                    <<std::setw(14)<<acfpot[h]<<" "
                    <<std::setw(14)<<acfkin[h]<<" "
                    <<std::setw(14)<<acftot[h]<<" "
                    <<std::setw(14)<<acftz[h]<<"\n";
        acout.close();
    }
} 

