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
#include "interpol.hpp"
#include <fstream>
/*******************************************************************
  INTERNAL UNITS ARE ATOMIC UNITS
  length=1 bohr=5.29177 10^-11 m
*******************************************************************/


using namespace std;
using namespace toolbox;
enum Potential { PHarmonic, PLennardJones, PQuasiHarmo, PFlexiWell, PDoubleWell, PSquareWell, PMorse, PFree, PGauss, PTable };
enum Forcing { FNone, FSine, FTable, FEnvelope };
class odops 
{
public:
    FMatrix<double> A, C;
    UnitConv units;
    std::string prefix;
    bool f_print_traj, f_histogram; 
    unsigned long h_nwin;
    long seed;
    double dt;
    unsigned long nstep, stride;
    unsigned long maxlag;
    unsigned long nseries;
    unsigned long drop;
    unsigned long nm_test;
    Potential pot;
    double temp, mass;
    double f_alpha, f_beta;
    Forcing force;
    std::valarray<double> f_pars;
    unsigned long autobins; double autoinit;
    std::valarray<double> pars;
};

std::ostream& operator<< (std::ostream& ostr, const odops& oo)
{
    std::vector<string> ppars, pforce;
    
    switch (oo.force)
    {
    case FNone:
        pforce.push_back("none");
        break;
    case FSine:
        if (oo.f_pars.size()<2)
            ERROR("Wrong number of parameters for sinusoidal forcing.");
        pforce.push_back("sine");
        pforce.push_back(float2str(oo.f_pars[0])); pforce.push_back(float2str(oo.f_pars[1]));
        break;
    case FTable:
        if (oo.f_pars.size()%2==1)
            ERROR("Wrong number of parameters for tabulated forcing.");
        pforce.push_back("table");
        for (int i=0; i<oo.f_pars.size(); ++i)
            pforce.push_back(float2str(oo.f_pars[i]));
        break;
    case FEnvelope:
        if (oo.f_pars.size()%2==1)
            ERROR("Wrong number of parameters for sinusoidal forcing w/tabulated envelope.");
        pforce.push_back("envelope");
        for (int i=0; i<oo.f_pars.size(); ++i)
            pforce.push_back(float2str(oo.f_pars[i]));
        break;
    }

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
    case PGauss:
        if (oo.pars.size()<3 || (oo.pars.size())%3 != 0)
            ERROR("Wrong number of parameters for 'gauss' well potential.");
        ppars.push_back("gauss");
        for (int i=0; i<oo.pars.size(); ++i)
          ppars.push_back(float2str(oo.pars[i]));
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
    case PTable:
        if (oo.pars.size()%2==1)
            ERROR("Wrong number of parameters for tabulated potential.");
        ppars.push_back("table");
        for (int i=0; i<oo.pars.size(); ++i)
            ppars.push_back(float2str(oo.pars[i]));
        break;
    }
    
    IOMap iom; 
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
    iom.insert(oo.maxlag,"max_lag");
    iom.insert(oo.nseries,"max_data");
    iom.insert(oo.prefix,"prefix");
    iom.insert(oo.seed,"seed");
    iom.insert(oo.f_print_traj,"print_traj");
    iom.insert(oo.f_histogram,"histogram");
    iom.insert(oo.h_nwin,"histo-bins");
    iom.insert(oo.nm_test,"non-mark-test");
    iom.insert(pforce,"forcing");
    iom.insert(oo.f_alpha,"alpha");
    iom.insert(oo.f_beta,"beta");
    
    iom.insert(oo.autoinit,"autoinit");
    iom.insert(oo.autobins,"autobins");
    ostr<<iom;
    
    return ostr;
}

std::istream& operator>> (std::istream& istr, odops& oo)
{
    IOMap iom;
    double ascale, cscale;
    std::vector<string> ppars, pforce;
    FMatrix<double> Z(0,0);
    iom.insert(ascale,"a_scale",1.0);
    iom.insert(cscale,"c_scale",1.0);
    iom.insert(oo.A,"a_matrix");
    iom.insert(oo.C,"c_matrix");
    iom.insert(oo.dt,"timestep");
    iom.insert(oo.nstep,"steps");
    iom.insert(oo.stride,"stride",(unsigned long) 1);
    iom.insert(ppars,"potential");
    iom.insert(oo.temp,"temperature",0.);
    iom.insert(oo.mass,"mass",1.);
    iom.insert(oo.maxlag,"max_lag");
    iom.insert(oo.nseries,"max_data");
    iom.insert(oo.drop,"discard",(unsigned long) 0);
    iom.insert(oo.prefix,"prefix",std::string("oned"));
    iom.insert(oo.seed,"seed",(long)1234321);
    iom.insert(oo.f_print_traj,"print_traj",true);
    iom.insert(oo.f_histogram,"histogram",true);
    iom.insert(oo.h_nwin,"histo-bins",(unsigned long) 0);
    iom.insert(oo.nm_test,"non-mark-test",(unsigned long) 0);
    iom.insert(oo.units, "units",UnitConv());
    iom.insert(oo.f_alpha,"alpha",0.0);
    iom.insert(oo.f_beta,"beta",0.0);
    iom.insert(pforce,"forcing");
    
    iom.insert(oo.autoinit,"autoinit",0.0);
    iom.insert(oo.autobins,"autobins",(unsigned long) 1000);

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
    else if (ppars[0]=="gauss")
    {
        oo.pot=PGauss;
        if (ppars.size()<4 || (ppars.size()-1)%3!=0)
            ERROR("Wrong number of parameters for Gauss potential.");
        oo.pars.resize(ppars.size()-1); 
        for (int i=0; i<ppars.size()-1; ++i) oo.pars[i]=str2float(ppars[i+1]);
        std::cerr<<ppars<<" read gauss parameters\n"<<oo.pars;
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
    else if (ppars[0]=="table")
    {
        oo.pot=PTable;
        if (ppars.size()<5 || (ppars.size()-1)%2!=0)
            ERROR("Wrong number of parameters for Gauss potential.");
        oo.pars.resize(ppars.size()-1); 
        for (int i=0; i<ppars.size()-1; ++i) oo.pars[i]=str2float(ppars[i+1]);
        std::cerr<<ppars<<" read tabular parameters\n"<<oo.pars;
    }
    
    oo.force=FNone;
    if (pforce.size()>0) {
    if (pforce[0]=="sine")
    {
        oo.force=FSine;
        if (pforce.size()<3)
            ERROR("Wrong number of parameters for sinusoidal forcing.");
        oo.f_pars.resize(2); oo.f_pars[0]=str2float(pforce[1]); oo.f_pars[1]=str2float(pforce[2]);
    }
    else if (pforce[0]=="table")
    {
        oo.force=FTable;
        if (pforce.size()<5 || (pforce.size()-1)%2!=0)
            ERROR("Wrong number of parameters for tabulated forcing.");
        oo.f_pars.resize(pforce.size()-1); 
        for (int i=0; i<pforce.size()-1; ++i) oo.f_pars[i]=str2float(pforce[i+1]);
        std::cerr<<pforce<<" read tabular parameters\n"<<oo.f_pars;
    }
    else if (pforce[0]=="envelope")
    {
        oo.force=FEnvelope;
        if (pforce.size()<7 || (pforce.size()-1)%2!=0)
            ERROR("Wrong number of parameters for sinusoidal forcing with tabulated envelope.");
        oo.f_pars.resize(pforce.size()-1); 
        for (int i=0; i<pforce.size()-1; ++i) oo.f_pars[i]=str2float(pforce[i+1]);
        std::cerr<<pforce<<" read tabular parameters\n"<<oo.f_pars;
    }
    }
    
    oo.A*=ascale;
    if (oo.C.rows()==0)
    {
        if (oo.temp<=0.) ERROR("C matrix or temperature should be given"); 
        //classical matrix for the target temperature
        oo.C.resize(oo.A.rows(),oo.A.cols());
        for (unsigned int i=0; i<oo.A.rows(); ++i) oo.C(i,i)=oo.temp;
    }
    else oo.C*=cscale;
    
    return istr;
}

double (*pf) (double);
std::valarray<double> fl;
InterpolateSpline ftable;
double fTable(double t)
{
    return ftable(t);
}

double fEnv(double t)
{
    return fl[0]*sin(t*fl[1])*ftable(t);
}

double fSine(double t)
{
    return fl[0]*sin(t*fl[1]);
}

double (*pma) (double), (*pe)(double);
std::valarray<double> pl;

InterpolateSpline ptable;
double aTable(double mq)
{
    return -ptable.d(mq);
}

double eTable(double mq)
{
    return ptable(mq);
}


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

/* 
  GAUSSIAN SUM POTENTIAL:
  E=-sum h exp-((x-d)^2*w^2/2)
  pars[3i]=h    (well depth)
  pars[3i+1]=d  (well position)
  pars[3i+2]=w^2/h  (curvature)
*/
double aGAUSS(double mq)
{
    double s=0.;
    for (int i=0; i<pl.size()/3; ++i)
        s+=pl[i*3]*pl[i*3+2]*(mq-pl[i*3+1])*
                exp(-(mq-pl[i*3+1])*(mq-pl[i*3+1])*0.5*pl[i*3+2]);
    //std::cerr<<mq<<","<<s<<pl<<"\n";
    return -s;
}
double eGAUSS(double mq)
{
    double s=0.;
    
    for (int i=0; i<pl.size()/3; ++i)
        s+=pl[i*3]*exp(-(mq-pl[i*3+1])*(mq-pl[i*3+1])*0.5*pl[i*3+2]);
    return -s;
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


void verlet(double& mq, double& mp, double& dt)
{
    static double olda=666e100;
    if (olda==666e100) olda=pma(mq);
    mp+=0.5*dt*olda;
    mq+=dt*mp;
    olda=pma(mq);
    mp+=0.5*dt*olda;
}

void sw_verlet(double& mq, double& mp, double& dt, double& L)
{
    //std::cerr<<" debug: "<<mq<<" ("<<L<<") "<<dt*mp<<"\n";
    mq+=dt*mp;
    if (fabs(mq)>L*0.5) 
    {
        mp=-mp;
        if (mq>0) mq=L-mq; else mq=-(L+mq);
    }
}

int main(int argc, char **argv)
{
    RndGaussian<double, MTRndUniform> rgauss(RGPars<double>(0.,1.));
    
    FMatrix<double>  S, T, rC; 
    odops ops;
    std::cin>>ops;
    std::cerr<<"SEED "<<ops.seed<<"\n";
    rgauss.RUGenerator().seed(ops.seed);
    
    //we convert the options in internal units
    ops.C*=ops.units.u2i(UnitConv::UEnergy,1.); 
    ops.A*=1./ops.units.u2i(UnitConv::UTime,1.); 
    std::cerr<<"QH: "<<ops.temp<<","<<ops.mass<<std::endl;
    ops.temp=ops.units.u2i(UnitConv::UEnergy,ops.temp); 
    ops.mass=ops.units.u2i(UnitConv::UMass,ops.mass); 
    std::cerr<<"QH: "<<ops.temp<<","<<ops.mass<<std::endl;
    ops.dt=ops.units.u2i(UnitConv::UTime,ops.dt); 
    std::cerr<<"QH: "<<ops.dt<<std::endl;
    unsigned long n=ops.A.rows();
    if (ops.temp>0.) 
    {
        //temperature control!
        std::cerr<<"RESCALING MATRICES FOR TEMPERATURE CONTROL!"<<std::endl;
        ops.A*=ops.temp/ops.C(0,0);
        ops.C*=ops.temp/ops.C(0,0);
    }
    StabCholesky(ops.C,rC);
    
    double q=0;     pl.resize(ops.pars.size());
    std::valarray<double> tx(ops.pars.size()/2), ty(ops.pars.size()/2);
    switch(ops.pot)
    {
        case PHarmonic:
            //pars[0] is the oscillator FREQUENCY (nu)
            ops.pars[0]=ops.units.u2i(UnitConv::UFrequency,ops.pars[0]); 
            pl[0]=ops.pars[0]*ops.pars[0];
            pma=&aHarmonic; pe=&eHarmonic;
            q=rgauss()/ops.pars[0]*rC(0,0);
            std::cerr<<q<<","<<q*q/2*pl[0]<<","<<pe(q)<<"\n";
            std::cerr<<"conv q: "<<ops.units.i2u(UnitConv::ULength,q)<<"\n";
            std::cerr<<"conv q: "<<ops.units.i2u(UnitConv::ULength,q/sqrt(ops.mass))<<"\n";
            break;
        case PQuasiHarmo:
            std::cerr<<"QH: "<<ops.pars[0]<<","<<ops.pars[1]<<std::endl;
            //pars[0] is the oscillator FREQUENCY (nu)
            ops.pars[0]=ops.units.u2i(UnitConv::UFrequency,ops.pars[0]); 
            //pars[1] is an inverse length, and we should also consider the mass-scaling
            ops.pars[1]*=1./(ops.units.u2i(UnitConv::ULength,1.)*sqrt(ops.mass)); 
            std::cerr<<"QH: "<<ops.pars[0]<<","<<ops.pars[1]<<std::endl;
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
            q=ops.pars[0]/sqrt(ops.mass); //it will decorrelate afterwards, we start in the bottom of a well
            std::cerr<<"initializing double well into "<<q<<std::endl;
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
        case PGauss:   //input parameters are well depths, positions and harmonic frequency
            pma=&aGAUSS; pe=&eGAUSS;
            for (int i=0; i<ops.pars.size()/3; ++i)
            {
                ops.pars[3*i]=ops.units.u2i(UnitConv::UEnergy,ops.pars[3*i]);
                ops.pars[3*i+1]=ops.units.u2i(UnitConv::ULength,ops.pars[3*i+1])*sqrt(ops.mass);
                ops.pars[3*i+2]=ops.units.u2i(UnitConv::UFrequency,ops.pars[3*i+2]); 
                pl[3*i]=ops.pars[3*i];
                pl[3*i+1]=ops.pars[3*i+1];
                pl[3*i+2]=ops.pars[3*i+2]*ops.pars[3*i+2]/ops.pars[3*i];
            }
            std::cerr<<ops.pars<<pl<<"gauss para\n";
            q=0.;
            break;
        case PSquareWell:
            pma=&aSW; pe=&eSW;
            ops.pars[0]=ops.units.u2i(UnitConv::ULength,ops.pars[0])*sqrt(ops.mass);
            std::cerr<<"UNIT LENGTH "<<ops.pars[0]<<"\n";
            pl[0]=ops.pars[0];
            q=0.;
            break;
        case PTable:
            pma=&aTable; pe=&eTable;
            for (int i=0; i<ops.pars.size()/2; ++i)
            {
                tx[i]=ops.units.u2i(UnitConv::ULength,ops.pars[2*i])*sqrt(ops.mass);
                ty[i]=ops.units.u2i(UnitConv::UEnergy,ops.pars[2*i+1]);
            }
            ptable.set_table(tx,ty);
            q=0.;
            break;
        case PFree:
            pma=&aFree; pe=&eFree;
            q=0.;
            break;
        default:
            ERROR("Potential not yet implemented!\n");
    }
    
    if (ops.autoinit>0.0) 
    {
        ops.autoinit=ops.units.u2i(UnitConv::ULength,ops.autoinit)*sqrt(ops.mass);
        std::valarray<double> autov(ops.autobins);
        double minv=pe(-ops.autoinit), dx=2*ops.autoinit/(ops.autobins-1), totp=0.0, pchoice;
        for (unsigned long i=0; i<ops.autobins; i++) { autov[i]=pe(-ops.autoinit+i*dx); if (minv>autov[i]) minv=autov[i]; }
        for (unsigned long i=0; i<ops.autobins; i++) { autov[i]=exp(-(autov[i]-minv)/ops.temp); totp+=autov[i]; }
        pchoice=totp*rgauss.RUGenerator()();
        std::cerr<<minv<<" of "<<totp<<" CHOSEN: "<<pchoice<<"\n";
        unsigned long ib=0;
        totp=0.0; while (pchoice>totp) { totp+=autov[ib]; ib++; }
        q=(ib-0.5)*dx-ops.autoinit;
    }

    fl.resize(ops.f_pars.size());
    
    switch(ops.force)
    {
    case FSine:
        fl[0]=ops.f_pars[0];  
        fl[1]=ops.units.u2i(UnitConv::UFrequency,ops.f_pars[1]);
        pf=fSine;
        break;
    case FTable:
        pf=fTable; 
        tx.resize(ops.f_pars.size()/2); ty.resize(ops.f_pars.size()/2);
        for (int i=0; i<ops.f_pars.size()/2; ++i)
        {
            tx[i]=ops.units.u2i(UnitConv::UTime,ops.f_pars[2*i]);
            ty[i]=ops.f_pars[2*i+1];
        }
        fl[0]=tx[0];
        fl[1]=tx[ops.f_pars.size()/2-1];
        ftable.set_table(tx,ty);
        break;
    case FEnvelope:
        pf=fEnv; 
        fl[0]=ops.f_pars[0];  
        fl[1]=ops.units.u2i(UnitConv::UFrequency,ops.f_pars[1]);
        tx.resize(ops.f_pars.size()/2-2); ty.resize(ops.f_pars.size()/2-2);
        for (int i=2; i<ops.f_pars.size()/2; ++i)
        {
            tx[i-2]=ops.units.u2i(UnitConv::UTime,ops.f_pars[2*i]);
            ty[i-2]=ops.f_pars[2*i+1];
        }
        fl[2]=tx[0];
        fl[3]=tx[ops.f_pars.size()/2-1-2];
        ftable.set_table(tx,ty);
        break;
    }

    std::ofstream acout((ops.prefix+std::string(".acf")).c_str());
    std::ofstream tout;
    if (ops.f_print_traj)  tout.open((ops.prefix+std::string(".traj")).c_str());
    acout.precision(6); tout.precision(6);
    std::cerr<<ops.dt*0.5<<" TIMESTEP IN INT. UNITS\n";
    get_TS(ops.A, ops.C, ops.dt*0.5, T, S);
    
    std::cerr<<"################# MATRICES ###################\n";
    std::cerr<<"****************** A *************************\n"<<ops.A<<"\n";
    std::cerr<<"****************** C *************************\n"<<ops.C<<"\n";
    std::cerr<<"****************** T *************************\n"<<T<<"\n";
    std::cerr<<"****************** S *************************\n"<<S<<"\n";
    
    FMatrix<double> chk,Tt,St,tmp;
    transpose(T,Tt); transpose(S,St);
    mult(ops.C,Tt,chk); mult(T,chk,Tt); chk=Tt; chk-=ops.C;
    mult(S,St,tmp);
    chk+=tmp;
    std::cerr<<"******************chk*************************\n"<<chk<<"\n";
    
    std::ofstream mout;
    mout.setf(std::ios::scientific); mout.width(14); mout.precision(6);
    mout.open((ops.prefix+std::string(".A")).c_str());
    mout <<ops.A;
    mout.close();
    mout.open((ops.prefix+std::string(".C")).c_str());
    mout <<ops.C;
    mout.close();
    
    std::valarray<double> pR(n), npR(n), rnd(n);
    unsigned long i, j;
    for (i=0; i<n; ++i) rnd[i]=rgauss();
    //init 'noise' terms
    pR=0;
    for (i=0; i<n; ++i)
        for (j=0; j<n; ++j)
            pR[i]+=rC(i,j)*rnd[j];
    
    //check of non-markovian correspondance
    std::valarray<double> K(ops.nm_test), p_hist(ops.nm_test);
    if (ops.nm_test>0) A2Kt(ops.A,ops.dt,K);
    unsigned long is=0;
    
    AutoCorrelation<double> acfp(ops.maxlag), acfq(ops.maxlag), acfpot(ops.maxlag), acfkin(ops.maxlag),
                          acftot(ops.maxlag), acftz(ops.maxlag), acfnmZ(ops.nm_test);
    ACOptions<AutoCorrelation<double> > aco;
    acfp.get_options(aco); aco.timestep=ops.units.i2u(UnitConv::UTime,ops.dt); 
    acfp.set_options(aco); acfq.set_options(aco); acfpot.set_options(aco); acfkin.set_options(aco);
    acftot.set_options(aco); acftz.set_options(aco);  acfnmZ.set_options(aco);
    
    HGOptions<Histogram<double> > hgo;
    hgo.window=HGWTriangle;
    
    hgo.boundaries.resize(ops.h_nwin+1);
    Histogram<double> hisp2(hgo), hisq2(hgo), hisp(hgo), hisq(hgo);
    double hpa=pR[0],hpb=pR[0],hqa=q,hqb=q,hp2=pR[0]*pR[0],hq2=pe(q)*2.;
    double dp, op;
    double dvdq;
    double ecns, eth=0;
    std::cerr<<"starting\n";
    for (is=0; is< ops.nstep; ++is)
    {
        //thermostat
        npR=dp=0; 
        op=pR[0]; //for markov-check!
        eth+=pR[0]*pR[0]/2.;
        for (i=0; i<n; ++i) rnd[i]=rgauss();
        for (i=0; i<n; ++i)
            for (j=0; j<n; ++j)
                npR[i]+=T(i,j)*pR[j]+S(i,j)*rnd[j];
        dp+=npR[0]-pR[0];
        pR=npR;
        eth-=pR[0]*pR[0]/2.;

        if (ops.pot==PSquareWell) sw_verlet(q,pR[0],ops.dt,pl[0]);
        else verlet(q,pR[0],ops.dt);
        
        if (ops.force!=FNone && (ops.f_alpha>0 || ops.f_beta>0))
        {
            double fnow=pf(is*ops.dt);
            //std::cerr<<fnow<<" forcing\n";
            pR[0]+=fnow*(ops.f_alpha+ops.f_beta*fnow)*ops.dt;
        }

        //std::cerr<<" TEST " <<q<<" "<<eTable(q)<<" "<<aTable(q)<<"\n";
        eth+=pR[0]*pR[0]/2.;
        npR=0;
        for (i=0; i<n; ++i) rnd[i]=rgauss();
        for (i=0; i<n; ++i)   
            for (j=0; j<n; ++j)
                npR[i]+=T(i,j)*pR[j]+S(i,j)*rnd[j];
        dp+=npR[0]-pR[0];
        pR=npR;
        eth-=pR[0]*pR[0]/2.;
    
        double epot, ekin, etz;
        ekin=ops.units.i2u(UnitConv::UEnergy,pR[0]*pR[0]/2.);
        epot=ops.units.i2u(UnitConv::UEnergy,pe(q));
        etz=ekin+epot;
        ecns=etz+ops.units.i2u(UnitConv::UEnergy,eth);
        for (i=1; i<n; ++i) etz+=ops.units.i2u(UnitConv::UEnergy,pR[i]*pR[i]/2.);
        
        if(ops.nm_test>0)
        {
            //non-markovian test
            p_hist[is%ops.nm_test]=op;
            double kint;
            dp*=1./ops.dt;
            if(is>ops.drop && is>ops.nm_test)
            {
                kint=op*0.5*K[0]; 
                for (int s=1; s<ops.nm_test; ++s)
                    kint+=p_hist[(is-s)%ops.nm_test]*K[s];
                kint*=ops.dt;
                //std::cout<<kint<<"  "<<dp<<"  "<<op<<"\n";
                //std::cerr<<op*(0.5*K[0]*ops.dt-ops.A(0,0))+ops.dt*kint<<"\n";
                dp+=kint; //this is the noise!
                acfnmZ<<dp;
            }
            else dp=0;
            dvdq=pma(q);
        }
        if (ops.f_histogram && (is<ops.drop || is<100)) //minimum number of data to accumulate before histogram
        {
            double nq=ops.units.i2u(UnitConv::ULength,q)/sqrt(ops.mass);
            double np=ops.units.i2u(UnitConv::ULength,pR[0])*sqrt(ops.mass)*ops.units.i2u(UnitConv::UMass,1.)/ops.units.i2u(UnitConv::UTime,1.);
            hpa=(hpa<np?hpa:np);
            hpb=(hpb>np?hpb:np);
            hqa=(hqa<nq?hqa:nq);
            hqb=(hqb>nq?hqb:nq);
            hp2=(hp2>ekin?hp2:ekin);
            hq2=(hq2>epot?hq2:epot);
        }
        if (ops.f_histogram && is==(ops.drop<100?100:ops.drop))
        {
            std::cerr<<"histogram boundaries "<<hqa<<","<<hqb<<" : "<<hpa<<","<<hpb<<" : "<<hq2<<" : "<<hp2<<"\n";
            for (int ij=0;ij<=ops.h_nwin;++ij)
            {
                hgo.boundaries[ij]=(hqa+hqb)/2.-(hqb-hqa)*2.*(1.-(2.*ij)/ops.h_nwin);
                std::cerr<<hgo.boundaries[ij]<<" ";}
            hgo.window_width=(hqb-hqa)/(ops.h_nwin*8.);
            hisq.set_options(hgo);
            for (int ij=0;ij<=ops.h_nwin;++ij)
                hgo.boundaries[ij]=(hpa+hpb)/2.-(hpb-hpa)*2.*(1.-(2.*ij)/ops.h_nwin);
            hgo.window_width=(hpb-hpa)/(ops.h_nwin*8.);
            hisp.set_options(hgo);
            for (int ij=0;ij<=ops.h_nwin;++ij)
                hgo.boundaries[ij]=(hp2)*(3.*ij)/ops.h_nwin;
            hgo.window_width=hp2/(ops.h_nwin*5.);
            hisp2.set_options(hgo);
            for (int ij=0;ij<=ops.h_nwin;++ij)
                hgo.boundaries[ij]=(hq2)*(3.*ij)/ops.h_nwin;
            hgo.window_width=hq2/(ops.h_nwin*5.);
            hisq2.set_options(hgo);
        }
        if (is>ops.drop+ops.nm_test)
        {
            if ((is-ops.drop) % ops.nseries==0)
            {
                std::cerr<<is<<" steps done.\n";
                //after the first series, we collapse on series zero
                if (ops.maxlag>0) 
                {
                    if (is-ops.drop>ops.nseries) 
                    {
                        acfp.series_collapse();
                        acfq.series_collapse();
                        acfpot.series_collapse();
                        acfkin.series_collapse();
                        acftot.series_collapse();
                        acftz.series_collapse();
                    }
                    ++acfp;
                    ++acfq;
                    ++acfkin;
                    ++acfpot;
                    ++acftot;
                    ++acftz;
                    if(ops.nm_test>0 && is>ops.nm_test)
                    {
                        if (is-ops.drop>ops.nseries) 
                            acfnmZ.series_collapse();
                        ++acfnmZ;
                    }
                }
            }
                
            
            if (ops.maxlag>0)
            {
                acfp<<ops.units.i2u(UnitConv::ULength,pR[0])*sqrt(ops.mass)*ops.units.i2u(UnitConv::UMass,1.)/ops.units.i2u(UnitConv::UTime,1.);
                acfq<<ops.units.i2u(UnitConv::ULength,q)/sqrt(ops.mass);
                acfpot<<epot;
                acfkin<<ekin;
                acftot<<(ekin+epot);
                acftz<<etz;
            }
            
            if (ops.f_histogram && is>100)
            {
                hisp<<ops.units.i2u(UnitConv::ULength,pR[0])*sqrt(ops.mass)*ops.units.i2u(UnitConv::UMass,1.)/ops.units.i2u(UnitConv::UTime,1.);
                hisq<<ops.units.i2u(UnitConv::ULength,q)/sqrt(ops.mass);
                hisp2<<ekin;
                hisq2<<epot;
            }

            if (is%ops.stride==0 && ops.f_print_traj)
            { 
                tout<<setw(12)<<std::scientific<<is*ops.units.i2u(UnitConv::UTime,ops.dt)<<" "<<
                        setw(12)<<ops.units.i2u(UnitConv::ULength,q/sqrt(ops.mass))<<" "<<
                        setw(12)<<ops.units.i2u(UnitConv::ULength,pR[0]*sqrt(ops.mass))*ops.units.i2u(UnitConv::UMass,1.)/ops.units.i2u(UnitConv::UTime,1.)<<" "<<
                        setw(12)<<epot<<" " <<
                        setw(12)<< ekin<<" "<<
                        setw(12)<< ecns<<" ";
                //for (i=0; i<n; ++i) tout<<setw(12)<<pR[i]<<" ";
                if (ops.nm_test>0)
                    tout<<setw(12)<<ops.units.i2u(UnitConv::ULength,dp*sqrt(ops.mass))*ops.units.i2u(UnitConv::UMass,1.)/pow(ops.units.i2u(UnitConv::UTime,1.),2)<<" "<<setw(12)<<ops.units.i2u(UnitConv::UEnergy,dvdq)/ops.units.i2u(UnitConv::ULength,1.)<<" ";
                tout<<std::endl;
            }
        }
    }
    if (ops.nm_test>0)
    {
        std::ofstream onmt;
        onmt.open((ops.prefix+std::string(".nmt")).c_str());
        for (int s=0; s<ops.nm_test; ++s)
            onmt<<s*ops.units.i2u(UnitConv::UTime,ops.dt)<<"  "<< K[s]/pow(ops.units.i2u(UnitConv::UTime,1.),0)<<"  "<<acfnmZ[s]*pow(acfnmZ.sigma(),2.)<<"\n";
        onmt.close();
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
    acout
            <<"#######################################################################\n"
            <<"# OUTPUT DATA:                                                        #\n"
            <<"# total trajectory time: "<<std::setw(14)<<ops.units.i2u(UnitConv::UTime,(ops.nstep-ops.drop)*ops.dt)<<"                               #\n"
            <<"# <Q> :  "<<std::setw(14)<<acfq.mean()
            <<" s(Q) :  "<<std::setw(14)<<acfq.sigma()
            <<" t(Q) :  "<<std::setw(14)<<acfq.actime()<<" #\n"
            <<"# <P> :  "<<std::setw(14)<<acfp.mean()
            <<" s(P) :  "<<std::setw(14)<<acfp.sigma()
            <<" t(P) :  "<<std::setw(14)<<acfp.actime()<<" #\n"
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
                <<std::setw(14)<<acfq [h]<<" "
                <<std::setw(14)<<acfp [h]<<" "
                <<std::setw(14)<<acfpot[h]<<" "
                <<std::setw(14)<<acfkin[h]<<" "
                <<std::setw(14)<<acftot[h]<<" "
                <<std::setw(14)<<acftz[h]<<"\n";
    acout.close();
} 

