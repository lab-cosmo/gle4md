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
enum Potential { PHarmonic, PTable };
enum Forcing { FNone, FSine, FTable, FEnvelope };
class odops 
{
public:
    FMatrix<double> A, C;
    UnitConv units;
    std::string prefix;
    bool f_print_traj;
    long seed;
    double dt;
    unsigned long nstep, stride;
    unsigned long drop;
    Potential pot;
    double temp, mass;
    double f_alpha, f_beta;
    Forcing forcex, forcey;
    std::valarray<double> f_parsx, f_parsy;
    unsigned long autobins; double autoinitx, autoinity;
    std::valarray<double> pars;
};

std::ostream& operator<< (std::ostream& ostr, const odops& oo)
{
    std::vector<string> ppars, pforcex, pforcey;
    
    switch (oo.forcex)
    {
    case FNone:
        pforcex.push_back("none");
        break;
    case FSine:
        if (oo.f_parsx.size()<2)
            ERROR("Wrong number of parameters for sinusoidal forcing.");
        pforcex.push_back("sine");
        pforcex.push_back(float2str(oo.f_parsx[0])); pforcex.push_back(float2str(oo.f_parsx[1]));
        break;
    case FTable:
        if (oo.f_parsx.size()%2==1)
            ERROR("Wrong number of parameters for tabulated forcing.");
        pforcex.push_back("table");
        for (int i=0; i<oo.f_parsx.size(); ++i)
            pforcex.push_back(float2str(oo.f_parsx[i]));
        break;
    case FEnvelope:
        if (oo.f_parsx.size()%2==1)
            ERROR("Wrong number of parameters for sinusoidal forcing w/tabulated envelope.");
        pforcex.push_back("envelope");
        for (int i=0; i<oo.f_parsx.size(); ++i)
            pforcex.push_back(float2str(oo.f_parsx[i]));
        break;
    }

    switch (oo.forcey)
    {
        case FNone:
            pforcey.push_back("none");
            break;
        case FSine:
            if (oo.f_parsy.size()<2)
                ERROR("Wrong number of parameters for sinusoidal forcing.");
            pforcey.push_back("sine");
            pforcey.push_back(float2str(oo.f_parsy[0])); pforcey.push_back(float2str(oo.f_parsy[1]));
            break;
        case FTable:
            if (oo.f_parsy.size()%2==1)
                ERROR("Wrong number of parameters for tabulated forcing.");
            pforcey.push_back("table");
            for (int i=0; i<oo.f_parsy.size(); ++i)
                pforcey.push_back(float2str(oo.f_parsy[i]));
            break;
        case FEnvelope:
            if (oo.f_parsy.size()%2==1)
                ERROR("Wrong number of parameters for sinusoidal forcing w/tabulated envelope.");
            pforcey.push_back("envelope");
            for (int i=0; i<oo.f_parsy.size(); ++i)
                pforcey.push_back(float2str(oo.f_parsy[i]));
            break;
    }
    
    switch (oo.pot)
    {
    case PHarmonic:
        if (oo.pars.size()<3)
            ERROR("Wrong number of parameters for harmonic potential.");
        ppars.push_back("harmonic");
        for (unsigned long i=0; i<3; ++i) ppars.push_back(float2str(oo.pars[i]));
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
    iom.insert(oo.prefix,"prefix");
    iom.insert(oo.seed,"seed");
    iom.insert(oo.f_print_traj,"print_traj");
    iom.insert(pforcex,"forcing-x");
    iom.insert(pforcey,"forcing-y");
    iom.insert(oo.f_alpha,"alpha");
    iom.insert(oo.f_beta,"beta");
    
    iom.insert(oo.autoinitx,"autoinit-x");
    iom.insert(oo.autoinity,"autoinit-y");
    iom.insert(oo.autobins,"autobins");
    ostr<<iom;
    
    return ostr;
}

std::istream& operator>> (std::istream& istr, odops& oo)
{
    IOMap iom;
    double ascale, cscale;
    std::vector<string> ppars, pforcex, pforcey;
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
    iom.insert(oo.drop,"discard",(unsigned long) 0);
    iom.insert(oo.prefix,"prefix",std::string("oned"));
    iom.insert(oo.seed,"seed",(long)1234321);
    iom.insert(oo.f_print_traj,"print_traj",true);
    iom.insert(oo.units, "units",UnitConv());
    iom.insert(oo.f_alpha,"alpha",0.0);
    iom.insert(oo.f_beta,"beta",0.0);
    iom.insert(pforcex,"forcing-x");
    iom.insert(pforcey,"forcing-y");
    
    iom.insert(oo.autoinitx,"autoinit-x",0.0); iom.insert(oo.autoinity,"autoinit-y",0.0);
    iom.insert(oo.autobins,"autobins",(unsigned long) 100);

    istr>>iom;
    //sets INTERNAL units
    oo.units.ITime=UnitConv::AtomicTime; oo.units.IEnergy=UnitConv::Hartree;
    oo.units.IFrequency=UnitConv::AtomicFrequency; oo.units.IMass=UnitConv::AtomicMass;

    if (ppars[0]=="harmonic")
    {
        oo.pot=PHarmonic;
        if (ppars.size()<4)
            ERROR("Wrong number of parameters for harmonic potential.");
        oo.pars.resize(3); oo.pars[0]=str2float(ppars[1]); oo.pars[1]=str2float(ppars[2]); oo.pars[2]=str2float(ppars[3]);
    }
    else if (ppars[0]=="table")
    {
        oo.pot=PTable;
        std::cerr<<"Table pots " <<ppars<<"\n";
        if (ppars.size()<2) ERROR("Wrong number of parameters for tabulated potential.");
        
        std::ifstream spfile(ppars[1].c_str());
        std::vector<double> vpars; double vn;
        while (spfile.good())  { spfile>>vn; vpars.push_back(vn); }
        oo.pars.resize(vpars.size()); 
        for (int i=0; i<vpars.size(); ++i) oo.pars[i]=vpars[i];
        //std::cerr<<ppars<<" read tabular parameters\n"<<oo.pars;
    }
    
    oo.forcex=FNone;
    if (pforcex.size()>0) {
    if (pforcex[0]=="sine")
    {
        oo.forcex=FSine;
        if (pforcex.size()<3)
            ERROR("Wrong number of parameters for sinusoidal forcing.");
        oo.f_parsx.resize(2); oo.f_parsx[0]=str2float(pforcex[1]); oo.f_parsx[1]=str2float(pforcex[2]);
    }
    else if (pforcex[0]=="table")
    {
        oo.forcex=FTable;
        if (pforcex.size()<5 || (pforcex.size()-1)%2!=0)
            ERROR("Wrong number of parameters for tabulated forcing.");
        oo.f_parsx.resize(pforcex.size()-1); 
        for (int i=0; i<pforcex.size()-1; ++i) oo.f_parsx[i]=str2float(pforcex[i+1]);
    }
    else if (pforcex[0]=="envelope")
    {
        oo.forcex=FEnvelope;
        if (pforcex.size()<7 || (pforcex.size()-1)%2!=0)
            ERROR("Wrong number of parameters for sinusoidal forcing with tabulated envelope.");
        oo.f_parsx.resize(pforcex.size()-1); 
        for (int i=0; i<pforcex.size()-1; ++i) oo.f_parsx[i]=str2float(pforcex[i+1]);
    }
    }
    
    oo.forcey=FNone;
    if (pforcey.size()>0) {
    if (pforcey[0]=="sine")
    {
        oo.forcey=FSine;
        if (pforcey.size()<3)
            ERROR("Wrong number of parameters for sinusoidal forcing.");
        oo.f_parsy.resize(2); oo.f_parsy[0]=str2float(pforcey[1]); oo.f_parsy[1]=str2float(pforcey[2]);
    }
    else if (pforcey[0]=="table")
    {
        oo.forcey=FTable;
        if (pforcey.size()<5 || (pforcey.size()-1)%2!=0)
            ERROR("Wrong number of parameters for tabulated forcing.");
        oo.f_parsy.resize(pforcey.size()-1); 
        for (int i=0; i<pforcey.size()-1; ++i) oo.f_parsy[i]=str2float(pforcey[i+1]);
    }
    else if (pforcey[0]=="envelope")
    {
        oo.forcey=FEnvelope;
        if (pforcey.size()<7 || (pforcey.size()-1)%2!=0)
            ERROR("Wrong number of parameters for sinusoidal forcing with tabulated envelope.");
        oo.f_parsy.resize(pforcey.size()-1); 
        for (int i=0; i<pforcey.size()-1; ++i) oo.f_parsy[i]=str2float(pforcey[i+1]);
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

double (*pfx) (double);
std::valarray<double> flx;
InterpolateSpline ftablex;
double fTablex(double t)
{
    return ftablex(t);
}

double fEnvx(double t)
{
    return flx[0]*sin(t*flx[1])*ftablex(t);
}

double fSinex(double t)
{
    return flx[0]*sin(t*flx[1]);
}

double (*pfy) (double);
std::valarray<double> fly;
InterpolateSpline ftabley;
double fTabley(double t)
{
    return ftabley(t);
}

double fEnvy(double t)
{
    return fly[0]*sin(t*fly[1])*ftabley(t);
}

double fSiney(double t)
{
    return fly[0]*sin(t*fly[1]);
}

void (*pea) (const std::valarray<double>&, double&, std::valarray<double>&);
std::valarray<double> pl;

InterpolateBicubic ptable;
void eaTable(const std::valarray<double>& q, double& e, std::valarray<double>& a)
{
    toolbox::fixarray<double,2> x,dy;
    x[0]=q[0]; x[1]=q[1];
    ptable.get_ydy(x,e,dy);
    a[0]=-dy[0]; a[1]=-dy[1];
}

/* 
  HARMONIC POTENTIAL:
  E=1/2w^2 x^2
  pars[0]=wx (frequency x)
  pars[1]=wy (frequency y)
  pars[2]=theta (angle for eigenvecs)
  pl[0]=Hxx
  pl[1]=Hyy
  pl[1]=Hxy
*/
void eaHarmonic(const std::valarray<double>& q, double& e, std::valarray<double>& a)
{
    a[0]=-(q[0]*pl[0]+q[1]*pl[2]);
    a[1]=-(q[0]*pl[2]+q[1]*pl[1]);
    e=-0.5*(a[0]*q[0]+a[1]*q[1]);
}

double verlet(std::valarray<double>& mq, std::valarray<double>& mp, const double& dt)
{
    static valarray<double> olda(666e100,2); double ne;
    if (olda[0]==666e100) { pea(mq,ne,olda); olda*=0.5*dt; }
    mp+=olda;
    olda=mp; olda*=dt; mq+=olda;
    pea(mq,ne,olda); olda*=0.5*dt;
    mp+=olda;
    return ne;
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
    
    std::valarray<double> q(0.,2);     pl.resize(ops.pars.size());
    unsigned long nx,ny,k;
    std::valarray<double> tx, ty; FMatrix<double> tz;
    switch(ops.pot)
    {
        case PHarmonic:
            //pars[0] is the oscillator FREQUENCY (nu)
            ops.pars[0]=ops.units.u2i(UnitConv::UFrequency,ops.pars[0]); 
            ops.pars[1]=ops.units.u2i(UnitConv::UFrequency,ops.pars[1]); 
            
            pl[0]=ops.pars[0]*ops.pars[0]*cos(ops.pars[2])*cos(ops.pars[2])+ops.pars[1]*ops.pars[1]*sin(ops.pars[2])*sin(ops.pars[2]);
            pl[1]=ops.pars[1]*ops.pars[1]*cos(ops.pars[2])*cos(ops.pars[2])+ops.pars[0]*ops.pars[0]*sin(ops.pars[2])*sin(ops.pars[2]);
            pl[2]=(ops.pars[1]*ops.pars[1]-ops.pars[0]*ops.pars[0])* cos(ops.pars[2])*sin(ops.pars[2]);
            
            pea=&eaHarmonic; 
            q=0.; //rgauss()/ops.pars[0]*rC(0,0);
            break;
        case PTable:
            pea=&eaTable; 
            nx=ops.pars[0], ny=ops.pars[1]; k=2;
            tx.resize(nx); ty.resize(ny); tz.resize(nx,ny);
            std::cerr<<"n of pars "<<nx<<","<<ny<<" >> "<<ops.pars.size()<<"\n";
            for (int i=0; i<nx; ++i)
                tx[i]=ops.units.u2i(UnitConv::ULength,ops.pars[k++])*sqrt(ops.mass);
            for (int i=0; i<ny; ++i)
                ty[i]=ops.units.u2i(UnitConv::ULength,ops.pars[k++])*sqrt(ops.mass);
            for (int i=0; i<nx; ++i) for (int j=0; j<ny; ++j)
                    tz(i,j)=ops.units.u2i(UnitConv::UEnergy,ops.pars[k++]);
            ptable.set_table(tx,ty,tz,FTensor<double,3>());
            q=0.;
            break;
        default:
            ERROR("Potential not yet implemented!\n");
    }
    
    std::valarray<double> qa(2);
    if (ops.autoinitx>0.0) 
    {
        ops.autoinitx=ops.units.u2i(UnitConv::ULength,ops.autoinitx)*sqrt(ops.mass);
        ops.autoinity=ops.units.u2i(UnitConv::ULength,ops.autoinity)*sqrt(ops.mass);
        
        std::valarray<double> autov(ops.autobins*ops.autobins);
        q[0]=-ops.autoinitx; q[1]=-ops.autoinity; 
        double minv, dx, dy, totp, pchoice;
        pea(q, minv, qa); dx=2*ops.autoinitx/(ops.autobins-1); dy=2*ops.autoinity/(ops.autobins-1);
        totp=0.0; k=0;
        std::ofstream pdbg("pot-debug");
        //nullstream pdbg;
        
        for (unsigned long i=0; i<ops.autobins; i++) 
        {
            q[0]=-ops.autoinitx+i*dx;
            for (unsigned long j=0; j<ops.autobins; j++) 
            {
                q[1]=-ops.autoinity+j*dy;
                pea(q,autov[k],qa); 
                pdbg<<ops.units.i2u(UnitConv::ULength,q[0])/sqrt(ops.mass)<<" "
                    <<ops.units.i2u(UnitConv::ULength,q[1])/sqrt(ops.mass)<<" "
                    <<ops.units.i2u(UnitConv::UEnergy,autov[k])<<"\n";
                if (minv>autov[k]) minv=autov[k]; 
                k++;
            }
            pdbg<<"\n";
        }
        for (unsigned long i=0; i<autov.size(); i++) 
        { autov[i]=exp(-(autov[i]-minv)/ops.temp); totp+=autov[i]; }
        pchoice=totp*rgauss.RUGenerator()();
        std::cerr<<minv<<" of "<<totp<<" CHOSEN: "<<pchoice<<"\n";
        
        unsigned long ib=0;
        totp=0.0;
        for (unsigned long i=0; i<ops.autobins; i++) 
        {
            q[0]=-ops.autoinitx+i*dx;
            for (unsigned long j=0; j<ops.autobins; j++) 
            {
                q[1]=-ops.autoinity+j*dy;
                totp+=autov[ib]; ib++;
                if (pchoice<=totp) { i=ops.autobins; break; }
            }
        }
    }

    flx.resize(ops.f_parsx.size());
    switch(ops.forcex)
    {
    case FSine:
        flx[0]=ops.f_parsx[0];  
        flx[1]=ops.units.u2i(UnitConv::UFrequency,ops.f_parsx[1]);
        pfx=fSinex;
        break;
    case FTable:
        pfx=fTablex; 
        tx.resize(ops.f_parsx.size()/2); ty.resize(ops.f_parsx.size()/2);
        for (int i=0; i<ops.f_parsx.size()/2; ++i)
        {
            tx[i]=ops.units.u2i(UnitConv::UTime,ops.f_parsx[2*i]);
            ty[i]=ops.f_parsx[2*i+1];
        }
        flx[0]=tx[0];
        flx[1]=tx[ops.f_parsx.size()/2-1];
        ftablex.set_table(tx,ty);
        break;
    case FEnvelope:
        pfx=fEnvx; 
        flx[0]=ops.f_parsx[0];  
        flx[1]=ops.units.u2i(UnitConv::UFrequency,ops.f_parsx[1]);
        tx.resize(ops.f_parsx.size()/2-2); ty.resize(ops.f_parsx.size()/2-2);
        for (int i=2; i<ops.f_parsx.size()/2; ++i)
        {
            tx[i-2]=ops.units.u2i(UnitConv::UTime,ops.f_parsx[2*i]);
            ty[i-2]=ops.f_parsx[2*i+1];
        }
        flx[2]=tx[0];
        flx[3]=tx[ops.f_parsx.size()/2-1-2];
        ftablex.set_table(tx,ty);
        break;
    }
    
    fly.resize(ops.f_parsy.size());
    switch(ops.forcey)
    {
        case FSine:
            fly[0]=ops.f_parsy[0];  
            fly[1]=ops.units.u2i(UnitConv::UFrequency,ops.f_parsy[1]);
            pfy=fSiney;
            break;
        case FTable:
            pfy=fTabley; 
            tx.resize(ops.f_parsy.size()/2); ty.resize(ops.f_parsy.size()/2);
            for (int i=0; i<ops.f_parsy.size()/2; ++i)
            {
                tx[i]=ops.units.u2i(UnitConv::UTime,ops.f_parsy[2*i]);
                ty[i]=ops.f_parsy[2*i+1];
            }
            fly[0]=tx[0];
            fly[1]=tx[ops.f_parsy.size()/2-1];
            ftabley.set_table(tx,ty);
            break;
        case FEnvelope:
            pfy=fEnvy; 
            fly[0]=ops.f_parsy[0];  
            fly[1]=ops.units.u2i(UnitConv::UFrequency,ops.f_parsy[1]);
            tx.resize(ops.f_parsy.size()/2-2); ty.resize(ops.f_parsy.size()/2-2);
            for (int i=2; i<ops.f_parsy.size()/2; ++i)
            {
                tx[i-2]=ops.units.u2i(UnitConv::UTime,ops.f_parsy[2*i]);
                ty[i-2]=ops.f_parsy[2*i+1];
            }
            fly[2]=tx[0];
            fly[3]=tx[ops.f_parsy.size()/2-1-2];
            ftabley.set_table(tx,ty);
            break;
    }

    
    std::ofstream tout;
    if (ops.f_print_traj)  tout.open((ops.prefix+std::string(".traj")).c_str());
    tout.precision(6);
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
    
    std::valarray<std::valarray<double> > pR(std::valarray<double>(n),2), npR(std::valarray<double>(n),2); 
    std::valarray<double> rnd(n);
    unsigned long i, j;
    
    for (k=0; k<2; k++)
    {
        for (i=0; i<n; ++i) rnd[i]=rgauss();
        //init 'noise' terms
        pR[k]=0.0;
        for (i=0; i<n; ++i)
            for (j=0; j<n; ++j)
                pR[k][i]+=rC(i,j)*rnd[j];
    }
    
    unsigned long is=0;
    
    
    double dvdq;
    double ecns, eth=0;
    std::cerr<<"starting\n";
    std::valarray<double> p(2);
    for (is=0; is< ops.nstep; ++is)
    {
        double epot, ekin, etz;
        
        //thermostat
        for (k=0; k<2; k++)
        {
            eth+=pR[k][0]*pR[k][0]/2.;
            for (i=0; i<n; ++i) rnd[i]=rgauss(); npR[k]=0.0;
            for (i=0; i<n; ++i)
                for (j=0; j<n; ++j)
                    npR[k][i]+=T(i,j)*pR[k][j]+S(i,j)*rnd[j];
            pR[k]=npR[k];
            eth-=pR[k][0]*pR[k][0]/2.;
        }
        
        for (k=0; k<2; k++) p[k]=pR[k][0];
        epot=verlet(q,p,ops.dt);
        for (k=0; k<2; k++) pR[k][0]=p[k];
        
        if (ops.f_alpha>0 || ops.f_beta>0) {
        if (ops.forcex!=FNone)
        {
            double fnow=pfx(is*ops.dt);
            pR[0][0]+=fnow*(ops.f_alpha+ops.f_beta*fnow)*ops.dt;
        }
        if (ops.forcey!=FNone)
        {
            double fnow=pfy(is*ops.dt);
            pR[1][0]+=fnow*(ops.f_alpha+ops.f_beta*fnow)*ops.dt;
        }
        }
        
        //std::cerr<<" TEST " <<q<<" "<<eTable(q)<<" "<<aTable(q)<<"\n";
        for (k=0; k<2; k++)
        {
            eth+=pR[k][0]*pR[k][0]/2.;
            for (i=0; i<n; ++i) rnd[i]=rgauss(); npR[k]=0.0;
            for (i=0; i<n; ++i)
                for (j=0; j<n; ++j)
                    npR[k][i]+=T(i,j)*pR[k][j]+S(i,j)*rnd[j];
            pR[k]=npR[k];
            eth-=pR[k][0]*pR[k][0]/2.;
        }
    
        ekin=ops.units.i2u(UnitConv::UEnergy,(pR[0][0]*pR[0][0]+pR[1][0]*pR[1][0])*0.5);
        epot=ops.units.i2u(UnitConv::UEnergy,epot);
        etz=ekin+epot;
        ecns=etz+ops.units.i2u(UnitConv::UEnergy,eth);
        for (k=0; k<2; k++) for (i=1; i<n; ++i) etz+=ops.units.i2u(UnitConv::UEnergy,pR[k][i]*pR[k][i]/2.);
        
        if (ops.f_print_traj && is%ops.stride==0 && is>ops.drop)
        { 
            tout<<setw(12)<<std::scientific<<is*ops.units.i2u(UnitConv::UTime,ops.dt)<<" "<<
                    setw(12)<<ops.units.i2u(UnitConv::ULength,q[0]/sqrt(ops.mass))<<" "<<
                    setw(12)<<ops.units.i2u(UnitConv::ULength,q[1]/sqrt(ops.mass))<<" "<<
                    setw(12)<<ops.units.i2u(UnitConv::ULength,pR[0][0]*sqrt(ops.mass))*ops.units.i2u(UnitConv::UMass,1.)/ops.units.i2u(UnitConv::UTime,1.)<<" "<<
                    setw(12)<<ops.units.i2u(UnitConv::ULength,pR[1][0]*sqrt(ops.mass))*ops.units.i2u(UnitConv::UMass,1.)/ops.units.i2u(UnitConv::UTime,1.)<<" "<<
                    setw(12)<<epot<<" " <<
                    setw(12)<< ekin<<" "<<
                    setw(12)<< ecns<<" ";
            tout<<std::endl;
        }
    }
    if (ops.f_print_traj) tout.close();
} 

