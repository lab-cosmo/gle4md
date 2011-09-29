#define TB_FFTAC 1
 
#include "color.hpp"
#include "clparser.hpp"
#include "tools-autocorr.hpp"
#include "conv-units.hpp"
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

namespace toolbox {
    //specialization for doubles (just to have it non-inlined for benchmarking)
    template<>
    inline double RndGaussian<double, MTRndUniform >::extract()
    {
        if (this->pstate.fg2ready) 
        {
            this->pstate.fg2ready=false;
            return this->pstate.gauss2*this->pstate.pars.sigma+this->pstate.pars.mean;
        }

        double x1,x2,w;

        do
        {
            x1=2.0*rugen()-1;
            x2=2.0*rugen()-1;
            w=x1*x1+x2*x2;
        } while (w>=1.0);

        w=sqrt(-2.0*log(w)/w);
        this->pstate.fg2ready=true;
        this->pstate.gauss2=x2*w;

        return x1*w*this->pstate.pars.sigma+this->pstate.pars.mean;
    }    
};

using namespace std;
using namespace toolbox;
class NDOps;

class Potential { 
public:
    friend std::istream& operator>> (std::istream& istr, NDOps& oo);
    enum PType {PHarmonic, PChain, PDChain, PUniform, PFree, PTable};
private:
    //harmonic potential
    FMatrix<double> H;
    void harm(const std::valarray<double>& x, double &e, std::valarray<double>& f);
    void chain(const std::valarray<double>& x, double &e, std::valarray<double>& f);
    void diagonal(const std::valarray<double>& x, double &e, std::valarray<double>& f);
    
    PType ptype;
    std::vector<double> pars;
    unsigned long xsize;
public:
    void set_pot(unsigned long nsz, const PType& ntyp, const std::vector<double>& npar);
    void get_pot(unsigned long& nsz, PType& ntyp, std::vector<double>& npar) const { nsz=xsize; ntyp=ptype, npar=pars; }
    void ef(const std::valarray<double>& x, double &e, std::valarray<double>& f)
    {
        switch(ptype)
        {
            case PHarmonic:
                harm(x,e,f);
                break;
            case PChain:
                chain(x,e,f);
                break;
            case PUniform:
                diagonal(x,e,f);
                break;
            case PDChain:
                diagonal(x,e,f);
                break;
            case PFree:
                f.resize(xsize); f=0.; e=0.;
                break;
        }
    }
    void verlet(std::valarray<double>& x, std::valarray<double>& v,const double& dt , const std::valarray<double>& m,  double &e, std::valarray<double>& f);
};

void Potential::verlet(std::valarray<double>& x, std::valarray<double>& v, const double& dt, const std::valarray<double>& m, double &e, std::valarray<double>& f)
{
    static std::valarray<double> olda;
    if (olda.size()==0)
    {
        olda.resize(xsize);
        ef(x,e,f);  for (int i=0; i<xsize; ++i) olda[i]=f[i]/m[i];
    }
    //std::cerr<<H<<x<<olda<<" data\n";
    v+=0.5*dt*olda;
    x+=dt*v;
    ef(x,e,f);  for (int i=0; i<xsize; ++i) olda[i]=f[i]/m[i];
    v+=0.5*dt*olda;
}

void Potential::harm(const std::valarray<double>& x, double &e, std::valarray<double>& f)
{
    f.resize(xsize); f=0.;
    for (unsigned int i=0; i<xsize;++i)
        for (unsigned int j=0; j<xsize;++j) 
            f[i]+=H(i,j)*x[j];
    e=0.;
    for (unsigned int i=0; i<xsize;++i) e+=f[i]*x[i];
    f*=-1.; e*=0.5;
}

void Potential::chain(const std::valarray<double>& x, double &e, std::valarray<double>& f)
{
    f.resize(xsize); f=0.;
    for (unsigned int i=1; i<xsize-1;++i)
        f[i]=H(0,0)*x[i]-H(0,1)*(x[i+1]+x[i-1]);
    f[0]=H(0,0)*x[0]-H(0,1)*(x[1]+x[xsize-1]);
    f[xsize-1]=H(0,0)*x[xsize-1]-H(0,1)*(x[0]+x[xsize-2]);
    e=0.;
    for (unsigned int i=0; i<xsize;++i) e+=f[i]*x[i];
    f*=-1.; e*=0.5;
}

void Potential::diagonal(const std::valarray<double>& x, double &e, std::valarray<double>& f)
{
    f.resize(xsize); f=0.;
    for (unsigned int i=0; i<xsize;++i)
        f[i]+=H(0,i)*x[i];
    e=0.;
    for (unsigned int i=0; i<xsize;++i) e+=f[i]*x[i];
    f*=-1.; e*=0.5;
}

void Potential::set_pot(unsigned long nsz, const PType& ntyp, const std::vector<double>& npar)
{
    xsize=nsz; ptype=ntyp; pars=npar;
    switch(ntyp)
    {
    case PHarmonic:
        if (npar.size()!=nsz*nsz) ERROR("Parameters don't fit with size");
        H.resize(nsz,nsz); H*=0.;
        for (unsigned int i=0; i<nsz;++i)
            for (unsigned int j=0; j<nsz;++j) 
                H(i,j)=npar[i*nsz+j];
        break;
    case PChain:
        if (npar.size()!=2) ERROR("Wrong parameters for linear chain with PBC");
        H.resize(1,2); H(0,0)=npar[0]; H(0,1)=npar[1]; 
        break;
    case PDChain:
        if (npar.size()!=2) ERROR("Wrong parameters for linear chain with PBC");
        H.resize(1,nsz); 
        for (unsigned int i=0; i<nsz;++i) H(0,i)=npar[0]-2*npar[1]*cos(2*toolbox::constant::pi*double(i)/double(nsz));
        break;
    case PUniform:
        if (npar.size()!=2) ERROR("Wrong parameters for uniform");
        H.resize(1,nsz); for (unsigned int i=0; i<nsz;++i) H(0,i)=npar[0]+(npar[1]-npar[0])*i/double(nsz-1.); 
        break;
    case PFree:
        if (npar.size()!=0) ERROR("Wrong parameters for free particle");
        H.resize(0,0); 
        break;
    default:
        ERROR("Unsupported potential code");
    }
}

class NDOps 
{
public:
    FMatrix<double> A, C;
    UnitConv units;
    std::string prefix;
    bool f_print_traj, f_init; 
    long seed;
    double dt;
    unsigned long nstep, stride;
    unsigned long maxlag;
    unsigned long nseries;
    unsigned long drop;
    
    Potential pot;
    std::valarray<double> mass;
    std::valarray<double> rx;
    std::valarray<double> rv;
};

std::ostream& operator<< (std::ostream& ostr, const Potential& oo)
{
    Potential::PType optype; std::vector<double> opars; unsigned long oxsize;
    oo.get_pot(oxsize, optype, opars);
    std::string spot;
    switch (optype)
    {
        case Potential::PHarmonic:
            spot="harmonic";
            break;
        case Potential::PDChain:
            spot="d-chain";
            break;
        case Potential::PChain:
            spot="chain";
            break;
        case Potential::PUniform:
            spot="uniform";
            break;
        case Potential::PFree:
            spot="free";
            break;
        default:
            ERROR("Unsupported potential type.");
    }
    
    IOMap iom;
    iom.insert(spot,"type"); 
    iom.insert(oxsize,"size");
    iom.insert(opars,"pars");
    
    return ostr<<iom;
}
std::istream& operator>> (std::istream& istr, Potential& oo)
{
    IOMap iom;
    Potential::PType optype; std::vector<double> opars; unsigned long oxsize;
    std::string spot;
    iom.insert(spot,"type"); 
    iom.insert(oxsize,"size");
    iom.insert(opars,"pars");
    std::cerr<<"READING POT\n";
    istr>>iom;
    std::cerr<<"READ POT\n";
    if (spot=="harmonic") optype=Potential::PHarmonic;
    else if (spot=="chain") optype=Potential::PChain;
    else if (spot=="d-chain") optype=Potential::PDChain;
    else if (spot=="uniform") optype=Potential::PUniform;
    else if (spot=="free") optype=Potential::PFree;
    else  ERROR("Unsupported potential type.");
    
    oo.set_pot(oxsize, optype, opars);
    return istr;
}

namespace toolbox {
    __MK_IT_IOFIELD(Potential);
    __MK_IT_IOFIELD(NDOps);
};

std::ostream& operator<< (std::ostream& ostr, const NDOps& oo)
{
    IOMap iom;
    iom.insert(oo.A,"a_matrix"); 
    iom.insert(oo.C,"c_matrix");
    iom.insert(oo.units, "units");
    iom.insert(oo.dt,"timestep");
    iom.insert(oo.drop,"discard");
    iom.insert(oo.pot,"potential");
    iom.insert(oo.mass,"mass");
    iom.insert(oo.nstep,"steps");
    iom.insert(oo.stride,"stride");
    iom.insert(oo.maxlag,"max_lag");
    iom.insert(oo.nseries,"max_data");
    iom.insert(oo.prefix,"prefix");
    iom.insert(oo.seed,"seed");
    iom.insert(oo.f_print_traj,"print_traj");
    iom.insert(oo.f_init,"init_s");
    iom.insert(oo.rx,"restart-x");
    iom.insert(oo.rv,"restart-v");
    iom.insert(oo.mass,"masses");
    ostr<<iom;
    
    return ostr;
}

std::istream& operator>> (std::istream& istr, NDOps& oo)
{
    IOMap iom;
    
    FMatrix<double> Z(0,0); double ascale, cscale,temp;
    iom.insert(ascale,"ascale",1.);
    iom.insert(cscale,"cscale",1.);
    iom.insert(oo.A,"a_matrix");
    iom.insert(oo.C,"c_matrix");
    iom.insert(oo.dt,"timestep");
    iom.insert(oo.nstep,"steps");
    iom.insert(oo.stride,"stride",(unsigned long)1);
    iom.insert(temp,"temperature",0.);
    iom.insert(oo.maxlag,"max_lag");
    iom.insert(oo.nseries,"max_data");
    iom.insert(oo.drop,"discard",(unsigned long) 0);
    iom.insert(oo.prefix,"prefix",std::string("oned"));
    iom.insert(oo.seed,"seed",(long)1234321);
    iom.insert(oo.f_print_traj,"print_traj",true);
    iom.insert(oo.f_init,"init_s",true);
    iom.insert(oo.units, "units",UnitConv());
    iom.insert(oo.pot, "potential");
    iom.insert(oo.rx,"restart-x",std::valarray<double>(0));
    iom.insert(oo.rv,"restart-v",std::valarray<double>(0));
    iom.insert(oo.mass,"masses");
    
    istr>>iom;
    oo.units.ITime=UnitConv::AtomicTime; oo.units.IEnergy=UnitConv::Hartree;
    oo.units.IFrequency=UnitConv::AtomicFrequency; oo.units.IMass=UnitConv::AtomicMass;
    
    if (oo.C.rows()==0)
    {
        if (temp<=0.) ERROR("C matrix or temperature should be given"); 
        //classical matrix for the target temperature
        oo.C.resize(oo.A.rows(),oo.A.cols());
        for (int i=0; i<oo.A.rows(); ++i) oo.C(i,i)=temp;
    }
    if (oo.mass.size()==0) { oo.mass.resize(oo.pot.xsize); oo.mass=1./oo.units.u2i(UnitConv::UMass,1.);; }
    else if (oo.mass.size()==1) { double tmpm=oo.mass[0]; oo.mass.resize(oo.pot.xsize); oo.mass=tmpm; }
    
    oo.A*=ascale; oo.C*=cscale;
    
    //we convert the options in internal units
    oo.C*=oo.units.u2i(UnitConv::UEnergy,1.); 
    oo.A*=1./oo.units.u2i(UnitConv::UTime,1.); 
    oo.dt=oo.units.u2i(UnitConv::UTime,oo.dt); 
    for (int i=0; i<oo.mass.size();++i) oo.mass[i]=oo.units.u2i(UnitConv::UMass,oo.mass[i]);
    for (int i=0; i<oo.rx.size();++i) oo.rx[i]=oo.units.u2i(UnitConv::ULength,oo.rx[i]);
    for (int i=0; i<oo.rv.size();++i) oo.rv[i]=oo.units.u2i(UnitConv::ULength,oo.rv[i])/oo.units.u2i(UnitConv::UTime,1.0);
    std::cerr<<"conversion "<<oo.units.u2i(UnitConv::UMass,1.)/pow(oo.units.u2i(UnitConv::UTime,1.),2.)<<"\n";
    
    //std::cerr<<oo.pot.H(0,0)<<" >>> ";
    oo.pot.H*=oo.units.u2i(UnitConv::UMass,1.)/pow(oo.units.u2i(UnitConv::UTime,1.),2.);
    /*
    std::cerr<<oo.pot.H(0,0)<<" H ELEMENT \n";
    std::cerr<<"should be 1/100!!!"<<1./oo.units.i2u(UnitConv::UTime,1./(sqrt(oo.pot.H(0,0)/oo.mass[0])))<<"\n";
    std::cerr<<"expected frequencies between "<<sqrt((oo.pot.H(0,0)-2*oo.pot.H(0,1))/oo.mass[0])*1./oo.units.i2u(UnitConv::UTime,1.)
            <<" and "<<sqrt((oo.pot.H(0,0)+2*oo.pot.H(0,1))/oo.mass[0])*1./oo.units.i2u(UnitConv::UTime,1.)<<
            "\n";
    std::cerr<<"uniform frequencies between "<<sqrt(oo.pot.H(0,0)/oo.mass[0])*1./oo.units.i2u(UnitConv::UTime,1.)<<
            " and "<<sqrt(oo.pot.H(0,oo.pot.H.cols()-1)/oo.mass[0])*1./oo.units.i2u(UnitConv::UTime,1.)<<"\n";
    */
    return istr;
}


int main(int argc, char **argv)
{
    RndGaussian<double, MTRndUniform> rgauss(RGPars<double>(0.,1.));
    unsigned long i, j, k;
    FMatrix<double>  S, T, rC;
    NDOps ops;
    std::cin>>ops;
    std::cerr<<"SEED "<<ops.seed<<"\n";
    rgauss.RUGenerator().seed(ops.seed);
    
    unsigned long n=ops.A.rows();
    StabCholesky(ops.C,rC);
    
    /*!*******Initialization */
    unsigned long sz; std::vector<double> ppar; Potential::PType ptype;
    ops.pot.get_pot(sz,ptype,ppar);
    std::cerr<<sz<<" SIZE\n";
    std::valarray<double> x(sz);
    FMatrix<double>gp(sz,n); gp*=0.; 

    std::valarray<double> rg(n);
    std::valarray<double> smass(sz);
    
    for (i=0; i<sz; ++i)
    {
        std::cerr<<i<<":"<<ops.mass[i]<<"  ";
        smass[i]=sqrt(ops.mass[i]);
        for (j=0; j<n; ++j) rg[j]=rgauss();
        if (ops.f_init) for (j=0; j<n; ++j) for (k=0; k<n; ++k) gp(i,j)+=rC(j,k)*rg[k]*smass[i];
    }
    
    x=0.; if (ops.rx.size()==sz) x=ops.rx;
    if (ops.rv.size()==sz) for (i=0; i<sz; ++i) gp(i,0)=ops.rv[i]*smass[i];
    
    //std::ofstream acout((ops.prefix+std::string(".acf")).c_str());
    std::ofstream tout, vout, xout;
    std::cerr<<ops.f_print_traj<< " true?\n";
    std::cerr<<ops.prefix<< " opening\n";
    tout.open((ops.prefix+std::string(".traj")).c_str());
    if (ops.f_print_traj) 
    {
        xout.open((ops.prefix+std::string(".xtraj")).c_str());
        vout.open((ops.prefix+std::string(".vtraj")).c_str());
    }
    //acout.precision(6); 
    tout.precision(6); vout.precision(6); xout.precision(6);
    
    get_TS(ops.A, ops.C, ops.dt*0.5, T, S);
    
    
    
    std::cerr<<"################# MATRICES ###################\n";
    std::cerr<<"****************** A *************************\n"<<ops.A<<"\n";
    std::cerr<<"****************** C *************************\n"<<ops.C<<"\n";
    std::cerr<<"****************** T *************************\n"<<T<<"\n";
    std::cerr<<"****************** S *************************\n"<<S<<"\n";
    
    FMatrix<double> chk,tmp,lt,ls;
    transpose(T,lt); transpose(S,ls);
    mult(ops.C,lt,chk); mult(T,chk,lt); chk=lt; chk-=ops.C;
    mult(S,ls,tmp);
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
    
    
    //double dp, op;
    FMatrix<double> ngp(gp), vrg(gp);
    FMatrix<double> Tt,St; transpose(T,Tt); transpose(S,St);
    std::valarray<double> nv(sz);
    std::valarray<double> rf;
    double epot, ekin, etz, egle=0.,ecns; unsigned long is;
    for (is=0; is< ops.nstep; ++is)
    {
        //thermostat
        //npR=dp=0; 
        //op=pR[0]; //for markov-check!
        
        
        for (k=0;k<sz;++k) { gp(k,0)=nv[k]*smass[k]; egle+=gp(k,0)*gp(k,0); }

        for (k=0;k<sz;++k) for (i=0; i<n; ++i) vrg(k,i)=rgauss();
        mult(gp,Tt,ngp);
        mult(vrg,St,gp);
        gp+=ngp;
        
        for (k=0;k<sz;++k) { nv[k]=gp(k,0)/smass[k]; egle-=gp(k,0)*gp(k,0); }
        
        /*
        ngp=0.;
        for (k=0;k<sz;++k)
        {
            for (i=0; i<n; ++i) rg[i]=rgauss();
            
            for (i=0; i<n; ++i)
                for (j=0; j<n; ++j)
                    ngp(k,i)+=T(i,j)*gp(k,j)+S(i,j)*rg[j];
            nv[k]=ngp(k,0)/smass[k];
        }
        gp=ngp;
        */
        //dp+=npR[0]-pR[0];
        
        //pR=npR;
        //if (ops.pot==PSquareWell) sw_verlet(q,pR[0],ops.dt,pl[0]);
        //else verlet(q,pR[0],ops.dt);
        ops.pot.verlet(x,nv,ops.dt,ops.mass,epot,rf);
        //npR=0;
        /*
        ngp=0.;
        for (k=0;k<sz;++k)
        {
            gp(k,0)=nv[k]*smass[k];
            for (i=0; i<n; ++i) rg[i]=rgauss();
            for (i=0; i<n; ++i)
                for (j=0; j<n; ++j)
                    ngp(k,i)+=T(i,j)*gp(k,j)+S(i,j)*rg[j];
            nv[k]=ngp(k,0)/smass[k];
            
        }
        gp=ngp;
        */
        
        for (k=0;k<sz;++k) { gp(k,0)=nv[k]*smass[k]; egle+=gp(k,0)*gp(k,0); }
        
        for (k=0;k<sz;++k) for (i=0; i<n; ++i) vrg(k,i)=rgauss();
        mult(gp,Tt,ngp);
        mult(vrg,St,gp);
        gp+=ngp;
        
        for (k=0;k<sz;++k) { nv[k]=gp(k,0)/smass[k]; egle-=gp(k,0)*gp(k,0); }
        
        //dp+=npR[0]-pR[0];
        //pR=npR;
        ekin=etz=0.;
        for (k=0;k<sz;++k)
        {
            ekin+=gp(k,0)*gp(k,0);
            for (i=0; i<n; ++i) etz+=gp(k,i)*gp(k,i);
        }
        ecns=ops.units.i2u(UnitConv::UEnergy,epot+0.5*(ekin+egle));
        ekin=ops.units.i2u(UnitConv::UEnergy,ekin*0.5); etz=ops.units.i2u(UnitConv::UEnergy,etz*0.5);
        epot=ops.units.i2u(UnitConv::UEnergy,epot); 
        etz+=epot;
        
        if (is>ops.drop)
        {
            if (is%ops.stride==0)
            { 
                tout<<setw(12)<<is*ops.units.i2u(UnitConv::UTime,ops.dt)<<" "<<
                        setw(12)<<epot<<" "<<setw(12)<< ekin<<" "
                        <<setw(12)<<ecns<<" "<<setw(12)<< etz<<std::endl;
                if (ops.f_print_traj)
                {
                    for (i=0; i<sz; ++i) xout<<setw(6)<<ops.units.i2u(UnitConv::ULength,x[i])<<" ";
                    xout<<std::endl;
                    for (i=0; i<sz; ++i) vout<<setw(6)<<ops.units.i2u(UnitConv::ULength,nv[i])/ops.units.i2u(UnitConv::UTime,1.)<<" ";
                    vout<<std::endl;
                }
            }
        }
    }
    tout.close();  if (ops.f_print_traj) { vout.close(); xout.close(); }
}

