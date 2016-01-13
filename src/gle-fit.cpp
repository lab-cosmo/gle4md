#include "glefit.hpp"
#include "clparser.hpp"
#include "matrix-io.hpp"
#include <fstream>
using namespace std;
using namespace toolbox;

/*******************************************************************
 in matrices and vectors the order of variables is 
 q p r1 r2 .... rn
*******************************************************************/

void banner() 
{
    std::cerr
            << " USAGE: gle-fit < input                                                        \n"
            << " compute best-fit parameters for the autocorrelation time of a set of           \n"
            << " thermostats, and for fluctuations as a function of w.                          \n"
            << "                                                                                \n";
}
class GLEFRestart {
public:
    std::valarray<double> pars;
    FMatrix<double> A, C, B, D;
    GLEFRestart(): pars(0), A(0,0), C(0,0) {}
};
std::ostream& operator<<(std::ostream& os, const GLEFRestart& gr)
{
    IOMap iom;
    iom.insert(gr.pars, "parameters");
    iom.insert(gr.A, "a_matrix");
    iom.insert(gr.C, "c_matrix");
    iom.insert(gr.B, "b_matrix");
    iom.insert(gr.D, "d_matrix");
    os << " {\n"<<iom<<"\n}\n";
    return os;
}  

std::istream& operator>>(std::istream& is, GLEFRestart& gr)
{
    IOMap iom;
    iom.insert(gr.A, "a_matrix", FMatrix<double>(0));
    iom.insert(gr.C, "c_matrix", FMatrix<double>(0));
    iom.insert(gr.B, "b_matrix", FMatrix<double>(0));
    iom.insert(gr.D, "d_matrix", FMatrix<double>(0));
    iom.insert(gr.pars, "parameters", std::valarray<double>(0));
    is>>iom;
    return is;
}

namespace toolbox{
    __MK_IT_IOFIELD(GLEFRestart);
}
int main(int argc, char **argv)
{
    /*************************************************************************
     The code tries to find the best fit to a required target autocorrelation
     time and, if required, fluctuation-dissipation theorem (in frequency space).
     Since this starts getting complex, I rather use an input file than
     command-line options.
     Therefore, we have separate options for initial conditions, fitting method,
     relative weight of response and FDT constraints, etc. 
     We expose both 'simple' and 'advanced' fitting options, so that the code is
     both simple to use and to customize.
    **************************************************************************/
    std::ofstream fmat;
    fmat.setf(std::ios::scientific);
    fmat.precision(9);
    fmat.width(16);
    IOMap iom;
    GLEFError ercls;
    GLEFFitOptions ofit;
    GLEFParOptions opar;
    GLEFSearchOptions osrc;
    GLEFRestart gres;
    std::string prefix;
    
    //std::valarray<double> ;
    //int n_res;
    
    //reads input
    iom.insert(ofit,"fitting");
    iom.insert(opar,"parameters");
    iom.insert(osrc,"search");
    iom.insert(prefix,"prefix", std::string("gle"));

    iom.insert(gres, "restart", GLEFRestart());
    std::cerr<<"gonna read"<<std::endl;
    std::cin>>iom;
    std::cerr<<"fkd read"<<std::endl;
    std::cerr<<"pars"<<opar<<std::endl;
    fmat.open((prefix+std::string(".tgt")).c_str());
    fmat<<ofit;
    fmat.close();

    ercls.slog=new std::ofstream((prefix+std::string(".log")).c_str());
    ercls.spars=new std::ofstream((prefix+std::string(".spars")).c_str());
    ercls.seva=new std::ofstream((prefix+std::string(".seva")).c_str());
    
    ercls.set_ops(ofit,opar);
    
    //implements restart
    bool fres_ok=false;
    std::valarray<double> ip; init_pars(opar,ip); ercls.set_vars(ip);
    double ttemp=ercls.C(0,0);
    std::cerr<<"READING RESTART\n";
    if (gres.pars.size()!=0) 
    {  
        std::cerr<<"PARAMETERS GIVEN\n";
        //parameters is the preferred restart method
        if (gres.pars.size()!=npars(opar))
            ERROR("Parameter size "<<gres.pars.size()<<" mismatches with required n. "<<npars(opar));
        ercls.set_vars(gres.pars);
        fres_ok=true;
    }
    if (gres.A.rows()!=0 && ! fres_ok)
    {
        if (gres.C.rows()==0) 
        {
            if (gres.D.rows()==0)
            {
                gres.C.resize(gres.A.rows(),gres.A.cols());
                for (int i=0; i<gres.A.rows(); ++i) gres.C(i,i)=ttemp;
            }
            else
            {
                GLEABC myabc;
                myabc.set_A(gres.A); myabc.set_BBT(gres.D); 
                myabc.get_C(gres.C);
            }
        }
            
        std::cerr<<"MATRIX GIVEN\n";
        if (gres.A.rows()!=gres.A.cols() || gres.A.rows()!=gres.C.rows() || gres.C.rows()!=gres.C.cols())
            ERROR("Wrong dimension for restart matrices.");
        int mn=(ercls.A.rows()<gres.A.rows()?ercls.A.rows():gres.A.rows());
        for (int i=0; i<mn; ++i)      // a certain freedom for resizing matrices
            for (int j=0; j<mn; ++j)
        { ercls.A(i,j)=gres.A(i,j); ercls.C(i,j)=gres.C(i,j); }
        /*if (ercls.C(0,0)!=ttemp)
        {
            ercls.A*=ttemp/ercls.C(0,0); ercls.C*=ttemp/ercls.C(0,0); 
        }*/
        std::cerr<<ercls.A<<"and"<<ercls.C<<"\n";
        ercls.AC2pars();   //backfits pars to given matrices!
        ercls.get_vars(ip);
        fres_ok=true;
    }
    if (!fres_ok) std::cerr<<"STARTING FROM SCRATCH\n";
    
    double ierr; ercls.get_vars(ip); ercls.get_value(ierr); 
    std::cerr<<" INITIAL ERROR: "<< ierr<<std::endl;
    std::cerr<<" INITIAL PARS: "<<ip<<"\n";
    ercls.set_vars(ip); ercls.get_value(ierr); 
    std::cerr<<" PARS ERROR: "<< ierr<<std::endl;
    /*if (fsimplex)
    {
        //simplex
        IterOptions<double,2> iops=IterOptions<double,2>(
                nsteps,
                fixarray<double,2>(1e-5, 1e-5),
                fixarray<double,2>(0.,0.),
                fixarray<double,2>(ichk_change, ichk_default));
        min_simplex(chi2min,make_simplex(pars,sap1,1+sap2),pars,ferr, iops);
    }
    else
    {*/
    double ferr;
    std::cerr<<osrc<<"\n"; 
    std::cerr<<opar<<"\n";
    AnnealingOptions ao;
    IterOptions<double,2> iops;
    switch(osrc.mode)
    {
        case Annealing:
            ao.temp_init=osrc.ti*ierr; 
            ao.temp_final=osrc.tf*ierr;
            ao.steps=osrc.steps; ao.mc_step=osrc.step; 
            
            ao.adapt=osrc.adapt_mult;
            ao.fpowell=true; ao.drnd=osrc.rand;
            sim_annealing(ercls,ip,ip,ferr,ao);
            break;
        case Simplex:
            iops=IterOptions<double,2>( 
                        osrc.steps,
                        fixarray<double,2>(osrc.tol, osrc.tol),fixarray<double,2>(0.,0.),
                        fixarray<double,2>(ichk_change, ichk_default));
            min_simplex(ercls,make_simplex(ip,osrc.step,1.),ip,ferr, iops);
            break;
        case Powell:
            std::valarray<std::valarray<double> > u;
            PowellOpts po; 
            po.fadjstep=(osrc.adapt_mult!=1.); po.tol=osrc.tol;
            po.linesearch.lstol=osrc.tol; po.drnd=osrc.rand; 
            po.linesearch.maxiter=osrc.nlinesrc; po.linesearch.istep=osrc.step;
            po.maxiter=osrc.steps;
            std::cerr<<"CALLING POWELL\n"; 
            min_powell(ercls,ip,ip,ferr,po,u);
    }
    std::cerr<<"final error is "<<ferr<<"\n";
    
    ercls.get_vars(ip);
    fmat.open((prefix+std::string(".pars")).c_str());
    fmat<<ip;
    fmat.close();
    std::cerr<<ercls.A.size()<<" "<<ercls.C.size()<<" SIZES\n";
    fmat.open((prefix+std::string(".A")).c_str());
    fmat<<ercls.A;
    fmat.close();
    fmat.open((prefix+std::string(".C")).c_str());
    FMatrix<double> Ct; transpose(ercls.C,Ct); 
    fmat<<0.5*(ercls.C+Ct);
    fmat.close();
    FMatrix<double> oB;
    ercls.abc.get_BBT(oB);
    fmat.open((prefix+std::string(".BBT")).c_str());
    fmat<<oB;
    fmat.close();
    
    
    *ercls.slog<<std::endl;
}
