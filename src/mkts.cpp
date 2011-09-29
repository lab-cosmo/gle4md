#include "color.hpp"
#include "clparser.hpp"
#include "ioparser.hpp"
#include "matrix-io.hpp"
#include "minsearch.hpp"
#include <fstream>
using namespace std;
using namespace toolbox;

int main(int argc, char **argv)
{
    std::valarray<ThermoOptions> thermo;
    IField<std::valarray<ThermoOptions > > ifthermo(thermo,"");
    std::cin>>ifthermo;
    
    unsigned long n=thermo.size()*2+1;
    FMatrix<double> A(n,n), B(n,n), T(n,n), S(n,n);
    /* builds A and B matrices only for the momentum-stochastic part */
    
    A*=0;
    for (unsigned long i=0; i<thermo.size(); ++i)
    {
        double tt=1./(thermo[i].tauf*thermo[i].tauf)+1./(thermo[i].tauw*thermo[i].tauw);
        A(0,2*i+1)=-(A(2*i+1,0)=sqrt(thermo[i].tauf*tt/thermo[i].taut));
        A(2*i+1,2*i+1)=A(2*i+2,2*i+2)=1./thermo[i].tauf;
        A(2*i+1,2*i+2)=-(A(2*i+2,2*i+1)=1./thermo[i].tauw);
        B(2*i+1,2*i+1)=B(2*i+2,2*i+2)=sqrt(2./(thops[i].tauf*thops[i].beta));
    }
    
    build_A(chi2min.thermo,A);
    build_B(chi2min.thermo,B);
    
    return 1;
}