#include "glefit.hpp"
#include "color.hpp"
#include "clparser.hpp"
#include "ioparser.hpp"
#include "matrix-io.hpp"
#include <fstream>
using namespace std;
using namespace toolbox;

std::string expn(double x);

void banner()
{
    std::cerr
            << " USAGE: par-test -a a-file [(-b b-file | -c c-file | -d d-file)]                \n"
            << "      -as astyle   -cs cstyle                                                   \n"
            << " simply tests the matrix parameterization by loading the matrix, converting to  \n"
            << " parametric form given by the style options and converting back to matrices.    \n";
}

int main(int argc, char **argv)
{
    CLParser clp(argc, argv);
    bool fhelp;
    std::string amat, bmat, cmat, dmat, astyle, cstyle;
    long np;
    bool fok=
            clp.getoption(amat,"a") &&
            clp.getoption(bmat,"b",std::string("")) &&
            clp.getoption(cmat,"c",std::string("")) &&
            clp.getoption(dmat,"d",std::string("")) &&
            clp.getoption(astyle,"as") &&
            clp.getoption(cstyle,"cs");

    if (!fok || fhelp) {banner(); return 0;}

    GLEABC abc;
    FMatrix<double> iA, iC;
    std::ifstream ifile;

    IField<FMatrix<double> > ifa(iA,"");
    ifile.open(amat.c_str());
    if (ifile.bad())
        ERROR("Unable to open file for A!");
    ifile >> ifa;
    ifile.close();
    unsigned long n=iA.rows()+1;
    abc.set_A(iA);

    //now we get somehow the C matrix
    if (cmat!="")
    {
        ifile.open(cmat.c_str());
        if (ifile.bad())
            ERROR("Unable to open file for C!");
        if (bmat!="") WARNING("Both B and C have been supplied. B will be ignored.");
        IField<FMatrix<double> > ifc(iC,"");
        ifile >> ifc;
        if (iC.rows()!=iC.cols() || iC.rows()!=n-1)
            ERROR("Invalid size or format of matrix C.");
        abc.set_C(iC);
    }
    else if (bmat!="")
    {
        //read B and get C
        ifile.open(bmat.c_str());
        if (ifile.bad())
            ERROR("Unable to open file for B!");
        FMatrix<double> iB, iBt, iBBt;
        IField<FMatrix<double> > ifb(iB,"");
        ifile >> ifb;
        if (iB.rows()!=iB.cols() || iB.rows()!=n-1)
            ERROR("Invalid size or format of matrix B.");
        transpose(iB,iBt); mult(iB,iBt,iBBt);
        abc.set_BBT(iBBt);
    }
    else if (dmat!="")
    {
        //read B and get C
        ifile.open(dmat.c_str());
        if (ifile.bad())
            ERROR("Unable to open file for D!");
        FMatrix<double> iD;
        IField<FMatrix<double> > ifd(iD,"");
        ifile >> ifd;
        if (iD.rows()!=iD.cols() || iD.rows()!=n-1)
            ERROR("Invalid size or format of matrix B.");
        abc.set_BBT(iD);
        abc.get_C(iD);
    }
    else
    {
        //unit matrix
        iC.resize(n-1,n-1); iC*=0.;
        for (unsigned long i=0; i<n-1; ++i) iC(i,i)=1.0;
        abc.set_C(iC);
    }


    FMatrix<double> iBBT; abc.get_BBT(iBBT); abc.get_C(iC);
    std::cout.precision(10);
    std::cout.setf(std::ios::left);
    std::cout.setf(std::ios::scientific);

    //computes analytics
    GLEFError ercls;
    GLEFParOptions opar;
    GLEFFitOptions ofit;

    opar.ns=iA.rows()-1;
    std::stringstream istr;
    istr.str(astyle); istr >> opar.pstyleA;
    istr.clear();
    istr.str(cstyle); istr >> opar.pstyleC;
    ercls.set_ops(ofit,opar);
    std::cerr<<opar;
    ercls.A=iA; ercls.C=iC;
    ercls.AC2pars();
    ercls.pars2AC();

    std::cout<<"# Input A matrix: \n";
    std::cout<< iA;
    std::cout<<"# Computed A matrix: \n";
    std::cout<<ercls.A;
    std::cout<<"# Input C matrix: \n";
    std::cout<< iC;
    std::cout<<"# Computed C matrix: \n";
    std::cout<<ercls.C;
}

