#include "color.hpp"
#include "glefit.hpp"
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
            << " USAGE: gle-scaler  -a a-file [(-b b-file | -c c-file | -d d-file)] [-sa scale] \n" 
            << "     [-sc scale] [-oa|-ob|-oc|-od] [ -oraw | -ostd ]                            \n"
            << " performs appropriate scaling of gle matrices                                   \n";
}

enum OMFormat { Standard, Raw }; 
void omatrix(std::ostream& ostr, const FMatrix<double>& mat, OMFormat fmt)
{
    ostr.setf(std::ios::scientific);
    ostr.precision(20);
    switch (fmt)
    {
    case Standard: 
        ostr<<mat; break;
    case Raw: 
        for (int i=0; i<mat.rows(); i++) { for (int j=0; j<mat.rows(); j++) ostr<<std::setw(20)<<mat(i,j)<<" ";  ostr<<std::endl; } break;
    }
}

int main(int argc, char **argv)
{
    CLParser clp(argc, argv);
    std::string amat, bmat, cmat, dmat, oraw, ostd;
    bool foa, fob, foc, fod;
    double shifta, shiftc, shiftkh;
    bool fok=clp.getoption(amat,"a") &&
            clp.getoption(bmat,"b",std::string("")) &&
            clp.getoption(cmat,"c",std::string("")) &&
            clp.getoption(dmat,"d",std::string("")) &&
            clp.getoption(oraw,"oraw",std::string("")) &&
            clp.getoption(ostd,"ostd",std::string("")) &&
            clp.getoption(foa,"oa",false) &&
            clp.getoption(fob,"ob",false) &&
            clp.getoption(foc,"oc",false) &&
            clp.getoption(fod,"od",false) &&
            clp.getoption(shiftkh,"skh",1.) &&
            clp.getoption(shifta,"sa",1.) &&
            clp.getoption(shiftc,"sc",1.);

    GLEABC abc;
    FMatrix<double> iA, iC, iD;
    std::ifstream ifile;
    
    IField<FMatrix<double> > ifa(iA,"");
    ifile.open(amat.c_str());
    if (ifile.bad())
        ERROR("Unable to open file for A!");
    ifile >> ifa;
    ifile.close();
    unsigned long n=iA.rows()-1; iA*=shifta;
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
        if (iC.rows()!=iC.cols() || iC.rows()!=n+1)
            ERROR("Invalid size or format of matrix C.");
        iC*=shiftc;
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
        if (iB.rows()!=iB.cols() || iB.rows()!=n+1)
            ERROR("Invalid size or format of matrix B.");
        transpose(iB,iBt); mult(iB,iBt,iBBt);
        iBBt*=shifta*shiftc;
        abc.set_BBT(iBBt);
    }
    else if (dmat!="")
    {
        //read B and get C
        ifile.open(dmat.c_str());
        if (ifile.bad())
            ERROR("Unable to open file for D!");
        IField<FMatrix<double> > ifd(iD,"");
        ifile >> ifd;
        if (iD.rows()!=iD.cols() || iD.rows()!=n+1)
            ERROR("Invalid size or format of matrix B.");
        iD*=shifta*shiftc;
        abc.set_BBT(iD);
        abc.get_C(iD);
    }
    else
    {
        //unit matrix
        iC.resize(n+1,n+1); iC*=0.;
        for (unsigned long i=0; i<n+1; ++i) iC(i,i)=shiftc;
        abc.set_C(iC);
    }
    
    if (shiftkh!=1.)
    {
        abc.get_A(iA); abc.get_BBT(iD); 
        iA(0,0)*=shiftkh; iD(0,0)*=shiftkh;
        double sskh=sqrt(shiftkh);
        for (unsigned long i=1; i<n+1; ++i) { iA(0,i)*=sskh; iA(i,0)*=sskh; iD(0,i)*=sskh; iD(i,0)*=sskh; }
        abc.set_A(iA); abc.set_BBT(iD);
    }
    
    FMatrix<double> oA, oB, oC, oD;
    abc.get_A(oA); abc.get_BBT(oD); abc.get_C(oC); 
    Cholesky(oD,oB);
    
    std::ofstream fmat;
    if (oraw!=std::string(""))
    {
        if (foa) { fmat.open((oraw+std::string("-A")).c_str()); omatrix(fmat,oA,Raw); fmat.close(); }
        if (fob) { fmat.open((oraw+std::string("-B")).c_str()); omatrix(fmat,oB,Raw); fmat.close(); }
        if (foc) { fmat.open((oraw+std::string("-C")).c_str()); omatrix(fmat,oC,Raw); fmat.close(); }
        if (fod) { fmat.open((oraw+std::string("-D")).c_str()); omatrix(fmat,oD,Raw); fmat.close(); }
    }
    if (ostd!=std::string(""))
    {
        if (foa) { fmat.open((ostd+std::string(".A")).c_str()); omatrix(fmat,oA,Standard); fmat.close(); }
        if (fob) { fmat.open((ostd+std::string(".B")).c_str()); omatrix(fmat,oB,Standard); fmat.close(); }
        if (foc) { fmat.open((ostd+std::string(".C")).c_str()); omatrix(fmat,oC,Standard); fmat.close(); }
        if (fod) { fmat.open((ostd+std::string(".D")).c_str()); omatrix(fmat,oD,Standard); fmat.close(); }
    }
}
