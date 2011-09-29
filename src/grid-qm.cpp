#include <iostream>
#include <fstream>
#include <valarray>
#include "clparser.hpp"
#include "ioparser.hpp"

extern "C" {
   void dstebz_(char * range, char* order, int * n, 
           double *vl, double* vu, int *il, int* iu, double *abstol,
            double* d, double* e, int* m, int*nsplit,
            double* w, int* iblock, int* isplit,
            double *work, int* iwork, int* info);
   void dstevx_(char * jobz, char * range, int * n, double* d, double* e, 
                double *vl, double* vu, int *il, int* iu, double *abstol,
                int* m, double* w, double* z, int *ldz, 
                double *work, int* iwork, int* ifail, int* info);
    
};

int main(int argc, char**argv)
{
    toolbox::CLParser clp(argc,argv);
    std::string gfile;
    int gsz, nev; double T,hbar,mass;
    std::string wff; int nwf; bool do_wf=false;
    bool fok=clp.getoption(gfile,"grid",std::string("grid")) &&
            clp.getoption(wff,"wf",std::string("")) &&
            clp.getoption(nev,"n",20) &&
            clp.getoption(hbar,"hbar",1.) &&
            clp.getoption(nwf,"nwf",10);
    
    std::ofstream wfstream;
    std::vector<std::valarray<double> > wfdata;
    if (wff!="") 
    {
        do_wf=true;
        wfstream.open(wff.c_str());
        if (wfstream.bad()) ERROR("Unable to open file for wavefunction output");
    }

    std::ifstream gstream(gfile.c_str());
    double dx;
    gstream >> dx;
    std::valarray<double> gv;
    std::cerr<<"dx "<<dx<<"  nev "<<nev<<"\n";
    toolbox::IField<std::valarray<double> > ifi(gv,"");
    gstream >> ifi;
    
    gsz=gv.size();
    std::cerr<<"grid size " << gsz<<"\n";
    //std::cerr<<"Running for w= " << f<<std::endl;
    std::valarray<double> d(gsz), e(gsz-1);
    
    //sets grid points for x and diagonal entries
    double dx2=1./(dx*dx);
    for(int i=0; i<gsz; ++i) d[i]=gv[i]+dx2*hbar*hbar;
    e=-0.5*dx2*hbar*hbar;
    
    char range='I', jobz='V'; 
    int n=gsz; double vl=0, vu=0.; 
    double abstol=0.; int m;
    int il=1, iu=nev, ldz=n; 
    std::valarray<double> w(gsz), work(5*n);
    std::valarray<double> z(gsz*nev);
    std::valarray<int> iwork(5*n), iblock(n), isplit(n);
    int info; std::valarray<int> ifail(n);
    int iwf=0;
    if (iu>nwf) iu=nwf;
    while(il<=n && iu<=nwf)
    {
        std::cerr<<"looking eigv between "<<il<<" "<<iu<<"\n";
        if (iu>n) iu=n; 
        dstevx_(&jobz, &range, &n, &d[0], &e[0], 
            &vl, &vu, &il, &iu, &abstol,
            &m, &w[0], &z[0], &ldz,
            &work[0], &iwork[0], &ifail[0], &info);
        if (info!=0) ERROR("dstevx failed with exit code "<<info);
        il=iu+1; iu=iu+nev;
        if (iu>nwf) iu=nwf;
        //now, w[i] contains the eigenvalue, while z contains the eigenvectors,
        //we might use to compute observables such as <V> and <K>
        int iz=0;
        for(int i=0; i<m; ++i)
        {
            double ak=0, av=0;
            for(int j=0; j<n; ++j) 
            { 
                av+=gv[j]*z[iz]*z[iz]; 
                if (j>0&&j<n-1) ak+=z[iz]*(z[iz-1]+z[iz+1]-2.*z[iz]);
                iz++; 
            }
            ak*=-hbar*hbar/2.*dx2;
            std::cout<<w[i]<<" "<<av<<" "<<ak<<"\n";
            if (do_wf && (nwf==0 || iwf<nwf)) 
            {
                std::valarray<double> nw(gsz);
                iz-=n;
                for(int j=0; j<n; ++j) nw[j]=z[iz++];
                wfdata.push_back(nw);
                iwf++;
            }
        }
        if (il>=iu) break;
        //if (w[m-1]>=vu) break;
    }
    if (do_wf) 
    {
        wfstream.precision(6);
        wfstream.width(12);
        wfstream<<dx<<"\n";
        for(int j=0; j<gsz; ++j)
        {
            //wfstream<<gx[j]<<"  ";
            for(int k=0; k<wfdata.size(); ++k) wfstream<<wfdata[k][j]<<"  ";
            wfstream<<std::endl;
        }
        wfstream.close();
    }
            
    std::cout.width(16); std::cout.precision(8); 
    
}