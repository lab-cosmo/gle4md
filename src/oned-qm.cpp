#include <iostream>
#include <fstream>
#include <valarray>
#include "clparser.hpp"

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

double pw, pk, ph, pd;
double vqh(double x)
{ return pw*pw/(pk+pk)*x*(1-exp(-pk*x)); }

double vdw(double x)
{ 
    double y=x/pd; y*=y;
    return ph*(y*(y-2)+1);
}

int main(int argc, char**argv)
{
    toolbox::CLParser clp(argc,argv);
    int gsz, nev; double T,hbar,mass,fbound;
    std::string wff; int nwf; bool do_wf=false;
    bool fok=clp.getoption(gsz,"grid",1000) &&
            clp.getoption(fbound,"bound",0.) &&
            clp.getoption(pw,"w",1.) &&
            clp.getoption(ph,"h",-1.) &&
            clp.getoption(pd,"d",0.) &&
            clp.getoption(pk,"k",1.) &&
            clp.getoption(nev,"n",20) &&
            clp.getoption(hbar,"hbar",1.) &&
            clp.getoption(mass,"mass",1.) &&
            clp.getoption(wff,"wf",std::string("")) &&
            clp.getoption(nwf,"nwf",0) &&
            clp.getoption(T,"t",10.);
    ph/=mass; //symetry with the anharmonic oscillator case
    std::ofstream wfstream;
    std::vector<std::valarray<double> > wfdata;
    if (wff!="") 
    {
        do_wf=true;
        wfstream.open(wff.c_str());
        if (wfstream.bad()) ERROR("Unable to open file for wavefunction output");
    }
    
    //std::cerr<<"Running for w= " << f<<std::endl;
    
    double (* v) (double x); bool fdw=false;
    if (ph<0) v=vqh; else { v=vdw; fdw=true; }
    std::valarray<double> gx(gsz);
    std::valarray<double> d(gsz), e(gsz-1);
    double tva, tvb, tvc, tfa, tfb, tfc;
    std::cerr<<pd<<"::"<<ph<<" params\n";
    double ax,bx; 
    if (fbound>0.) { ax=fbound; bx=-ax; } else {
    if (fdw) ax=-pd+pow(T,0.25); else ax=-1./pk*fabs(log((T+hbar*pw)*2*pk/(pw*pw)));
    tva=ax; tfa=mass*v(tva); 
    while (tfa<T) { tva*=2.; tfa=mass*v(tva);}
    tvb=0.; tvc=0.5*(tva+tvb); tfb=mass*v(tvb); tfc=mass*v(tvc);
    if(tfa<T) ax*=T/tfa;
    while (fabs((tfc-T)/T)>1.e-5)
    {
        if (tfc>T)
        { tva=tvc; tfa=tfc; }
        else
        { tvb=tvc; tfb=tfc; }
        tvc=0.5*(tva+tvb); tfc=mass*v(tvc);
    }
    ax=tvc;
    
    
    //gets approximate boundaries for the grid
    if (fdw) bx=-ax; else bx=(T+pw)*2.*pk/(pw*pw);
    tva=bx;  tfa=mass*v(tva);
    while (tfa<T) { tva*=2.; tfa=mass*v(tva);}
    tvb=0.; tvc=0.5*(tva+tvb); tfb=mass*v(tvb); tfc=mass*v(tvc);
    while (fabs((tfc-T)/T)>1.e-5)
    {
        if (tfc>T)
        { tva=tvc; tfa=tfc; }
        else
        { tvb=tvc; tfb=tfc; }
        tvc=0.5*(tva+tvb); tfc=mass*v(tvc);
    }
    bx=tvc;
    }
    
    std::cerr<<"Boundaries found: "<<ax<<","<<bx<< " >> "<<mass*v(ax)<<","<<mass*v(bx)<<std::endl;
    //sets grid points for x and diagonal entries
    double dx=(bx-ax)/(gsz-1), dx2=1./(dx*dx);
    for(int i=0; i<gsz; ++i) gx[i]=ax+dx*i;
    for(int i=0; i<gsz; ++i) d[i]=mass*v(gx[i])+dx2*hbar*hbar/mass;
    e=-0.5*dx2*hbar*hbar/mass;
    
    char range='I', jobz='V'; 
    int n=gsz; double vl=0, vu; 
    if (fdw) vu=(T+hbar*sqrt(8.*ph/(pd*pd)))*.5; else vu=(T+hbar*pw)*0.5; 
    double abstol=0.; int m;
    int il=1, iu=nev, ldz=n; 
    std::valarray<double> w(gsz), work(5*n);
    std::valarray<double> z(gsz*nev);
    std::valarray<int> iwork(5*n), iblock(n), isplit(n);
    int info; std::valarray<int> ifail(n);
    int iwf=0;
    std::cout<<" # WF.index  H  <V>   <T>   <V^2>   <V^3>\n";
    while(il<=n)
    {
        if (iu>n) iu=n; 
        dstevx_(&jobz, &range, &n, &d[0], &e[0], 
            &vl, &vu, &il, &iu, &abstol,
            &m, &w[0], &z[0], &ldz,
            &work[0], &iwork[0], &ifail[0], &info);
        if (info!=0) ERROR("dstevx failed with exit code "<<info);
        il=iu+1; iu=iu+nev;
        
        //now, w[i] contains the eigenvalue, while z contains the eigenvectors,
        //we might use to compute observables such as <V> and <K>
        int iz=0;
        for(int i=0; i<m; ++i)
        {
            double ak=0, av=0, av2=0, av3=0;
            for(int j=0; j<n; ++j) 
            { 
                av+=mass*v(gx[j])*z[iz]*z[iz]; 
                av2+=mass*v(gx[j])*mass*v(gx[j])*z[iz]*z[iz]; 
                av3+=mass*v(gx[j])*mass*v(gx[j])*mass*v(gx[j])*z[iz]*z[iz]; 
                if (j>0&&j<n-1) ak+=z[iz]*(z[iz-1]+z[iz+1]-2.*z[iz]);
                iz++; 
            }
            ak*=-hbar*hbar/(2.*mass)*dx2;
            std::cout<<i<<" "<<w[i]<<" "<<av<<" "<<ak<<" "<<av2<<" "<<av3<<"\n";
            if (do_wf && (nwf==0 || iwf<nwf)) 
            {
                std::valarray<double> nw(gsz);
                iz-=n;
                for(int j=0; j<n; ++j) nw[j]=z[iz++];
                wfdata.push_back(nw);
                iwf++;
            }
        }
        if (w[m-1]>=vu) break;
    }
    if (do_wf) 
    {
        for(int j=0; j<gsz; ++j)
        {
            wfstream.precision(12);
            wfstream.width(20);
        
            wfstream<<gx[j]<<"  ";
            wfstream.precision(6);
            wfstream.width(12);
            for(int k=0; k<wfdata.size(); ++k) wfstream<<wfdata[k][j]<<"  ";
            wfstream<<std::endl;
        }
        wfstream.close();
    }
            
    //for(int i=0; i<n*nev; ++i) std::cout<<z[i]<<"\n";
    /*dstebz_(&range, &order, &n, &vl, &vu, &il, &iu,
             &abstol, &d[0], &e[0], &m, &nsplit, &w[0],
             &iblock[0], & isplit[0], &work[0], &iwork[0], &info);
    */
    std::cout.width(16); std::cout.precision(8); 
    
}