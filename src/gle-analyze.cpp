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
            << " USAGE: gle-analyze -a a-file [(-b b-file | -c c-file | -d d-file)] [-sa scale] \n" 
            << "     [-sc scale] [-wi w_i] [-wf w_f] [-np np] [-ww wacf] [-tex]                 \n"
            << "     [-pk delta] [-dt dt] [-tfree maxt] [-wrp wrpmd] [-dwrp rpalpha]            \n"
            << " performs in-depth analytical study of the generalized Langevin equation        \n"
            << "                            dx=-A x dt + B dW                                   \n"
            << " A matrix must be provided, wherease for the diffusion term one can provide     \n"
            << " B, D=BB^T, or C (the static covariance, AC+CA^T=BB^T). If none is given, C=1   \n"
            << " is assumed. Matrices must be provided in files with the syntax                 \n"
            << " cols <n+1>                                                                     \n"
            << " rows <n+1>                                                                     \n"
            << " data <(n+1)^2>                                                                 \n"
            << " m00 m01 m02 .... m0n                                                           \n"
            << " m10 m11 ...      m1n   ....                                                    \n"
            << " For convenience, matrices can then be scaled according to the time (-sa) and   \n"
            << " energy (-sc) scaling factors.                                                  \n"
            << " If -ww is provided, diagnostics are printed for the harmonic oscillator at the \n"
            << " frequency value supplied. If -wi and -wf options are given, integral properties\n"
            << " are plotted as a function of the frequency. np is the number of data points to \n"
            << " be plotted, in both cases, on a log scale.                                     \n"
            << " -w0 specifies the frequency of a harmonic mode for which to compute the GLE    \n"
            << " modified vibrational spectrum.                                                 \n"
            << " If -tfree is provided, MSD and related properties are produced for a free      \n"
            << " particle, from time zero up to the specified time maxt.                        \n"
            << " -tex outputs data in a form which can be post-processed with latex to give     \n"
            << " a 'facts sheet' about the matrices provided.                                   \n"
            << " -pk also performs a 'disturbance analysis' with relative window delta.         \n"
            << " -dt computes also fluctuations which would result from finite-timestep verlet. \n"
            << " -wrp also computes the disturbance introduced on the dynamics from a rp mode.  \n"
            << " -dwrp indicates the strength of coupling with the rpmd mode.                   \n";
}

int main(int argc, char **argv)
{
    CLParser clp(argc, argv);
    bool fhelp, fsingle, ffree, ftex, funit=false;
    double wi, wf, ww, w0, shifta, shiftc,tfree,dpeak,deltat, wrpmd, rpalpha;
    std::string amat, bmat, cmat, dmat;
    long np;
    bool fok=
            clp.getoption(wi,"wi",1e-2) &&
            clp.getoption(wf,"wf",1e2) &&
            clp.getoption(ww,"ww",-1.) &&
            clp.getoption(w0,"w0",1.0) &&
            clp.getoption(tfree,"tfree",0.) &&
            clp.getoption(amat,"a") &&
            clp.getoption(bmat,"b",std::string("")) &&
            clp.getoption(cmat,"c",std::string("")) &&
            clp.getoption(dmat,"d",std::string("")) &&
            clp.getoption(shifta,"sa",1.) &&
            clp.getoption(shiftc,"sc",1.) &&
            clp.getoption(dpeak,"pk",0.) &&
            clp.getoption(deltat,"dt",0.) &&
            clp.getoption(wrpmd,"wrp",0.) &&
            clp.getoption(rpalpha,"dwrp",0.) &&
            clp.getoption(np,"np",(long)1000) &&
            clp.getoption(ftex,"tex",false) &&
            clp.getoption(fhelp,"h",false);

    if (!fok || fhelp) {banner(); return 0;}
    fsingle=(ww>0.); ffree=(tfree>0.); 
         
    GLEABC abc;
    FMatrix<double> iA, iC;
    std::ifstream ifile;
    
    IField<FMatrix<double> > ifa(iA,"");
    ifile.open(amat.c_str());
    if (ifile.bad())
        ERROR("Unable to open file for A!");
    ifile >> ifa;
    ifile.close();
    unsigned long n=iA.rows()+1; iA*=shifta;
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
        if (iB.rows()!=iB.cols() || iB.rows()!=n-1)
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
        FMatrix<double> iD;
        IField<FMatrix<double> > ifd(iD,"");
        ifile >> ifd;
        if (iD.rows()!=iD.cols() || iD.rows()!=n-1)
            ERROR("Invalid size or format of matrix B.");
        iD*=shifta*shiftc;
        abc.set_BBT(iD);
        abc.get_C(iD);
    }
    else
    {
        //unit matrix
        funit=true;
        iC.resize(n-1,n-1); iC*=0.;
        for (unsigned long i=0; i<n-1; ++i) iC(i,i)=shiftc;
        abc.set_C(iC);
    }
    
    FMatrix<double> iBBT; abc.get_BBT(iBBT); abc.get_C(iC);
    
    std::cout.precision(10);
    std::cout.setf(std::ios::left);
    std::cout.setf(std::ios::scientific);
    
    //computes analytics

    double diff;
    abc.get_diff(diff);
    if (ffree)
    {
        double kt, ht, vvac, msd, dmsd;
        std::cout<<"# Free-particle analysis\n# Diffusion coefficient: "<<diff<< std::endl;
        std::cout<<"# time  K(t)  H(t)  <v(t)v(0)>  <x^2(t)>/t   d<x^2>/dt\n";
        for (double it=0; it<tfree; it+=tfree/np)
        {
            abc.get_KHt(it,kt,ht);
            abc.get_acf(0,0,it,vvac);
            abc.get_msd(0,it,msd,dmsd);
            std::cout<<it<<" "<<kt<<" "<<ht<< " "<<vvac<< " "<<msd<< " "<<dmsd<<std::endl;
        }
    }
    else if (fsingle)
    {
        GLEABC abcww;
        toolbox::FMatrix<double> Aw(n,n), BBtw(n,n), Cw;
        double tq2, tp2, th, q2, p2, pq, dwq, dwp,lfp, dummy;
        
        make_harm_abc(iA, iBBT, ww, abcww);
        harm_check(abcww,tq2,tp2,th,q2,p2,pq,lfp);
        
        abcww.get_A(Aw); abcww.get_BBT(BBtw);   abcww.get_C(Cw);
    
        //writes out detailed information on a single frequency.
        std::cout<<"# GLE analysis for frequency  "<<ww << std::endl;
        std::cout<<"# Covariance matrix: "<<std::endl;
        for (unsigned long i=0; i<n; i++)
        {
            std::cout<<"# ";
            for (unsigned long j=0; j<n; j++) std::cout<<Cw(i,j)<<" ";
            std::cout<<std::endl;
        }
        std::cout<<"# Power spectrum peak width: Dwq, Dwp"<<dwq<<" ,  "<<dwp<<std::endl;
        std::cout<<"# Correlation functions of: q^2, p^2"<<std::endl;
        double dt=0.1/ww;
        double q20, p20; abcww.get_acf(0,0,0,0,0,q20); abcww.get_acf(1,1,1,1,0,p20);
        for (unsigned long ip=0; ip<np; ip++)
        {
            double wt=ip*dt;
            double acq2, acp2;
            abcww.get_acf(0,0,0,0,wt,acq2); abcww.get_acf(1,1,1,1,wt,acp2);
            std::cout<<wt<<" "<<acq2/q20<<" "<<acp2/p20<<std::endl;
        }
        
    }
    else
    {           
        std::valarray<double> w(np), kw(np), hw(np), tq2(np), tp2(np), th(np), q2(np), p2(np), pq(np), hdist(np), lfp(np),
                 rew(np), imw(np), qw(np), pw(np), 
                 q2dt(np), p2dt(np), pqdt(np), PWw0q(np), PWgq(np), PWshapeq(np), PWw0p(np), PWgp(np), PWshapep(np), 
                 sqq(np), spp(np), rp_PWw0q(np), rp_PWgq(np), rp_PWshapeq(np) ;
                
        for (unsigned long ip=0; ip<np; ip++) w[ip]=pow(wi,(np-ip-1.)/(np-1.))*pow(wf,(1.*ip)/(np-1.));        
        //harm_spectrum(iA, iBBT, w0,w, sqq, spp);
        double dummy;
        GLEABC abcip;
        for (unsigned long ip=0; ip<np; ip++)
        {
            abc.get_KH(w[ip], kw[ip], hw[ip]);
            make_harm_abc(iA, iBBT, w[ip], abcip);
        
            harm_check(abcip,tq2[ip],tp2[ip],th[ip],q2[ip],p2[ip],pq[ip],lfp[ip]);
            harm_shape(abcip, PWw0q[ip], PWgq[ip], PWshapeq[ip],0);
            harm_shape(abcip, PWw0p[ip], PWgp[ip], PWshapep[ip],1);
            if (dpeak>0.) harm_peak(iA,iBBT,w[ip],dpeak,hdist[ip]);
            if (deltat>0) verlet_check(iA,iC,w[ip],deltat,q2dt[ip],p2dt[ip],pqdt[ip]);
            if (wrpmd>0) rp_check(iA, iBBT, w[ip], wrpmd, rpalpha, rp_PWw0q[ip], rp_PWgq[ip], rp_PWshapeq[ip]);
        }
        if (!ftex)
        {
            std::cout<<"# D kT/m = "<<diff<<"\n";
            std::cout<<"# omega  1/tau_h  1/tau_q2  1/tau_p2  K(omega)  H(omega)  <q2>(omega) <p2>(omega) <pq>(omega) lFP(omega)  Cqq["<<w0<<"](w) Cpp["<<w0<<"](w)"<<" repeak  impeak  qpeakw  ppeakw  PWw0q PWgq PWshapeq" <<
                    (deltat>0.?" <q2>,<p2>,<pq>(dt=":"")<<(deltat>0.?float2str(deltat):std::string(""))<<(deltat>0.?")   ":"")<<
                    (dpeak>0.?" peak_dist(":"")<<(dpeak>0.?float2str(dpeak):std::string(""))<<(dpeak>0.?")":"")<<
                    (wrpmd>0.?" rpmd(":"")<<(wrpmd>0.?float2str(wrpmd):std::string(""))<<(wrpmd>0.?"): PWw0q PWgq PWshapeq":"")
                    <<"\n";
            for (unsigned long ip=0; ip<np; ip++)
            {
                std::cout<<w[ip]<<"  "<<1./th[ip] //1./((pppp+qqqq+ppqq+qqpp)/(pppp0+qqqq0+ppqq0+qqpp0))
                        <<"  "<<1./tq2[ip] //1./(qqqq/qqqq0)
                        <<"  "<<1./tp2[ip] //1./(pppp/pppp0)
                        <<"  "<<kw[ip]
                        <<"  "<<hw[ip]
                        <<"  "<<q2[ip]
                        <<"  "<<p2[ip]
                        <<"  "<<pq[ip]
                        <<"  "<<lfp[ip]
                        <<"  "<<sqq[ip]
                        <<"  "<<spp[ip]                        
                        <<"  "<<rew[ip]
                        <<"  "<<imw[ip]
                        <<"  "<<qw[ip]
                        <<"  "<<pw[ip]
                        <<"  "<<PWw0q[ip]
                        <<"  "<<PWgq[ip]
                        <<"  "<<PWshapeq[ip]
                        <<"  "<<PWw0p[ip]
                        <<"  "<<PWgp[ip]
                        <<"  "<<PWshapep[ip]                           
                        <<"  "<<(deltat>0.?float2str(q2dt[ip]):"")
                        <<"  "<<(deltat>0.?float2str(p2dt[ip]):"")
                        <<"  "<<(deltat>0.?float2str(pqdt[ip]):"")
                        <<"  "<<(dpeak>0.?float2str(hdist[ip]):"")
                        <<"  "<<(wrpmd>0.?float2str(rp_PWw0q[ip]):"")
                        <<"  "<<(wrpmd>0.?float2str(rp_PWgq[ip]):"")
                        <<"  "<<(wrpmd>0.?float2str(rp_PWshapeq[ip]):"")
                        <<std::endl; 
            }
        }
        else
        {
            double bxf=ceil(log(wf)/log(10.)), bxi=floor(log(wi)/log(10.));
            double byi, byf; double psfw=8.0, psfh=5.0;
            
            std::cout.precision(4); std::cout.setf(std::ios::scientific);
            // PREAMBLE
            std::cout<<
                    " %          Output LaTeX generated by gle-analyze          %\n"
                    "\\documentclass{article}\n"
                    "\\usepackage{amsmath}\n\\usepackage{pstricks}\n\\usepackage{pst-plot,pstricks-add}\n"
                    "\\usepackage[a4paper,top=2cm,bottom=2cm,left=1.5cm,right=1.5cm]{geometry}\n"
                    "\\pagestyle{empty}\\providecommand{\\e}[1]{\\ensuremath{\\times 10^{#1}}}"
                    "\\begin{document}\n"
                    "\\begin{center}{\\Large\\bf\\sc Generalized Langevin Equation analytics}\\end{center}\n"
                    "\\begin{itemize}\n";
            //A MATRIX
            std::cout<<
                    "\\begin{samepage}\n"
                    "\\item {\\large\\bf Drift matrix $\\mathbf{A}_p$:} \\newline\n"
                    "\\begin{equation*}\n{{\\left(\\begin{array}{";
            for(int i=0; i<n-1; ++i) std::cout<<"r"; std::cout<<"}\n";
            for(int i=0; i<n-1; ++i) for(int j=0; j<n-1; ++j) std::cout<<"\\scriptstyle{\\mathtt{"<<expn(iA(i,j))<<(j==n-2?"}}\\\\\n":"}} & ");
            std::cout<<"\\end{array}\\right)}}\\end{equation*}\n\\end{samepage}\n";
            byf=byi=1.0;
            for (int i=1;i<np; ++i) if (byi>kw[i]/kw[0]) byi=kw[i]/kw[0];
            for (int i=0;i<np; ++i) if (byi>hw[i]/hw[0]) byi=hw[i]/hw[0];
            for (int i=1;i<np; ++i) if (byf<kw[i]/kw[0]) byf=kw[i]/kw[0];
            for (int i=0;i<np; ++i) if (byf<hw[i]/hw[0]) byf=hw[i]/hw[0];
            byi=floor(log(byi)/log(10.)); byf=log(byf)/log(10.);

            if (funit) 
            {
                //NO NEED TO WRITE C
                std::cout<<
                    "\\begin{samepage}\n"
                    "\\item {\\large\\bf Fluctuation-Dissipation theorem is enforced, $\\mathbf{C}_p=k_B T$}\n";
                
                //WRITES K(w)=H(w)
                std::cout.setf(std::ios::fixed);
                std::cout<<
                        "\\item {\\large\\bf Memory kernel FT, $\\color{red}{K(\\omega)/K(0)=H(\\omega)/H(0)}$}\\newline\n"
                        "\\begin{center}\\begin{pspicture}("<<psfw<<","<<psfh<<")\n"
                        "\\psset{xunit="<<psfw/(bxf-bxi)<<",yunit="<<psfh/(byf-byi)<<"}\n"
                        "\\psaxes[axesstyle=frame,xlogBase=10,ylogBase=10,"
                        "Oy="<<byi<<",Ox="<<bxi<<"]"
                        "("<<0<<","<<0<<")("<<bxf-bxi<<","<<byf-byi<<")\n";
                
                std::cout.setf(std::ios::scientific);
                std::cout<<"\\listplot[plotstyle=curve,linecolor=red,linewidth=1.0pt]{\n";
                for(int i=0; i<np; ++i) std::cout<<log(w[i])/log(10.)-bxi<<" "<<log(kw[i]/kw[0])/log(10.)-byi<<" ";
                std::cout<<"}\n";
                std::cout.setf(std::ios::fixed);
                std::cout<<
                        "\\rput("<<(bxf-bxi)/2<<","<<-(byf-byi)*0.15<<"){$\\omega$}\n"
                        "\\rput{90}("<<-(bxf-bxi)*0.15<<","<<(byf-byi)/2<<"){$K(\\omega), H(\\omega)$}";
                std::cout<<"\\end{pspicture}\\end{center}\\vspace{0.5cm}\n"
                        "\\end{samepage}\n";
            }
            else
            {
                std::cout<<
                        "\\begin{samepage}\n"
                        "\\item {\\large\\bf Free-particle covariance matrix $\\mathbf{C}_p$:} \\newline\n"
                        "\\begin{equation*}\n{{\\left(\\begin{array}{";
                for(int i=0; i<n-1; ++i) std::cout<<"r"; std::cout<<"}\n";
                for(int i=0; i<n-1; ++i) for(int j=0; j<n-1; ++j) std::cout<<"\\scriptstyle{\\mathtt{"<<expn(iC(i,j))<<(j==n-2?"}}\\\\\n":"}} & ");
                std::cout<<"\\end{array}\\right)}}\\end{equation*}\n";
                //WRITES K(w) and H(w)
                std::cout.setf(std::ios::fixed);
                std::cout<<
                        "\\item {\\large\\bf Fourier-transform of memory kernels $\\color{red}{K(\\omega)/K(0)}$ and $\\color{blue}{H(\\omega)/H(0)}$}\\newline\n"
                        "\\begin{center}\\begin{pspicture}("<<psfw<<","<<psfh<<")\n"
                        "\\psset{xunit="<<psfw/(bxf-bxi)<<",yunit="<<psfh/(byf-byi)<<"}\n"
                        "\\psaxes[axesstyle=frame,xlogBase=10,ylogBase=10,"
                        "Oy="<<byi<<",Ox="<<bxi<<"]"
                        "("<<0<<","<<0<<")("<<bxf-bxi<<","<<byf-byi<<")\n";
                
                std::cout.setf(std::ios::scientific);
                std::cout<<"\\listplot[plotstyle=curve,linecolor=red,linewidth=1.0pt]{\n";
                for(int i=0; i<np; ++i) std::cout<<log(w[i])/log(10.)-bxi<<" "<<log(kw[i]/kw[0])/log(10.)-byi<<" ";
                std::cout<<"}\n";
                std::cout.setf(std::ios::scientific);
                std::cout<<"\\listplot[plotstyle=curve,linecolor=blue,linewidth=1.0pt]{\n";
                for(int i=0; i<np; ++i) std::cout<<log(w[i])/log(10.)-bxi<<" "<<log(hw[i]/hw[0])/log(10.)-byi<<" ";
                std::cout<<"}\n";
                std::cout.setf(std::ios::fixed);
                std::cout<<
                        "\\rput("<<(bxf-bxi)/2<<","<<-(byf-byi)*0.15<<"){$\\omega$}\n"
                        "\\rput{90}("<<-(bxf-bxi)*0.15<<","<<(byf-byi)/2<<"){$K(\\omega), H(\\omega)$}";
                std::cout<<"\\end{pspicture}\\end{center}\\vspace{0.5cm}\n"
                        "\\end{samepage}\n";
            }
            
            
            byf=log(2.)/log(10.);  byi=1./(tq2[0]*w[0]); 
            for (int i=1;i<np; ++i) if (byi>1./(tq2[i]*w[i])) byi=1./(tq2[i]*w[i]); 
            for (int i=0;i<np; ++i) if (byi>1./(th[i]*w[i])) byi=1./(th[i]*w[i]); 
            byi=floor(log(byi)/log(10.));
            
            std::cout.setf(std::ios::fixed);
            std::cout<<
                "\\begin{samepage}\n"
                "\\item {\\large\\bf Sampling efficiency, for $\\color{red}{q^2}$ and $\\color{blue}{p^2+\\omega^2 q^2}$:} \\newline\n"
                "\\begin{center}\\begin{pspicture}("<<psfw<<","<<psfh<<")\n"
                "\\psset{xunit="<<psfw/(bxf-bxi)<<",yunit="<<psfh/(byf-byi)<<"}\n"
                "\\psaxes[axesstyle=frame,xlogBase=10,ylogBase=10,"
                "Oy="<<byi<<",Ox="<<bxi<<"]"
                "("<<0<<","<<0<<")("<<bxf-bxi<<","<<byf-byi<<")\n";
            std::cout<<
                    "\\psline[linecolor=black,linewidth=0.5pt](0,"<<-byi<<")("<<bxf-bxi<<","<<-byi<<")\n";
            
            std::cout.setf(std::ios::scientific);
            std::cout<<"\\listplot[plotstyle=curve,linecolor=red,linewidth=1.0pt]{\n";
            for(int i=0; i<np; ++i) std::cout<<log(w[i])/log(10.)-bxi<<" "<<-log((tq2[i]*w[i]))/log(10.)-byi<<" ";
            std::cout<<"}\n";
            std::cout<<"\\listplot[plotstyle=curve,linecolor=blue,linewidth=1.0pt]{\n";
            std::cout.setf(std::ios::scientific);
            for(int i=0; i<np; ++i) std::cout<<log(w[i])/log(10.)-bxi<<" "<<-log((th[i]*w[i]))/log(10.)-byi<<" ";
            std::cout<<"}\n";
            
            std::cout.setf(std::ios::fixed);
            std::cout<<
                    "\\rput("<<(bxf-bxi)/2<<","<<-(byf-byi)*0.15<<"){$\\omega$}\n"
                    "\\rput{90}("<<-(bxf-bxi)*0.15<<","<<(byf-byi)/2<<"){$\\kappa(\\omega)$}";
            std::cout<<"\\end{pspicture}\\end{center}\\vspace{0.5cm}\n"
                    "\\end{samepage}\n";
            std::cout<<
                    "\\item {\\large\\bf Free-particle diffusion coeff. ($mD/k_B T$): $"<<expn(diff)<<"$}\n";
        if (!funit) 
        {
            byf=byi=q2[0];
            for (int i=1;i<np; ++i) if (byi>q2[i]) byi=q2[i];
            for (int i=0;i<np; ++i) if (byi>p2[i]) byi=p2[i];
            for (int i=1;i<np; ++i) if (byf<q2[i]) byf=q2[i];
            for (int i=0;i<np; ++i) if (byf<p2[i]) byf=p2[i];
            byi=floor(log(byi)/log(10.)); byf=log(byf)/log(10.);

            std::cout.setf(std::ios::fixed);
            std::cout<<
                    "\\begin{samepage}\n"
                    "\\item {\\large\\bf $\\omega$-dependent fluctuations, $\\color{red}{\\omega^2\\left<q^2\\right>}$ and $\\color{blue}{\\left<p^2\\right>}$}\\newline\n"
                    "\\begin{center}\\begin{pspicture}("<<psfw<<","<<psfh<<")\n"
                    "\\psset{xunit="<<psfw/(bxf-bxi)<<",yunit="<<psfh/(byf-byi)<<"}\n"
                    "\\psaxes[axesstyle=frame,xlogBase=10,ylogBase=10,"
                    "Oy="<<byi<<",Ox="<<bxi<<"]"
                    "("<<0<<","<<0<<")("<<bxf-bxi<<","<<byf-byi<<")\n";
                
            std::cout.setf(std::ios::scientific);
            std::cout<<"\\listplot[plotstyle=curve,linecolor=red,linewidth=1.0pt]{\n";
            for(int i=0; i<np; ++i) std::cout<<log(w[i])/log(10.)-bxi<<" "<<log(q2[i])/log(10.)-byi<<" ";
            std::cout<<"}\n";
            std::cout.setf(std::ios::scientific);
            std::cout<<"\\listplot[plotstyle=curve,linecolor=blue,linewidth=1.0pt]{\n";
            for(int i=0; i<np; ++i) std::cout<<log(w[i])/log(10.)-bxi<<" "<<log(p2[i])/log(10.)-byi<<" ";
            std::cout<<"}\n";
            std::cout.setf(std::ios::fixed);
            std::cout<<
                    "\\rput("<<(bxf-bxi)/2<<","<<-(byf-byi)*0.15<<"){$\\omega$}\n"
                    "\\rput{90}("<<-(bxf-bxi)*0.15<<","<<(byf-byi)/2<<"){$\\omega^2\\left<q^2\\right>, \\left<p^2\\right>$}";
            std::cout<<"\\end{pspicture}\\end{center}\\vspace{0.5cm}\n"
                    "\\end{samepage}\n";
            
        }
            
            std::cout<<
                    "\\end{itemize}\n"
                    "\\end{document}\n"<<std::endl;
        }
    }
}

std::string expn(double x)
{
    std::stringstream sstr;
    sstr.precision(4); sstr.setf(std::ios::scientific);
    sstr << x;
    std::string s; sstr >> s;
    s=s.replace(s.find('e',0),1,"\\e{")+"}";
    return s;
}
