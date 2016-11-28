#include "color.hpp"
#include "matrix-fun.hpp"
#include "matrix-io.hpp"

#include <iostream>
#include <fstream>

using namespace toolbox;
using namespace tblapack;

void GLEABC::set_A(const DMatrix& nA)
{
    if (n==0) set_size(nA.rows());
    if (BBT.rows()>0) fr_init=true;
    if (n!=nA.rows()) ERROR("Use set_size if you want to change matrix size.");
    if (n!=nA.cols()) ERROR("Square matrix expected.");
    A=nA; fr_eva=fr_hk=fr_spec=false;
    if (fr_init && fr_c)
    {
        //update BBT, it's so cheap we can do it routinely;
        DMatrix T; mult(A, C, AC); transpose(AC, T);  BBT=AC+T; 
    }
}

void GLEABC::set_BBT(const DMatrix& nBBT)
{
    if (n==0) set_size(nBBT.rows());
    if (A.rows()>0) fr_init=true;
    if (n!=nBBT.rows()) ERROR("Use set_size if you want to change matrix size.");
    if (n!=nBBT.cols()) ERROR("Square matrix expected.");
    BBT=nBBT; fr_c=fr_spec=false;
}

void GLEABC::set_C(const DMatrix& nC)
{
    if (n==0) set_size(nC.rows());
    if (A.rows()>0) fr_init=true;
    if (n!=nC.rows()) ERROR("Use set_size if you want to change matrix size.");
    if (n!=nC.cols()) ERROR("Square matrix expected.");
    C=nC; fr_c=true; fr_spec=false;
    if (fr_init && fr_c)
    {
        //update BBT, it's so cheap we can do it routinely;
        DMatrix T; mult(A, C, AC); transpose(AC, T);  BBT=AC+T; 
    }
}

void GLEABC::get_A(DMatrix& rA)
{
    if (A.rows()>0) rA=A;
    else ERROR("Matrix is not initialized yet");
}

void GLEABC::get_BBT(DMatrix& rBBT)
{
    if (BBT.rows()>0) rBBT=BBT;
    else ERROR("Matrix is not initialized yet");
}

void GLEABC::get_C(DMatrix& rC)
{
    if (!fr_c) prepare_C();
    rC=C;
}

void GLEABC::get_cov(unsigned long i, unsigned long j, double& cov)
{
    if (!fr_c) prepare_C();
    cov=C(i,j);
}

//returns diffusion coefficient for a free particle (better, D m/kT)
void GLEABC::get_diff(double& d)
{
    if (!fr_c) prepare_C();
    if (!fr_hk) prepare_hk();

    DMatrix T1,T2; MatrixInverse(A,T1);
    mult(T1,C,T2);
    d=T2(0,0)/C(0,0);
}

void GLEABC::get_evA(std::valarray<tblapack::complex>& ra)
{
    if (!fr_eva)
    {
        EigenDecomposition(A, O, O1, a); fr_eva=true;
    }
    
    ra.resize(a.size()); ra=a;
}

void GLEABC::get_esA(std::valarray<tblapack::complex>& ra, FMatrix<tblapack::complex>& u, FMatrix<tblapack::complex>& u1)
{
    if (!fr_eva)
    {
        EigenDecomposition(A, O, O1, a); fr_eva=true;
    }
    
    ra.resize(a.size()); ra=a;
    // MR: do I need to resize below in any way?
    u=O; u1=O1;
}

void GLEABC::get_KH(double w, double& kw, double& hw)
{
    if (!fr_c) prepare_C();
    if (!fr_hk) prepare_hk();
    
    double w2=w*w;
    unsigned long ln=n-1;
    DMatrix T1, AZ;
    //computes AZ=A/(A^2+w^2)=(A+w2 A^-1)^-1
    if (ln>0) 
    {
        T1=lA1*w2;
        T1+=lA;
        MatrixInverse(T1,AZ);
    }
    std::valarray<double> AWa(n);
    AWa=0.;
    for (int i=0; i<ln; ++i)
        for (int j=0; j<ln; ++j)
    { AWa[j]+=A(0,i+1)*AZ(i,j); }
    
    kw=A(0,0);
    for (int j=0; j<ln; ++j) kw-=AWa[j]*A(j+1,0);
    kw*=2;
    
    //compute matrices needed in the calculation of H(w)
    DMatrix Z,AZA;
    if (ln>0) { mult(AZ,lA1,Z);  mult(lA,AZ,AZA); }
    
    double t1, t2, t3;
    t1=t2=0.;
    for (int i=0; i<ln; ++i) for (int j=0; j<ln; ++j)
    { t3=A(0,i+1)*C(0,j+1); t1+=t3*AZ(i,j); t2+= t3*Z(i,j); }
    hw=kw*(C(0,0)-t1);
    
    t3=0.;
    for (int i=0; i<ln; ++i) for (int j=0; j<ln; ++j) t3+=A(0,i+1)*A(j+1,0)*Z(i,j);
    hw+=2.*w2*t2*(1.+t3);
}


void GLEABC::get_KHt(double& t, double& kt, double& ht)
{
    if (!fr_c) prepare_C();
    if (!fr_hk) prepare_hk();
    
    kt=ht=0.; 
    if (t==0.) { kt+=2.*A(0,0); ht+=BBT(0,0); }
    
    unsigned long ln=n-1;
    if (ln==0) return;
    
    CMatrix T1(ln,ln,0.), T2; DMatrix T(ln,ln);
    for (unsigned long i=0; i<ln; ++i) T1(i,i)=exp(-t*la[i]); mult(lAO,T1,T2); mult(T2,lAO1,T1);
    for (unsigned long i=0; i<ln; ++i) for (unsigned long j=0; j<ln; ++j) T(i,j)=real(T1(i,j));
    //now, T holds exp(-t lA)
    
    for (unsigned long i=0; i<ln; ++i) for (unsigned long j=0; j<ln; ++j) kt-=A(0,i+1)*T(i,j)*A(j+1,0);
    for (unsigned long i=0; i<ln; ++i) for (unsigned long j=0; j<ln; ++j) ht+=A(0,i+1)*T(i,j)*(lZap[j]-BBT(0,j+1));
}

void GLEABC::get_acf(unsigned long i, unsigned long j, double t, double& acf)
{
    if (!fr_c) prepare_C();
    FMatrix<double> T, TC;
    toolbox::exp(A*(-fabs(t)),T,1e-20);
    toolbox::mult(T,C,TC);
    acf=TC(i,j);
}

// here we plot quantities which are proportional to the mean square displacement and its derivative,
// i.e.
// MSD  = m/kT <r^2(t)>/2t
// DMSD = m/kT d/dt<r^2(t)>/2
// proportionality is defined such that for t->infty they should converge to the diffusion coefficient
void GLEABC::get_msd(unsigned long i, double t, double& msd, double& dmsd)
{
    if (!fr_c) prepare_C();
    DMatrix T, T1, iA, ACF, IACF;
    toolbox::exp(A*(-fabs(t)),T,1e-20);
    MatrixInverse(A,iA);
    toolbox::mult(T,C,ACF);
    T1=-T; for (int k=0; k<T.rows(); ++k) T1(k,k)+=1.;
    mult(iA,T1,T); mult(T,C,IACF);
    
    msd=IACF(i,i)/C(i,i);
    dmsd=(IACF(i,i)+t*ACF(i,i))/C(i,i);
}

void GLEABC::get_tau2(unsigned long i, unsigned long j, unsigned long k, unsigned long l, double& tau)
{
    tau=0.5*(x(i,j,k,l)+x(i,j,l,k)+x(k,l,i,j)+x(l,k,i,j));
}

void GLEABC::get_acf(unsigned long i, unsigned long j, unsigned long k, unsigned long l, double t, double& acf)
{
    if (!fr_c) prepare_C();
    FMatrix<double> T, TC;
    toolbox::exp(A*(-fabs(t)),T,1e-20);
    toolbox::mult(T,C,TC);
    acf=TC(i,k)*TC(j,l)+TC(i,l)*TC(j,k);
}

double GLEABC::x(unsigned long i, unsigned long j, unsigned long k, unsigned long l)
{
    if (!fr_c) prepare_C();
    
    double x=0.; complex inx;
    std::valarray<complex> qil(n), qjk(n);
    
    for (unsigned long m=0; m<n; ++m) { qil[m]=O(i,m)*O1CT(l,m); }
    if (i==j && l==k) //handles special case quickly
    {
        for (unsigned long m=0; m<n; ++m)
        {
            x+=0.5*real(qil[m]*qil[m]/a[m]);
            inx=0.;
            for (unsigned long p=0; p<m; ++p) inx+=qil[p]/(a[m]+a[p]);
            x+=2.*real(qil[m]*inx);
        }
    }
    else 
    {
        for (unsigned long m=0; m<n; ++m) { qjk[m]=O(j,m)*O1CT(k,m); }
        for (unsigned long m=0; m<n; ++m)
        {
            inx=0.;
            for (unsigned long p=0; p<n; ++p) inx+=qjk[p]/(a[m]+a[p]);
            x+=real(qil[m]*inx);
        }
    }
    return 0.5*x;
}

void GLEABC::prepare_hk()
{
    if (!fr_init) ERROR("Object is not initialized yet");
    
    unsigned long ln=n-1; lA.resize(ln,ln); lAO.resize(ln,ln); lAO1.resize(ln,ln); la.resize(ln); lZap.resize(ln);
    if (ln==0) return;
    for (int i=0; i<ln; ++i)
        for (int j=0; j<ln; ++j)
    { lA(i,j)=A(i+1,j+1); }

    MatrixInverse(lA,lA1);
    EigenDecomposition(lA,lAO,lAO1,la);

    //prepare the lZ matrix
    CMatrix T1,T2,T3;
    transpose(lAO1,T3);
    FMatrix<double> lBBT(ln,ln);
    for (int i=0; i<ln; ++i) for (int j=0; j<ln; ++j) lBBT(i,j)=BBT(i+1,j+1);
    mult(lAO1,lBBT,T1); mult(T1,T3,T2); 
    
    for (int i=0; i<ln; ++i) for (int j=0; j<ln; ++j) T2(i,j)*=(1./(la[i]+la[j]));
    
    transpose(lAO,T3);
    mult(lAO,T2,T1); mult(T1,T3,T2);
    lZap=0.; 
    for (int i=0; i<ln; ++i) for (int j=0; j<ln; ++j) lZap[i]+=real(T2(i,j))*A(0,i+1);
    
    fr_hk=true;
}

void GLEABC::prepare_C()
{
    //this computes C and (at no or little extra cost) gets matrix ready for tau!
    if (!fr_init) ERROR("Object is not initialized yet");
    if (!fr_eva)
    {
        EigenDecomposition(A, O, O1, a); fr_eva=true;
    }
    
    CMatrix T1, T2, G;
    mult(O1,BBT,T1);
    transpose(O1,T2); 
    mult(T1,T2,G);
    
    tblapack::complex aij;
    for (int i=0; i<n; ++i) for (int j=0; j<n; ++j) G(i,j)/=(a[i]+a[j]); 
    transpose(O,T2);
    mult(G,T2,T1);
    /*!*/ transpose(T1,O1CT); //BIG BIG SAVER!
    mult(O,T1,G);
    C.resize(n,n);
    for (int i=0; i<n; ++i) for (int j=0; j<n; ++j) C(i,j)=real(G(i,j));
    fr_c=true;
    mult(A,C,AC); mult(O1,AC,O1AC);
    fr_spec=true;
}

void rfsum(const std::valarray<tblapack::complex>& p1, const std::valarray<tblapack::complex>& q1,
           const std::valarray<tblapack::complex>& p2, const std::valarray<tblapack::complex>& q2,
           std::valarray<tblapack::complex>& p, std::valarray<tblapack::complex>& q)
{
    //given the coefficients of the polynomials in P1/Q1+P2/Q2 computes the coefficients
    //in the polynomials in P/Q (the reduced sum). NO SIMPLIFICATIONS ARE DONE (OF COURSE)
    q.resize(q1.size()+q2.size()-1); q=0.;
    for (int i=0; i<q1.size(); i++) for (int j=0; j<q2.size(); j++) q[i+j]+=q1[i]*q2[j];
            
    if (p1.size()==0 || (p2.size()>0 && q1.size()+p2.size()>p1.size()+q2.size())) p.resize(q1.size()+p2.size()-1);
    else if (p1.size()>0) p.resize(p1.size()+q2.size()-1);
    else p.resize(0);
    for (int i=0; i<q1.size(); i++) for (int j=0; j<p2.size(); j++) p[i+j]+=q1[i]*p2[j];
    for (int i=0; i<q2.size(); i++) for (int j=0; j<p1.size(); j++) p[i+j]+=q2[i]*p1[j];
}

//gets the minimum value of K(w) for real w.
void GLEABC::get_kw_min(double& rmin)
{
    //first, we need to construct the simple rational functions corresponding to 
    //blocks in ap A/(A^2+w^2) bap. here w^2=x.
    if (!fr_c) prepare_C();
    if (!fr_hk) prepare_hk();
    
    std::valarray<tblapack::complex> v(n-1), bv(n-1), a;
    FMatrix<tblapack::complex> S, S1; 
    EigenDecomposition(lA,S,S1,a);
    //builds the projected vectors (these are needed only here)
    v=0.; bv=0.;

    for (int i=0; i<n-1; ++i) for (int j=0; j<n-1; ++j) v[i]+=A(0,j+1)*S(j,i);
    for (int i=0; i<n-1; ++i) for (int j=0; j<n-1; ++j) bv[i]+=A(j+1,0)*S1(i,j);
    
    std::valarray<tblapack::complex> p(0),q(1),np,nq,vp(1),vq(3);
    q[0]=1.;
    for (int i=0; i<n-1; ++i)
    { 
        vp[0]=-a[i]*v[i]*bv[i]; vq[1]=a[i]*a[i]; vq[0]=vq[1]*vq[1]; vq[1]*=2.; vq[2]=1.; 
        if (i<n-2 && abs(real(a[i])-real(a[i+1]))<abs(imag(a[i])-imag(a[i+1]))) //complex root!
        {
            //for extra accuracy, I compute complex roots in pairs...
            std::valarray<complex> cvp(1), cvq(3), cnp, cnq;
            cvp[0]=conj(vp[0]); cvq[1]=conj(vq[1]); cvq[0]=conj(vq[0]); cvq[2]=1.;
            rfsum(vp,vq,cvp,cvq,cnp,cnq);
            rfsum(p,q,cnp,cnq,np,nq);
            i++;
        }
        else 
        {
            rfsum(p,q,vp,vq,np,nq);
        }
        p.resize(np.size()); p=np; q.resize(nq.size()); q=nq; 
    }
    
    //now builds companion matrix
    int m=p.size()-1;
    while(abs(real(p[m]))<1e-200) m--;
    FMatrix<tblapack::complex> CM(m,m), LQ, RQ; CM*=0.;
    std::valarray<tblapack::complex> z;
     
    if (m>=1)
    {
        for (int i=0; i<m-1; ++i)  CM(i+1,i)=1.;
        for (int i=0; i<m; ++i) CM(0,i)=-p[m-i-1]/p[m];

        EigenSolver(CM,LQ,RQ,z);
    }
    
    //now the eigenvalues in z correspond to extremal points in K(w). to be overly sure, 
    //I compute the value of K(w) for the modulus of ALL the ev, even when they're complex and so lies off the real axis.
    double mink, ck, ch;
    mink=A(0,0); //this is the value at infinity
    get_KH(0.,ck,ch); if (ck<mink) mink=ck; //this is the value at zero
    for (int i=0; i<m; ++i)
    { get_KH(abs(z[i]),ck,ch); if (ck<mink) mink=ck; }
    
    rmin=mink;
}

void verlet_check(const DMatrix& A, const DMatrix& C, double w, double dt, double& q2, double& p2, double& pq)
{
    DMatrix Tp,Sp,C1,C2,C3,W,U,T1,T2,T3;
    get_TS(A, C, dt*0.5, Tp, Sp);
    unsigned long n=A.rows();
    C1.resize(n+1,n+1); C1*=0.; for (int i=0; i<n+1;i++) C1(i,i)=1.0; C2=C1; C2(0,0)=0.0; C3=C1;
    for (int i=0; i<n;i++) for (int j=0; j<n;j++) { C1(i+1,j+1)=Tp(i,j); C2(i+1,j+1)=Sp(i,j); }
    double phi=w*w*dt*0.5;
    C3(0,0)=1.0-phi*dt;  C3(0,1)=dt; C3(1,0)=phi*(phi*dt-2.0); C3(1,1)=(1.0-phi*dt);
    
    mult(C1,C3,T1); mult(T1,C1,W); 
    mult(C1,C3,T1); mult(T1,C2,T2); transpose(T2,T1); mult(T2,T1,U);
    transpose(C2,T1); mult(C2,T1,T2); U+=T2;
    
    std::valarray<std::complex<double> >wv; toolbox::FMatrix<std::complex<double> > WO, WO1, R1, R2, R3;
    EigenDecomposition(W,WO,WO1,wv);
    mult(WO1,U,R1); transpose(WO1,R2); mult(R1,R2,R3); 
    for (int i=0; i<n+1;i++) for (int j=0; j<n+1;j++) R3(i,j)*=std::complex<double>(1.0,0.)/(std::complex<double>(1.0,0.)-wv[i]*wv[j]);
    mult(WO,R3,R1); transpose(WO,R2); mult(R1,R2,R3);
    q2=R3(0,0).real()*w*w; p2=R3(1,1).real(); pq=R3(0,1).real();
}

void spectral_analysis(GLEABC& abc, double& repole, double& impole, double& qres, double& pres)
{
    toolbox::FMatrix<double> xA, xC;
    abc.get_A(xA); abc.get_C(xC);
    unsigned int n=xA.rows();
    std::valarray<std::complex<double> >eva, resq(n), resp(n); 
    CMatrix evec, evec1, xvecC;
    abc.get_esA(eva, evec, evec1); 

    //get poles
    std::valarray<std::complex<double> > poles(n); poles=eva*std::complex<double>(0,1);

    // get residues or weights
    mult(evec1,xC,xvecC); // U-1 . C ... hoping evec1 is the inverse of evec        
    for (int k=0; k<(n);++k) {
        resq[k]=(evec(0, k)*xvecC(k, 0)/xC(0, 0));
        resp[k]=(evec(1, k)*xvecC(k, 1)/xC(1, 1));
        // std::cerr<< "SPEC " << k <<" " << std::real(poles[k])<< " "<<std::imag(poles[k])<<" "<<
        // std::real(resq[k])<<" "<< std::real(resp[k])<<"\n";
    }

    // this is all obsolete now
    // now here we are calculating a bunch of stuff. real and imaginary averages of the poles, their spread, as well as the
    // pole which has the highest weight and picking it for further characterization. 
    double mxw=0, kw, diffq, diffp;    
    // wavgq=wavgp=wspreadq=wspreadp=wimgq=wimgp=wkurtq=wkurtp=0;
    for (int k=0; k<(n);++k){
       //wavgq += std::abs(std::real(poles[k]))*std::real(resq[k]);
       //wavgp += std::abs(std::real(poles[k]))*std::real(resp[k]);
       //wimgq += std::abs(std::imag(poles[k]))*std::real(resq[k]);
       //wimgp += std::abs(std::imag(poles[k]))*std::real(resp[k]);       
       // the eig 
       kw  = std::real(resq[k]); // selects based only on q to avoid ambiguity
       if (kw > mxw){
          mxw = kw;
          repole = std::abs(std::real(poles[k]));
          impole = std::abs(std::imag(poles[k]));
          qres = 2.0*std::real(resq[k]);
          pres = 2.0*std::real(resp[k]);          
        }
    }
//
//    for (int k=0; k<(n);++k){
//        diffq=(std::abs(std::real(poles[k]))-wavgq);
//        diffp=(std::abs(std::real(poles[k]))-wavgp);
//        wspreadq += diffq*diffq*std::real(resq[k]);
//        wspreadp += diffp*diffp*std::real(resp[k]);
//        wkurtq += pow(diffq,4)*std::real(resq[k]);
//        wkurtp += pow(diffp,4)*std::real(resp[k]);
//    }
//    wkurtq=wkurtq/pow(wspreadq, 2);
//    wkurtp=wkurtp/pow(wspreadp, 2);
//    wspreadq=std::sqrt(wspreadq);
//    wspreadq=std::sqrt(wspreadq);    
}

tblapack::complex atan(tblapack::complex c)
{
    const tblapack::complex i(0.0,1.0);
    const tblapack::complex one(1.0,0.0);
    const tblapack::complex r=log( (one + i * c) / (one - i * c)) / tblapack::complex(0.0,2.0);
    return log( (one + i * c) / (one - i * c)) / tblapack::complex(0.0,2.0);
}
tblapack::complex ataninv(tblapack::complex x) { return atan(1./x); }

double lormodel(double& a, double& b, double& w){
    return 2*a*(a*a + b*b + w*w)/(toolbox::constant::pi*((a*a + b*b)*(a*a + b*b) + 2*(a - b)*(a + b)*w*w + w*w*w*w));
}

double ppspectrum(FMatrix<double>& xA, FMatrix<double>& xC, double& w){
    toolbox::FMatrix<double> xAC(xA), xIn(xA), xAux(xA);
    toolbox::mult(xA, xA, xAux);    
    for( int i = 0; i < xA.rows(); ++i )
    {    xAux(i,i)+=w*w;   }
    toolbox::MatrixInverse(xAux, xIn);
    toolbox::mult(xA, xC, xAC);
    toolbox::mult(xIn, xAC, xAux);
    return xAux(1,1)/xC(1,1)*2.0/toolbox::constant::pi; // pp part
}


double l2norm(double& a, double& b){
    return (a-b)*(a-b);
}

double adaptiveintegration(FMatrix<double>& xA, FMatrix<double>& xC, double& alor, double& blor, double a, double b, double epsilon, int bottom) {                 
    double c = (a + b)/2, h = b - a;                                                                  
    double d = (a + c)/2, e = (c + b)/2;
    double fa, fb, fc, fd, fe; 
    double spec, lmodel;

    spec=ppspectrum(xA, xC, a);
    lmodel=lormodel(alor, blor, a);
    fa=l2norm(spec, lmodel);
    spec=ppspectrum(xA, xC, b);
    lmodel=lormodel(alor, blor, b);
    fb=l2norm(spec, lmodel);
    spec=ppspectrum(xA, xC, c);
    lmodel=lormodel(alor, blor, c);
    fc=l2norm(spec, lmodel);
    spec=ppspectrum(xA, xC, d);
    lmodel=lormodel(alor, blor, d);
    fd=l2norm(spec, lmodel);
    spec=ppspectrum(xA, xC, e);
    lmodel=lormodel(alor, blor, e);
    fe=l2norm(spec, lmodel);                                                                
 
    double Q1 = h/6  * (fa + 4*fc + fb);
    double Q2 = h/12 * (fa + 4*fd + 2*fc + 4*fe + fb);
                                                                    
    if (bottom <= 0 || std::fabs(Q2- Q1) <= epsilon)                                         
      return Q2 + (Q2 - Q1)/15;            
    else                                                                    
    return adaptiveintegration(xA, xC, alor, blor, a, c, epsilon/2, bottom-1) +                    
           adaptiveintegration(xA, xC, alor, blor, c, b, epsilon/2, bottom-1);                     
}

// gets the normalized power spectrum for 
// a GLE dynamics described by A and C. selects frequency w and picks the index-th component
double corr_fun(FMatrix<double>& xA, FMatrix<double>& xC, double w, int index)
{
    FMatrix<double> A(xA), A2;
    mult(A,A,A2);
    for (int i=0; i<=A.rows(); ++i) A2(i,i)+=w*w;
    MatrixInverse(A2, A);
    mult(A,xA,A2);
    mult(A2,xC,A);
    return 2/toolbox::constant::pi*A(index,index)/xC(index, index); 
}

double GLEABC::get_pwspec(unsigned long i, unsigned long j, double w)
{
    if (!fr_eva)
    {
        EigenDecomposition(A, O, O1, a); fr_eva=true;
    }
    if (!fr_spec)
    {        
        mult(A,C,AC); mult(O1,AC,O1AC);
        fr_spec=true;
    }
    
    //COMPUTES [A/(A^2+w^2)C]_{i,j}
    std::valarray<tblapack::complex> ra(a); tblapack::complex one;
    one=1.0;  ra*=ra; ra+=w*w; ra=one/ra;
    ra*=O.row(i); ra*=O1AC.col(j);     
    return 2/toolbox::constant::pi*ra.sum().real()/C(i,j);
}

// gets the cumulative distribution of the normalized power spectrum for 
// a GLE dynamics described by A and C. Integrates up to frequency L and picks the index-th component
double corr_cdf(FMatrix<double>& xA, FMatrix<double>& xC, double L, int index=1)
{
    FMatrix<double> AWL(xA), atAWL;
    AWL*=(1./L);
    MatrixFunction(AWL,&ataninv,atAWL);
    mult(atAWL,xC,AWL); 
    return 2/toolbox::constant::pi*AWL(index,index)/xC(index, index); 
}

#define BISEC_ACCURACY 1e-8
void corr_cdf_bisect(FMatrix<double>& xA, FMatrix<double>& xC, double target, double low, double flow, double high, double fhigh, double &ret, int index=1)
{
    double mid=0.5*(low+high), fmid=corr_cdf(xA, xC, mid, index);
    //std::cerr << mid <<" "<< fmid<<"bisec\n";
    if (fabs((high-low)/mid)<BISEC_ACCURACY) { ret=mid;  return; }
    if (fmid>target) corr_cdf_bisect(xA, xC, target, low, flow, mid, fmid, ret, index);
    else corr_cdf_bisect(xA, xC, target, mid, fmid, high, fhigh, ret, index);
}

// L2 norm of a lorenzian with median w and interquartile distance g
double overlap_lorlor(double g, double w)
{
    //return (2*g*g + w*w)/(2*g*toolbox::constant::pi*(g*g +  w*w));
    return (2*g + w*w/g)/(2*toolbox::constant::pi*(g*g +  w*w));
}

// L2 norm of a the peak originating from an eigenvalue a+Ib with weight c+Id
double overlap_spec(double a, double b, double c, double d)
{
    //return (2*a*a*c*c + b*b*c*c + 2*a*b*c*d + b*b*d*d)/(2*a*toolbox::constant::pi*(a*a + b*b));
    return (2*a*c*c + b*b*(c*c+d*d)/a + 2*b*c*d)/(2*toolbox::constant::pi*(a*a + b*b));
}

double overlap_specspec(double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2)
{
    /*return (2*((a1 + a2)*((a1 + a2)*(a1 + a2) + (b1)*(b1) + (b2)*(b2))*c1*c2 + 
b1*((a1 + a2)*(a1 + a2) + (b1 - b2)*(b1 + b2))*c2*d1 + b2*(((a1 + 
a2)*(a1 + a2) - (b1)*(b1) + (b2)*(b2))*c1 + 2*(a1 + 
a2)*b1*d1)*d2))/(((a1 + a2)*(a1 + a2) + (b1)*(b1))*((a1 + a2)*(a1 + 
a2) + (b1)*(b1)) + 2*(a1 + a2 - b1)*(a1 + a2 + b1)*(b2)*(b2) + 
(b2)*(b2)*(b2)*(b2))/toolbox::constant::pi;*/
    double bbar=0.5*(fabs(b1)+fabs(b2));
    if (bbar > 1e-10)
    {
        b1=b1/bbar; b2=b2/bbar;
        return ((2*c2*(((a1 + a2)*(a1 + a2) + 2*(bbar)*(bbar))*c1 + (a1 + a2)*bbar*d1) + 2*bbar*((a1 + a2)*c1 + 2*bbar*d1)*d2)/((a1 + a2)*((a1 + a2)*(a1 + a2) + 4*(bbar)*(bbar)))+ (-2*bbar*((a1)*(a1)*(a1)*(a1)*(a1)*(c2*(d1 - b1*d1) - (-1 + b2)*c1*d2) + (a2)*(a2)*(a2)*(a2)*(a2)*(c2*(d1 - b1*d1) - (-1 + b2)*c1*d2) + 2*((b1)*(b1) - (b2)*(b2))*((b1)*(b1) - (b2)*(b2))*(bbar)*(bbar)*(bbar)*(bbar)*(bbar)*(c1*c2 + d1*d2) + (a2)*(a2)*(b1 - b2)*(b1 - b2)*(bbar)*(bbar)*(bbar)*((b1)*(b1)*c1*c2 + 2*b1*b2*c1*c2 + (b2)*(b2)*c1*c2 + 4*d1*d2) + (a2)*(a2)*(a2)*(a2)*bbar*((-2 + (b1)*(b1) + (b2)*(b2))*c1*c2 - 2*(-1 + b1*b2)*d1*d2) + a2*((b1)*(b1) - (b2)*(b2))*(bbar)*(bbar)*(bbar)*(bbar)*(-4*b1*c2*d1 + (b1)*(b1)*(c2*d1 + c1*d2) - b2*(b2*c2*d1 - 4*c1*d2 + b2*c1*d2)) + (a1)*(a1)*(a1)*(a1)*(5*a2*(c2*(d1 - b1*d1) - (-1 + b2)*c1*d2) + bbar*((-2 + (b1)*(b1) + (b2)*(b2))*c1*c2 - 2*(-1 + b1*b2)*d1*d2)) + (a2)*(a2)*(a2)*(bbar)*(bbar)*(-((b1)*(b1)*(b1)*c2*d1) + b1*(-4 + (b2)*(b2))*c2*d1 + (b1)*(b1)*(2*c2*d1 + (2 + b2)*c1*d2) + b2*(-4*c1*d2 - (b2)*(b2)*c1*d2 + 2*b2*(c2*d1 + c1*d2))) + (a1)*(a1)*(-10*(a2)*(a2)*(a2)*((-1 + b1)*c2*d1 + (-1 + b2)*c1*d2) + (b1 - b2)*(b1 - b2)*(bbar)*(bbar)*(bbar)*((b1)*(b1)*c1*c2 + 2*b1*b2*c1*c2 + (b2)*(b2)*c1*c2 + 4*d1*d2) + 6*(a2)*(a2)*bbar*((-2 + (b1)*(b1) + (b2)*(b2))*c1*c2 - 2*(-1 + b1*b2)*d1*d2) - 3*a2*(bbar)*(bbar)*((b1)*(b1)*(b1)*c2*d1 - b1*(-4 + (b2)*(b2))*c2*d1 - (b1)*(b1)*(2*c2*d1 + (2 + b2)*c1*d2) + b2*(4*c1*d2 + (b2)*(b2)*c1*d2 - 2*b2*(c2*d1 + c1*d2)))) + a1*(-5*(a2)*(a2)*(a2)*(a2)*((-1 + b1)*c2*d1 + (-1 + b2)*c1*d2) + 2*a2*(b1 - b2)*(b1 - b2)*(bbar)*(bbar)*(bbar)*((b1)*(b1)*c1*c2 + 2*b1*b2*c1*c2 + (b2)*(b2)*c1*c2 + 4*d1*d2) + 4*(a2)*(a2)*(a2)*bbar*((-2 + (b1)*(b1) + (b2)*(b2))*c1*c2 - 2*(-1 + b1*b2)*d1*d2) + ((b1)*(b1) - (b2)*(b2))*(bbar)*(bbar)*(bbar)*(bbar)*(-4*b1*c2*d1 + (b1)*(b1)*(c2*d1 + c1*d2) - b2*(b2*c2*d1 - 4*c1*d2 + b2*c1*d2)) - 3*(a2)*(a2)*(bbar)*(bbar)*((b1)*(b1)*(b1)*c2*d1 - b1*(-4 + (b2)*(b2))*c2*d1 - (b1)*(b1)*(2*c2*d1 + (2 + b2)*c1*d2) + b2*(4*c1*d2 + (b2)*(b2)*c1*d2 - 2*b2*(c2*d1 + c1*d2)))) + (a1)*(a1)*(a1)*(-10*(a2)*(a2)*((-1 + b1)*c2*d1 + (-1 + b2)*c1*d2) + 4*a2*bbar*((-2 + (b1)*(b1) + (b2)*(b2))*c1*c2 - 2*(-1 + b1*b2)*d1*d2) + (bbar)*(bbar)*(-((b1)*(b1)*(b1)*c2*d1) + b1*(-4 + (b2)*(b2))*c2*d1 + (b1)*(b1)*(2*c2*d1 + (2 + b2)*c1*d2) + b2*(-4*c1*d2 - (b2)*(b2)*c1*d2 + 2*b2*(c2*d1 + c1*d2))))))/((a1 + a2)*((a1)*(a1) + 2*a1*a2 + (a2)*(a2) + 4*(bbar)*(bbar))*((a1)*(a1) + 2*a1*a2 + (a2)*(a2) + (b1 - b2)*(b1 - b2)*(bbar)*(bbar))*((a1)*(a1) + 2*a1*a2 + (a2)*(a2) + (b1 + b2)*(b1 + b2)*(bbar)*(bbar))))/toolbox::constant::pi;
    }
    else
    {
        return (2*((a1 + a2)*((a1 + a2)*(a1 + a2) + (b1)*(b1) + (b2)*(b2))*c1*c2 + b1*((a1 + a2)*(a1 + a2) + (b1 - b2)*(b1 + b2))*c2*d1 + b2*(((a1 + a2)*(a1 + a2) - (b1)*(b1) + (b2)*(b2))*c1 + 2*(a1 + a2)*b1*d1)*d2))/(((a1 + a2)*(a1 + a2) + (b1)*(b1))*((a1 + a2)*(a1 + a2) + (b1)*(b1)) + 2*(a1 + a2 - b1)*(a1 + a2 + b1)*(b2)*(b2) + (b2)*(b2)*(b2)*(b2))/toolbox::constant::pi;
    }
}

// overlap of a peak and a lorentzian
double overlap_lorspec(double g, double w, double a, double b, double c, double d)
{
  /*  return (2*((b*d + c*(a + g))*((b)*(b) + (a + g)*(a + \
g))*toolbox::constant::pi + (-(b*d) + c*(a + \
g))*toolbox::constant::pi*(w)*(w)))/(((b)*(b) + (a + g)*(a + \
g))*((b)*(b) + (a + g)*(a + g)) + 2*(a - b + g)*(a + b + g)*(w)*(w) + \
(w)*(w)*(w)*(w))/(toolbox::constant::pi*toolbox::constant::pi);
 */
    b=b/w;  // more complicated but more stable expression
    return (2*c*(g)*(g) + c*(w)*(w))/(2*(g)*(g)*(g) + 2*g*(w)*(w))/toolbox::constant::pi
    + (-2*c*(a - g)*(g)*(g)*(a + g)*(a + g)*(a + g) + 4*b*d*(g)*(g)*(g)*(a 
+ g)*(a + g)*w - c*(a + g)*((a)*(a)*(a) - (a)*(a)*g + a*(-1 + 
4*(b)*(b))*(g)*(g) - 3*(g)*(g)*(g))*(w)*(w) + 4*b*d*g*((a)*(a) + 
2*a*g + (b)*(b)*(g)*(g))*(w)*(w)*(w) - 2*c*((a)*(a)*(1 + (b)*(b)) + 
(b)*(b)*(-3 + (b)*(b))*(g)*(g))*(w)*(w)*(w)*(w) + 4*b*(-1 + 
(b)*(b))*d*g*(w)*(w)*(w)*(w)*(w) - (-1 + b*b)*(-1 + b*b)
*c*(w)*(w)*(w)*(w)*(w)*(w))/(2.*g*((g)*(g) + (w)*(w))*((a 
+ g)*(a + g) + (-1 + b)*(-1 + b)*(w)*(w))*((a + g)*(a + g) + (1 + 
b)*(1 + b)*(w)*(w)))/toolbox::constant::pi;
}

//analyzes the shape of a peak in the harmonic distribution 
void harm_shape(GLEABC& abc, double w, double& specdiff, double& median, double& interq, int index=0)
{
    
    toolbox::FMatrix<double> xA, xC;
    abc.get_A(xA); abc.get_C(xC);
    unsigned int n=xA.rows();
    
    //find 1,2,3 quartiles
    //the total integral under the peak is 1
    
    double Lhigh = w, cdfhigh=corr_cdf(xA, xC, Lhigh, index);
    double Llow=w, cdflow=cdfhigh;
    double Lmid=w, cdfmid=cdfhigh;
    // computes the integral of Cpp from 0 to L
    while (cdfhigh<0.75) // upper bracket
    {   Lhigh*=2; cdfhigh= corr_cdf(xA, xC, Lhigh, index);  }
    while (cdflow>0.25) //lower bracket
    {   Llow/=2;  cdflow=corr_cdf(xA, xC, Llow, index);     }
    // std::cerr<<"CDFl "<<Llow<<" "<<cdflow<<std::endl;
    // std::cerr<<"CDFh "<<Lhigh<<" "<<cdfhigh<<std::endl;
    double LD50, LD25, LD75;
    corr_cdf_bisect(xA, xC, 0.5, Llow, cdflow, Lhigh, cdfhigh, LD50, index);
    corr_cdf_bisect(xA, xC, 0.75, LD50, 0.5, Lhigh, cdfhigh, LD75, index);
    corr_cdf_bisect(xA, xC, 0.25, Llow, cdflow, LD50, 0.5, LD25, index);
        
    median=LD50;
    interq=(LD75-LD25)*0.5;
    //std::cerr<<"LDs "<< LD25<<" "<< LD50 <<" "<< LD75 << std::endl;

    // integrate square distance between both from 0 to wf with an adaptive grid...
//    double wi=0, wf=LD50+100*(LD75-LD25);
 //   double alor, blor; 


    // Here define a fake Lorentzian with these parameters
    
   // alor=(LD75-LD25)*0.5;
   // blor=LD50*LD50-alor*alor;
    
   // if (blor>0.0) { blor=std::sqrt(blor);}
   // else{ blor=0.0; }
    
//    std::cerr<<"CDFMED "<<corr_cdf(xA, xC, LD50)<<std::endl;

    // std::cerr<<"A and B lorentizan "<< alor <<" "<< blor <<" \n";
 
    //specdiff=adaptiveintegration(xA, xC, alor, blor, wi, wf, 0.0001, 500);


    // now we set off to compute something tricky. we want a measure of how
    // much the spectrum deviates from a lorentzian. we estimate the Lorentzian
    // parameters from the median and interquartile distance, then compute the 2-norm
    // of the difference between such lorentzian and the actual power spectrum
    std::valarray<tblapack::complex> ra; CMatrix U, U1, U1C;
    abc.get_esA(ra, U, U1);
    mult(U1,xC,U1C);
    U1C*=1.0/xC(index,index);  // normalize for the evaluation of the desired (index,index) component
    
    double lw=LD50, lg=(LD75-LD25)*0.5; // parameters of the reference Lorentzian 
    double ai, bi, ci, di, aj, bj, cj, dj;
    double i11=0; int na=ra.size();
    double ri1, di1, ri2, di2;
#define PEAK_SMOOTH 1e-8    
    lg+=lw*PEAK_SMOOTH;    // stabilizes expressions by smoothening too-sharp peaks
    //std::cerr<<w<<" REFLOR "<<lw<<" "<<lg<<std::endl;
    for (int i=0; i<na; ++i)
    {
        
        ai = ra[i].real(); bi= ra[i].imag(); ci=(U(index,i) * U1C(i,index)).real(); di=(U(index,i) * U1C(i,index)).imag();   
        ai += lw*PEAK_SMOOTH;
        //std::cerr<<"PEAK DATA " << ai <<" "     << bi <<" "<< ci <<" "<< di <<" \n";
        //std::cerr<<"OVERLAP "<<overlap_lorspec(lg,lw,ai,bi,ci,di)/ci<< " SELF "<< overlap_spec(ai,bi,ci,di)/(ci*ci)<<" ref "<<overlap_lorlor(lg,lw)<<std::endl;
        //overlap_lorspec(lg,lw,ai,bi,ci,di,ri1,di1);
        //std::cerr<<"RESOLVE? "<< ri1/ci<<"  "<<di1/ci<<" = "<<(ri1+di1)/ci<<std::endl;
        //std::cerr<<"STABILITY CHECK \n";
        //overlap_lorlor(lg, lw, ri1, di1);
        //overlap_spec(ai,bi,ci,di, ri2, di2);
        //std::cerr<<"LOR-SPEC: div: "<< di1/na-di2<< "  reg: " << ri1/na-ri2<<std::endl ;
        //std::cerr<<"SEP-FULL: div: "<< (di1+ri1)-overlap_lorlor(lg,lw)<<std::endl ;
        i11 -= 2*overlap_lorspec(lg,lw,ai,bi,ci,di);
            
        for (int j=0; j<na; ++j)
        {
                        
            aj = ra[j].real(); bj= ra[j].imag(); cj=(U(index,j) * U1C(j,index)).real(); dj=(U(index,j) * U1C(j,index)).imag();
            aj += lw*PEAK_SMOOTH;
            if (i==j)  i11 += overlap_spec(ai,bi,ci,di);
            else  i11 += overlap_specspec(ai,bi,ci,di,aj,bj,cj,dj);
            //std::cerr<<"TEST SPECSPEC "<<i<<" " <<j<<" >> "<<overlap_specspec(ai,bi,ci,di,aj,bj,cj,dj)/(ci*cj)<<std::endl;
            
            //i11 -= overlap_lorspec(lg,lw,aj,bj,cj,dj)/(na*xC(1,1));            
            //std::cerr<<li<<" "<<lj<<" "<<symlor_sp(li, lj)<<"  "<<symlor_sp(lj, li)<<"\n";
            //i11 += U(1,i) * U1C(i,1) * U(1,j) * U1C(j,1) * symlor_sp(li, lj);
            //i11 += symlor_sp(llor, llor)/(na*na);
            //i11 -= U(1,i) * U1C(i,1) * (symlor_sp(li, llor) /na) ;
            //ei11 -= U(1,j) * U1C(j,1) * ( symlor_sp(llor, lj) / na);
        }
    }
    i11 += overlap_lorlor(lg,lw);       
    //std::cerr<<ra.size()<<" NUMERICAL: "<<specdiff<<"  ANALYTICAL:  "<<i11<<std::endl;
    specdiff = sqrt(fabs(i11)/overlap_lorlor(lg,lw));
    // std::cerr<<"LDs "<< LD25<<" "<< LD50 <<" "<< LD75 << std::endl;
    // std::cerr<<"sqdiff: "<<specdiff<<std::endl;


}


void harm_check(const DMatrix& A, const DMatrix& BBT, double w, double &tq2, double &tp2, double& th, double& q2, double& p2, double& pq, double& lambdafp, double& repole, double& impole, double& qres, double& pres, double& median, double& interq, double& specdiff)
{
    unsigned long n=A.rows();
    double w2=w*w, w4=w2*w2; 
    
    //prepares extended matrices
    toolbox::FMatrix<double> xA(n+1,n+1), xBBT(n+1,n+1), xC;
    xA*=0.; xBBT=xA;
    for (int i=0; i<n;++i)for (int j=0; j<n;++j)
    { xA(i+1,j+1)=A(i,j);  xBBT(i+1,j+1)=BBT(i,j); }
    xA(0,1)=-1; xA(1,0)=w2;   //sets the harmonic hamiltonian part
    GLEABC abc; abc.set_A(xA); abc.set_BBT(xBBT);
    
    std::valarray<std::complex<double> >eva; abc.get_evA(eva);
    lambdafp=eva[0].real(); 
    for (int i=0; i<n+1; i++) lambdafp=(eva[i].real()<lambdafp?eva[i].real():lambdafp);
    
//    std::cerr<<" ---  C in tauw "<<w<<" ----\n"<<xC<<" ---------- \n";
    abc.get_C(xC);
    double c00=xC(0,0), c01=xC(0,1), c11=xC(1,1); 
    
    q2=c00*w2;
    p2=c11;
    pq=c01*w;

    double t00, t01, t11;
    abc.get_tau2(0,0,0,0,t00); abc.get_tau2(0,0,1,1,t01); abc.get_tau2(1,1,1,1,t11);

    double qqqq=t00*w4;
    double qqpp=t01*w2;
    double ppqq=qqpp; //come on, the simmetry is evident. let's save time!
    double pppp=t11;

//    std::cerr<<"NEW: "<<qqqq<<","<<qqpp<<","<<ppqq<<","<<pppp<<"\n";
    double qqqq0=q2*q2;
    double qqpp0=pq*pq;
    double ppqq0=qqpp0;
    double pppp0=p2*p2;

    th=(pppp+qqqq+ppqq+qqpp)/(pppp0+qqqq0+ppqq0+qqpp0);
    tp2=pppp/pppp0;
    tq2=qqqq/qqqq0;
    
    //get power spectrum peak widths
    spectral_analysis(abc, repole, impole, qres, pres);

    // get difference in shapes of spectra
    harm_shape(abc, w, specdiff, median, interq);
}

void rp_check(const DMatrix& A, const DMatrix& BBT, double w, double wrp, double alpha, double& repole, double& impole, double& qres, double& pres, double& median, double& interq, double& specdiff)
{
    unsigned long n=A.rows();
    double w2=w*w, wrp2=wrp*wrp, dw; 

    // build a model of two coupled oscillators of frequencies w and wrp, coupled by a constant dw
    dw=alpha*w*wrp/(1+(std::abs(w-wrp)/w));
    
    toolbox::FMatrix<double> xA(n+3,n+3), xBBT(n+3,n+3), xC;
    xA*=0.; xBBT=xA;
    for (int i=0; i<n;++i)for (int j=0; j<n;++j)
    { xA(i+3,j+3)=A(i,j);  xBBT(i+3,j+3)=BBT(i,j); }
    xA(0,1)=-1; xA(1,0)=w2; xA(1,2)=dw; xA(2,3)=-1; xA(3,0)=dw; xA(3,2)=wrp2;   //sets the two coupled harmonic oscilators hamiltonian part
    GLEABC abc; abc.set_A(xA); abc.set_BBT(xBBT);
    abc.get_C(xC);

    spectral_analysis(abc, repole, impole, qres, pres);

    harm_shape(abc, w, specdiff, median, interq);
    
}
/*
//dirty, dirty way of making a recursive integration function...
toolbox::FMatrix<double> __xA, __xC;
double __vvac(double w)
{
    toolbox::FMatrix<double> t1, t2, xDELTA;
    mult(__xA,__xA,t1);                         //A^2
    for (int i=0; i<__xA.rows();++i) t1(i,i)+=w*w;   //A^2+w^2
    MatrixInverse(t1,t2);                   //1/(A^2+w^2)
    mult(__xA,t2,t1);                         //A/(A^2+w^2)
    mult(t1,__xC,xDELTA);                     //A/(A^2+w^2)C
    return xDELTA(0,0)/__xC(0,0);
}

#define __PW_INTACCU 1e-3
double __intme(double xa, double xb, double fa, double fb)
{
    double xc=0.5*(xb+xa), fc=__vvac(xc), 
       io=((xb-xa)*(fa+fb)*0.5), in=(xb-xc)*(fb+fc)*0.5+(xc-xa)*(fa+fc)*0.5;
   if (fabs((in-io)/in)<__PW_INTACCU) return in; else return __intme(xa,xc,fa,fc)+__intme(xc,xb,fc,fb);
}
*/


//integrates the peak of the velocity-velocity correlation function from w*(1-d) to w*(1+d)
void harm_peak(const DMatrix& A, const DMatrix& BBT, double w, double d, double &pi)
{
    unsigned long n=A.rows();
    
    //prepares extended matrices
    toolbox::FMatrix<double> xA(n+1,n+1), xBBT(n+1,n+1), xC;
    xA*=0.; xBBT=xA;
    for (int i=0; i<n;++i)for (int j=0; j<n;++j)
    { xA(i+1,j+1)=A(i,j);  xBBT(i+1,j+1)=BBT(i,j); }
    xA(0,1)=-1; xA(1,0)=w*w;   //sets the harmonic hamiltonian part
    GLEABC abc; abc.set_A(xA); abc.set_BBT(xBBT);
    
    //std::cerr<<" ---  C in tauw "<<w<<" ----\n"<<xC<<" ---------- \n";
    abc.get_C(xC);
    
    //get power spectrum peak intensity (should program a better way to integrate....)
    //the total integral under the peak is Pi/2
    FMatrix<double> AW1(xA), AW2(xA), atAW1, atAW2;
    
    AW1*=(1./(w*(1.-d/2.)));
    AW2*=(1./(w*(1.+d/2.)));
    
    MatrixFunction(AW1,&atan,atAW1);
    MatrixFunction(AW2,&atan,atAW2);
    mult(atAW1,xC,AW1); mult(atAW2,xC,AW2);
    
    pi=2./toolbox::constant::pi*(AW1(1,1)-AW2(1,1))/xC(1,1);
}

void harm_spectrum(const DMatrix& A, const DMatrix& BBT, double w, const std::valarray<double>& wl, std::valarray<double>& cqq, std::valarray<double>& cpp )
{
    unsigned long n=A.rows();
    
    //prepares extended matrices
    toolbox::FMatrix<double> xA(n+1,n+1), xBBT(n+1,n+1), xC;
    xA*=0.; xBBT=xA;
    for (int i=0; i<n;++i)for (int j=0; j<n;++j)
    { xA(i+1,j+1)=A(i,j);  xBBT(i+1,j+1)=BBT(i,j); }
    xA(0,1)=-1; xA(1,0)=w*w;   //sets the harmonic hamiltonian part
    GLEABC abc; abc.set_A(xA); abc.set_BBT(xBBT);
    
    abc.get_C(xC);
    
    //get power spectrum peak intensity (should program a better way to integrate....)
    //the total integral under the peak is Pi/2
    FMatrix<double> xA2(xA), xAC, A2w(xA), ixA2w(xA), Cww;
    
    mult(xA,xA,xA2); mult(xA,xC,xAC);
    
    cqq.resize(wl.size());
    cpp.resize(wl.size());
    for (int i=0; i<wl.size(); ++i)
    {
        A2w=xA2; for (int j=0; j<n+1; ++j) A2w(j,j)+=(wl[i]*wl[i]);
        MatrixInverse(A2w, ixA2w);
        mult(ixA2w, xAC, Cww);
        cqq[i]=Cww(0,0)/xC(0,0);
        cpp[i]=Cww(1,1)/xC(1,1);
    }    
}


void get_TS(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& C, const double& t, 
            toolbox::FMatrix<double>& T, toolbox::FMatrix<double>& S)
{
    FMatrix<double> tmp, lC=C;
    toolbox::exp(A*(-t),T,1e-20);
    toolbox::transpose(T,S);
    toolbox::mult(lC,S,tmp);
    toolbox::mult(T,tmp,S);
    lC-=S;
    //std::cerr<<"CHOLESKY \n"<<lC<<"\n\n";
    //now C holds in fact C-TCT^T, so we must cholesky-factor C in order to get S
    StabCholesky(lC,S);
}

void A2Kt(const toolbox::FMatrix<double>& A, double dt, std::valarray<double>& K)
{
    int n=A.rows();
    toolbox::FMatrix<double> la(n-1,n-1), le, len(n-1,n-1);
    len*=0.;
    for(int i=1; i<n;++i)
    {
        len(i-1,i-1)=1.;
        for(int j=1; j<n;++j)
            la(i-1,j-1)=A(i,j);
    }
    la*=(-dt);
    toolbox::exp(la,le,1e-20);
    K=0; K[0]=A(0,0)/dt; 
    for (int s=0;s<K.size();++s)
    {
        for(int i=1; i<n;++i)
            for(int j=1; j<n;++j)
        {
            K[s]-=A(0,i)*A(j,0)*len(i-1,j-1);
        }
        mult(len,le,la); len=la;
    }
    K[0]*=2.;
}
