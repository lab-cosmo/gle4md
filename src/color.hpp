#ifndef __COLOR_HPP
#define __COLOR_HPP 1

#define GLE_DEBUG 1
#define  __NOPEAK 1

#include "tbdefs.hpp"
#include "matrix-full-blas.hpp"
//#include "matrix-full.hpp"
#include "linalg.hpp"
#include <iostream>

typedef toolbox::FMatrix<double> DMatrix;
typedef toolbox::FMatrix<tblapack::complex> CMatrix;

class GLEABC {
private:
    unsigned long n;
    DMatrix A, BBT, C;
    bool fr_eva, fr_c, fr_init, fr_hk;
    CMatrix O, O1; std::valarray<tblapack::complex> a;
    DMatrix lA, lA1; std::valarray<double> lZap; CMatrix lAO, lAO1; std::valarray<tblapack::complex> la;
    CMatrix O1CT;

    void prepare_C();
    void prepare_hk();
    double x(unsigned long i, unsigned long j, unsigned long k, unsigned long l);

public:
    GLEABC(): n(0), A(0,0), BBT(0,0), C(0,0), fr_eva(false), fr_c(false), fr_init(false), fr_hk(false) {}
    void set_size(unsigned long nn) { fr_init=false; n=nn; }
    void set_A(const DMatrix& rA);
    void set_C(const DMatrix& rC);
    void set_BBT(const DMatrix& rBBT);

    void get_A(DMatrix& rA);
    void get_C(DMatrix& rC);
    void get_BBT(DMatrix& rBBT);
    void get_evA(std::valarray<tblapack::complex>& ra);

    void get_KH(double w, double& kw, double& hw);
    void get_tau2(unsigned long i, unsigned long j, unsigned long k, unsigned long l, double& tau);
    void get_acf(unsigned long i, unsigned long j, unsigned long k, unsigned long l, double t, double& acf);
    void get_cov(unsigned long i, unsigned long j, double& cov);
    void get_kw_min(double& rk);
    void get_diff(double& d);
    void get_KHt(double& t, double& kt, double& ht);
    void get_acf(unsigned long i, unsigned long j, double t, double& acf);
    void get_msd(unsigned long i, double t, double& msd, double& dmsd);
};

void harm_check(const DMatrix& A, const DMatrix& BBT, double w, double &tq2, double &tp2, double& th, double& q2, double& p2, double& pq, double& dwq, double& dwp, double& lambdafp);
void verlet_check(const DMatrix& A, const DMatrix& C, double w, double dt, double& q2, double& p2, double& pq);

void harm_peak(const DMatrix& A, const DMatrix& BBT, double w, double d, double &pi);

//options for a colored complex thermostat.
class ThermoOptions {
public: double beta,   // thermostat inverse temperature
        taut,      // thermostat friction
        tauf,      // real part of ac. decay
        tauw;      // imag part of ac. decay
};

std::ostream& operator<<(std::ostream& ostr, const ThermoOptions& to);
std::istream& operator>>(std::istream& istr, ThermoOptions& to);

inline void ij2q(unsigned long i, unsigned long j, unsigned long& q)
{
    if (i>=j) q=j+((i+1)*i)/2; else ij2q(j,i,q);
}

inline void q2ij(unsigned long q, unsigned long& i, unsigned long& j)
{
    i=(int) floor((sqrt(8.*q+1.)-1.)*0.5);
    j=q-(i*(i+1))/2;
}

inline void ij2q(unsigned long n, unsigned long i, unsigned long j, unsigned long& q)
{  q=i*n+j; }

inline void q2ij(unsigned long n, unsigned long q, unsigned long& i, unsigned long& j)
{
    i=q/n;
    j=q-n*i;
}

struct MDeco {
    toolbox::FMatrix<tblapack::complex> O, O1;
    std::valarray<tblapack::complex> w;
};

void EVDeco(const toolbox::FMatrix<tblapack::complex>& A, MDeco& deco);

class XTau {
private:
    unsigned long sz;
    toolbox::FMatrix<tblapack::complex> Q, Q1CT, AMN;
    std::valarray<tblapack::complex> a;
public:
    XTau() : sz(0) {}
    XTau(const MDeco& Adec, const toolbox::FMatrix<double>& C) { prepare(Adec,C); }
    XTau(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& C) { prepare(A,C); }
    void prepare(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& C);
    void prepare(const MDeco& Adec, const toolbox::FMatrix<double>& C);

    double operator() (unsigned long i, unsigned long j, unsigned long k, unsigned long l);
    double tau(unsigned long i, unsigned long j, unsigned long k, unsigned long l)
    {return 0.5*((*this)(i,j,k,l)+(*this)(i,j,l,k)+(*this)(k,l,i,j)+(*this)(l,k,i,j));}
};

void AB2C(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& B, toolbox::FMatrix<double>& C);
void AD2C(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& D, toolbox::FMatrix<double>& C);
void AC2B(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& C, toolbox::FMatrix<double>& B);

void sqrt(const toolbox::FMatrix<double>& M2, toolbox::FMatrix<double>& M);
void build_B(const std::valarray<ThermoOptions>& thops, toolbox::FMatrix<double>& B);
void build_A(const std::valarray<ThermoOptions>& thops, toolbox::FMatrix<double>& A);
void get_TS(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& B, const double& t,
         toolbox::FMatrix<double>& T, toolbox::FMatrix<double>& S);
void A2Kt(const toolbox::FMatrix<double>& A, double dt, std::valarray<double>& K);

void time_acf(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& C, double t, toolbox::FMatrix<double>& Ct);
void int_acf(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& C, toolbox::FMatrix<double>& T);
void get_Xab(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& C, unsigned long a, unsigned long b, toolbox::FMatrix<double>& Xkl);
double int_2acf(const toolbox::FMatrix<double>& Xij, const toolbox::FMatrix<double>& Xkl,
                unsigned long i, unsigned long j, unsigned long k, unsigned long l);
double time_2acf(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& C, double t,
                 unsigned long i, unsigned long j, unsigned long k, unsigned long l);

void AD2tau_w(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& D, const double w,
              double& tau_q2, double& tau_p2, double& tau_h, double& q2, double &p2, double &pq);
void FTKH(const toolbox::FMatrix<double>& A, const toolbox::FMatrix<double>& C,  double w, toolbox::FMatrix<double>& iA, double &Kw, double& Hw);

#endif //ends #define __COLOR_HPP
