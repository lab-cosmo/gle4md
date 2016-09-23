#ifndef __GLEF_LIB_HPP
#define __GLEF_LIB_HPP 1
#include "color.hpp"
#include "matrix-fun.hpp"
#include "minsearch.hpp"
#include <vector>
#include <map>


enum GLEFDir    { Equal, Greater, Smaller };
enum GLEFMetric { Linear, Logarithmic, Exponential };
/* 4MR */
enum GLEFPointType { TauP2, TauQ2, TauH, KP2, KQ2, KH, Hw, Kw, HonK, DwQ, DwP, rDwQ, rDwP, LFP, PI2, PI10, PI100, RPImPole, RPRePole, RPQRes, RPPRes, RPQimv, RPQavg, RPQspr, RPPavg, RPPimv, RPPspr, CqqDT, CppDT, Cqq, Cpp, PTMin=TauP2, PTMax=Cpp+1 };
enum GLEFGlobType  { GZero, GInf, TZero, TInf, ARatio, ACondNum, CRatio, CCondNum, DCondNum, PQRatio, AEigvMax, AEigvSpread, AEigvCenter, AEigvWeight, AECondNum, DeltaSpread, DeltaWeight };
enum GLEFSearchMode { Annealing, Simplex, Powell };
enum GLEFParStyleA { ARealOnly, AComplex, AFull, APositive, APReal, APGeneral, ADelta };
enum GLEFParStyleC { COne, CPositive, CIndirect, CBDiagonal, CDelta };

typedef struct {
    GLEFDir dir;
    GLEFMetric metric;
    double y, w, e;
} GLEFValue;

typedef struct {
    double x;
    std::map<GLEFPointType,GLEFValue> values;
} GLEFPoint;

class GLEFFitOptions
{
    public:
    std::vector<GLEFPoint> points;
    std::map<GLEFPointType,double> pexp;
    std::map<GLEFPointType,double> pweights;
    std::map<GLEFGlobType,GLEFValue> globs;
};

class GLEFParOptions
{
  public:
    unsigned long ns;
    double deltat; double rpomega; double rpalpha; /* 4MR */
    GLEFParStyleA pstyleA; GLEFParStyleC pstyleC;
};

class GLEFSearchOptions
{
public:
    unsigned long steps, nlinesrc;
    GLEFSearchMode mode;

    //annealing options
    double step, ti, tf, tol, rand;
    double adapt_mult;
};

std::ostream& operator<< (std::ostream& ostr, const GLEFDir& p);
std::ostream& operator<< (std::ostream& ostr, const GLEFMetric& p);
std::ostream& operator<< (std::ostream& ostr, const GLEFPointType& p);
std::ostream& operator<< (std::ostream& ostr, const GLEFGlobType& p);
std::ostream& operator<< (std::ostream& ostr, const GLEFParStyleA& p);
std::ostream& operator<< (std::ostream& ostr, const GLEFParStyleC& p);
std::ostream& operator<< (std::ostream& ostr, const GLEFSearchMode& p);
std::istream& operator>> (std::istream& istr, GLEFDir& p);
std::istream& operator>> (std::istream& istr, GLEFMetric& p);
std::istream& operator>> (std::istream& istr, GLEFPointType& p);
std::istream& operator>> (std::istream& istr, GLEFGlobType& p);
std::istream& operator>> (std::istream& istr, GLEFParStyleA& p);
std::istream& operator>> (std::istream& istr, GLEFParStyleC& p);
std::istream& operator>> (std::istream& istr, GLEFSearchMode& p);

namespace toolbox {
__MK_IT_IOFIELD(GLEFDir);
__MK_IT_IOFIELD(GLEFMetric);
__MK_IT_IOFIELD(GLEFPointType);
__MK_IT_IOFIELD(GLEFGlobType);
__MK_IT_IOFIELD(GLEFParStyleA);
__MK_IT_IOFIELD(GLEFParStyleC);
__MK_IT_IOFIELD(GLEFSearchMode);
};

/*
namespace toolbox {
    template<> inline bool toolbox::IField<GLEFDir>::operator<<(std::istream& istr);
    template<> inline bool toolbox::IField<GLEFMetric>::operator<<(std::istream& istr);
    template<> inline bool toolbox::IField<GLEFPointType>::operator<<(std::istream& istr);
    template<> inline bool toolbox::IField<GLEFGlobType>::operator<<(std::istream& istr);
    template<> inline bool toolbox::IField<GLEFSearchMode>::operator<<(std::istream& istr);
    template<> inline bool toolbox::IField<GLEFParStyleA>::operator<<(std::istream& istr);
    template<> inline bool toolbox::IField<GLEFParStyleC>::operator<<(std::istream& istr);
}*/


std::istream& operator>> (std::istream& istr, GLEFValue& p);
std::ostream& operator<< (std::ostream& ostr, const GLEFValue& p);

std::istream& operator>> (std::istream& istr, GLEFPoint& p);
std::ostream& operator<< (std::ostream& ostr, const GLEFPoint& p);


std::istream& operator>> (std::istream& istr, GLEFSearchOptions& p);
std::ostream& operator<< (std::ostream& ostr, const GLEFSearchOptions& p);

std::ostream& operator<< (std::ostream& ostr, const GLEFFitOptions& op);
std::istream& operator>> (std::istream& istr, GLEFFitOptions& oo);

std::ostream& operator<< (std::ostream& ostr, const GLEFParOptions& op);
std::istream& operator>> (std::istream& istr, GLEFParOptions& oo);

#define MAX(a,b) ((a)<=(b)?(b):(a))
unsigned long npars(const GLEFParOptions& op);
void init_pars(const GLEFParOptions& op, std::valarray<double>& ip);
//void resize_pars(const GLEFFitOptions& fop, unsigned long osz, std::valarray<double>& ip);
class GLEFError
{
    private:
        unsigned long np;
    public:
        GLEABC abc;
        void compute_globs(std::map<GLEFGlobType, double>& lims);
        void compute_points(const std::vector<double>& xp, std::vector<std::map<GLEFPointType, double> >& val);

        toolbox::nullstream *ogarbage;
        std::ostream *slog, *seva, *spars;
        toolbox::FMatrix<double> A, C;
        std::valarray<double> p;

        GLEFFitOptions ofit;
        GLEFParOptions opar;

        GLEFError() : ogarbage(new toolbox::nullstream) {slog=ogarbage; spars=ogarbage; seva=ogarbage; }

        void set_ops(const GLEFFitOptions& nof, const GLEFParOptions& nop)
        {
            ofit=nof; opar=nop; p.resize(np=npars(opar)); p=0.;
            A.resize(opar.ns+1,opar.ns+1); A*=0.; C=A;
        }

        void pars2AC();
        void AC2pars();

        inline unsigned long size() { return opar.ns; }

        inline void set_vars(const std::valarray<double>& rv)
        {
            if (rv.size()!=np) ERROR("Wrong size for var list.");
            p=rv;  for (int i=0; i<p.size(); ++i) *spars<<p[i]<<" "; *spars<<std::endl;
            pars2AC();
        }

        inline void get_vars(std::valarray<double>& rv) const
        {
            if (rv.size()!=np) rv.resize(np); rv=p;
        }

        void get_value(double& rv);

        inline void get_gradient(std::valarray<double>& rv) const
        { ERROR("Minimization function class does not define a gradient function"); }
        inline void get_hessian(std::valarray<std::valarray<double> >& rv) const
        { ERROR("Minimization function class does not define a hessian function"); }
        bool bailout() { return false; }
};

#endif //ends #ifndef __GLED_LIB_HPP
