/* TwoPunctures:  File  "TwoPunctures.h"*/
#ifndef __TWO_PUNCTURES_MAIN__
#define __TWO_PUNCTURES_MAIN__

// Chombo includes
#include "MayDay.H"
#include "SPMD.H"
#include "parstream.H"

// Our includes
#include "TP_ConservedVectorNames.hpp"
#include "TP_Parameters.hpp"
#include "TP_Utilities.hpp"
#include <cstdio>

// Chombo namespace
#include "UsingNamespace.H"

namespace TP
{

/** We compile TwoPunctures-Standalone as C++, and C++ does not yet
    now about the C99 language feature "restrict". If you compile
    with GCC, you can however use the following. Else just leave
    it empty.

    In any case, restrict is a keyword to improve performance at
    array access.
**/
#define TP_RESTRICT __restrict__

/* Give MPI information about rank */
inline int TP_MyProc() { return procID(); }

#define TP_INFO(args...)                                                       \
    {                                                                          \
        pout() << "TP: ";                                                      \
        char buf[1024];                                                        \
        sprintf(buf, args);                                                    \
        pout() << buf << std::endl;                                            \
    }
#define TP_ERROR MayDay::Error
#define TP_WARN MayDay::Warning

/* end of quick adaptions */

typedef struct DERIVS
{
    double *d0, *d1, *d2, *d3, *d11, *d12, *d13, *d22, *d23, *d33;
} derivs;

enum GRID_SETUP_METHOD
{
    GSM_Taylor_expansion,
    GSM_evaluation
};

/*
Files of "TwoPunctures":
        TwoPunctures.c
        FuncAndJacobian.c
        CoordTransf.c
        Equations.c
        Newton.c
        utilities.c (see utilities.h)
**************************
*/

class TwoPunctures : public Parameters
{
  public:
    // Work variables are defined in this class.
    // Parameters go to the class TwoPunctures::Parameters.

    // shared variables between TP Setup and TP Interpolation
    GRID_SETUP_METHOD gsm;
    int nvar, n1, n2, n3;
    derivs u, v, cf_v;
    int antisymmetric_lapse, averaged_lapse, pmn_lapse, brownsville_lapse;

    /* Consructor */

    TwoPunctures()
    {
        nvar = 0;
        n1 = 0;
        n2 = 0;
        n3 = 0;
    }

    /* Routines in  "TwoPunctures.c"*/
    void set_initial_guess(derivs v);
    double TestSolution(double A, double B, double X, double R, double phi);
    void TestVector_w(double *par, int nvar, int n1, int n2, int n3, double *w);
    void Run();
    void Interpolate(const double *pos, double *Q) const;

    /* Routines in  "FuncAndJacobian.c"*/
    int Index(int ivar, int i, int j, int k, int nvar, int n1, int n2,
              int n3) const;
    void allocate_derivs(derivs *v, int n) const;
    void free_derivs(derivs *v, int n) const;
    void Derivatives_AB3(int nvar, int n1, int n2, int n3, derivs v);
    void F_of_v(int nvar, int n1, int n2, int n3, derivs v, double *F,
                derivs u);
    void J_times_dv(int nvar, int n1, int n2, int n3, derivs dv, double *Jdv,
                    derivs u);
    void JFD_times_dv(int i, int j, int k, int nvar, int n1, int n2, int n3,
                      derivs dv, derivs u, double *values);
    void SetMatrix_JFD(int nvar, int n1, int n2, int n3, derivs u, int *ncols,
                       int **cols, double **Matrix);
    double PunctEvalAtArbitPosition(double *v, int ivar, double A, double B,
                                    double phi, int nvar, int n1, int n2,
                                    int n3);
    void calculate_derivs(int i, int j, int k, int ivar, int nvar, int n1,
                          int n2, int n3, derivs v, derivs vv) const;
    double interpol(double a, double b, double c, derivs v) const;
    double PunctTaylorExpandAtArbitPosition(int ivar, int nvar, int n1, int n2,
                                            int n3, derivs v, double x,
                                            double y, double z) const;
    double PunctIntPolAtArbitPosition(int ivar, int nvar, int n1, int n2,
                                      int n3, derivs v, double x, double y,
                                      double z);
    void SpecCoef(int n1, int n2, int n3, int ivar, double *v, double *cf);
    double PunctEvalAtArbitPositionFast(double *v, int ivar, double A, double B,
                                        double phi, int nvar, int n1, int n2,
                                        int n3) const;
    double PunctIntPolAtArbitPositionFast(int ivar, int nvar, int n1, int n2,
                                          int n3, derivs v, double x, double y,
                                          double z) const;

    /* Routines in  "CoordTransf.c"*/
    void AB_To_XR(int nvar, double A, double B, double *X, double *R, derivs U);
    void C_To_c(int nvar, double X, double R, double *x, double *r, derivs U);
    void rx3_To_xyz(int nvar, double x, double r, double phi, double *y,
                    double *z, derivs U);

    /* Routines in  "Equations.c"*/
    double BY_KKofxyz(double x, double y, double z);
    void BY_Aijofxyz(double x, double y, double z, double Aij[3][3]) const;
    void NonLinEquations(double rho_adm, double A, double B, double X, double R,
                         double x, double r, double phi, double y, double z,
                         derivs U, double *values);
    void LinEquations(double A, double B, double X, double R, double x,
                      double r, double phi, double y, double z, derivs dU,
                      derivs U, double *values);

    /* Routines in  "Newton.c"*/
    void TestRelax(int nvar, int n1, int n2, int n3, derivs v, double *dv);
    void Newton(int nvar, int n1, int n2, int n3, derivs v, double tol,
                int itmax);

    /* Former static routines in "Newton.c"*/
    int bicgstab(int const nvar, int const n1, int const n2, int const n3,
                 derivs v, derivs dv, int const output, int const itmax,
                 double const tol, double *TP_RESTRICT const normres);
    double norm_inf(double const *TP_RESTRICT const F, int const ntotal);
    void relax(double *TP_RESTRICT const dv, int const nvar, int const n1,
               int const n2, int const n3, double const *TP_RESTRICT const rhs,
               int const *TP_RESTRICT const ncols,
               int const *TP_RESTRICT const *TP_RESTRICT const cols,
               double const *TP_RESTRICT const *TP_RESTRICT const JFD);
    void resid(double *TP_RESTRICT const res, int const ntotal,
               double const *TP_RESTRICT const dv,
               double const *TP_RESTRICT const rhs,
               int const *TP_RESTRICT const ncols,
               int const *TP_RESTRICT const *TP_RESTRICT const cols,
               double const *TP_RESTRICT const *TP_RESTRICT const JFD);
    void LineRelax_al(double *TP_RESTRICT const dv, int const j, int const k,
                      int const nvar, int const n1, int const n2, int const n3,
                      double const *TP_RESTRICT const rhs,
                      int const *TP_RESTRICT const ncols,
                      int const *TP_RESTRICT const *TP_RESTRICT const cols,
                      double const *TP_RESTRICT const *TP_RESTRICT const JFD);
    void LineRelax_be(double *TP_RESTRICT const dv, int const i, int const k,
                      int const nvar, int const n1, int const n2, int const n3,
                      double const *TP_RESTRICT const rhs,
                      int const *TP_RESTRICT const ncols,
                      int const *TP_RESTRICT const *TP_RESTRICT const cols,
                      double const *TP_RESTRICT const *TP_RESTRICT const JFD);

}; // class TwoPunctures

} // namespace TP

#endif /* __TWO_PUNCTURES_MAIN__ */
