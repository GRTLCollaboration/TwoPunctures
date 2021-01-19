#ifndef __TwoPunctures_Utilities__
#define __TwoPunctures_Utilities__

/* TwoPunctures:  File  "utilities.h"*/
/* these functions have no dependency on the TP parameters */

#include "TwoPunctures.hpp"
#include <math.h>

namespace TP
{
namespace Utilities
{

#define Pi 3.14159265358979323846264338328
#define Pih 1.57079632679489661923132169164 /* Pi/2*/
#define Piq 0.78539816339744830961566084582 /* Pi/4*/

#define TINY 1.0e-20
#define SWAP(a, b)                                                             \
    {                                                                          \
        double temp = (a);                                                     \
        (a) = (b);                                                             \
        (b) = temp;                                                            \
    }

void nrerror(char error_text[]);
int *ivector(long nl, long nh);
double *dvector(long nl, long nh);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_ivector(int *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh);

int minimum2(int i, int j);
int minimum3(int i, int j, int k);
int maximum2(int i, int j);
int maximum3(int i, int j, int k);
int pow_int(int mantisse, int exponent);

void chebft_Zeros(double u[], int n, int inv);
void chebft_Extremes(double u[], int n, int inv);
void chder(double *c, double *cder, int n);
double chebev(double a, double b, double c[], int m, double x);
void fourft(double *u, int N, int inv);
void fourder(double u[], double du[], int N);
void fourder2(double u[], double d2u[], int N);
double fourev(double *u, int N, double x);

double norm1(double *v, int n);
double norm2(double *v, int n);
double scalarproduct(double *v, double *w, int n);

/*
typedef nrerror TP_nrerror;
typedef ivector TP_ivector;
typedef dvector TP_dvector;
typedef imatrix TP_imatrix;
typedef dmatrix TP_dmatrix;
typedef d3tensor TP_d3tensor;
typedef free_ivector TP_free_ivector;
typedef free_dvector TP_free_dvector;
typedef free_imatrix TP_free_imatrix;
typedef free_dmatrix TP_free_dmatrix;
typedef free_d3tensor TP_free_d3tensor;
*/

} // namespace Utilities
} // namespace TP

#endif /* __TwoPunctures_Utilities__ */
