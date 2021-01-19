#ifndef __TP_CONSERVED_VECTOR_NAMES__
#define __TP_CONSERVED_VECTOR_NAMES__

/**
 * This file contains aliases for making access to the long state vector Q
 * as used eg. in ExaHyPE and at the functional interface of this code,
 * more readable.
 *
 * With these constants, you can write Q[g11], Q[K22] and Q[B33] instead of
 * Q[0], Q[9] and Q[34]. Your code is again meaningful in terms of physics.
 *
 * The list os converted from the PDENames function from Trentos Fortran code
 * with this expression:
 *
 * > perl -p -i -e "s/^\s+Name\((.*)\)\s* = '(.+)'/\"const int \$2 =
 * \".(\$1-1).';';/eg"  FileHoldingPDEDefs.tmp
 *
 * Methods have been implemented already in the ExaHyPE SRHD application to
 * deal even more naturally with physical expressions which are stored in
 * a long and unintuitive vector. See end of this file for comments.
 */

namespace TP
{
namespace Z4VectorShortcuts
{

const int g11 = 0;
const int g12 = 1;
const int g13 = 2;
const int g22 = 3;
const int g23 = 4;
const int g33 = 5;
const int K11 = 6;
const int K12 = 7;
const int K13 = 8;
const int K22 = 9;
const int K23 = 10;
const int K33 = 11;
const int Z1 = 12;
const int Z2 = 13;
const int Z3 = 14;
const int Theta = 15;
const int lapse = 16;
const int shift1 = 17;
const int shift2 = 18;
const int shift3 = 19;
const int b1 = 20;
const int b2 = 21;
const int b3 = 22;
const int A1 = 23;
const int A2 = 24;
const int A3 = 25;
const int B11 = 26;
const int B21 = 27;
const int B31 = 28;
const int B12 = 29;
const int B22 = 30;
const int B32 = 31;
const int B13 = 32;
const int B23 = 33;
const int B33 = 34;
const int D111 = 35;
const int D112 = 36;
const int D113 = 37;
const int D122 = 38;
const int D123 = 39;
const int D133 = 40;
const int D211 = 41;
const int D212 = 42;
const int D213 = 43;
const int D222 = 44;
const int D223 = 45;
const int D233 = 46;
const int D311 = 47;
const int D312 = 48;
const int D313 = 49;
const int D322 = 50;
const int D323 = 51;
const int D333 = 52;
const int K0 = 53;

// total length of conserved vector
const int Qlen = 53;

/************************** Syntactic sugar (example from SRHD)
 * ********************************/
/************************** Just for placing it somewhere, currently not used
 * ******************/

// Method 1: C++ namespaced constants, as present here.

// Method 2:
// Labels for indices into Q and V, used only in the local scope. However, the
// C++ namespace allows the same just with a cleaner notation.
///  #define DEFINE_LABELS \
///  	const int rho = 0; \
///  	const int vx = 1; \
///  	const int vy = 2; \
///  	const int vz = 3; \
///  	const int p = 4
/* Usage:
     function foo(double* Q, double* V) {
        DEFINE_LABELS;
        Q[rho] = V[p];
     }
   You can archieve the same when writing
     function foo(double* Q, double* V) {
        using namespace TP::Z4VectorShortcuts;
        Q[rho] = V[p];
     }
*/

// direct references deep into X  (might be named V, Q, S, whatever)
// Pro: More natural, shorter names for fields
// Contra: Can only be used for one state vector in a scope.
//         Will produce a lot of warnings when -Wunused-but-set-variable
//         warnings is active
///  #define SHORTHANDS(X)   double &rho=X[0], &vx=X[1], &vy=X[2], &vz=X[3],
///  &p=X[4]
/* Usage:
   function foo(double* Q) {
      SHORTHANDS(Q);
      rho = p*p + vx;
   }
*/

} // namespace Z4VectorShortcuts
} // namespace TP
#endif /* __TP_CONSERVED_VECTOR_NAMES__ */
