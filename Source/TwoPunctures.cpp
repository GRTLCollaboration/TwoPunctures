/* TwoPunctures:  File  "TwoPunctures.c"*/

#include "TwoPunctures.hpp"
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace TP
{
using namespace Utilities;
using namespace Z4VectorShortcuts;

void TwoPunctures::set_initial_guess(derivs v)
{
    int nvar = 1;

    double *s_x, *s_y, *s_z;
    double al, A, Am1, be, B, phi, R, r, X;
    int ivar, i, j, k, i3D, indx;
    derivs U;
    FILE *debug_file;

    if (solve_momentum_constraint)
        nvar = 4;

    s_x = new double[n1 * n2 * n3]();
    s_y = new double[n1 * n2 * n3]();
    s_z = new double[n1 * n2 * n3]();
    allocate_derivs(&U, nvar);
    for (ivar = 0; ivar < nvar; ivar++)
        for (i = 0; i < n1; i++)
            for (j = 0; j < n2; j++)
                for (k = 0; k < n3; k++)
                {
                    i3D = Index(ivar, i, j, k, 1, n1, n2, n3);

                    al = Pih * (2 * i + 1) / n1;
                    A = -cos(al);
                    be = Pih * (2 * j + 1) / n2;
                    B = -cos(be);
                    phi = 2. * Pi * k / n3;

                    /* Calculation of (X,R)*/
                    AB_To_XR(nvar, A, B, &X, &R, U);
                    /* Calculation of (x,r)*/
                    C_To_c(nvar, X, R, &(s_x[i3D]), &r, U);
                    /* Calculation of (y,z)*/
                    rx3_To_xyz(nvar, s_x[i3D], r, phi, &(s_y[i3D]), &(s_z[i3D]),
                               U);
                }
    // @TODO where is the function Set_Initial_Guess_for_u implemented!???
    // When looking in EinsteinInitialData, I find it at
    /*
  TOVSolver/interface.ccl
  35:CCTK_INT FUNCTION Set_Initial_Guess_for_u(  \
  43:PROVIDES FUNCTION Set_Initial_Guess_for_u \
  44:         WITH TOV_Set_Initial_Guess_for_u \
    */
    //// Set_Initial_Guess_for_u(n1*n2*n3, v.d0, s_x, s_y, s_z);
    for (ivar = 0; ivar < nvar; ivar++)
        for (i = 0; i < n1; i++)
            for (j = 0; j < n2; j++)
                for (k = 0; k < n3; k++)
                {
                    indx = Index(ivar, i, j, k, 1, n1, n2, n3);
                    v.d0[indx] /= (-cos(Pih * (2 * i + 1) / n1) - 1.0);
                }
    Derivatives_AB3(nvar, n1, n2, n3, v);
    if (do_initial_debug_output && TP_MyProc() == 0)
    {
        debug_file = fopen("initial.dat", "w");
        assert(debug_file);
        for (ivar = 0; ivar < nvar; ivar++)
            for (i = 0; i < n1; i++)
                for (j = 0; j < n2; j++)
                {
                    al = Pih * (2 * i + 1) / n1;
                    A = -cos(al);
                    Am1 = A - 1.0;
                    be = Pih * (2 * j + 1) / n2;
                    B = -cos(be);
                    phi = 0.0;
                    indx = Index(ivar, i, j, 0, 1, n1, n2, n3);
                    U.d0[0] = Am1 * v.d0[indx];                    /* U*/
                    U.d1[0] = v.d0[indx] + Am1 * v.d1[indx];       /* U_A*/
                    U.d2[0] = Am1 * v.d2[indx];                    /* U_B*/
                    U.d3[0] = Am1 * v.d3[indx];                    /* U_3*/
                    U.d11[0] = 2 * v.d1[indx] + Am1 * v.d11[indx]; /* U_AA*/
                    U.d12[0] = v.d2[indx] + Am1 * v.d12[indx];     /* U_AB*/
                    U.d13[0] = v.d3[indx] + Am1 * v.d13[indx];     /* U_AB*/
                    U.d22[0] = Am1 * v.d22[indx];                  /* U_BB*/
                    U.d23[0] = Am1 * v.d23[indx];                  /* U_B3*/
                    U.d33[0] = Am1 * v.d33[indx];                  /* U_33*/
                    /* Calculation of (X,R)*/
                    AB_To_XR(nvar, A, B, &X, &R, U);
                    /* Calculation of (x,r)*/
                    C_To_c(nvar, X, R, &(s_x[indx]), &r, U);
                    /* Calculation of (y,z)*/
                    rx3_To_xyz(nvar, s_x[i3D], r, phi, &(s_y[indx]),
                               &(s_z[indx]), U);
                    fprintf(debug_file,
                            "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g "
                            "%.16g %.16g "
                            "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g "
                            "%.16g\n",
                            (double)s_x[indx], (double)s_y[indx], (double)A,
                            (double)B, (double)U.d0[0],
                            (double)(-cos(Pih * (2 * i + 1) / n1) - 1.0),
                            (double)U.d1[0], (double)U.d2[0], (double)U.d3[0],
                            (double)U.d11[0], (double)U.d22[0],
                            (double)U.d33[0], (double)v.d0[indx],
                            (double)v.d1[indx], (double)v.d2[indx],
                            (double)v.d3[indx], (double)v.d11[indx],
                            (double)v.d22[indx], (double)v.d33[indx]);
                }
        fprintf(debug_file, "\n\n");
        for (i = n2 - 10; i < n2; i++)
        {
            double d;
            indx = Index(0, 0, i, 0, 1, n1, n2, n3);
            d = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, s_x[indx],
                                           0.0, 0.0);
            fprintf(debug_file, "%.16g %.16g\n", (double)s_x[indx], (double)d);
        }
        fprintf(debug_file, "\n\n");
        for (i = n2 - 10; i < n2 - 1; i++)
        {
            double d;
            int ip = Index(0, 0, i + 1, 0, 1, n1, n2, n3);
            indx = Index(0, 0, i, 0, 1, n1, n2, n3);
            for (j = -10; j < 10; j++)
            {
                d = PunctIntPolAtArbitPosition(
                    0, nvar, n1, n2, n3, v,
                    s_x[indx] + (s_x[ip] - s_x[indx]) * j / 10, 0.0, 0.0);
                fprintf(debug_file, "%.16g %.16g\n",
                        (double)(s_x[indx] + (s_x[ip] - s_x[indx]) * j / 10),
                        (double)d);
            }
        }
        fprintf(debug_file, "\n\n");
        for (i = 0; i < n1; i++)
            for (j = 0; j < n2; j++)
            {
                X = 2 * (2.0 * i / n1 - 1.0);
                R = 2 * (1.0 * j / n2);
                if (X * X + R * R > 1.0)
                {
                    C_To_c(nvar, X, R, &(s_x[indx]), &r, U);
                    rx3_To_xyz(nvar, s_x[i3D], r, phi, &(s_y[indx]),
                               &(s_z[indx]), U);
                    *U.d0 = s_x[indx] * s_x[indx];
                    *U.d1 = 2 * s_x[indx];
                    *U.d2 = 0.0;
                    *U.d3 = 0.0;
                    *U.d11 = 2.0;
                    *U.d22 = 0.0;
                    *U.d33 = *U.d12 = *U.d23 = *U.d13 = 0.0;
                    C_To_c(nvar, X, R, &(s_x[indx]), &r, U);
                    fprintf(debug_file,
                            "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g "
                            "%.16g %.16g %.16g\n",
                            (double)s_x[indx], (double)r, (double)X, (double)R,
                            (double)U.d0[0], (double)U.d1[0], (double)U.d2[0],
                            (double)U.d3[0], (double)U.d11[0], (double)U.d22[0],
                            (double)U.d33[0]);
                }
            }
        fclose(debug_file);
    }
    free(s_z);
    free(s_y);
    free(s_x);
    free_derivs(&U, nvar);
    /*exit(0);*/
}

/* -------------------------------------------------------------------*/
void TwoPunctures::Run()
{

    nvar = 1;
    n1 = npoints_A;
    n2 = npoints_B;
    n3 = npoints_phi;

    mp = par_m_plus;
    mm = par_m_minus;

    int imin[3], imax[3];
    int const ntotal = n1 * n2 * n3 * nvar;
#if 0
  int percent10 = 0;
#endif
    static double *F = NULL;
    double admMass;

    if (!F)
    {
        double up, um;
        /* Solve only when called for the first time */
        F = dvector(0, ntotal - 1);
        allocate_derivs(&u, ntotal);
        allocate_derivs(&v, ntotal);
        allocate_derivs(&cf_v, ntotal);

        if (use_sources)
        {
            TP_INFO("Solving puncture equation for BH-NS/NS-NS system");
        }
        else
        {
            TP_INFO("Solving puncture equation for BH-BH system");
        }
        TP_INFO("b = %g", par_b);

        /* initialise to 0 */
        for (int j = 0; j < ntotal; j++)
        {
            cf_v.d0[j] = 0.0;
            cf_v.d1[j] = 0.0;
            cf_v.d2[j] = 0.0;
            cf_v.d3[j] = 0.0;
            cf_v.d11[j] = 0.0;
            cf_v.d12[j] = 0.0;
            cf_v.d13[j] = 0.0;
            cf_v.d22[j] = 0.0;
            cf_v.d23[j] = 0.0;
            cf_v.d33[j] = 0.0;
            v.d0[j] = 0.0;
            v.d1[j] = 0.0;
            v.d2[j] = 0.0;
            v.d3[j] = 0.0;
            v.d11[j] = 0.0;
            v.d12[j] = 0.0;
            v.d13[j] = 0.0;
            v.d22[j] = 0.0;
            v.d23[j] = 0.0;
            v.d33[j] = 0.0;
        }
        /* call for external initial guess */
        if (use_external_initial_guess)
        {
            set_initial_guess(v);
        }

        /* If bare masses are not given, iteratively solve for them given the
           target ADM masses target_M_plus and target_M_minus and with initial
           guesses given by par_m_plus and par_m_minus. */
        if (!(give_bare_mass))
        {
            double tmp, mp_adm_err, mm_adm_err;
            char valbuf[100];

            double M_p = target_M_plus;
            double M_m = target_M_minus;

            TP_INFO("Attempting to find bare masses.");
            TP_INFO("Target ADM masses: M_p=%g and M_m=%g", (double)M_p,
                    (double)M_m);
            TP_INFO("ADM mass tolerance: %g", (double)adm_tol);

            /* Loop until both ADM masses are within adm_tol of their target */
            do
            {
                TP_INFO("Bare masses: mp=%.15g, mm=%.15g", (double)mp,
                        (double)mm);
                Newton(nvar, n1, n2, n3, v, Newton_tol, 1);

                F_of_v(nvar, n1, n2, n3, v, F, u);

                up = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, par_b,
                                                0., 0.);
                um = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, -par_b,
                                                0., 0.);

                /* Calculate the ADM masses from the current bare mass guess */
                mp_adm = (1 + up) * mp + mp * mm / (4. * par_b);
                mm_adm = (1 + um) * mm + mp * mm / (4. * par_b);

                /* Check how far the current ADM masses are from the target */
                mp_adm_err = fabs(M_p - mp_adm);
                mm_adm_err = fabs(M_m - mm_adm);
                TP_INFO("ADM mass error: M_p_err=%.15g, M_m_err=%.15g",
                        (double)mp_adm_err, (double)mm_adm_err);

                /* Invert the ADM mass equation and update the bare mass guess
                   so that it gives the correct target ADM masses */
                tmp =
                    -4 * par_b * (1 + um + up + um * up) +
                    sqrt(16 * par_b * M_m * (1 + um) * (1 + up) +
                         pow(-M_m + M_p + 4 * par_b * (1 + um) * (1 + up), 2));
                mp = (tmp + M_p - M_m) / (2. * (1 + up));
                mm = (tmp - M_p + M_m) / (2. * (1 + um));

                /* Set the par_m_plus and par_m_minus parameters */
                par_m_plus = mp;
                par_m_minus = mm;

            } while ((mp_adm_err > adm_tol) || (mm_adm_err > adm_tol));

            TP_INFO("Found bare masses.");
        }

        Newton(nvar, n1, n2, n3, v, Newton_tol, Newton_maxit);

        F_of_v(nvar, n1, n2, n3, v, F, u);

        SpecCoef(n1, n2, n3, 0, v.d0, cf_v.d0);

        TP_INFO("The two puncture masses are mp=%.17g and mm=%.17g", (double)mp,
                (double)mm);

        up = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, par_b, 0., 0.);
        um = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, -par_b, 0., 0.);

        /* Calculate the ADM masses from the current bare mass guess */
        mp_adm = (1 + up) * mp + mp * mm / (4. * par_b);
        mm_adm = (1 + um) * mm + mp * mm / (4. * par_b);

        TP_INFO("Puncture 1 ADM mass is %.17g", (double)mp_adm);
        TP_INFO("Puncture 2 ADM mass is %.17g", (double)mm_adm);

        /* print out ADM mass, eq.: \Delta M_ADM=2*r*u=4*b*V for A=1,B=0,phi=0
         */
        admMass =
            (mp + mm -
             4 * par_b *
                 PunctEvalAtArbitPosition(v.d0, 0, 1, 0, 0, nvar, n1, n2, n3));
        TP_INFO("The total ADM mass is %.17g", (double)admMass);

        /*
          Run this in Mathematica (version 8 or later) with
            math -script <file>

          Needs["SymbolicC`"];
          co = Table["center_offset[" <> ToString[i] <> "]", {i, 0, 2}];
          r1 = co + {"par_b", 0, 0};
          r2 = co + {-"par_b", 0, 0};
          {p1, p2} = Table["par_P_" <> bh <> "[" <> ToString[i] <> "]", {bh,
          {"plus", "minus"}}, {i, 0, 2}]; {s1, s2} = Table["par_S_" <> bh <> "["
          <> ToString[i] <> "]", {bh, {"plus", "minus"}}, {i, 0, 2}];

          J = Cross[r1, p1] + Cross[r2, p2] + s1 + s2;

          JVar = Table["*J" <> ToString[i], {i, 1, 3}];
          Print[OutputForm@StringReplace[
            ToCCodeString@MapThread[CAssign[#1, CExpression[#2]] &, {JVar, J}],
            "\"" -> ""]];
         */

        J1 = -(center_offset[2] * par_P_minus[1]) +
             center_offset[1] * par_P_minus[2] -
             center_offset[2] * par_P_plus[1] +
             center_offset[1] * par_P_plus[2] + par_S_minus[0] + par_S_plus[0];
        J2 = center_offset[2] * par_P_minus[0] -
             center_offset[0] * par_P_minus[2] + par_b * par_P_minus[2] +
             center_offset[2] * par_P_plus[0] -
             center_offset[0] * par_P_plus[2] - par_b * par_P_plus[2] +
             par_S_minus[1] + par_S_plus[1];
        J3 = -(center_offset[1] * par_P_minus[0]) +
             center_offset[0] * par_P_minus[1] - par_b * par_P_minus[1] -
             center_offset[1] * par_P_plus[0] +
             center_offset[0] * par_P_plus[1] + par_b * par_P_plus[1] +
             par_S_minus[2] + par_S_plus[2];
    }

    if (grid_setup_method == "Taylor expansion")
    {
        gsm = GSM_Taylor_expansion;
    }
    else if (grid_setup_method == "evaluation")
    {
        gsm = GSM_evaluation;
    }
    else
    {
        TP_WARN("internal error");
    }

    if (initial_lapse == "twopunctures")
        TP_WARN("Please specify a lapse which we can use");
    antisymmetric_lapse = (initial_lapse == "twopunctures-antisymmetric");
    averaged_lapse = (initial_lapse == "twopunctures-averaged");
    pmn_lapse = (initial_lapse == "psi^n");
    if (pmn_lapse)
        TP_INFO("Setting initial lapse to psi^%f profile.",
                (double)initial_lapse_psi_exponent);
    brownsville_lapse = (initial_lapse == "brownsville");
    if (brownsville_lapse)
        TP_INFO("Setting initial lapse to a Brownsville-style profile "
                "with exp %f.",
                (double)initial_lapse_psi_exponent);

    TP_INFO("Preparing interpolation of result");
    if (metric_type == "static conformal")
    {
        if (conformal_storage == "factor")
        {
            conformal_state = 1;
        }
        else if (conformal_storage == "factor+derivs")
        {
            conformal_state = 2;
        }
        else if (conformal_storage == "factor+derivs+2nd derivs")
        {
            conformal_state = 3;
        }
    }
    else
    {
        conformal_state = 0;
    }
} /* End of TwoPointures_Setup() */

/**
 * Interpolation function for an external caller.
 **/
void TwoPunctures::Interpolate(const double *pos, double *Q) const
{
    // First, zero everything
    // @TODO: it should be checked if this represents Vaccuum
    for (int i = 0; i < Qlen; i++)
        Q[i] = 0;

    double x = pos[0], y = pos[1], z = pos[2];
    double xx, yy, zz;
    xx = x - center_offset[0];
    yy = y - center_offset[1];
    zz = z - center_offset[2];

    /* We implement swapping the x and z coordinates as follows.
       The bulk of the code that performs the actual calculations
       is unchanged.  This code looks only at local variables.
       Before the bulk --i.e., here-- we swap all x and z tensor
       components, and after the code --i.e., at the end of this
       main loop-- we swap everything back.  */
    if (swap_xz)
    {
        /* Swap the x and z coordinates */
        SWAP(xx, zz);
    }

    double r_plus = sqrt(pow(xx - par_b, 2) + pow(yy, 2) + pow(zz, 2));
    double r_minus = sqrt(pow(xx + par_b, 2) + pow(yy, 2) + pow(zz, 2));

    double U;
    switch (gsm)
    {
    case GSM_Taylor_expansion:
        U = PunctTaylorExpandAtArbitPosition(0, nvar, n1, n2, n3, v, xx, yy,
                                             zz);
        break;
    case GSM_evaluation:
        U = PunctIntPolAtArbitPositionFast(0, nvar, n1, n2, n3, cf_v, xx, yy,
                                           zz);
        break;
    default:
        assert(0);
    }
    r_plus = pow(pow(r_plus, 4) + pow(TP_epsilon, 4), 0.25);
    r_minus = pow(pow(r_minus, 4) + pow(TP_epsilon, 4), 0.25);
    if (r_plus < TP_Tiny)
        r_plus = TP_Tiny;
    if (r_minus < TP_Tiny)
        r_minus = TP_Tiny;
    double psi1 = 1 + 0.5 * mp / r_plus + 0.5 * mm / r_minus + U;
#define EXTEND(M, r)                                                           \
    (M * (3. / 8 * pow(r, 4) / pow(TP_Extend_Radius, 5) -                      \
          5. / 4 * pow(r, 2) / pow(TP_Extend_Radius, 3) +                      \
          15. / 8 / TP_Extend_Radius))
    if (r_plus < TP_Extend_Radius)
    {
        psi1 = 1 + 0.5 * EXTEND(mp, r_plus) + 0.5 * mm / r_minus + U;
    }
    if (r_minus < TP_Extend_Radius)
    {
        psi1 = 1 + 0.5 * EXTEND(mm, r_minus) + 0.5 * mp / r_plus + U;
    }
    double static_psi = 1;

    double Aij[3][3];
    BY_Aijofxyz(xx, yy, zz, Aij);

    double old_alp = 1.0;
    if (multiply_old_lapse)
        TP_WARN("@todo: Old Lapse is not be possible, Q as input vector is "
                "undefined.");

    if ((conformal_state > 0) || (pmn_lapse) || (brownsville_lapse))
    {

        double xp, yp, zp, rp, ir;
        double s1, s3, s5;
        double p, px, py, pz, pxx, pxy, pxz, pyy, pyz, pzz;
        p = 1.0;
        px = py = pz = 0.0;
        pxx = pxy = pxz = 0.0;
        pyy = pyz = pzz = 0.0;

        /* first puncture */
        xp = xx - par_b;
        yp = yy;
        zp = zz;
        rp = sqrt(xp * xp + yp * yp + zp * zp);
        rp = pow(pow(rp, 4) + pow(TP_epsilon, 4), 0.25);
        if (rp < TP_Tiny)
            rp = TP_Tiny;
        ir = 1.0 / rp;

        if (rp < TP_Extend_Radius)
        {
            ir = EXTEND(1., rp);
        }

        s1 = 0.5 * mp * ir;
        s3 = -s1 * ir * ir;
        s5 = -3.0 * s3 * ir * ir;

        p += s1;

        px += xp * s3;
        py += yp * s3;
        pz += zp * s3;

        pxx += xp * xp * s5 + s3;
        pxy += xp * yp * s5;
        pxz += xp * zp * s5;
        pyy += yp * yp * s5 + s3;
        pyz += yp * zp * s5;
        pzz += zp * zp * s5 + s3;

        /* second puncture */
        xp = xx + par_b;
        yp = yy;
        zp = zz;
        rp = sqrt(xp * xp + yp * yp + zp * zp);
        rp = pow(pow(rp, 4) + pow(TP_epsilon, 4), 0.25);
        if (rp < TP_Tiny)
            rp = TP_Tiny;
        ir = 1.0 / rp;

        if (rp < TP_Extend_Radius)
        {
            ir = EXTEND(1., rp);
        }

        s1 = 0.5 * mm * ir;
        s3 = -s1 * ir * ir;
        s5 = -3.0 * s3 * ir * ir;

        p += s1;

        px += xp * s3;
        py += yp * s3;
        pz += zp * s3;

        pxx += xp * xp * s5 + s3;
        pxy += xp * yp * s5;
        pxz += xp * zp * s5;
        pyy += yp * yp * s5 + s3;
        pyz += yp * zp * s5;
        pzz += zp * zp * s5 + s3;

        // @TODO: psix, psiy, psiz, psixx, psixy are not defined, currently.
        // They were stored by Thorn einsteinbase/StaticConformal

        if (conformal_state >= 1)
        {
            static_psi = p;
            /// psi[ind] = static_psi;
        }
        if (conformal_state >= 2)
        {
            /// psix[ind] = px / static_psi;
            /// psiy[ind] = py / static_psi;
            /// psiz[ind] = pz / static_psi;
        }
        if (conformal_state >= 3)
        {
            /// psixx[ind] = pxx / static_psi;
            /// psixy[ind] = pxy / static_psi;
            /// psixz[ind] = pxz / static_psi;
            /// psiyy[ind] = pyy / static_psi;
            /// psiyz[ind] = pyz / static_psi;
            /// psizz[ind] = pzz / static_psi;
        }

        if (pmn_lapse)
            Q[lapse] = pow(p, initial_lapse_psi_exponent);
        if (brownsville_lapse)
            Q[lapse] = 2.0 / (1.0 + pow(p, initial_lapse_psi_exponent));

    } /* if conformal-state > 0 */

    /// puncture_u[ind] = U; /// @TODO: Also no storage for this

    Q[g11] = pow(psi1 / static_psi, 4);
    Q[g12] = 0;
    Q[g13] = 0;
    Q[g22] = pow(psi1 / static_psi, 4);
    Q[g23] = 0;
    Q[g33] = pow(psi1 / static_psi, 4);

    Q[K11] = Aij[0][0] / pow(psi1, 2);
    Q[K12] = Aij[0][1] / pow(psi1, 2);
    Q[K13] = Aij[0][2] / pow(psi1, 2);
    Q[K22] = Aij[1][1] / pow(psi1, 2);
    Q[K23] = Aij[1][2] / pow(psi1, 2);
    Q[K33] = Aij[2][2] / pow(psi1, 2);

    if (antisymmetric_lapse || averaged_lapse)
    {
        Q[lapse] = ((1.0 - 0.5 * mp / r_plus - 0.5 * mm / r_minus) /
                    (1.0 + 0.5 * mp / r_plus + 0.5 * mm / r_minus));

        if (r_plus < TP_Extend_Radius)
        {
            Q[lapse] = ((1.0 - 0.5 * EXTEND(mp, r_plus) - 0.5 * mm / r_minus) /
                        (1.0 + 0.5 * EXTEND(mp, r_plus) + 0.5 * mm / r_minus));
        }
        if (r_minus < TP_Extend_Radius)
        {
            Q[lapse] = ((1.0 - 0.5 * EXTEND(mm, r_minus) - 0.5 * mp / r_plus) /
                        (1.0 + 0.5 * EXTEND(mp, r_minus) + 0.5 * mp / r_plus));
        }

        if (averaged_lapse)
        {
            Q[lapse] = 0.5 * (1.0 + Q[lapse]);
        }
    }
    if (multiply_old_lapse)
        Q[lapse] *= old_alp;

    if (swap_xz)
    {
        /* Swap the x and z components of all tensors */
        if (conformal_state >= 2)
        {
            /// SWAP (psix[ind], psiz[ind]);
        }
        if (conformal_state >= 3)
        {
            /// SWAP (psixx[ind], psizz[ind]);
            /// SWAP (psixy[ind], psiyz[ind]);
        }
        SWAP(Q[g11], Q[g33]);
        SWAP(Q[g12], Q[g23]);
        SWAP(Q[K11], Q[K33]);
        SWAP(Q[K12], Q[K23]);
    } /* if swap_xz */

    if (use_sources && rescale_sources)
    {
        TP_WARN("@TODO Rescale_Sources is in some other thorn and has not been "
                "copied");
#if 0 == 1
    Rescale_Sources(cctkGH,
                    cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2],
                    x, y, z,
                    (conformal_state > 0) ? psi : NULL,
                    gxx, gyy, gzz,
                    gxy, gxz, gyz);
#endif
    }

#if 0 == 1
    /* Cleanup not possible as we don't know when the grid queries are finished */
    /* Keep the result around for the next time */
    free_dvector (F, 0, ntotal - 1);
    free_derivs (&u, ntotal);
    free_derivs (&v, ntotal);
    free_derivs (&cf_v, ntotal);
#endif
}

} // namespace TP
