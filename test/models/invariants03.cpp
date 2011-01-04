/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

// Std Lib
#include <iostream>

// GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// MechSys
#include <mechsys/numerical/numdiff.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/models/anisoinvs.h>

using std::cout;
using std::endl;
using Util::_8s;
using Util::PI;
using Util::SQ2;
using Util::SQ3;
using Util::SQ6;
using Numerical::Diff;
#ifdef HAS_TENSORS
  using namespace TensorsLib;
#endif

class Problem
{
public:
    // Constructor & Destructor
    Problem () : test(2)
    {
        Sig0  .change_dim(6);
        Sig   .change_dim(6);
        dSigdt.change_dim(6);
        M     .change_dim(6,6);
        dMdt  .change_dim(6,6);
        //Sig0 = -1.5, -2.0, -0.5,  -0.1*SQ2, -0.5*SQ2, -0.8*SQ2;
        //Sig0 = -1.5, -2.0, -0.5,  0., 0., 0;
        Sig0 = -1.5, -2.0, -0.5,  -0.1*SQ2, -0.5*SQ2, 0.0;
        //Sig0 = -1.5, -2.0, -0.5,  -0.1*SQ2, -0.5*SQ2, -0.00000001*SQ2;

        Vec3_t a(0.2, 0.3, 1.0);
        AI = new AnisoInvs (0.5, 0.15, a, true); // b, alpha, a, obliq
        //AI = new AnisoInvs (0.5, 0.0, a, false); // b, alpha, a, obliq

        dvdt.Resize (3);
    }
    ~Problem () { delete AI; }

    // Methods
    void CalcState (double t)
    {
        if (test==1)
        {
            Sig    = (1.+t*t)*Sig0;
            dSigdt = (2.*t)*Sig0;
        }
        else
        {
            M      = t,    0.,      0.,   0.,    0.,      0.,
                     0.,  t*t,      0.,   0.,    0.,      0.,
                     0.,   0.,  sin(t),   0.,    0.,      0.,
                     0.,   0.,      0.,    t,    0.,      0.,
                     0.,   0.,      0.,   0.,   t*t,      0.,
                     0.,   0.,      0.,   0.,    0.,  cos(t);
            dMdt   = 1.,   0.,      0.,   0.,    0.,      0.,
                     0., 2.*t,      0.,   0.,    0.,      0.,
                     0.,   0.,  cos(t),   0.,    0.,      0.,
                     0.,   0.,      0.,   1.,    0.,      0.,
                     0.,   0.,      0.,   0.,  2.*t,      0.,
                     0.,   0.,      0.,   0.,    0., -sin(t);
            Sig    = 0.1*Sig0 + M*Sig0;
            dSigdt = dMdt*Sig0;
        }

        // stress rate
        Ten2_t dAdt;
        Ten2Tensor (dSigdt, dAdt);

        // calculate state
        AI->Calc (Sig, true);

        // eigenvalues
        dLdt[0] = AI->E0 % dAdt;
        dLdt[1] = AI->E1 % dAdt;
        dLdt[2] = AI->E2 % dAdt;

        // eigenvectors
        dvdt[0] = AI->dv0dSig % dAdt;
        dvdt[1] = AI->dv1dSig % dAdt;
        dvdt[2] = AI->dv2dSig % dAdt;

        // normal to SMP in principal system
        dgdt = AI->dN123dL * dLdt;

        // normal to SMP in laboratory system
        dNdt = AI->dNdSig % dAdt;

        // unit normal to SMP in laboratory system
        dNudt = AI->dNudSig % dAdt;

        // normal to AMP in laboratory system
        dndt = AI->dndSig % dAdt;

        // unit normal to AMP in laboratory system
        dnudt = AI->dnudSig % dAdt;

        // traction
        dtdt = AI->dtdSig % dAdt;
        //dtdt = (dAdt * AI->tnu) + (AI->tSig * dnudt);

        // normal projection
        dpdt = AI->dpdSig % dAdt;

        // on-plane projection
        dqdt = AI->dqdSig % dAdt;

        // invariants
        dspdt = dot (AI->dspdSig, dSigdt);
        dsqdt = dot (AI->dsqdSig, dSigdt);
    }

    // Functions
    double L0Fun  (double t) { CalcState(t); return AI->L(0); }
    double L1Fun  (double t) { CalcState(t); return AI->L(1); }
    double L2Fun  (double t) { CalcState(t); return AI->L(2); }

    double v00Fun (double t) { CalcState(t); return AI->v0(0); }
    double v01Fun (double t) { CalcState(t); return AI->v0(1); }
    double v02Fun (double t) { CalcState(t); return AI->v0(2); }

    double v10Fun (double t) { CalcState(t); return AI->v1(0); }
    double v11Fun (double t) { CalcState(t); return AI->v1(1); }
    double v12Fun (double t) { CalcState(t); return AI->v1(2); }

    double v20Fun (double t) { CalcState(t); return AI->v2(0); }
    double v21Fun (double t) { CalcState(t); return AI->v2(1); }
    double v22Fun (double t) { CalcState(t); return AI->v2(2); }

    double g0Fun  (double t) { CalcState(t); return AI->N123(0); }
    double g1Fun  (double t) { CalcState(t); return AI->N123(1); }
    double g2Fun  (double t) { CalcState(t); return AI->N123(2); }

    double N0Fun  (double t) { CalcState(t); return AI->N(0); }
    double N1Fun  (double t) { CalcState(t); return AI->N(1); }
    double N2Fun  (double t) { CalcState(t); return AI->N(2); }

    double Nu0Fun (double t) { CalcState(t); return AI->Nu(0); }
    double Nu1Fun (double t) { CalcState(t); return AI->Nu(1); }
    double Nu2Fun (double t) { CalcState(t); return AI->Nu(2); }

    double n0Fun  (double t) { CalcState(t); return AI->n(0);  }
    double n1Fun  (double t) { CalcState(t); return AI->n(1);  }
    double n2Fun  (double t) { CalcState(t); return AI->n(2);  }

    double nu0Fun (double t) { CalcState(t); return AI->nu(0); }
    double nu1Fun (double t) { CalcState(t); return AI->nu(1); }
    double nu2Fun (double t) { CalcState(t); return AI->nu(2); }

    double t0Fun  (double t) { CalcState(t); return AI->t(0); }
    double t1Fun  (double t) { CalcState(t); return AI->t(1); }
    double t2Fun  (double t) { CalcState(t); return AI->t(2); }

    double p0Fun  (double t) { CalcState(t); return AI->p(0); }
    double p1Fun  (double t) { CalcState(t); return AI->p(1); }
    double p2Fun  (double t) { CalcState(t); return AI->p(2); }

    double q0Fun  (double t) { CalcState(t); return AI->q(0); }
    double q1Fun  (double t) { CalcState(t); return AI->q(1); }
    double q2Fun  (double t) { CalcState(t); return AI->q(2); }

    double spFun  (double t) { CalcState(t); return AI->sp; }
    double sqFun  (double t) { CalcState(t); return AI->sq; }

    // Data
    int           test;              // test number
    Vec_t         Sig0, Sig, dSigdt; // stress
    Mat_t         M, dMdt;           // multiplier
    AnisoInvs   * AI;                // invariants
    Ten1_t        dLdt;              // eigenvalues
    Array<Ten1_t> dvdt;              // eigenvectors
    Ten1_t        dgdt;              // normal to SMP in principal system
    Ten1_t        dNdt, dNudt;       // normal to SMP in laboratory system
    Ten1_t        dndt, dnudt;       // normal to AMP in laboratory system
    Ten1_t        dtdt;              // traction
    Ten1_t        dpdt;              // normal projection
    Ten1_t        dqdt;              // on-plane projection
    double        dspdt, dsqdt;      // invariants
};

typedef double (Problem::*pFun) (double t);

int main(int argc, char **argv) try
{
    // initialize problem
    Problem       prob;
    Diff<Problem> nd(&prob);
    bool   verbose = false;
    size_t ndiv    = 20;
    if (argc>1) prob.test = atoi(argv[1]);
    if (argc>2) verbose   = atoi(argv[2]);
    if (argc>3) ndiv      = atoi(argv[3]);



    // L0, L1, L2 and derivatives
    if (verbose)
    {
        printf("\n");
        printf("%6s %12s %12s %12s %16s  %12s %12s %12s %16s  %12s %12s %12s %16s\n",
                "t","L0","dL0dt_num","dL0dt","err(dL0dt)",
                    "L1","dL1dt_num","dL1dt","err(dL1dt)",
                    "L2","dL2dt_num","dL2dt","err(dL2dt)");
    }
    double max_err_dL0dt = 0.0;
    double max_err_dL1dt = 0.0;
    double max_err_dL2dt = 0.0;
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t         = (double)i/(double)ndiv;
        double dL0dt_num = nd.DyDx (&Problem::L0Fun, t);
        double dL1dt_num = nd.DyDx (&Problem::L1Fun, t);
        double dL2dt_num = nd.DyDx (&Problem::L2Fun, t);
        prob.CalcState (t);
        double err_dL0dt = fabs(dL0dt_num - prob.dLdt[0]);
        double err_dL1dt = fabs(dL1dt_num - prob.dLdt[1]);
        double err_dL2dt = fabs(dL2dt_num - prob.dLdt[2]);
        if (err_dL0dt > max_err_dL0dt) max_err_dL0dt = err_dL0dt;
        if (err_dL1dt > max_err_dL1dt) max_err_dL1dt = err_dL1dt;
        if (err_dL2dt > max_err_dL2dt) max_err_dL2dt = err_dL2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.AI->L(0), dL0dt_num, prob.dLdt[0], err_dL0dt,
                       prob.AI->L(1), dL1dt_num, prob.dLdt[1], err_dL1dt,
                       prob.AI->L(2), dL2dt_num, prob.dLdt[2], err_dL2dt);
        }
    }
    double tol_dL0dt = 1.0e-6;
    double tol_dL1dt = 1.0e-6;
    double tol_dL2dt = 1.0e-6;
    printf("  max_err_dL0dt  = %s%16.8e%s\n",(max_err_dL0dt >tol_dL0dt ?TERM_RED:TERM_GREEN),max_err_dL0dt, TERM_RST);
    printf("  max_err_dL1dt  = %s%16.8e%s\n",(max_err_dL1dt >tol_dL1dt ?TERM_RED:TERM_GREEN),max_err_dL1dt, TERM_RST);
    printf("  max_err_dL2dt  = %s%16.8e%s\n",(max_err_dL2dt >tol_dL2dt ?TERM_RED:TERM_GREEN),max_err_dL2dt, TERM_RST);



    // eigenvectors
    double max_err_dvdt[3][3] = {{0.,0.,0.},
                                 {0.,0.,0.},
                                 {0.,0.,0.}};
    pFun vfuncs[3][3] = {{&Problem::v00Fun, &Problem::v01Fun, &Problem::v02Fun},
                         {&Problem::v10Fun, &Problem::v11Fun, &Problem::v12Fun},
                         {&Problem::v20Fun, &Problem::v21Fun, &Problem::v22Fun}};
    for (size_t k=0; k<3; ++k)
    {
        if (verbose)
        {
            if (verbose) printf("\n%6s","t");
            for (size_t j=0; j<3; ++j)
            {
                char str0[32];
                char str1[32];
                char str2[32];
                char str3[32];
                sprintf(str0,"v%zd%zd",          k,j);
                sprintf(str1,"dv%zd%zddt_num",   k,j);
                sprintf(str2,"dv%zd%zddt",       k,j);
                sprintf(str3,"error(dv%zd%zddt)",k,j);
                printf("%12s %12s %12s %16s  ",str0,str1,str2,str3);
            }
            printf("\n");
        }
        for (size_t i=0; i<ndiv+1; ++i)
        {
            double t = (double)i/(double)ndiv;
            prob.CalcState (t);
            if (verbose) printf("%6.3f",t);
            for (size_t j=0; j<3; ++j)
            {
                double dvdt_num = nd.DyDx (vfuncs[k][j], t);
                double err      = fabs(dvdt_num - prob.dvdt[k][j]);
                if (err > max_err_dvdt[k][j]) max_err_dvdt[k][j] = err;
                Vec3_t const & vv = (k==0 ? prob.AI->v0 : (k==1 ? prob.AI->v1 : prob.AI->v2));
                if (verbose) printf("%12.8f %12.8f %12.8f %16.8e  ", vv(j), dvdt_num, prob.dvdt[k][j], err);
            }
            if (verbose) printf("\n");
        }
    }
    double tol_dvdt[3][3]= {{1.0e-5, 1.0e-6, 1.0e-6},
                            {1.0e-6, 1.0e-5, 1.0e-6},
                            {1.0e-6, 1.0e-6, 1.0e-6}};
    for (size_t k=0; k<3; ++k)
    for (size_t i=0; i<3; ++i)
        printf("  max_err_dv%zd%zddt = %s%16.8e%s\n",k,i,(max_err_dvdt[k][i]>tol_dvdt[k][i]?TERM_RED:TERM_GREEN),max_err_dvdt[k][i],TERM_RST);



    // normal to the SMP in principal system : g
    if (verbose)
    {
        printf("\n");
        printf("%6s %12s %12s %12s %16s  %12s %12s %12s %16s  %12s %12s %12s %16s\n",
                "t","g0","dg0dt_num","dg0dt","err(dg0dt)",
                    "g1","dg1dt_num","dg1dt","err(dg1dt)",
                    "g2","dg2dt_num","dg2dt","err(dg2dt)");
    }
    double max_err_dg0dt = 0.0;
    double max_err_dg1dt = 0.0;
    double max_err_dg2dt = 0.0;
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t         = (double)i/(double)ndiv;
        double dg0dt_num = nd.DyDx (&Problem::g0Fun, t);
        double dg1dt_num = nd.DyDx (&Problem::g1Fun, t);
        double dg2dt_num = nd.DyDx (&Problem::g2Fun, t);
        prob.CalcState (t);
        double err_dg0dt = fabs(dg0dt_num - prob.dgdt[0]);
        double err_dg1dt = fabs(dg1dt_num - prob.dgdt[1]);
        double err_dg2dt = fabs(dg2dt_num - prob.dgdt[2]);
        if (err_dg0dt > max_err_dg0dt) max_err_dg0dt = err_dg0dt;
        if (err_dg1dt > max_err_dg1dt) max_err_dg1dt = err_dg1dt;
        if (err_dg2dt > max_err_dg2dt) max_err_dg2dt = err_dg2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.AI->t(0), dg0dt_num, prob.dgdt[0], err_dg0dt,
                       prob.AI->t(1), dg1dt_num, prob.dgdt[1], err_dg1dt,
                       prob.AI->t(2), dg2dt_num, prob.dgdt[2], err_dg2dt);
        }
    }
    double tol_dg0dt = 1.0e-5;
    double tol_dg1dt = 1.0e-6;
    double tol_dg2dt = 1.0e-6;
    printf("  max_err_dg0dt  = %s%16.8e%s\n",(max_err_dg0dt >tol_dg0dt?TERM_RED:TERM_GREEN),max_err_dg0dt, TERM_RST);
    printf("  max_err_dg1dt  = %s%16.8e%s\n",(max_err_dg1dt >tol_dg1dt?TERM_RED:TERM_GREEN),max_err_dg1dt, TERM_RST);
    printf("  max_err_dg2dt  = %s%16.8e%s\n",(max_err_dg2dt >tol_dg2dt?TERM_RED:TERM_GREEN),max_err_dg2dt, TERM_RST);



    // normal to the SMP in laboratory system : N
    if (verbose)
    {
        printf("\n");
        printf("%6s %12s %12s %12s %16s  %12s %12s %12s %16s  %12s %12s %12s %16s\n",
                "t","N0","dN0dt_num","dN0dt","err(dN0dt)",
                    "N1","dN1dt_num","dN1dt","err(dN1dt)",
                    "N2","dN2dt_num","dN2dt","err(dN2dt)");
    }
    double max_err_dN0dt = 0.0;
    double max_err_dN1dt = 0.0;
    double max_err_dN2dt = 0.0;
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t         = (double)i/(double)ndiv;
        double dN0dt_num = nd.DyDx (&Problem::N0Fun, t);
        double dN1dt_num = nd.DyDx (&Problem::N1Fun, t);
        double dN2dt_num = nd.DyDx (&Problem::N2Fun, t);
        prob.CalcState (t);
        double err_dN0dt = fabs(dN0dt_num - prob.dNdt[0]);
        double err_dN1dt = fabs(dN1dt_num - prob.dNdt[1]);
        double err_dN2dt = fabs(dN2dt_num - prob.dNdt[2]);
        if (err_dN0dt > max_err_dN0dt) max_err_dN0dt = err_dN0dt;
        if (err_dN1dt > max_err_dN1dt) max_err_dN1dt = err_dN1dt;
        if (err_dN2dt > max_err_dN2dt) max_err_dN2dt = err_dN2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.AI->t(0), dN0dt_num, prob.dNdt[0], err_dN0dt,
                       prob.AI->t(1), dN1dt_num, prob.dNdt[1], err_dN1dt,
                       prob.AI->t(2), dN2dt_num, prob.dNdt[2], err_dN2dt);
        }
    }
    double tol_dN0dt = 1.0e-5;
    double tol_dN1dt = 1.0e-6;
    double tol_dN2dt = 1.0e-5;
    printf("  max_err_dN0dt  = %s%16.8e%s\n",(max_err_dN0dt >tol_dN0dt?TERM_RED:TERM_GREEN),max_err_dN0dt, TERM_RST);
    printf("  max_err_dN1dt  = %s%16.8e%s\n",(max_err_dN1dt >tol_dN1dt?TERM_RED:TERM_GREEN),max_err_dN1dt, TERM_RST);
    printf("  max_err_dN2dt  = %s%16.8e%s\n",(max_err_dN2dt >tol_dN2dt?TERM_RED:TERM_GREEN),max_err_dN2dt, TERM_RST);



    // unit normal to the SMP in laboratory system : Nu
    if (verbose)
    {
        printf("\n");
        printf("%6s %12s %12s %12s %16s  %12s %12s %12s %16s  %12s %12s %12s %16s\n",
                "t","Nu0","dNu0dt_num","dNu0dt","err(dNu0dt)",
                    "Nu1","dNu1dt_num","dNu1dt","err(dNu1dt)",
                    "Nu2","dNu2dt_num","dNu2dt","err(dNu2dt)");
    }
    double max_err_dNu0dt = 0.0;
    double max_err_dNu1dt = 0.0;
    double max_err_dNu2dt = 0.0;
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t         = (double)i/(double)ndiv;
        double dNu0dt_num = nd.DyDx (&Problem::Nu0Fun, t);
        double dNu1dt_num = nd.DyDx (&Problem::Nu1Fun, t);
        double dNu2dt_num = nd.DyDx (&Problem::Nu2Fun, t);
        prob.CalcState (t);
        double err_dNu0dt = fabs(dNu0dt_num - prob.dNudt[0]);
        double err_dNu1dt = fabs(dNu1dt_num - prob.dNudt[1]);
        double err_dNu2dt = fabs(dNu2dt_num - prob.dNudt[2]);
        if (err_dNu0dt > max_err_dNu0dt) max_err_dNu0dt = err_dNu0dt;
        if (err_dNu1dt > max_err_dNu1dt) max_err_dNu1dt = err_dNu1dt;
        if (err_dNu2dt > max_err_dNu2dt) max_err_dNu2dt = err_dNu2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.AI->t(0), dNu0dt_num, prob.dNudt[0], err_dNu0dt,
                       prob.AI->t(1), dNu1dt_num, prob.dNudt[1], err_dNu1dt,
                       prob.AI->t(2), dNu2dt_num, prob.dNudt[2], err_dNu2dt);
        }
    }
    double tol_dNu0dt = 1.0e-6;
    double tol_dNu1dt = 1.0e-6;
    double tol_dNu2dt = 1.0e-7;
    printf("  max_err_dNu0dt = %s%16.8e%s\n",(max_err_dNu0dt >tol_dNu0dt?TERM_RED:TERM_GREEN),max_err_dNu0dt, TERM_RST);
    printf("  max_err_dNu1dt = %s%16.8e%s\n",(max_err_dNu1dt >tol_dNu1dt?TERM_RED:TERM_GREEN),max_err_dNu1dt, TERM_RST);
    printf("  max_err_dNu2dt = %s%16.8e%s\n",(max_err_dNu2dt >tol_dNu2dt?TERM_RED:TERM_GREEN),max_err_dNu2dt, TERM_RST);



    // unit normal to the SMP in laboratory system : n
    if (verbose)
    {
        printf("\n");
        printf("%6s %12s %12s %12s %16s  %12s %12s %12s %16s  %12s %12s %12s %16s\n",
                "t","n0","dn0dt_num","dn0dt","err(dn0dt)",
                    "n1","dn1dt_num","dn1dt","err(dn1dt)",
                    "n2","dn2dt_num","dn2dt","err(dn2dt)");
    }
    double max_err_dn0dt = 0.0;
    double max_err_dn1dt = 0.0;
    double max_err_dn2dt = 0.0;
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t         = (double)i/(double)ndiv;
        double dn0dt_num = nd.DyDx (&Problem::n0Fun, t);
        double dn1dt_num = nd.DyDx (&Problem::n1Fun, t);
        double dn2dt_num = nd.DyDx (&Problem::n2Fun, t);
        prob.CalcState (t);
        double err_dn0dt = fabs(dn0dt_num - prob.dndt[0]);
        double err_dn1dt = fabs(dn1dt_num - prob.dndt[1]);
        double err_dn2dt = fabs(dn2dt_num - prob.dndt[2]);
        if (err_dn0dt > max_err_dn0dt) max_err_dn0dt = err_dn0dt;
        if (err_dn1dt > max_err_dn1dt) max_err_dn1dt = err_dn1dt;
        if (err_dn2dt > max_err_dn2dt) max_err_dn2dt = err_dn2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.AI->t(0), dn0dt_num, prob.dndt[0], err_dn0dt,
                       prob.AI->t(1), dn1dt_num, prob.dndt[1], err_dn1dt,
                       prob.AI->t(2), dn2dt_num, prob.dndt[2], err_dn2dt);
        }
    }
    double tol_dn0dt = 1.0e-6;
    double tol_dn1dt = 1.0e-6;
    double tol_dn2dt = 1.0e-7;
    printf("  max_err_dn0dt  = %s%16.8e%s\n",(max_err_dn0dt >tol_dn0dt?TERM_RED:TERM_GREEN),max_err_dn0dt, TERM_RST);
    printf("  max_err_dn1dt  = %s%16.8e%s\n",(max_err_dn1dt >tol_dn1dt?TERM_RED:TERM_GREEN),max_err_dn1dt, TERM_RST);
    printf("  max_err_dn2dt  = %s%16.8e%s\n",(max_err_dn2dt >tol_dn2dt?TERM_RED:TERM_GREEN),max_err_dn2dt, TERM_RST);



    // unit normal to the AMP in laboratory system : nu
    if (verbose)
    {
        printf("\n");
        printf("%6s %12s %12s %12s %16s  %12s %12s %12s %16s  %12s %12s %12s %16s\n",
                "t","nu0","dnu0dt_num","dnu0dt","err(dnu0dt)",
                    "nu1","dnu1dt_num","dnu1dt","err(dnu1dt)",
                    "nu2","dnu2dt_num","dnu2dt","err(dnu2dt)");
    }
    double max_err_dnu0dt = 0.0;
    double max_err_dnu1dt = 0.0;
    double max_err_dnu2dt = 0.0;
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t         = (double)i/(double)ndiv;
        double dnu0dt_num = nd.DyDx (&Problem::nu0Fun, t);
        double dnu1dt_num = nd.DyDx (&Problem::nu1Fun, t);
        double dnu2dt_num = nd.DyDx (&Problem::nu2Fun, t);
        prob.CalcState (t);
        double err_dnu0dt = fabs(dnu0dt_num - prob.dnudt[0]);
        double err_dnu1dt = fabs(dnu1dt_num - prob.dnudt[1]);
        double err_dnu2dt = fabs(dnu2dt_num - prob.dnudt[2]);
        if (err_dnu0dt > max_err_dnu0dt) max_err_dnu0dt = err_dnu0dt;
        if (err_dnu1dt > max_err_dnu1dt) max_err_dnu1dt = err_dnu1dt;
        if (err_dnu2dt > max_err_dnu2dt) max_err_dnu2dt = err_dnu2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.AI->t(0), dnu0dt_num, prob.dnudt[0], err_dnu0dt,
                       prob.AI->t(1), dnu1dt_num, prob.dnudt[1], err_dnu1dt,
                       prob.AI->t(2), dnu2dt_num, prob.dnudt[2], err_dnu2dt);
        }
    }
    double tol_dnu0dt = 1.0e-6;
    double tol_dnu1dt = 1.0e-7;
    double tol_dnu2dt = 1.0e-7;
    printf("  max_err_dnu0dt = %s%16.8e%s\n",(max_err_dnu0dt >tol_dnu0dt?TERM_RED:TERM_GREEN),max_err_dnu0dt, TERM_RST);
    printf("  max_err_dnu1dt = %s%16.8e%s\n",(max_err_dnu1dt >tol_dnu1dt?TERM_RED:TERM_GREEN),max_err_dnu1dt, TERM_RST);
    printf("  max_err_dnu2dt = %s%16.8e%s\n",(max_err_dnu2dt >tol_dnu2dt?TERM_RED:TERM_GREEN),max_err_dnu2dt, TERM_RST);



    // traction
    if (verbose)
    {
        printf("\n");
        printf("%6s %12s %12s %12s %16s  %12s %12s %12s %16s  %12s %12s %12s %16s\n",
                "t","t0","dt0dt_num","dt0dt","err(dt0dt)",
                    "t1","dt1dt_num","dt1dt","err(dt1dt)",
                    "t2","dt2dt_num","dt2dt","err(dt2dt)");
    }
    double max_err_dt0dt = 0.0;
    double max_err_dt1dt = 0.0;
    double max_err_dt2dt = 0.0;
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t         = (double)i/(double)ndiv;
        double dt0dt_num = nd.DyDx (&Problem::t0Fun, t);
        double dt1dt_num = nd.DyDx (&Problem::t1Fun, t);
        double dt2dt_num = nd.DyDx (&Problem::t2Fun, t);
        prob.CalcState (t);
        double err_dt0dt = fabs(dt0dt_num - prob.dtdt[0]);
        double err_dt1dt = fabs(dt1dt_num - prob.dtdt[1]);
        double err_dt2dt = fabs(dt2dt_num - prob.dtdt[2]);
        if (err_dt0dt > max_err_dt0dt) max_err_dt0dt = err_dt0dt;
        if (err_dt1dt > max_err_dt1dt) max_err_dt1dt = err_dt1dt;
        if (err_dt2dt > max_err_dt2dt) max_err_dt2dt = err_dt2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.AI->t(0), dt0dt_num, prob.dtdt[0], err_dt0dt,
                       prob.AI->t(1), dt1dt_num, prob.dtdt[1], err_dt1dt,
                       prob.AI->t(2), dt2dt_num, prob.dtdt[2], err_dt2dt);
        }
    }
    double tol_dt0dt = 1.0e-6;
    double tol_dt1dt = 1.0e-6;
    double tol_dt2dt = 1.0e-6;
    printf("  max_err_dt0dt  = %s%16.8e%s\n",(max_err_dt0dt >tol_dt0dt?TERM_RED:TERM_GREEN),max_err_dt0dt, TERM_RST);
    printf("  max_err_dt1dt  = %s%16.8e%s\n",(max_err_dt1dt >tol_dt1dt?TERM_RED:TERM_GREEN),max_err_dt1dt, TERM_RST);
    printf("  max_err_dt2dt  = %s%16.8e%s\n",(max_err_dt2dt >tol_dt2dt?TERM_RED:TERM_GREEN),max_err_dt2dt, TERM_RST);



    // normal projection
    if (verbose)
    {
        printf("\n");
        printf("%6s %12s %12s %12s %16s  %12s %12s %12s %16s  %12s %12s %12s %16s\n",
                "t","p0","dp0dt_num","dp0dt","err(dp0dt)",
                    "p1","dp1dt_num","dp1dt","err(dp1dt)",
                    "p2","dp2dt_num","dp2dt","err(dp2dt)");
    }
    double max_err_dp0dt = 0.0;
    double max_err_dp1dt = 0.0;
    double max_err_dp2dt = 0.0;
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t         = (double)i/(double)ndiv;
        double dp0dt_num = nd.DyDx (&Problem::p0Fun, t);
        double dp1dt_num = nd.DyDx (&Problem::p1Fun, t);
        double dp2dt_num = nd.DyDx (&Problem::p2Fun, t);
        prob.CalcState (t);
        double err_dp0dt = fabs(dp0dt_num - prob.dpdt[0]);
        double err_dp1dt = fabs(dp1dt_num - prob.dpdt[1]);
        double err_dp2dt = fabs(dp2dt_num - prob.dpdt[2]);
        if (err_dp0dt > max_err_dp0dt) max_err_dp0dt = err_dp0dt;
        if (err_dp1dt > max_err_dp1dt) max_err_dp1dt = err_dp1dt;
        if (err_dp2dt > max_err_dp2dt) max_err_dp2dt = err_dp2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.AI->p(0), dp0dt_num, prob.dpdt[0], err_dp0dt,
                       prob.AI->p(1), dp1dt_num, prob.dpdt[1], err_dp1dt,
                       prob.AI->p(2), dp2dt_num, prob.dpdt[2], err_dp2dt);
        }
    }
    double tol_dp0dt = 1.0e-6;
    double tol_dp1dt = 1.0e-6;
    double tol_dp2dt = 1.0e-6;
    printf("  max_err_dp0dt  = %s%16.8e%s\n",(max_err_dp0dt >tol_dp0dt?TERM_RED:TERM_GREEN),max_err_dp0dt, TERM_RST);
    printf("  max_err_dp1dt  = %s%16.8e%s\n",(max_err_dp1dt >tol_dp1dt?TERM_RED:TERM_GREEN),max_err_dp1dt, TERM_RST);
    printf("  max_err_dp2dt  = %s%16.8e%s\n",(max_err_dp2dt >tol_dp2dt?TERM_RED:TERM_GREEN),max_err_dp2dt, TERM_RST);



    // on-plane projection
    if (verbose)
    {
        printf("\n");
        printf("%6s %12s %12s %12s %16s  %12s %12s %12s %16s  %12s %12s %12s %16s\n",
                "t","q0","dq0dt_num","dq0dt","err(dq0dt)",
                    "q1","dq1dt_num","dq1dt","err(dq1dt)",
                    "q2","dq2dt_num","dq2dt","err(dq2dt)");
    }
    double max_err_dq0dt = 0.0;
    double max_err_dq1dt = 0.0;
    double max_err_dq2dt = 0.0;
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t         = (double)i/(double)ndiv;
        double dq0dt_num = nd.DyDx (&Problem::q0Fun, t);
        double dq1dt_num = nd.DyDx (&Problem::q1Fun, t);
        double dq2dt_num = nd.DyDx (&Problem::q2Fun, t);
        prob.CalcState (t);
        double err_dq0dt = fabs(dq0dt_num - prob.dqdt[0]);
        double err_dq1dt = fabs(dq1dt_num - prob.dqdt[1]);
        double err_dq2dt = fabs(dq2dt_num - prob.dqdt[2]);
        if (err_dq0dt > max_err_dq0dt) max_err_dq0dt = err_dq0dt;
        if (err_dq1dt > max_err_dq1dt) max_err_dq1dt = err_dq1dt;
        if (err_dq2dt > max_err_dq2dt) max_err_dq2dt = err_dq2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.AI->p(0), dq0dt_num, prob.dqdt[0], err_dq0dt,
                       prob.AI->p(1), dq1dt_num, prob.dqdt[1], err_dq1dt,
                       prob.AI->p(2), dq2dt_num, prob.dqdt[2], err_dq2dt);
        }
    }
    double tol_dq0dt = 1.0e-6;
    double tol_dq1dt = 1.0e-6;
    double tol_dq2dt = 1.0e-6;
    printf("  max_err_dq0dt  = %s%16.8e%s\n",(max_err_dq0dt >tol_dq0dt?TERM_RED:TERM_GREEN),max_err_dq0dt, TERM_RST);
    printf("  max_err_dq1dt  = %s%16.8e%s\n",(max_err_dq1dt >tol_dq1dt?TERM_RED:TERM_GREEN),max_err_dq1dt, TERM_RST);
    printf("  max_err_dq2dt  = %s%16.8e%s\n",(max_err_dq2dt >tol_dq2dt?TERM_RED:TERM_GREEN),max_err_dq2dt, TERM_RST);



    // invariants
    double max_err_dspdt = 0.0;
    double max_err_dsqdt = 0.0;
    if (verbose)
    {
        printf("\n%6s","t");
        printf("%12s %12s %16s  %12s %12s %16s\n", "dspdt_num","dspdt","err(dspdt)", "dsqdt_num","dsqdt","err(dsqdt)");
    }
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t = (double)i/(double)ndiv;
        prob.CalcState (t);
        if (verbose) printf("%6.3f",t);

        // dspdt
        double dspdt_num = nd.DyDx (&Problem::spFun, t);
        double err_dspdt = fabs(dspdt_num - prob.dspdt);
        if (err_dspdt > max_err_dspdt) max_err_dspdt = err_dspdt;

        // dsqdt
        double dsqdt_num = nd.DyDx (&Problem::sqFun, t);
        double err_dsqdt = fabs(dsqdt_num - prob.dsqdt);
        if (err_dsqdt > max_err_dsqdt) max_err_dsqdt = err_dsqdt;

        if (verbose) printf("%12.8f %12.8f %16.8e  %12.8f %12.8f %16.8e\n", dspdt_num,prob.dspdt,err_dspdt, dsqdt_num,prob.dsqdt,err_dsqdt);
    }
    double tol_dspdt = 1.0e-7;
    double tol_dsqdt = 1.0e-7;
    printf("  max_err_dspdt  = %s%16.8e%s\n",(max_err_dspdt>tol_dspdt?TERM_RED:TERM_GREEN),max_err_dspdt,TERM_RST);
    printf("  max_err_dsqdt  = %s%16.8e%s\n",(max_err_dsqdt>tol_dsqdt?TERM_RED:TERM_GREEN),max_err_dsqdt,TERM_RST);
    printf("\n");



    // end
    if (max_err_dL0dt > tol_dL0dt) return 1;
    if (max_err_dL1dt > tol_dL1dt) return 1;
    if (max_err_dL2dt > tol_dL2dt) return 1;
    for (size_t k=0; k<3; ++k)
    for (size_t i=0; i<3; ++i) if (max_err_dvdt[k][i] > tol_dvdt[k][i]) return 1;
    if (max_err_dg0dt  > tol_dg0dt ) return 1;
    if (max_err_dg1dt  > tol_dg1dt ) return 1;
    if (max_err_dg2dt  > tol_dg2dt ) return 1;
    if (max_err_dN0dt  > tol_dN0dt ) return 1;
    if (max_err_dN1dt  > tol_dN1dt ) return 1;
    if (max_err_dN2dt  > tol_dN2dt ) return 1;
    if (max_err_dNu0dt > tol_dNu0dt) return 1;
    if (max_err_dNu1dt > tol_dNu1dt) return 1;
    if (max_err_dNu2dt > tol_dNu2dt) return 1;
    if (max_err_dn0dt  > tol_dn0dt ) return 1;
    if (max_err_dn1dt  > tol_dn1dt ) return 1;
    if (max_err_dn2dt  > tol_dn2dt ) return 1;
    if (max_err_dnu0dt > tol_dnu0dt) return 1;
    if (max_err_dnu1dt > tol_dnu1dt) return 1;
    if (max_err_dnu2dt > tol_dnu2dt) return 1;
    if (max_err_dt0dt  > tol_dt0dt ) return 1;
    if (max_err_dt1dt  > tol_dt1dt ) return 1;
    if (max_err_dt2dt  > tol_dt2dt ) return 1;
    if (max_err_dp0dt  > tol_dp0dt ) return 1;
    if (max_err_dp1dt  > tol_dp1dt ) return 1;
    if (max_err_dp2dt  > tol_dp2dt ) return 1;
    if (max_err_dq0dt  > tol_dq0dt ) return 1;
    if (max_err_dq1dt  > tol_dq1dt ) return 1;
    if (max_err_dq2dt  > tol_dq2dt ) return 1;
    if (max_err_dspdt  > tol_dspdt ) return 1;
    if (max_err_dsqdt  > tol_dsqdt ) return 1;
    return 0;
}
MECHSYS_CATCH
