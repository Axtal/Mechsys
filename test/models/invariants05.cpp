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
#include <mechsys/models/smpinvs.h>

using std::cout;
using std::endl;
using Util::_8s;
using Util::PI;
using Util::SQ2;
using Util::SQ3;
using Util::SQ6;
using Numerical::Diff;

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
        //Sig0 = -1.5, -2.0, -0.5,  -0.1*SQ2, -0.5*SQ2, 0.0;
        Sig0 = -1.5, -2.0, -0.5,  -0.1*SQ2, -0.5*SQ2, -0.00000001*SQ2;

        SI.b = 0.5;
    }

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

        // calculate state
        SI.Calc (Sig, true);

        // eigenvalues
        dLdt(0) = dot(SI.E0, dSigdt);
        dLdt(1) = dot(SI.E1, dSigdt);
        dLdt(2) = dot(SI.E2, dSigdt);

        // normal to SMP
        dNdt = product(SI.dNdL, dLdt);

        // unit normal to SMP
        dNudt = product(SI.dNudN, dNdt);

        // traction
        dtdt = product(SI.dtdL, dLdt);

        // normal projection
        dpdt = product(SI.dpdL, dLdt);

        // on-plane projection
        dqdt = product(SI.dqdL, dLdt);

        // invariants
        dspdt = dot (SI.dspdL, dLdt);
        dsqdt = dot (SI.dsqdL, dLdt);
        double dspdt_tmp = dot (SI.dspdSig, dSigdt);
        double dsqdt_tmp = dot (SI.dsqdSig, dSigdt);
        double err_dspdt = fabs(dspdt_tmp - dspdt);
        double err_dsqdt = fabs(dsqdt_tmp - dsqdt);
        //printf("dspdt = %16.8f = %16.8f : error = %16.8f       dsqdt = %16.8f = %16.8f : error = %16.8f\n", dspdt,dspdt_tmp,err_dspdt, dsqdt,dsqdt_tmp,err_dsqdt);
        if (err_dspdt>1.0e-15) throw new Fatal("Error on dspdt = %g",err_dspdt);
        if (err_dsqdt>1.0e-15) throw new Fatal("Error on dsqdt = %g",err_dsqdt);
    }

    // Functions
    double L0Fun  (double t) { CalcState(t); return SI.L(0); }
    double L1Fun  (double t) { CalcState(t); return SI.L(1); }
    double L2Fun  (double t) { CalcState(t); return SI.L(2); }

    double v00Fun (double t) { CalcState(t); return SI.v0(0); }
    double v01Fun (double t) { CalcState(t); return SI.v0(1); }
    double v02Fun (double t) { CalcState(t); return SI.v0(2); }

    double v10Fun (double t) { CalcState(t); return SI.v1(0); }
    double v11Fun (double t) { CalcState(t); return SI.v1(1); }
    double v12Fun (double t) { CalcState(t); return SI.v1(2); }

    double v20Fun (double t) { CalcState(t); return SI.v2(0); }
    double v21Fun (double t) { CalcState(t); return SI.v2(1); }
    double v22Fun (double t) { CalcState(t); return SI.v2(2); }

    double N0Fun  (double t) { CalcState(t); return SI.N(0); }
    double N1Fun  (double t) { CalcState(t); return SI.N(1); }
    double N2Fun  (double t) { CalcState(t); return SI.N(2); }

    double Nu0Fun (double t) { CalcState(t); return SI.Nu(0); }
    double Nu1Fun (double t) { CalcState(t); return SI.Nu(1); }
    double Nu2Fun (double t) { CalcState(t); return SI.Nu(2); }

    double t0Fun  (double t) { CalcState(t); return SI.t(0); }
    double t1Fun  (double t) { CalcState(t); return SI.t(1); }
    double t2Fun  (double t) { CalcState(t); return SI.t(2); }

    double p0Fun  (double t) { CalcState(t); return SI.p(0); }
    double p1Fun  (double t) { CalcState(t); return SI.p(1); }
    double p2Fun  (double t) { CalcState(t); return SI.p(2); }

    double q0Fun  (double t) { CalcState(t); return SI.q(0); }
    double q1Fun  (double t) { CalcState(t); return SI.q(1); }
    double q2Fun  (double t) { CalcState(t); return SI.q(2); }

    double spFun  (double t) { CalcState(t); return SI.sp; }
    double sqFun  (double t) { CalcState(t); return SI.sq; }

    // Data
    int      test;              // test number
    Vec_t    Sig0, Sig, dSigdt; // stress
    Mat_t    M, dMdt;           // multiplier
    SMPInvs  SI;                // invariants
    Vec3_t   dLdt;              // eigenvalues
    Vec3_t   dNdt, dNudt;       // normal to SMP
    Vec3_t   dtdt;              // traction
    Vec3_t   dpdt;              // normal projection
    Vec3_t   dqdt;              // on-plane projection
    double   dspdt, dsqdt;      // invariants
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
        double err_dL0dt = fabs(dL0dt_num - prob.dLdt(0));
        double err_dL1dt = fabs(dL1dt_num - prob.dLdt(1));
        double err_dL2dt = fabs(dL2dt_num - prob.dLdt(2));
        if (err_dL0dt > max_err_dL0dt) max_err_dL0dt = err_dL0dt;
        if (err_dL1dt > max_err_dL1dt) max_err_dL1dt = err_dL1dt;
        if (err_dL2dt > max_err_dL2dt) max_err_dL2dt = err_dL2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.SI.L(0), dL0dt_num, prob.dLdt(0), err_dL0dt,
                       prob.SI.L(1), dL1dt_num, prob.dLdt(1), err_dL1dt,
                       prob.SI.L(2), dL2dt_num, prob.dLdt(2), err_dL2dt);
        }
    }
    double tol_dL0dt = 1.0e-6;
    double tol_dL1dt = 1.0e-6;
    double tol_dL2dt = 1.0e-6;
    printf("  max_err_dL0dt  = %s%16.8e%s\n",(max_err_dL0dt >tol_dL0dt ?TERM_RED:TERM_GREEN),max_err_dL0dt, TERM_RST);
    printf("  max_err_dL1dt  = %s%16.8e%s\n",(max_err_dL1dt >tol_dL1dt ?TERM_RED:TERM_GREEN),max_err_dL1dt, TERM_RST);
    printf("  max_err_dL2dt  = %s%16.8e%s\n",(max_err_dL2dt >tol_dL2dt ?TERM_RED:TERM_GREEN),max_err_dL2dt, TERM_RST);



    // normal to the SMP
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
        double err_dN0dt = fabs(dN0dt_num - prob.dNdt(0));
        double err_dN1dt = fabs(dN1dt_num - prob.dNdt(1));
        double err_dN2dt = fabs(dN2dt_num - prob.dNdt(2));
        if (err_dN0dt > max_err_dN0dt) max_err_dN0dt = err_dN0dt;
        if (err_dN1dt > max_err_dN1dt) max_err_dN1dt = err_dN1dt;
        if (err_dN2dt > max_err_dN2dt) max_err_dN2dt = err_dN2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.SI.t(0), dN0dt_num, prob.dNdt(0), err_dN0dt,
                       prob.SI.t(1), dN1dt_num, prob.dNdt(1), err_dN1dt,
                       prob.SI.t(2), dN2dt_num, prob.dNdt(2), err_dN2dt);
        }
    }
    double tol_dN0dt = 1.0e-5;
    double tol_dN1dt = 1.0e-6;
    double tol_dN2dt = 1.0e-5;
    printf("  max_err_dN0dt  = %s%16.8e%s\n",(max_err_dN0dt >tol_dN0dt?TERM_RED:TERM_GREEN),max_err_dN0dt, TERM_RST);
    printf("  max_err_dN1dt  = %s%16.8e%s\n",(max_err_dN1dt >tol_dN1dt?TERM_RED:TERM_GREEN),max_err_dN1dt, TERM_RST);
    printf("  max_err_dN2dt  = %s%16.8e%s\n",(max_err_dN2dt >tol_dN2dt?TERM_RED:TERM_GREEN),max_err_dN2dt, TERM_RST);



    // unit normal to the SMP
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
        double err_dNu0dt = fabs(dNu0dt_num - prob.dNudt(0));
        double err_dNu1dt = fabs(dNu1dt_num - prob.dNudt(1));
        double err_dNu2dt = fabs(dNu2dt_num - prob.dNudt(2));
        if (err_dNu0dt > max_err_dNu0dt) max_err_dNu0dt = err_dNu0dt;
        if (err_dNu1dt > max_err_dNu1dt) max_err_dNu1dt = err_dNu1dt;
        if (err_dNu2dt > max_err_dNu2dt) max_err_dNu2dt = err_dNu2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.SI.t(0), dNu0dt_num, prob.dNudt(0), err_dNu0dt,
                       prob.SI.t(1), dNu1dt_num, prob.dNudt(1), err_dNu1dt,
                       prob.SI.t(2), dNu2dt_num, prob.dNudt(2), err_dNu2dt);
        }
    }
    double tol_dNu0dt = 1.0e-6;
    double tol_dNu1dt = 1.0e-6;
    double tol_dNu2dt = 1.0e-7;
    printf("  max_err_dNu0dt = %s%16.8e%s\n",(max_err_dNu0dt >tol_dNu0dt?TERM_RED:TERM_GREEN),max_err_dNu0dt, TERM_RST);
    printf("  max_err_dNu1dt = %s%16.8e%s\n",(max_err_dNu1dt >tol_dNu1dt?TERM_RED:TERM_GREEN),max_err_dNu1dt, TERM_RST);
    printf("  max_err_dNu2dt = %s%16.8e%s\n",(max_err_dNu2dt >tol_dNu2dt?TERM_RED:TERM_GREEN),max_err_dNu2dt, TERM_RST);



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
        double err_dt0dt = fabs(dt0dt_num - prob.dtdt(0));
        double err_dt1dt = fabs(dt1dt_num - prob.dtdt(1));
        double err_dt2dt = fabs(dt2dt_num - prob.dtdt(2));
        if (err_dt0dt > max_err_dt0dt) max_err_dt0dt = err_dt0dt;
        if (err_dt1dt > max_err_dt1dt) max_err_dt1dt = err_dt1dt;
        if (err_dt2dt > max_err_dt2dt) max_err_dt2dt = err_dt2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.SI.t(0), dt0dt_num, prob.dtdt(0), err_dt0dt,
                       prob.SI.t(1), dt1dt_num, prob.dtdt(1), err_dt1dt,
                       prob.SI.t(2), dt2dt_num, prob.dtdt(2), err_dt2dt);
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
        double err_dp0dt = fabs(dp0dt_num - prob.dpdt(0));
        double err_dp1dt = fabs(dp1dt_num - prob.dpdt(1));
        double err_dp2dt = fabs(dp2dt_num - prob.dpdt(2));
        if (err_dp0dt > max_err_dp0dt) max_err_dp0dt = err_dp0dt;
        if (err_dp1dt > max_err_dp1dt) max_err_dp1dt = err_dp1dt;
        if (err_dp2dt > max_err_dp2dt) max_err_dp2dt = err_dp2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.SI.p(0), dp0dt_num, prob.dpdt(0), err_dp0dt,
                       prob.SI.p(1), dp1dt_num, prob.dpdt(1), err_dp1dt,
                       prob.SI.p(2), dp2dt_num, prob.dpdt(2), err_dp2dt);
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
        double err_dq0dt = fabs(dq0dt_num - prob.dqdt(0));
        double err_dq1dt = fabs(dq1dt_num - prob.dqdt(1));
        double err_dq2dt = fabs(dq2dt_num - prob.dqdt(2));
        if (err_dq0dt > max_err_dq0dt) max_err_dq0dt = err_dq0dt;
        if (err_dq1dt > max_err_dq1dt) max_err_dq1dt = err_dq1dt;
        if (err_dq2dt > max_err_dq2dt) max_err_dq2dt = err_dq2dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.SI.p(0), dq0dt_num, prob.dqdt(0), err_dq0dt,
                       prob.SI.p(1), dq1dt_num, prob.dqdt(1), err_dq1dt,
                       prob.SI.p(2), dq2dt_num, prob.dqdt(2), err_dq2dt);
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
    double tol_dspdt = 1.0e-6;
    double tol_dsqdt = 1.0e-6;
    printf("  max_err_dspdt  = %s%16.8e%s\n",(max_err_dspdt>tol_dspdt?TERM_RED:TERM_GREEN),max_err_dspdt,TERM_RST);
    printf("  max_err_dsqdt  = %s%16.8e%s\n",(max_err_dsqdt>tol_dsqdt?TERM_RED:TERM_GREEN),max_err_dsqdt,TERM_RST);
    printf("\n");



    // end
    if (max_err_dL0dt  > tol_dL0dt ) return 1;
    if (max_err_dL1dt  > tol_dL1dt ) return 1;
    if (max_err_dL2dt  > tol_dL2dt ) return 1;
    if (max_err_dN0dt  > tol_dN0dt ) return 1;
    if (max_err_dN1dt  > tol_dN1dt ) return 1;
    if (max_err_dN2dt  > tol_dN2dt ) return 1;
    if (max_err_dNu0dt > tol_dNu0dt) return 1;
    if (max_err_dNu1dt > tol_dNu1dt) return 1;
    if (max_err_dNu2dt > tol_dNu2dt) return 1;
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
