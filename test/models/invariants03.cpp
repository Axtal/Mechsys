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
        Sig0 = -1.5, -2.0, -0.5,  -0.1*SQ2, -0.5*SQ2, -0.8*SQ2;
        //Sig0 = -1.5, -2.0, -0.5,  0., 0., 0;

        Vec_t a(3);
        a  = 0., 0., 1.;
        AI = new AnisoInvs (0.5, 0.1, a, true); // b, alpha, a, obliq

        dsdt.change_dim (3);
        dNdt.change_dim (3);
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

        AI->Calc (Sig, true);
        dsdt(0) = dot(AI->P0, dSigdt);
        dsdt(1) = dot(AI->P1, dSigdt);
        dsdt(2) = dot(AI->P2, dSigdt);
        dNdt    = AI->dNds  * dsdt;
        dNudt   = AI->dNuds * dsdt;
        dndt    = AI->dnds  * dsdt;
        dnudt   = AI->dnuds * dsdt;
    }

    // Functions
    double L0Fun  (double t) { CalcState(t); return AI->L(0);  }
    double L1Fun  (double t) { CalcState(t); return AI->L(1);  }
    double L2Fun  (double t) { CalcState(t); return AI->L(2);  }
    double spFun  (double t) { CalcState(t); return AI->sp;    }
    double sqFun  (double t) { CalcState(t); return AI->sq;    }
    double N0Fun  (double t) { CalcState(t); return AI->N(0);  }
    double N1Fun  (double t) { CalcState(t); return AI->N(1);  }
    double N2Fun  (double t) { CalcState(t); return AI->N(2);  }
    double Nu0Fun (double t) { CalcState(t); return AI->Nu(0); }
    double Nu1Fun (double t) { CalcState(t); return AI->Nu(1); }
    double Nu2Fun (double t) { CalcState(t); return AI->Nu(2); }
    double n0Fun  (double t) { CalcState(t); return AI->n(0);  }
    double n1Fun  (double t) { CalcState(t); return AI->n(1);  }
    double n2Fun  (double t) { CalcState(t); return AI->n(2);  }
    double nu0Fun (double t) { CalcState(t); return AI->nu(0); }
    double nu1Fun (double t) { CalcState(t); return AI->nu(1); }
    double nu2Fun (double t) { CalcState(t); return AI->nu(2); }

    // Data
    int         test;              // test number
    Vec_t       Sig0, Sig, dSigdt; // stress
    Mat_t       M, dMdt;           // multiplier
    AnisoInvs * AI;                // invariants
    double      dspdt, dsqdt;      // derivs
    Vec_t       dsdt;              // derivs
    Vec_t       dNdt, dNudt;       // derivs
    Vec_t       dndt, dnudt;       // derivs
};

typedef double (Problem::*pFun) (double t);

int main(int argc, char **argv) try
{
    // initialize problem
    Problem       prob;
    Diff<Problem> nd(&prob);
    bool   verbose = true;
    size_t ndiv    = 20;
    if (argc>1) prob.test = atoi(argv[1]);
    if (argc>2) verbose   = atoi(argv[2]);
    if (argc>3) ndiv      = atoi(argv[3]);



    // L0, L1, L2 and derivatives
    if (verbose)
    {
        printf("\n");
        printf("%6s %12s %12s %12s %16s  %12s %12s %12s %16s  %12s %12s %12s %16s\n",
                "t","L0","dL0dt_num","dL0dt1","err(dL0dt1)",
                    "L1","dL1dt_num","dL1dt1","err(dL1dt1)",
                    "L2","dL2dt_num","dL2dt1","err(dL2dt1)");
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
        double err_dL0dt1 = fabs(dL0dt_num    - prob.dsdt(0));
        double err_dL1dt1 = fabs(dL1dt_num    - prob.dsdt(1));
        double err_dL2dt1 = fabs(dL2dt_num    - prob.dsdt(2));
        if (err_dL0dt1 > max_err_dL0dt) max_err_dL0dt = err_dL0dt1;
        if (err_dL1dt1 > max_err_dL1dt) max_err_dL1dt = err_dL1dt1;
        if (err_dL2dt1 > max_err_dL2dt) max_err_dL2dt = err_dL2dt1;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.AI->L(0), dL0dt_num, prob.dsdt(0), err_dL0dt1,
                       prob.AI->L(1), dL1dt_num, prob.dsdt(1), err_dL1dt1,
                       prob.AI->L(2), dL2dt_num, prob.dsdt(2), err_dL2dt1);
        }
    }
    double tol_dL0dt = 1.0e-6;
    double tol_dL1dt = 1.0e-6;
    double tol_dL2dt = 1.0e-6;
    printf("  max_err_dL0dt  = %s%16.8e%s\n",(max_err_dL0dt >tol_dL0dt ?TERM_RED:TERM_GREEN),max_err_dL0dt, TERM_RST);
    printf("  max_err_dL1dt  = %s%16.8e%s\n",(max_err_dL1dt >tol_dL1dt ?TERM_RED:TERM_GREEN),max_err_dL1dt, TERM_RST);
    printf("  max_err_dL2dt  = %s%16.8e%s\n",(max_err_dL2dt >tol_dL2dt ?TERM_RED:TERM_GREEN),max_err_dL2dt, TERM_RST);



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
    printf("  max_err_dspdt = %s%16.8e%s\n",(max_err_dspdt>tol_dspdt?TERM_RED:TERM_GREEN),max_err_dspdt,TERM_RST);
    printf("  max_err_dsqdt = %s%16.8e%s\n",(max_err_dsqdt>tol_dsqdt?TERM_RED:TERM_GREEN),max_err_dsqdt,TERM_RST);
    printf("\n");



    // dNds
    double max_err_dNdt[3] = {0., 0., 0.};
    pFun   NFun[3] = {&Problem::N0Fun, &Problem::N1Fun, &Problem::N2Fun};
    if (verbose)
    {
        printf("\n%6s","t");
        printf("%12s %12s %16s  %12s %12s %16s  %12s %12s %16s\n", "dN0dt_num","dN0dt","err(dN0dt)", "dN1dt_num","dN1dt","err(dN1dt)", "dN2dt_num","dN2dt","err(dN2dt)");
    }
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t = (double)i/(double)ndiv;
        prob.CalcState (t);
        if (verbose) printf("%6.3f",t);
        for (size_t k=0; k<3; ++k)
        {
            double dNdt_num = nd.DyDx (NFun[k], t);
            double err_dNdt = fabs(dNdt_num - prob.dNdt(k));
            if (err_dNdt > max_err_dNdt[k]) max_err_dNdt[k] = err_dNdt;
            if (verbose) printf("%12.8f %12.8f %16.8e  ", dNdt_num,prob.dNdt(k),err_dNdt);
        }
        if (verbose) printf("\n");
    }
    double tol_dNdt[3] = {1.0e-3, 1.0e-6, 1.0e-6};
    for (size_t k=0; k<3; ++k) printf("  max_err_dN%zddt  = %s%16.8e%s\n",k,(max_err_dNdt[k]>tol_dNdt[k]?TERM_RED:TERM_GREEN),max_err_dNdt[k], TERM_RST);



    // dNuds
    double max_err_dNudt[3] = {0., 0., 0.};
    pFun   NuFun[3] = {&Problem::Nu0Fun, &Problem::Nu1Fun, &Problem::Nu2Fun};
    if (verbose)
    {
        printf("\n%6s","t");
        printf("%12s %12s %16s  %12s %12s %16s  %12s %12s %16s\n", "dNu0dt_num","dNu0dt","err(dNu0dt)", "dNu1dt_num","dNu1dt","err(dNu1dt)", "dNu2dt_num","dNu2dt","err(dNu2dt)");
    }
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t = (double)i/(double)ndiv;
        prob.CalcState (t);
        if (verbose) printf("%6.3f",t);
        for (size_t k=0; k<3; ++k)
        {
            double dNudt_num = nd.DyDx (NuFun[k], t);
            double err_dNudt = fabs(dNudt_num - prob.dNudt(k));
            if (err_dNudt > max_err_dNudt[k]) max_err_dNudt[k] = err_dNudt;
            if (verbose) printf("%12.8f %12.8f %16.8e  ", dNudt_num,prob.dNudt(k),err_dNudt);
        }
        if (verbose) printf("\n");
    }
    double tol_dNudt[3] = {1.0e-6, 1.0e-6, 1.0e-6};
    for (size_t k=0; k<3; ++k) printf("  max_err_dNu%zddt  = %s%16.8e%s\n",k,(max_err_dNudt[k]>tol_dNudt[k]?TERM_RED:TERM_GREEN),max_err_dNudt[k], TERM_RST);



    // dnds
    double max_err_dndt[3] = {0., 0., 0.};
    pFun   nFun[3] = {&Problem::n0Fun, &Problem::n1Fun, &Problem::n2Fun};
    if (verbose)
    {
        printf("\n%6s","t");
        printf("%12s %12s %16s  %12s %12s %16s  %12s %12s %16s\n", "dn0dt_num","dn0dt","err(dn0dt)", "dn1dt_num","dn1dt","err(dn1dt)", "dn2dt_num","dn2dt","err(dn2dt)");
    }
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t = (double)i/(double)ndiv;
        prob.CalcState (t);
        if (verbose) printf("%6.3f",t);
        for (size_t k=0; k<3; ++k)
        {
            double dndt_num = nd.DyDx (nFun[k], t);
            double err_dndt = fabs(dndt_num - prob.dndt(k));
            if (err_dndt > max_err_dndt[k]) max_err_dndt[k] = err_dndt;
            if (verbose) printf("%12.8f %12.8f %16.8e  ", dndt_num,prob.dndt(k),err_dndt);
        }
        if (verbose) printf("\n");
    }
    double tol_dndt[3] = {1.0e-6, 1.0e-6, 1.0e-6};
    for (size_t k=0; k<3; ++k) printf("  max_err_dn%zddt  = %s%16.8e%s\n",k,(max_err_dndt[k]>tol_dndt[k]?TERM_RED:TERM_GREEN),max_err_dndt[k], TERM_RST);



    // dnuds
    double max_err_dnudt[3] = {0., 0., 0.};
    pFun   nuFun[3] = {&Problem::nu0Fun, &Problem::nu1Fun, &Problem::nu2Fun};
    if (verbose)
    {
        printf("\n%6s","t");
        printf("%12s %12s %16s  %12s %12s %16s  %12s %12s %16s\n", "dnu0dt_num","dnu0dt","err(dnu0dt)", "dnu1dt_num","dnu1dt","err(dnu1dt)", "dnu2dt_num","dnu2dt","err(dnu2dt)");
    }
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t = (double)i/(double)ndiv;
        prob.CalcState (t);
        if (verbose) printf("%6.3f",t);
        for (size_t k=0; k<3; ++k)
        {
            double dnudt_num = nd.DyDx (nuFun[k], t);
            double err_dnudt = fabs(dnudt_num - prob.dnudt(k));
            if (err_dnudt > max_err_dnudt[k]) max_err_dnudt[k] = err_dnudt;
            if (verbose) printf("%12.8f %12.8f %16.8e  ", dnudt_num,prob.dnudt(k),err_dnudt);
        }
        if (verbose) printf("\n");
    }
    double tol_dnudt[3] = {1.0e-6, 1.0e-6, 1.0e-6};
    for (size_t k=0; k<3; ++k) printf("  max_err_dnu%zddt  = %s%16.8e%s\n",k,(max_err_dnudt[k]>tol_dnudt[k]?TERM_RED:TERM_GREEN),max_err_dnudt[k], TERM_RST);



    // end
    if (max_err_dL0dt > tol_dL0dt) return 1;
    if (max_err_dL1dt > tol_dL1dt) return 1;
    if (max_err_dL2dt > tol_dL2dt) return 1;
    if (max_err_dspdt > tol_dspdt) return 1;
    if (max_err_dsqdt > tol_dsqdt) return 1;
    return 0;
}
MECHSYS_CATCH
