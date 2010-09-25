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
    Problem ()
    {
        test = 2;
        I   .change_dim(6);
        Sig0.change_dim(6);
        I    =  1.0,  1.0,  1.0,   0.0,      0.0,      0.0;
        Sig0 = -1.5, -2.0, -0.5,  -0.1*SQ2, -0.5*SQ2, -0.8*SQ2;
        Sig   .change_dim(6);
        S     .change_dim(6);
        dSigdt.change_dim(6);
        dSdt  .change_dim(6);
        M     .change_dim(6,6);
        dMdt  .change_dim(6,6);
    }
    void CalcState (double t) const
    {
        if (test==1)
        {
            Sig    = (1.+t*t)*Sig0;
            dSigdt = (2.*t)*Sig0;
        }
        else
        {
            M      = t,    0.,      0.,   0.,   0.,      0.,
                     0.,  t*t,      0.,   0.,   0.,      0.,
                     0.,   0.,  sin(t),   0.,   0.,      0.,
                     0.,   0.,      0.,    t,   0.,      0.,
                     0.,   0.,      0.,   0.,  t*t,      0.,
                     0.,   0.,      0.,   0.,    0., cos(t);
            dMdt   = 1.,   0.,      0.,   0.,   0.,      0.,
                     0., 2.*t,      0.,   0.,   0.,      0.,
                     0.,   0., cos(t),    0.,   0.,      0.,
                     0.,   0.,      0.,   1.,   0.,      0.,
                     0.,   0.,      0.,   0., 2.*t,      0.,
                     0.,   0.,      0.,   0.,   0., -sin(t);
            Sig    = M*Sig0;
            dSigdt = dMdt*Sig0;
        }
        EigenProj (Sig,    L,     P0,  P1,  P2, /*sort*/true);
        //EigenProj (dSigdt, dLdt, dP0, dP1, dP2);
        dLdt(0) = dot(P0,dSigdt);
        dLdt(1) = dot(P1,dSigdt);
        dLdt(2) = dot(P2,dSigdt);
        Dev          (Sig,    S);
        Dev          (dSigdt, dSdt);
        OctInvs      (L, p1, q1, T1, dpdL, dqdL, dTdL);
        OctDerivs    (L, p2, q2, T2, dpqthdL);
        InvOctDerivs (L, p3, q3, T3, dLdpqth);
        Pow2         (S,  SS);
        Dev          (SS, devSS);
        p = Calc_poct (Sig);
        q = Calc_qoct (Sig);
        if (fabs(p-p1) >1.0e-14) throw new Fatal("Error with p=%g, p1=%g",p,p1);
        if (fabs(q-q1) >1.0e-14) throw new Fatal("Error with q=%g, q1=%g",q,q1);
        if (fabs(p-p2) >1.0e-14) throw new Fatal("Error with p=%g, p2=%g",p,p1);
        if (fabs(q-q2) >1.0e-14) throw new Fatal("Error with q=%g, q2=%g",q,q1);
        if (fabs(T1-T2)>1.0e-14) throw new Fatal("Error with T1=%g, T2=%g",T1,T2);
        if (fabs(p-p3) >1.0e-14) throw new Fatal("Error with p=%g, p3=%g",p,p3);
        if (fabs(q-q3) >1.0e-14) throw new Fatal("Error with q=%g, q3=%g",q,q3);
        if (fabs(T1-T3)>1.0e-14) throw new Fatal("Error with T1=%g, T3=%g",T1,T3);
        th      = asin(T1)/3.0;
        c       = -1.0/(pow(q,3.)*cos(3.0*th));
        dthdL   = (1.0/(3.0*cos(3.0*th)))*dTdL;
        dpdt1   = (-1.0/SQ3)*dot(I,dSigdt);
        dpdt2   = dot(dpdL,dLdt);
        dqdt1   = (1.0/q)*dot(S,dSigdt);
        dqdt2   = dot(dqdL,dLdt);
        dthdSig = c*SQ6*devSS + c*q*T1*S;
        dthdt1  = dot(dthdSig,dSigdt);
        dthdt2  = dot(dthdL,dLdt);
        dL0dt2  = dLdpqth(0,0)*dpdt1 + dLdpqth(0,1)*dqdt1 + dLdpqth(0,2)*dthdt1;
        dL1dt2  = dLdpqth(1,0)*dpdt1 + dLdpqth(1,1)*dqdt1 + dLdpqth(1,2)*dthdt1;
        dL2dt2  = dLdpqth(2,0)*dpdt1 + dLdpqth(2,1)*dqdt1 + dLdpqth(2,2)*dthdt1;
        pqth2L (p, q, th, L1);
        for (size_t i=0; i<3; ++i)
        {
            if (fabs(dpdL (i)-dpqthdL(0,i))>1.0e-15) throw new Fatal("Error with dpdL(%d)=%g, dpqthdL(0,%d)=%g", i,dpdL(i),i,dpqthdL(0,i));
            if (fabs(dqdL (i)-dpqthdL(1,i))>1.0e-15) throw new Fatal("Error with dqdL(%d)=%g, dpqthdL(1,%d)=%g", i,dqdL(i),i,dpqthdL(1,i));
            if (fabs(dthdL(i)-dpqthdL(2,i))>1.0e-15) throw new Fatal("Error with dthdL(%d)=%g, dpqthdL(2,%d)=%g",i,dthdL(i),i,dpqthdL(2,i));
        }
    }
    double pFun  (double t) const { CalcState(t); return p;    }
    double qFun  (double t) const { CalcState(t); return q;    }
    double thFun (double t) const { CalcState(t); return th;   }
    double L0Fun (double t) const { CalcState(t); return L(0); }
    double L1Fun (double t) const { CalcState(t); return L(1); }
    double L2Fun (double t) const { CalcState(t); return L(2); }

    // Data
    int    test;
    Vec_t  Sig0, I;
    mutable Mat_t M, dMdt;
    mutable Mat3_t dpqthdL, dLdpqth;
    mutable Vec3_t L, L1, dLdt, dpdL, dqdL, dTdL, dthdL;
    mutable Vec_t  Sig, S, SS, devSS, P0, P1, P2, dSigdt, dSdt, dthdSig;
    //mutable Vec_t dP0, dP1, dP2;
    mutable double c, p, q, th, dpdt1, dqdt1, dqdt2, dpdt2, dthdt1, dthdt2, dL0dt2, dL1dt2, dL2dt2;
    mutable double p1,q1,T1, p2,q2,T2, p3,q3,T3;
};

int main(int argc, char **argv) try
{
    // initialize problem
    Problem       prob;
    Diff<Problem> nd(&prob);
    size_t ndiv = 20;
    if (argc>1) prob.test = atoi(argv[1]);
    if (argc>2) ndiv      = atoi(argv[2]);

    // p, q, t and derivatives
    printf("%6s %12s %12s %12s %12s %16s %16s  %12s %12s %12s %12s %16s %16s  %12s %12s %12s %12s %16s %16s\n",
            "t","p", "dpdt_num", "dpdt1", "dpdt2", "err(dpdt1)", "err(dpdt2)",
                "q", "dqdt_num", "dqdt1", "dqdt2", "err(dqdt1)", "err(dqdt2)",
                "th","dthdt_num","dthdt1","dthdt2","err(dthdt1)","err(dthdt2)");
    double max_err_dpdt  = 0.0;
    double max_err_dqdt  = 0.0;
    double max_err_dthdt = 0.0;
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t         = (double)i/(double)ndiv;
        double dpdt_num  = nd.DyDx (&Problem::pFun,  t);
        double dqdt_num  = nd.DyDx (&Problem::qFun,  t);
        double dthdt_num = nd.DyDx (&Problem::thFun, t);
        prob.CalcState (t);
        double err_dpdt1  = fabs(dpdt_num     - prob.dpdt1);
        double err_dpdt2  = fabs(prob.dpdt1   - prob.dpdt2);
        double err_dqdt1  = fabs(dqdt_num     - prob.dqdt1);
        double err_dqdt2  = fabs(prob.dqdt1   - prob.dqdt2);
        double err_dthdt1 = fabs(dthdt_num    - prob.dthdt1);
        double err_dthdt2 = fabs(prob.dthdt1  - prob.dthdt2);
        if (err_dpdt1  > max_err_dpdt ) max_err_dpdt  = err_dpdt1;
        if (err_dpdt2  > max_err_dpdt ) max_err_dpdt  = err_dpdt2;
        if (err_dqdt1  > max_err_dqdt ) max_err_dqdt  = err_dqdt1;
        if (err_dqdt2  > max_err_dqdt ) max_err_dqdt  = err_dqdt2;
        if (err_dthdt1 > max_err_dthdt) max_err_dthdt = err_dthdt1;
        if (err_dthdt2 > max_err_dthdt) max_err_dthdt = err_dthdt2;
        printf("%6.3f %12.8f %12.8f %12.8f %12.8f %16.8e %16.8e  %12.8f %12.8f %12.8f %12.8f %16.8e %16.8e  %12.8f %12.8f %12.8f %12.8f %16.8e %16.8e\n",
                t, prob.p,    dpdt_num,  prob.dpdt1,   prob.dpdt2,  err_dpdt1,  err_dpdt2,
                   prob.q,    dqdt_num,  prob.dqdt1,   prob.dqdt2,  err_dqdt1,  err_dqdt2,
                   prob.th,   dthdt_num, prob.dthdt1,  prob.dthdt2, err_dthdt1, err_dthdt2);
    }

    // L0, L1, L2 and derivatives
    printf("\n");
    printf("%6s %12s %12s %12s %12s %12s %16s %16s  %12s %12s %12s %12s %12s %16s %16s  %12s %12s %12s %12s %12s %16s %16s\n",
            "t","L0","L0","dL0dt_num","dL0dt1","dL0dt2","err(dL0dt1)","err(dL0dt2)",
                "L1","L1","dL1dt_num","dL1dt1","dL1dt2","err(dL1dt1)","err(dL1dt2)",
                "L2","L2","dL2dt_num","dL2dt1","dL2dt2","err(dL2dt1)","err(dL2dt2)");
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
        double err_dL0dt1 = fabs(dL0dt_num    - prob.dLdt(0));
        double err_dL1dt1 = fabs(dL1dt_num    - prob.dLdt(1));
        double err_dL2dt1 = fabs(dL2dt_num    - prob.dLdt(2));
        double err_dL0dt2 = fabs(prob.dLdt(0) - prob.dL0dt2);
        double err_dL1dt2 = fabs(prob.dLdt(1) - prob.dL1dt2);
        double err_dL2dt2 = fabs(prob.dLdt(2) - prob.dL2dt2);
        if (err_dL0dt1 > max_err_dL0dt) max_err_dL0dt = err_dL0dt1;
        if (err_dL0dt2 > max_err_dL0dt) max_err_dL0dt = err_dL0dt2;
        if (err_dL1dt1 > max_err_dL1dt) max_err_dL1dt = err_dL1dt1;
        if (err_dL1dt2 > max_err_dL1dt) max_err_dL1dt = err_dL1dt2;
        if (err_dL2dt1 > max_err_dL2dt) max_err_dL2dt = err_dL2dt1;
        if (err_dL2dt2 > max_err_dL2dt) max_err_dL2dt = err_dL2dt2;
        printf("%6.3f %12.8f %12.8f %12.8f %12.8f %12.8f %16.8e %16.8e  %12.8f %12.8f %12.8f %12.8f %12.8f %16.8e %16.8e  %12.8f %12.8f %12.8f %12.8f %12.8f %16.8e %16.8e\n",
                t, prob.L(0), prob.L1(0), dL0dt_num, prob.dLdt(0), prob.dL0dt2, err_dL0dt1, err_dL0dt2,
                   prob.L(1), prob.L1(1), dL1dt_num, prob.dLdt(1), prob.dL1dt2, err_dL1dt1, err_dL1dt2,
                   prob.L(2), prob.L1(2), dL2dt_num, prob.dLdt(2), prob.dL2dt2, err_dL2dt1, err_dL2dt2);
    }

    // error
    double tol_dpdt  = 1.0e-6;
    double tol_dqdt  = 1.0e-7;
    double tol_dthdt = 1.0e-7;
    double tol_dL0dt = 1.0e-6;
    double tol_dL1dt = 1.0e-6;
    double tol_dL2dt = 1.0e-6;
    printf("\n");
    printf("  max_err_dpdt  = %s%16.8e%s\n",(max_err_dpdt >tol_dpdt ?TERM_RED:TERM_GREEN),max_err_dpdt, TERM_RST);
    printf("  max_err_dqdt  = %s%16.8e%s\n",(max_err_dqdt >tol_dqdt ?TERM_RED:TERM_GREEN),max_err_dqdt, TERM_RST);
    printf("  max_err_dthdt = %s%16.8e%s\n",(max_err_dthdt>tol_dthdt?TERM_RED:TERM_GREEN),max_err_dthdt,TERM_RST);
    printf("  max_err_dL0dt = %s%16.8e%s\n",(max_err_dL0dt>tol_dL0dt?TERM_RED:TERM_GREEN),max_err_dL0dt,TERM_RST);
    printf("  max_err_dL1dt = %s%16.8e%s\n",(max_err_dL1dt>tol_dL1dt?TERM_RED:TERM_GREEN),max_err_dL1dt,TERM_RST);
    printf("  max_err_dL2dt = %s%16.8e%s\n",(max_err_dL2dt>tol_dL2dt?TERM_RED:TERM_GREEN),max_err_dL2dt,TERM_RST);

    // end
    if (max_err_dpdt  > tol_dpdt)  return 1;
    if (max_err_dqdt  > tol_dqdt)  return 1;
    if (max_err_dthdt > tol_dthdt) return 1;
    if (max_err_dL0dt > tol_dL0dt) return 1;
    if (max_err_dL1dt > tol_dL1dt) return 1;
    if (max_err_dL2dt > tol_dL2dt) return 1;
    return 0;
}
MECHSYS_CATCH
