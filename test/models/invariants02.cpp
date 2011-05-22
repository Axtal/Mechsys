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
#ifdef HAS_TENSORS
  using namespace TensorsLib;
#endif

class Problem
{
public:
    Problem ()
        : test    (2),
          use_iod (false),
          sort    (false)
    {
        I   .change_dim(6);
        Sig0.change_dim(6);
        I    =  1.0,  1.0,  1.0,   0.0,      0.0,      0.0;
        Sig0 = -1.5, -2.0, -0.5,  -0.1*SQ2, -0.5*SQ2, -0.8*SQ2;
        //Sig0 = -1.5, -2.0, -0.5,  -0.1*SQ2, -0.5*SQ2, -0.00000001*SQ2;
        Sig   .change_dim(6);
        S     .change_dim(6);
        dSigdt.change_dim(6);
        dSdt  .change_dim(6);
        invSig.change_dim(6);
        diSdt .change_dim(6);
        M     .change_dim(6,6);
        dMdt  .change_dim(6,6);
        diSig_dSig.change_dim(6,6);
        dPdt.Resize(3);
#ifdef HAS_TENSORS
        dvdt.Resize(3);
#endif
    }
    void CalcState (double t)
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
            Sig    = 0.1*Sig0 + M*Sig0;
            dSigdt = dMdt*Sig0;
        }

        // characteristic invariants
        CharInvs (Sig, I1,I2,I3, dI1dsig,dI2dsig,dI3dsig);
        dI1dt = dot (dI1dsig , dSigdt);
        dI2dt = dot (dI2dsig , dSigdt);
        dI3dt = dot (dI3dsig , dSigdt);

        EigenProjDerivs (Sig, L, v0,v1,v2, P0,P1,P2, dP0dSig, dP1dSig, dP2dSig, sort);
        //EigenProj (Sig,    L,     P0,  P1,  P2, sort);
        //EigenProj (dSigdt, dLdt, dP0, dP1, dP2);
        dLdt(0) = dot(P0,dSigdt);
        dLdt(1) = dot(P1,dSigdt);
        dLdt(2) = dot(P2,dSigdt);
        Dev          (Sig,    S);
        Dev          (dSigdt, dSdt);
        OctInvs      (L, p1, q1, T1, dpdL, dqdL, dTdL);
        OctInvs      (Sig, p4, q4, T4, th4, devSig4, 1.0e-8, &dthdSig4);
        OctDerivs    (L, p2, q2, T2, dpqthdL);
        if (use_iod) InvOctDerivs (L, p3, q3, T3, dLdpqth); // L need to be sorted
        else         Inv          (dpqthdL, dLdpqth);
        Pow2         (S,  SS);
        Dev          (SS, devSS);
        Inv          (Sig, invSig);
        DerivInv     (Sig, invSig, diSig_dSig);
        p = Calc_poct (Sig);
        q = Calc_qoct (Sig);
        Vec_t dif_devSig4(S - devSig4);
        double err_devSig4 = Norm(dif_devSig4);
        if (err_devSig4>1.0e-14) throw new Fatal("Error with devSig4 = %g",err_devSig4);
        if (fabs(p-p1) >1.0e-14) throw new Fatal("Error with p=%g, p1=%g",p,p1);
        if (fabs(q-q1) >1.0e-14) throw new Fatal("Error with q=%g, q1=%g",q,q1);
        if (fabs(p-p2) >1.0e-14) throw new Fatal("Error with p=%g, p2=%g",p,p1);
        if (fabs(q-q2) >1.0e-14) throw new Fatal("Error with q=%g, q2=%g",q,q1);
        if (fabs(p-p4) >1.0e-14) throw new Fatal("Error with p=%g, p4=%g",p,p4);
        if (fabs(q-q4) >1.0e-14) throw new Fatal("Error with q=%g, q4=%g",q,q4);
        if (fabs(T1-T2)>1.0e-14) throw new Fatal("Error with T1=%g, T2=%g",T1,T2);
        if (fabs(T1-T4)>1.0e-14) throw new Fatal("Error with T1=%g, T4=%g",T1,T4);
        if (use_iod)
        {
            if (fabs(p-p3) >1.0e-14) throw new Fatal("Error with p=%g, p3=%g",p,p3);
            if (fabs(q-q3) >1.0e-14) throw new Fatal("Error with q=%g, q3=%g",q,q3);
            if (fabs(T1-T3)>1.0e-14) throw new Fatal("Error with T1=%g, T3=%g",T1,T3);
        }
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
        diSdt   = diSig_dSig * dSigdt;
        pqth2L (p, q, th, L1);
        Vec_t dif_dthdSig(dthdSig - dthdSig4);
        double err_dthdSig = Norm(dif_dthdSig);
        if (err_dthdSig>1.0e-14) throw new Fatal("Error with dthdSig4 = %g",err_dthdSig);
        for (size_t i=0; i<3; ++i)
        {
            if (fabs(dpdL (i)-dpqthdL(0,i))>1.0e-15) throw new Fatal("Error with dpdL(%d)=%g, dpqthdL(0,%d)=%g", i,dpdL(i),i,dpqthdL(0,i));
            if (fabs(dqdL (i)-dpqthdL(1,i))>1.0e-15) throw new Fatal("Error with dqdL(%d)=%g, dpqthdL(1,%d)=%g", i,dqdL(i),i,dpqthdL(1,i));
            if (fabs(dthdL(i)-dpqthdL(2,i))>1.0e-13) throw new Fatal("Error with dthdL(%d)=%g, dpqthdL(2,%d)=%g",i,dthdL(i),i,dpqthdL(2,i));
        }

#ifdef HAS_TENSORS
        // derivative of inverse w.r.t. sig
        Tensor2<double,3> A, iA, dAdt, dInvAdt;
        Tensor4<double,3> dInvAdA;
        Ten2Tensor (Sig,    A);
        Ten2Tensor (dSigdt, dAdt);
        Inv        (A,      iA);
        dInvAdA = ((iA^iA) + (iA|iA)) * (-0.5);
        dInvAdt = (dInvAdA % dAdt);
        Tensor2Ten (dInvAdt, diAdt, /*ncp*/6);
        Mat_t dInvAdAmat;
        Tensor2Ten (dInvAdA,dInvAdAmat,6);
        for (size_t i=0; i<6; ++i)
        for (size_t j=0; j<6; ++j)
        {
            double diff = fabs(diSig_dSig(i,j) - dInvAdAmat(i,j));
            if (diff>1.0e-9) throw new Fatal("diSig_dSig is not equal to dInvAdAmat. Error=%17.9e",diff);
        }

        // derivative of eigenvectors
        Ten1_t V0, V1, V2;
        Ten3_t dv0dSig, dv1dSig, dv2dSig;
        Vec2Tensor     (v0,V0);
        Vec2Tensor     (v1,V1);
        Vec2Tensor     (v2,V2);
        EigenVecDerivs (L, V0,V1,V2, dv0dSig,dv1dSig,dv2dSig);
        dvdt[0] = dv0dSig % dAdt;
        dvdt[1] = dv1dSig % dAdt;
        dvdt[2] = dv2dSig % dAdt;
#endif
        // derivative of eigenprojectors
        dPdt[0] = dP0dSig * dSigdt;
        dPdt[1] = dP1dSig * dSigdt;
        dPdt[2] = dP2dSig * dSigdt;
    }

    double I1Fun  (double t) { CalcState(t); return I1;        }
    double I2Fun  (double t) { CalcState(t); return I2;        }
    double I3Fun  (double t) { CalcState(t); return I3;        }

    double pFun   (double t) { CalcState(t); return p;         }
    double qFun   (double t) { CalcState(t); return q;         }
    double thFun  (double t) { CalcState(t); return th;        }
    double L0Fun  (double t) { CalcState(t); return L(0);      }
    double L1Fun  (double t) { CalcState(t); return L(1);      }
    double L2Fun  (double t) { CalcState(t); return L(2);      }
    double iS0Fun (double t) { CalcState(t); return invSig(0); }
    double iS1Fun (double t) { CalcState(t); return invSig(1); }
    double iS2Fun (double t) { CalcState(t); return invSig(2); }
    double iS3Fun (double t) { CalcState(t); return invSig(3); }
    double iS4Fun (double t) { CalcState(t); return invSig(4); }
    double iS5Fun (double t) { CalcState(t); return invSig(5); }

    double P00Fun (double t) { CalcState(t); return P0(0);     }
    double P01Fun (double t) { CalcState(t); return P0(1);     }
    double P02Fun (double t) { CalcState(t); return P0(2);     }
    double P03Fun (double t) { CalcState(t); return P0(3);     }
    double P04Fun (double t) { CalcState(t); return P0(4);     }
    double P05Fun (double t) { CalcState(t); return P0(5);     }

    double P10Fun (double t) { CalcState(t); return P1(0);     }
    double P11Fun (double t) { CalcState(t); return P1(1);     }
    double P12Fun (double t) { CalcState(t); return P1(2);     }
    double P13Fun (double t) { CalcState(t); return P1(3);     }
    double P14Fun (double t) { CalcState(t); return P1(4);     }
    double P15Fun (double t) { CalcState(t); return P1(5);     }

    double P20Fun (double t) { CalcState(t); return P2(0);     }
    double P21Fun (double t) { CalcState(t); return P2(1);     }
    double P22Fun (double t) { CalcState(t); return P2(2);     }
    double P23Fun (double t) { CalcState(t); return P2(3);     }
    double P24Fun (double t) { CalcState(t); return P2(4);     }
    double P25Fun (double t) { CalcState(t); return P2(5);     }

    double v00Fun (double t) { CalcState(t); return v0(0);     }
    double v01Fun (double t) { CalcState(t); return v0(1);     }
    double v02Fun (double t) { CalcState(t); return v0(2);     }

    double v10Fun (double t) { CalcState(t); return v1(0);     }
    double v11Fun (double t) { CalcState(t); return v1(1);     }
    double v12Fun (double t) { CalcState(t); return v1(2);     }

    double v20Fun (double t) { CalcState(t); return v2(0);     }
    double v21Fun (double t) { CalcState(t); return v2(1);     }
    double v22Fun (double t) { CalcState(t); return v2(2);     }

    // Data
    int    test;
    bool   use_iod; // Use: InvOctDerivs ?
    bool   sort;    // sort eigenvalues ?
    Vec_t  Sig0, I;
    Vec_t  invSig, diSdt, diAdt;
    Mat_t  M, dMdt, diSig_dSig;
    Mat3_t dpqthdL, dLdpqth;
    Vec3_t L, L1, dLdt, dpdL, dqdL, dTdL, dthdL;
    Vec3_t v0,v1,v2;
    Vec_t  Sig, S, SS, devSS, P0, P1, P2, dSigdt, dSdt, dthdSig;
    double c, p, q, th, dpdt1, dqdt1, dqdt2, dpdt2, dthdt1, dthdt2, dL0dt2, dL1dt2, dL2dt2;
    double p1,q1,T1, p2,q2,T2, p3,q3,T3;
    Mat_t  dP0dSig, dP1dSig, dP2dSig;
    double p4, q4, T4, th4;
    Vec_t  devSig4, dthdSig4;
    Vec_t  dI1dsig,dI2dsig,dI3dsig; ///< derivative of characteristic invariants
    double I1, I2, I3, dI1dt, dI2dt, dI3dt; ///< Invariants

    Array<Vec_t> dPdt;
#ifdef HAS_TENSORS
    Array<Ten1_t> dvdt;
#endif
};

typedef double (Problem::*pFun) (double t);

int main(int argc, char **argv) try
{
    // initialize problem
    Problem       prob;
    Diff<Problem> nd(&prob);
    bool   verbose = false;
    size_t ndiv    = 20;
    if (argc>1) prob.test    = atoi(argv[1]);
    if (argc>2) prob.use_iod = atoi(argv[2]);
    if (argc>3) prob.sort    = atoi(argv[3]);
    if (argc>4) verbose      = atoi(argv[4]);
    if (argc>5) ndiv         = atoi(argv[5]);

    if (prob.use_iod) printf("  . . . Using InvOctDerivs method\n");
    if (prob.sort)    printf("  . . . Sorting eigenvalues\n");

    // I1, I2, I3 and derivatives
    if (verbose)
    {
        printf("%6s %12s %12s %12s %16s  %12s %12s %12s %16s  %12s %12s %12s %16s\n",
                "t", "I1", "dI1dt_num", "dI1dt", "err(dI1dt)",
                     "I2", "dI2dt_num", "dI2dt", "err(dI2dt)",
                     "I3", "dI3dt_num", "dI3dt", "err(dI3dt)");
    }
    double max_err_dI1dt = 0.0;
    double max_err_dI2dt = 0.0;
    double max_err_dI3dt = 0.0;
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t         = (double)i/(double)ndiv;
        double dI1dt_num = nd.DyDx (&Problem::I1Fun, t);
        double dI2dt_num = nd.DyDx (&Problem::I2Fun, t);
        double dI3dt_num = nd.DyDx (&Problem::I3Fun, t);
        prob.CalcState (t);
        double err_dI1dt = fabs(dI1dt_num - prob.dI1dt);
        double err_dI2dt = fabs(dI2dt_num - prob.dI2dt);
        double err_dI3dt = fabs(dI3dt_num - prob.dI3dt);
        if (err_dI1dt > max_err_dI1dt) max_err_dI1dt = err_dI1dt;
        if (err_dI2dt > max_err_dI2dt) max_err_dI2dt = err_dI2dt;
        if (err_dI3dt > max_err_dI3dt) max_err_dI3dt = err_dI3dt;
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e  %12.8f %12.8f %12.8f %16.8e\n",
                    t, prob.I1, dI1dt_num, prob.dI1dt, err_dI1dt,
                       prob.I2, dI2dt_num, prob.dI2dt, err_dI2dt,
                       prob.I3, dI3dt_num, prob.dI3dt, err_dI3dt);
        }
    }

    // p, q, t and derivatives
    if (verbose)
    {
        printf("\n");
        printf("%6s %12s %12s %12s %12s %16s %16s  %12s %12s %12s %12s %16s %16s  %12s %12s %12s %12s %16s %16s\n",
                "t","p", "dpdt_num", "dpdt1", "dpdt2", "err(dpdt1)", "err(dpdt2)",
                    "q", "dqdt_num", "dqdt1", "dqdt2", "err(dqdt1)", "err(dqdt2)",
                    "th","dthdt_num","dthdt1","dthdt2","err(dthdt1)","err(dthdt2)");
    }
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
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %12.8f %16.8e %16.8e  %12.8f %12.8f %12.8f %12.8f %16.8e %16.8e  %12.8f %12.8f %12.8f %12.8f %16.8e %16.8e\n",
                    t, prob.p,    dpdt_num,  prob.dpdt1,   prob.dpdt2,  err_dpdt1,  err_dpdt2,
                       prob.q,    dqdt_num,  prob.dqdt1,   prob.dqdt2,  err_dqdt1,  err_dqdt2,
                       prob.th,   dthdt_num, prob.dthdt1,  prob.dthdt2, err_dthdt1, err_dthdt2);
        }
    }

    // L0, L1, L2 and derivatives
    if (verbose)
    {
        printf("\n");
        printf("%6s %12s %12s %12s %12s %12s %16s %16s  %12s %12s %12s %12s %12s %16s %16s  %12s %12s %12s %12s %12s %16s %16s\n",
                "t","L0","L0","dL0dt_num","dL0dt1","dL0dt2","err(dL0dt1)","err(dL0dt2)",
                    "L1","L1","dL1dt_num","dL1dt1","dL1dt2","err(dL1dt1)","err(dL1dt2)",
                    "L2","L2","dL2dt_num","dL2dt1","dL2dt2","err(dL2dt1)","err(dL2dt2)");
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
        if (verbose)
        {
            printf("%6.3f %12.8f %12.8f %12.8f %12.8f %12.8f %16.8e %16.8e  %12.8f %12.8f %12.8f %12.8f %12.8f %16.8e %16.8e  %12.8f %12.8f %12.8f %12.8f %12.8f %16.8e %16.8e\n",
                    t, prob.L(0), prob.L1(0), dL0dt_num, prob.dLdt(0), prob.dL0dt2, err_dL0dt1, err_dL0dt2,
                       prob.L(1), prob.L1(1), dL1dt_num, prob.dLdt(1), prob.dL1dt2, err_dL1dt1, err_dL1dt2,
                       prob.L(2), prob.L1(2), dL2dt_num, prob.dLdt(2), prob.dL2dt2, err_dL2dt1, err_dL2dt2);
        }
    }

    // inverse
    if (verbose)
    {
        printf("\n");
        printf("%6s %14s %14s %14s %16s  %14s %14s %16s  %14s %14s %16s  %14s %14s %16s  %14s %14s %16s  %14s %14s %16s\n",
                "t","diS0dt_num","diS0dt","diA0dt","err(diS0dt)",
                    "diS1dt_num","diS1dt","err(diS1dt)",
                    "diS2dt_num","diS2dt","err(diS2dt)",
                    "diS3dt_num","diS3dt","err(diS3dt)",
                    "diS4dt_num","diS4dt","err(diS4dt)",
                    "diS5dt_num","diS5dt","err(diS5dt)");
    }
    double max_err_diS0dt = 0.0;
    double max_err_diS1dt = 0.0;
    double max_err_diS2dt = 0.0;
    double max_err_diS3dt = 0.0;
    double max_err_diS4dt = 0.0;
    double max_err_diS5dt = 0.0;
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t          = (double)i/(double)ndiv;
        double diS0dt_num = nd.DyDx (&Problem::iS0Fun, t);
        double diS1dt_num = nd.DyDx (&Problem::iS1Fun, t);
        double diS2dt_num = nd.DyDx (&Problem::iS2Fun, t);
        double diS3dt_num = nd.DyDx (&Problem::iS3Fun, t);
        double diS4dt_num = nd.DyDx (&Problem::iS4Fun, t);
        double diS5dt_num = nd.DyDx (&Problem::iS5Fun, t);
        prob.CalcState (t);
        double err_diS0dt = fabs(diS0dt_num    - prob.diSdt(0));
        double err_diA0dt = fabs(prob.diSdt(0) - prob.diAdt(0));
        double err_diS1dt = fabs(diS1dt_num    - prob.diSdt(1));
        double err_diS2dt = fabs(diS2dt_num    - prob.diSdt(2));
        double err_diS3dt = fabs(diS3dt_num    - prob.diSdt(3));
        double err_diS4dt = fabs(diS4dt_num    - prob.diSdt(4));
        double err_diS5dt = fabs(diS5dt_num    - prob.diSdt(5));
        if (err_diS0dt > max_err_diS0dt) max_err_diS0dt = err_diS0dt;
        if (err_diA0dt > max_err_diS0dt) max_err_diS0dt = err_diA0dt;
        if (err_diS1dt > max_err_diS1dt) max_err_diS1dt = err_diS1dt;
        if (err_diS2dt > max_err_diS2dt) max_err_diS2dt = err_diS2dt;
        if (err_diS3dt > max_err_diS3dt) max_err_diS3dt = err_diS3dt;
        if (err_diS4dt > max_err_diS4dt) max_err_diS4dt = err_diS4dt;
        if (err_diS5dt > max_err_diS5dt) max_err_diS5dt = err_diS5dt;
        if (verbose)
        {
            printf("%6.3f %14.8f %14.8f %14.8f %16.8e  %14.8f %14.8f %16.8e  %14.8f %14.8f %16.8e  %14.8f %14.8f %16.8e  %14.8f %14.8f %16.8e  %14.8f %14.8f %16.8e\n",
                    t, diS0dt_num, prob.diSdt(0), prob.diAdt(0), err_diS0dt,
                       diS1dt_num, prob.diSdt(1),                err_diS1dt,
                       diS2dt_num, prob.diSdt(2),                err_diS2dt,
                       diS3dt_num, prob.diSdt(3),                err_diS3dt,
                       diS4dt_num, prob.diSdt(4),                err_diS4dt,
                       diS5dt_num, prob.diSdt(5),                err_diS5dt);
        }
    }

    // eigenprojectors
    double max_err_dPdt[3][6] = {{0.,0.,0.,0.,0.,0.},
                                 {0.,0.,0.,0.,0.,0.},
                                 {0.,0.,0.,0.,0.,0.}};
    pFun funcs[3][6] = {{&Problem::P00Fun, &Problem::P01Fun, &Problem::P02Fun, &Problem::P03Fun, &Problem::P04Fun, &Problem::P05Fun},
                        {&Problem::P10Fun, &Problem::P11Fun, &Problem::P12Fun, &Problem::P13Fun, &Problem::P14Fun, &Problem::P15Fun},
                        {&Problem::P20Fun, &Problem::P21Fun, &Problem::P22Fun, &Problem::P23Fun, &Problem::P24Fun, &Problem::P25Fun}};
    for (size_t k=0; k<3; ++k)
    {
        if (verbose)
        {
            printf("\n%6s","t");
            for (size_t j=0; j<6; ++j)
            {
                char str0[32];
                char str1[32];
                char str2[32];
                sprintf(str0,"dP%zd%zddt_num",   k,j);
                sprintf(str1,"dP%zd%zddt",       k,j);
                sprintf(str2,"error(dP%zd%zddt)",k,j);
                printf("%12s %12s %16s  ",str0,str1,str2);
            }
            printf("\n");
        }
        for (size_t i=0; i<ndiv+1; ++i)
        {
            double t = (double)i/(double)ndiv;
            prob.CalcState (t);
            if (verbose) printf("%6.3f",t);
            for (size_t j=0; j<6; ++j)
            {
                double dPdt_num = nd.DyDx (funcs[k][j], t);
                double err      = fabs(dPdt_num - prob.dPdt[k](j));
                if (err > max_err_dPdt[k][j]) max_err_dPdt[k][j] = err;
                if (verbose) printf("%12.8f %12.8f %16.8e  ", dPdt_num, prob.dPdt[k](j), err);
            }
            if (verbose) printf("\n");
        }
    }

#ifdef HAS_TENSORS
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
            printf("\n%6s","t");
            for (size_t j=0; j<3; ++j)
            {
                char str0[32];
                char str1[32];
                char str2[32];
                sprintf(str0,"dv%zd%zddt_num",   k,j);
                sprintf(str1,"dv%zd%zddt",       k,j);
                sprintf(str2,"error(dv%zd%zddt)",k,j);
                printf("%12s %12s %16s  ",str0,str1,str2);
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
                if (verbose) printf("%12.8f %12.8f %16.8e  ", dvdt_num, prob.dvdt[k][j], err);
            }
            if (verbose) printf("\n");
        }
    }
#endif

    // error
    double tol_dI1dt  = 1.0e-6;
    double tol_dI2dt  = 1.0e-6;
    double tol_dI3dt  = 1.0e-7;
    double tol_dpdt   = 1.0e-6;
    double tol_dqdt   = 1.0e-7;
    double tol_dthdt  = 1.0e-6;
    double tol_dL0dt  = 1.0e-6;
    double tol_dL1dt  = 1.0e-6;
    double tol_dL2dt  = 1.0e-6;
    double tol_diS0dt = 1.0e-5;
    double tol_diS1dt = 1.0e-5;
    double tol_diS2dt = 1.0e-4;
    double tol_diS3dt = 1.0e-5;
    double tol_diS4dt = 1.0e-4;
    double tol_diS5dt = 1.0e-4;
    double tol_dPdt[3][6]= {{1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-7, 1.0e-6},
                            {1.0e-6, 1.0e-6, 1.0e-7, 1.0e-6, 1.0e-6, 1.0e-6},
                            {1.0e-6, 1.0e-6, 1.0e-7, 1.0e-6, 1.0e-6, 1.0e-7}};
    printf("\n");
    printf("  max_err_dI1dt  = %s%16.8e%s\n",(max_err_dI1dt >tol_dI1dt ?TERM_RED:TERM_GREEN),max_err_dI1dt, TERM_RST);
    printf("  max_err_dI2dt  = %s%16.8e%s\n",(max_err_dI2dt >tol_dI2dt ?TERM_RED:TERM_GREEN),max_err_dI2dt, TERM_RST);
    printf("  max_err_dI3dt  = %s%16.8e%s\n",(max_err_dI3dt >tol_dI3dt ?TERM_RED:TERM_GREEN),max_err_dI3dt, TERM_RST);
    printf("  max_err_dpdt   = %s%16.8e%s\n",(max_err_dpdt  >tol_dpdt  ?TERM_RED:TERM_GREEN),max_err_dpdt,  TERM_RST);
    printf("  max_err_dqdt   = %s%16.8e%s\n",(max_err_dqdt  >tol_dqdt  ?TERM_RED:TERM_GREEN),max_err_dqdt,  TERM_RST);
    printf("  max_err_dthdt  = %s%16.8e%s\n",(max_err_dthdt >tol_dthdt ?TERM_RED:TERM_GREEN),max_err_dthdt, TERM_RST);
    printf("  max_err_dL0dt  = %s%16.8e%s\n",(max_err_dL0dt >tol_dL0dt ?TERM_RED:TERM_GREEN),max_err_dL0dt, TERM_RST);
    printf("  max_err_dL1dt  = %s%16.8e%s\n",(max_err_dL1dt >tol_dL1dt ?TERM_RED:TERM_GREEN),max_err_dL1dt, TERM_RST);
    printf("  max_err_dL2dt  = %s%16.8e%s\n",(max_err_dL2dt >tol_dL2dt ?TERM_RED:TERM_GREEN),max_err_dL2dt, TERM_RST);
    printf("  max_err_diS0dt = %s%16.8e%s\n",(max_err_diS0dt>tol_diS0dt?TERM_RED:TERM_GREEN),max_err_diS0dt,TERM_RST);
    printf("  max_err_diS1dt = %s%16.8e%s\n",(max_err_diS1dt>tol_diS1dt?TERM_RED:TERM_GREEN),max_err_diS1dt,TERM_RST);
    printf("  max_err_diS2dt = %s%16.8e%s\n",(max_err_diS2dt>tol_diS2dt?TERM_RED:TERM_GREEN),max_err_diS2dt,TERM_RST);
    printf("  max_err_diS3dt = %s%16.8e%s\n",(max_err_diS3dt>tol_diS3dt?TERM_RED:TERM_GREEN),max_err_diS3dt,TERM_RST);
    printf("  max_err_diS4dt = %s%16.8e%s\n",(max_err_diS4dt>tol_diS4dt?TERM_RED:TERM_GREEN),max_err_diS4dt,TERM_RST);
    printf("  max_err_diS5dt = %s%16.8e%s\n",(max_err_diS5dt>tol_diS5dt?TERM_RED:TERM_GREEN),max_err_diS5dt,TERM_RST);
    for (size_t k=0; k<3; ++k)
    for (size_t i=0; i<6; ++i)
        printf("  max_err_dP%zd%zddt = %s%16.8e%s\n",k,i,(max_err_dPdt[k][i]>tol_dPdt[k][i]?TERM_RED:TERM_GREEN),max_err_dPdt[k][i],TERM_RST);
#ifdef HAS_TENSORS
    double tol_dvdt[3][3]= {{1.0e-7, 1.0e-7, 1.0e-7},
                            {1.0e-6, 1.0e-6, 1.0e-7},
                            {1.0e-6, 1.0e-6, 1.0e-6}};
    for (size_t k=0; k<3; ++k)
    for (size_t i=0; i<3; ++i)
        printf("  max_err_dv%zd%zddt = %s%16.8e%s\n",k,i,(max_err_dvdt[k][i]>tol_dvdt[k][i]?TERM_RED:TERM_GREEN),max_err_dvdt[k][i],TERM_RST);
#endif
    printf("\n");

    // end
    if (max_err_dpdt   > tol_dpdt)   return 1;
    if (max_err_dqdt   > tol_dqdt)   return 1;
    if (max_err_dthdt  > tol_dthdt)  return 1;
    if (max_err_dL0dt  > tol_dL0dt)  return 1;
    if (max_err_dL1dt  > tol_dL1dt)  return 1;
    if (max_err_dL2dt  > tol_dL2dt)  return 1;
    if (max_err_diS0dt > tol_diS0dt) return 1;
    if (max_err_diS1dt > tol_diS1dt) return 1;
    if (max_err_diS2dt > tol_diS2dt) return 1;
    if (max_err_diS3dt > tol_diS3dt) return 1;
    if (max_err_diS4dt > tol_diS4dt) return 1;
    if (max_err_diS5dt > tol_diS5dt) return 1;
    for (size_t k=0; k<3; ++k)
    for (size_t i=0; i<6; ++i)
        if (max_err_dPdt[k][i] > tol_dPdt[k][i]) return 1;
#ifdef HAS_TENSORS
    for (size_t k=0; k<3; ++k)
    for (size_t i=0; i<3; ++i)
        if (max_err_dvdt[k][i] > tol_dvdt[k][i]) return 1;
#endif
    return 0;
}
MECHSYS_CATCH
