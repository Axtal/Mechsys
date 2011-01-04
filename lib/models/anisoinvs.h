/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

#ifndef MECHSYS_ANISOINVS_H
#define MECHSYS_ANISOINVS_H

// Std Lib
#include <iostream>
#include <fstream>

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>

class AnisoInvs
{
public:
    // Constructor & Destructor
     AnisoInvs (double b, double Alpha, Vec3_t const & a, bool Obliq=true);
    ~AnisoInvs () {}

    // Methods
    void Calc (Vec_t const & Sig, bool WithDerivs=false);

    // Constants
    double  b, Alpha;
    Vec3_t  a;
    bool    Obliq;
    double  Tol;
    Vec3_t  au;

    // Data
    static double  sp, sq;      // invariants
    static Vec3_t  N123, N, Nu; // SMP
    static Vec3_t  n, nu;       // AMP
    static Mat3_t  Q;           // rotation matrix
    static Mat3_t  mSig;        // matrix of Sig
    static Vec3_t  L;           // eigenvalues
    static Vec3_t  v0,v1,v2;    // eigenvectors
    static Vec3_t  t, p, q;     // traction and projections
    static Mat3_t  P;           // projector

    // Derivs
#ifdef HAS_TENSORS
    static Ten1_t V0, V1, V2;                // eigenvectors
    static Ten2_t E0, E1, E2;                // eigenprojectors
    static Ten3_t dv0dSig, dv1dSig, dv2dSig; // eigenvectors
    static Ten2_t dN123dL;                   // SMP(123)
    static Ten3_t dNdSig, dNudSig;           // SMP
    static Ten2_t dNudN;                     // SMP
    static Ten3_t dndSig, dnudSig;           // AMP
    static Ten2_t dnudn;                     // AMP
    static Ten2_t tSig, I;                   // tensor Sig and Identity
    static Ten1_t tN, tNu, tn, tnu;          // tensors: normal vectors
    static Ten3_t dtdSig;                    // traction
    static Ten1_t tt, tp, tq;                // tensor: traction and projections
    static Ten2_t tP;                        // tensor: projector
    static Ten2_t dsdSig;                    // auxiliary tensor for dpdSig
    static Ten3_t dpdSig, dqdSig;            // projections
    static Ten2_t tdspdSig, tdsqdSig;        // invariants (tensor)
    static Vec_t  dspdSig, dsqdSig;          // invariants (Mandel's basis)
#endif

private:
    bool  _allocated;
};

// Data
double AnisoInvs::sp;       double AnisoInvs::sq;
Vec3_t AnisoInvs::N123;     Vec3_t AnisoInvs::N;       Vec3_t AnisoInvs::Nu;
Vec3_t AnisoInvs::n;        Vec3_t AnisoInvs::nu;
Mat3_t AnisoInvs::Q;
Mat3_t AnisoInvs::mSig;
Vec3_t AnisoInvs::L;
Vec3_t AnisoInvs::v0;       Vec3_t AnisoInvs::v1;      Vec3_t AnisoInvs::v2;
Vec3_t AnisoInvs::t;        Vec3_t AnisoInvs::p;       Vec3_t AnisoInvs::q;
Mat3_t AnisoInvs::P;

// Derivs
#ifdef HAS_TENSORS
Ten1_t AnisoInvs::V0;        Ten1_t AnisoInvs::V1;        Ten1_t AnisoInvs::V2;
Ten3_t AnisoInvs::dv0dSig;   Ten3_t AnisoInvs::dv1dSig;   Ten3_t AnisoInvs::dv2dSig;
Ten2_t AnisoInvs::E0;        Ten2_t AnisoInvs::E1;        Ten2_t AnisoInvs::E2;
Ten2_t AnisoInvs::dN123dL;
Ten3_t AnisoInvs::dNdSig;    Ten3_t AnisoInvs::dNudSig;
Ten2_t AnisoInvs::dNudN;
Ten3_t AnisoInvs::dndSig;    Ten3_t AnisoInvs::dnudSig;
Ten2_t AnisoInvs::dnudn;
Ten2_t AnisoInvs::tSig;      Ten2_t AnisoInvs::I;
Ten1_t AnisoInvs::tN;        Ten1_t AnisoInvs::tNu;       Ten1_t AnisoInvs::tn;    Ten1_t AnisoInvs::tnu;
Ten3_t AnisoInvs::dtdSig;
Ten1_t AnisoInvs::tt;        Ten1_t AnisoInvs::tp;        Ten1_t AnisoInvs::tq;
Ten2_t AnisoInvs::tP;
Ten2_t AnisoInvs::dsdSig;
Ten3_t AnisoInvs::dpdSig;    Ten3_t AnisoInvs::dqdSig;
Ten2_t AnisoInvs::tdspdSig;  Ten2_t AnisoInvs::tdsqdSig;
Vec_t  AnisoInvs::dspdSig;   Vec_t  AnisoInvs::dsqdSig;
#endif


/////////////////////////////////////////////////////////////////////////////////////////// Implementation


inline AnisoInvs::AnisoInvs (double Theb, double TheAlpha, Vec3_t const & Thea, bool TheObliq)
    : b(Theb), Alpha(TheAlpha), a(Thea), Obliq(TheObliq), Tol(1.0e-8), _allocated(false)
{
    if (!_allocated)
    {
        I.SetDiagonal (1.0);
        _allocated = true;
    }

    // check
    double norm_a = Norm(a);
    if (norm_a<Tol) throw new Fatal("AnisoInvs::AnisoInvs 'a' vector must be non zero. a=[%g,%g,%g]. norm(a)=%g",a(0),a(1),a(2),norm_a);
    au = a / norm_a;
}

inline void AnisoInvs::Calc (Vec_t const & Sig, bool WithDerivs)
{
    // principal values and eigenprojectors
    Mat3_t  mSig;
    Ten2Mat (Sig, mSig);
    Eig     (mSig, L, v0, v1, v2);
    //if (L(0)>0.0 || L(1)>0.0 || L(2)>0.0) throw new Fatal("AnisoInvs::Calc: This method only works when all principal values are negative (compression octant). L=[%g,%g,%g]",L(0),L(1),L(2));

    // rotation matrix
    Q(0,0) = v0(0);   Q(0,1) = v1(0);   Q(0,2) = v2(0);
    Q(1,0) = v0(1);   Q(1,1) = v1(1);   Q(1,2) = v2(1);
    Q(2,0) = v0(2);   Q(2,1) = v1(2);   Q(2,2) = v2(2);

    // normal vectors
    bool zero = (fabs(L(0))<Tol || fabs(L(1))<Tol || fabs(L(2))<Tol);
    if (zero) N123 = 1.0, 1.0, 1.0;                                                       // SMP(123)
    else      N123 = 1.0/pow(fabs(L(0)),b), 1.0/pow(fabs(L(1)),b), 1.0/pow(fabs(L(2)),b); // SMP(123)
    N  = product (Q, N123); // SMP(xyz)
    Nu = N / Norm(N);       // SMP
    n  = Alpha*au + Nu;     // AMP
    nu = n / Norm(n);       // AMP

    // traction
    t = product (mSig, nu);

    // P projector
    double s = 0.0;
    if (Obliq)
    {
        s = 1.0 / dot(N, n);
        Dyad (s, N, n, P); // P = s*(N dy n)
    }
    else Dyad (nu, nu, P);

    // invariants
    p  = product (P, t);
    q  = t - p;
    sp = Norm(p);
    sq = Norm(q);

    // derivatives
    if (WithDerivs)
    {
#ifdef HAS_TENSORS
        // normal to SMP in principal system
        if (zero) dN123dL.SetDiagonal (1.0);
        else dN123dL =            -b/(L(0)*pow(fabs(L(0)),b)), 0.0, 0.0,
                        0.0,      -b/(L(1)*pow(fabs(L(1)),b)), 0.0,
                        0.0, 0.0, -b/(L(2)*pow(fabs(L(2)),b));

        // eigenvectors and eigenprojectors
        Vec2Tensor     (v0,V0);
        Vec2Tensor     (v1,V1);
        Vec2Tensor     (v2,V2);
        EigenVecDerivs (L, V0,V1,V2, dv0dSig,dv1dSig,dv2dSig);
        E0 = (V0 & V0);
        E1 = (V1 & V1);
        E2 = (V2 & V2);

        // normal to SMP in xyz system
        dNdSig = ((V0&E0)*dN123dL[0][0]) + ((V1&E1)*dN123dL[1][1]) + ((V2&E2)*dN123dL[2][2]) + (dv0dSig*N123(0)) + (dv1dSig*N123(1)) + (dv2dSig*N123(2));

        // unit normal to SMP in xyz system
        UnitVecDeriv (N, Nu, dNudN, Tol);
        dNudSig = dNudN * dNdSig;

        // normal to AMP in xyz system
        dndSig = dNudSig;

        // unit normal to AMP in xyz system
        UnitVecDeriv (n, nu, dnudn, Tol);
        dnudSig = dnudn * dndSig;

        // convert to tensors
        Mat2Tensor (P,    tP);
        Mat2Tensor (mSig, tSig);
        Vec2Tensor (t,    tt);
        Vec2Tensor (N,    tN);
        Vec2Tensor (Nu,   tNu);
        Vec2Tensor (n,    tn);
        Vec2Tensor (nu,   tnu);
        Vec2Tensor (p,    tp);
        Vec2Tensor (q,    tq);

        // traction
        dtdSig = (I & tnu) + (tSig * dnudSig);

        // normal projection
        if (Obliq)
        {
            double m = tn * tt;
            dsdSig = (tn*dNdSig + tN*dndSig)*(-s*s);
            dpdSig = (tN&dsdSig)*m + dNdSig*(s*m) + (tN&(tt*dndSig))*s + tP*dtdSig;
        }
        else
        {
            double m = tnu * tt;
            dpdSig = dnudSig*m + (tnu&(tt*dnudSig)) + tP*dtdSig;
        }

        // on-plane projection
        dqdSig = dtdSig - dqdSig;

        // invariants
        if (sp>Tol) tdspdSig = (tp*(1.0/sp)) * dpdSig;   else tdspdSig.SetDiagonal (1.0);
        if (sq>Tol) tdsqdSig = (tq*(1.0/sq)) * dqdSig;   else tdsqdSig.SetDiagonal (1.0);
        size_t ncp = size(Sig);
        //Tensor2Ten (tdspdSig, dspdSig, ncp);
        //Tensor2Ten (tdsqdSig, dsqdSig, ncp);

#else
        throw new Fatal("AnisoInvs::Calc: Tensors library is required in order to calculate derivatives");
#endif
    }
}

#endif // MECHSYS_ANISOINVS_H
