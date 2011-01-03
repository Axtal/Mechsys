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
     AnisoInvs (double b, double Alpha, Vec_t const & a, bool Obliq=true);
    ~AnisoInvs () {}

    // Methods
    void Calc (Vec_t const & Sig, bool WithDerivs=false);

    // Constants
    double b, Alpha;
    Vec_t  a;
    bool   Obliq;
    double Tol;
    Vec_t  au;

    // Data
    static double sp, sq;
    static Vec_t  N, Nu, n, nu;
    static Mat_t  Qt;
    static Vec3_t L;
    static Vec3_t v0,v1,v2;
    static Vec_t  P0,P1,P2;
    static Vec_t  t, p, q;
    static Mat_t  P;

    // Derivs
    static Mat_t  dNudN,   dNds,  dNuds; // SMP
    static Mat_t  dnudn,   dnds,  dnuds; // AMP
    static Mat_t  dtds,    dpds,  dqds;
    static double dspds,   dsqds;
    static Vec_t  dspdSig, dsqdSig;

private:
    bool  _allocated;
};

double AnisoInvs::sp;  double AnisoInvs::sq;
Vec_t  AnisoInvs::N;   Vec_t  AnisoInvs::Nu;  Vec_t AnisoInvs::n;  Vec_t AnisoInvs::nu;
Mat_t  AnisoInvs::Qt;
Vec3_t AnisoInvs::L;
Vec3_t AnisoInvs::v0;  Vec3_t AnisoInvs::v1;  Vec3_t AnisoInvs::v2;
Vec_t  AnisoInvs::P0;  Vec_t  AnisoInvs::P1;  Vec_t  AnisoInvs::P2;
Vec_t  AnisoInvs::t;   Vec_t  AnisoInvs::p;   Vec_t  AnisoInvs::q;
Mat_t  AnisoInvs::P;

Mat_t  AnisoInvs::dNudN;    Mat_t  AnisoInvs::dNds;    Mat_t AnisoInvs::dNuds;
Mat_t  AnisoInvs::dnudn;    Mat_t  AnisoInvs::dnds;    Mat_t AnisoInvs::dnuds;
Mat_t  AnisoInvs::dtds;     Mat_t  AnisoInvs::dpds;    Mat_t AnisoInvs::dqds;
double AnisoInvs::dspds;    double AnisoInvs::dsqds;
Vec_t  AnisoInvs::dspdSig;  Vec_t  AnisoInvs::dsqdSig;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation


inline AnisoInvs::AnisoInvs (double Theb, double TheAlpha, Vec_t const & Thea, bool TheObliq)
    : b(Theb), Alpha(TheAlpha), a(Thea), Obliq(TheObliq), Tol(1.0e-8), _allocated(false)
{
    if (!_allocated)
    {
        N .change_dim (3);
        Nu.change_dim (3);
        n .change_dim (3);
        nu.change_dim (3);
        Qt.change_dim (3,3);
        t .change_dim (3);
        p .change_dim (3);
        q .change_dim (3);
        P .change_dim (3,3);

        dNudN.change_dim (3,3);
        dNds .change_dim (3,3);
        dNuds.change_dim (3,3);

        dnudn.change_dim (3,3);
        dnds .change_dim (3,3);
        dnuds.change_dim (3,3);

        dtds.change_dim (3,3);
        dpds.change_dim (3,3);
        dqds.change_dim (3,3);
    }

    // check
    double norm_a = Norm(a);
    if (norm_a<Tol) throw new Fatal("AnisoInvs::AnisoInvs 'a' vector must be non zero. a=[%g,%g,%g]. norm(a)=%g",a(0),a(1),a(2),norm_a);
    au = a / norm_a;
}

inline void AnisoInvs::Calc (Vec_t const & Sig, bool WithDerivs)
{
    // principal values and eigenprojectors
    EigenProj (Sig, L, v0, v1, v2, P0, P1, P2);
    for (size_t j=0; j<3; ++j)
    {
        Qt(0,j) = v0(j);
        Qt(1,j) = v1(j);
        Qt(2,j) = v2(j);
    }
    //if (L(0)>0.0 || L(1)>0.0 || L(2)>0.0) throw new Fatal("AnisoInvs::Calc: This method only works when all principal values are negative (compression octant). L=[%g,%g,%g]",L(0),L(1),L(2));

    // normal vectors
    bool zero = (fabs(L(0))<Tol || fabs(L(1))<Tol || fabs(L(2))<Tol);
    if (zero) N = 1.0, 1.0, 1.0;                                                       // SMP
    else      N = 1.0/pow(fabs(L(0)),b), 1.0/pow(fabs(L(1)),b), 1.0/pow(fabs(L(2)),b); // SMP
    Nu = N / Norm(N);      // SMP
    n  = Alpha*Qt*au + Nu; // AMP
    //n  = Alpha*au + Nu; // AMP
    nu = n / Norm(n);      // AMP

    // traction
    t = L(0)*nu(0), L(1)*nu(1), L(2)*nu(2);

    // P projector
    double s = 0.0;
    if (Obliq)
    {
        s = 1.0 / dot(N, n);
        Dyad (N, n, P);
        P *= s;
    }
    else Dyad (nu, nu, P);

    // invariants
    p  = P * t;
    q  = t - p;
    sp = Norm(p);
    sq = Norm(q);

    // derivatives
    if (WithDerivs)
    {
        // normal vectors
        UnitVecDeriv (N, Nu, dNudN, Tol);
        if (zero) Identity (3, dNds);
        else dNds = -b/(L(0)*pow(fabs(L(0)),b)), 0.0, 0.0,
                    0.0, -b/(L(1)*pow(fabs(L(1)),b)), 0.0,
                    0.0, 0.0, -b/(L(2)*pow(fabs(L(2)),b)); // SMP
        dNuds = dNudN * dNds;             // SMP
        dnds  = dNuds;                    // AMP
        UnitVecDeriv (n, nu, dnudn, Tol); // AMP
        dnuds = dnudn * dnds;             // AMP
    }
}

#endif // MECHSYS_ANISOINVS_H
