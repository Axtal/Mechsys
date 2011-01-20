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

#ifndef MECHSYS_SMPINVS_H
#define MECHSYS_SMPINVS_H

// Std Lib
#include <iostream>
#include <fstream>

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/linalg/jacobirot.h>

class SMPInvs
{
public:
    SMPInvs () : b(0.5), Tol(1.0e-8) {}

    // Methods
    void Calc (Vec_t const & Sig, bool WithDerivs=false);

    // Variables
    double b;   ///< coefficient for normal to SMP
    double Tol; ///< tol for zero sp and sq

    // Data
    static double  sp, sq;   ///< invariants
    static Vec3_t  N, Nu;    ///< normal to SMP
    static Mat3_t  Q;        ///< rotation matrix
    static Vec3_t  L;        ///< eigenvalues
    static Vec3_t  v0,v1,v2; ///< eigenvectors
    static Vec_t   P0,P1,P2; ///< eigenprojectors
    static Vec3_t  t, p, q;  ///< traction and projections
    static Mat3_t  P;        ///< projector

    // Data for derivatives
    static Vec_t  E0, E1, E2;         ///< eigenprojectors
    static Mat3_t dNdL, dNudN, dNudL; ///< normal to SMP
    static Mat3_t mL, NN, dtdL;       ///< traction
    static Vec3_t vTmp;               ///< projections
    static Mat3_t pTmp;               ///< projections
    static Mat3_t mTmp, dpdL, dqdL;   ///< projections
    static Vec3_t dspdL, dsqdL;       ///< invariants
    static Vec_t  dspdSig, dsqdSig;   ///< invariants w.r.t Sig
};

// Data
double SMPInvs::sp;       double SMPInvs::sq;
Vec3_t SMPInvs::N;        Vec3_t SMPInvs::Nu;
Mat3_t SMPInvs::Q;
Vec3_t SMPInvs::L;
Vec3_t SMPInvs::v0;       Vec3_t SMPInvs::v1;       Vec3_t SMPInvs::v2;
Vec3_t SMPInvs::t;        Vec3_t SMPInvs::p;        Vec3_t SMPInvs::q;
Mat3_t SMPInvs::P;

// Derivs
Vec_t  SMPInvs::E0;       Vec_t  SMPInvs::E1;       Vec_t  SMPInvs::E2;
Mat3_t SMPInvs::dNdL;     Mat3_t SMPInvs::dNudN;    Mat3_t SMPInvs::dNudL;
Mat3_t SMPInvs::mL;       Mat3_t SMPInvs::dtdL;
Vec3_t SMPInvs::vTmp;     Mat3_t SMPInvs::pTmp;
Mat3_t SMPInvs::mTmp;     Mat3_t SMPInvs::dpdL;     Mat3_t SMPInvs::dqdL;
Vec3_t SMPInvs::dspdL;    Vec3_t SMPInvs::dsqdL;
Vec_t  SMPInvs::dspdSig;  Vec_t  SMPInvs::dsqdSig;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation


inline void SMPInvs::Calc (Vec_t const & Sig, bool WithDerivs)
{
    if (WithDerivs) EigenProj (Sig, L, v0, v1, v2, E0, E1, E2);
    else            Eig       (Sig, L, v0, v1, v2);
    //if (L(0)>0.0 || L(1)>0.0 || L(2)>0.0) throw new Fatal("SMPInvs::Calc: This method only works when all principal values are negative (compression octant). L=[%g,%g,%g]",L(0),L(1),L(2));

    // normal vector
    bool zero = (fabs(L(0))<Tol || fabs(L(1))<Tol || fabs(L(2))<Tol);
    if (zero) N = 1.0, 1.0, 1.0;
    else      N = 1.0/pow(fabs(L(0)),b), 1.0/pow(fabs(L(1)),b), 1.0/pow(fabs(L(2)),b);
    Nu = N / Norm(N);

    // traction
    t(0) = L(0)*Nu(0);
    t(1) = L(1)*Nu(1);
    t(2) = L(2)*Nu(2);

    // P projector
    Dyad (Nu, Nu, P);

    // invariants
    p  = product (P, t);
    q  = t - p;
    sp = Norm(p);
    sq = Norm(q);
    if ((L(0)+L(1)+L(2))>0.0) sp *= (-1.0);

    // derivatives
    if (WithDerivs)
    {
        // normal to SMP
        if (zero) Identity (dNdL);
        else dNdL =           -b/(L(0)*pow(fabs(L(0)),b)), 0.0, 0.0,
                    0.0,      -b/(L(1)*pow(fabs(L(1)),b)), 0.0,
                    0.0, 0.0, -b/(L(2)*pow(fabs(L(2)),b));
        UnitVecDeriv (N, Nu, dNudN, Tol);
        dNudL = product (dNudN, dNdL);

        // traction
        mL = L(0),   0.0,   0.0,
              0.0,  L(1),   0.0,
              0.0,   0.0,  L(2);
        dtdL = product (mL, dNudL);
        dtdL(0,0) += Nu(0);
        dtdL(1,1) += Nu(1);
        dtdL(2,2) += Nu(2);

        // projections
        double sc = dot(Nu, t);
        Mult (t, dNudL, vTmp);
        Dyad (Nu, vTmp, mTmp);
        pTmp = product (P, dtdL);
        for (size_t i=0; i<3; ++i)
        for (size_t j=0; j<3; ++j)
        {
            dpdL(i,j) = sc*dNudL(i,j) + mTmp(i,j) + pTmp(i,j);
            dqdL(i,j) = dtdL(i,j) - dpdL(i,j);
        }

        // invariants
        if (sp>Tol) { Mult (p, dpdL, dspdL);  dspdL /= sp; } else dspdL = 0.,0.,0.;
        if (sq>Tol) { Mult (q, dqdL, dsqdL);  dsqdL /= sq; } else dsqdL = 0.,0.,0.;

        // invariants w.r.t Sig
        dspdSig = dspdL(0)*E0 + dspdL(1)*E1 + dspdL(2)*E2;
        dsqdSig = dsqdL(0)*E0 + dsqdL(1)*E1 + dsqdL(2)*E2;
    }
}

#endif // MECHSYS_SMPINVS_H
