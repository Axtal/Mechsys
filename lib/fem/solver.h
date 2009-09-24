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

#ifndef MECHSYS_FEM_SOLVER_H
#define MECHSYS_FEM_SOLVER_H

// Std Lib
#include <cstring>  // for strcmp
#include <ctime>    // for clock
#include <iostream> // for cout

// Blitz++
#include <blitz/tinyvec-et.h>

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "fem/domain.h"
#include "linalg/sparse_triplet.h"
#include "linalg/sparse_matrix.h"
#include "linalg/umfpack.h"

namespace FEM
{

class Solver
{
public:
    // enum
    enum Scheme_t { FE_t, ME_t }; ///< Integration scheme

    // Constructor
    Solver (Domain const & Dom);

    // Methods
    void Solve         (size_t NDiv=1, double tf=1.0, bool Transient=false);
    void AssembleK     ();
    void AssembleKandM ();
    void TgIncs        (double t, double dt, Vec_t & dU, Vec_t & dF);
    void TransTgIncs   (double t, double dt, Vec_t & dU, Vec_t & dF);

    // Data
    Domain const & Dom;   ///< Domain
    double         Time;  ///< Current time (t)
    size_t         NEQ;   ///< Total number of equations (DOFs)
    Array<size_t>  pDOFs; ///< prescribed DOFs (known equations)

    // Triplets and sparse matrices
    Sparse::Triplet<double,int> K11,K12,K21,K22; ///< Stiffness matrices
    Sparse::Triplet<double,int> M11,M12,M21,M22; ///< Mass matrices
    Sparse::Matrix <double,int> k11;             ///< Augmented K11 in compressed-column format
    Sparse::Matrix <double,int> m11;             ///< Augmented M11 in compressed-column format

    // Vectors
    Vec_t U, F, F_int, W; // U, F, F_int, and Workspace

    // Constants for integration
    Scheme_t Scheme; ///< Scheme: FE_t (Forward-Euler), ME_t (Modified-Euler)
    size_t   nSS;    ///< FE: number of substeps
    double   STOL;   ///< ME:
    double   dTini;  ///< ME:
    double   mMin;   ///< ME:
    double   mMax;   ///< ME:
    size_t   maxSS;  ///< ME:
    bool     CteTg;  ///< Constant tangent matrices (linear problems) => K and M will be calculated once

    double   tsw; ///< t_switch

private:
    typedef void (Solver::*p2TgIncs) (double t, double dt, Vec_t & dU, Vec_t & dF);
    void _initialize (bool Transient=false);       ///< Initialize global matrices and vectors
    void _FE_update  (p2TgIncs TgIncs, double tf);
    void _ME_update  (p2TgIncs TgIncs, double tf);
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Solver::Solver (Domain const & TheDom)
    : Dom    (TheDom),
      Time   (0.0),
      Scheme (ME_t),
      nSS    (1),
      STOL   (1.0e-5),
      dTini  (1.0),
      mMin   (0.1),
      mMax   (10.0),
      maxSS  (2000),
      CteTg  (false),
      tsw    (0.1)
{
}

inline void Solver::Solve (size_t NDiv, double tf, bool Transient)
{
    // initialize global matrices and vectors
    _initialize (Transient);

    // residual
    Vec_t R(F-F_int);
    double norm_R = Norm(R);
    double tol_R  = 1.0e-9;

    // output initial state
    std::cout << "\n[1;37m--- Stage solution --------- (" << (Scheme==FE_t?"FE":"ME") << ") --------------------------------------------\n";
    std::cout << Util::_6_3 << "Time" <<                                         Util::_8s <<"Norm(R)" << "[0m\n";
    std::cout << Util::_6_3 <<  Time  << (norm_R>tol_R?"[1;31m":"[1;32m") << Util::_8s << norm_R   << "[0m\n";
    Dom.OutResults (Time);

    // pointer to function that calculates tangent increments
    p2TgIncs pfun = NULL;
    if (Transient) pfun = &Solver::TransTgIncs; // transient analysis
    else           pfun = &Solver::TgIncs;      // steady (default) analysis

    // solve
    if (!Transient) tf = Time + 1.0;
    double Dt   = tf - Time;
    double dt   = Dt/NDiv;
    double tout = Time + dt;
    for (size_t inc=0; inc<NDiv; ++inc)
    {
        // update Time and elements to tout
             if (Scheme==FE_t) _FE_update (pfun, tout);
        else if (Scheme==ME_t) _ME_update (pfun, tout);

        // update nodes to tout
        for (size_t i=0; i<Dom.Nods.Size(); ++i)
        {
            for (size_t j=0; j<Dom.Nods[i]->nDOF(); ++j)
            {
                long eq = Dom.Nods[i]->EQ[j];
                Dom.Nods[i]->U[j] = U(eq);
                Dom.Nods[i]->F[j] = F(eq);
            }
        }

        // residual
        R      = F - F_int;
        norm_R = Norm(R);

        // output
        std::cout << Util::_6_3 << Time << (norm_R>tol_R?"[1;31m":"[1;32m") << Util::_8s << norm_R << "[0m\n";
        Dom.OutResults (Time);

        // next tout
        tout = Time + dt;
    }
}

inline void Solver::AssembleK ()
{
    K11.ResetTop(); // reset top (position to insert new values) => clear triplet
    K12.ResetTop();
    K21.ResetTop();
    K22.ResetTop();
    for (size_t k=0; k<Dom.Eles.Size(); ++k)
    {
        Mat_t         K;   // K matrix
        Array<size_t> loc; // location array
        Array<bool>   pre; // prescribed U ?
        Dom.Eles[k]->CalcK  (K);
        Dom.Eles[k]->GetLoc (loc, pre);
        for (size_t i=0; i<loc.Size(); ++i)
        {
            for (size_t j=0; j<loc.Size(); ++j)
            {
                     if (!pre[i] && !pre[j]) K11.PushEntry (loc[i], loc[j], K(i,j));
                else if (!pre[i] &&  pre[j]) K12.PushEntry (loc[i], loc[j], K(i,j));
                else if ( pre[i] && !pre[j]) K21.PushEntry (loc[i], loc[j], K(i,j));
                else if ( pre[i] &&  pre[j]) K22.PushEntry (loc[i], loc[j], K(i,j));
            }
        }
    }
    // augment K11 and set k11
    for (size_t i=0; i<pDOFs.Size(); ++i) K11.PushEntry (pDOFs[i],pDOFs[i], 1.0);
    k11.Set (K11);
    //Mat_t D11;  k11.GetDense (D11);
    //std::cout << "K11 =\n" << PrintMatrix(D11);
}

inline void Solver::AssembleKandM ()
{
    K11.ResetTop(); // reset top (position to insert new values) => clear triplet
    K12.ResetTop();
    K21.ResetTop();
    K22.ResetTop();
    M11.ResetTop();
    M12.ResetTop();
    M21.ResetTop();
    M22.ResetTop();
    for (size_t k=0; k<Dom.Eles.Size(); ++k)
    {
        Mat_t         K, M; // matrices
        Array<size_t> loc;  // location array
        Array<bool>   pre;  // prescribed U ?
        Dom.Eles[k]->CalcK  (K);
        Dom.Eles[k]->CalcM  (M);
        Dom.Eles[k]->GetLoc (loc, pre);
        for (size_t i=0; i<loc.Size(); ++i)
        {
            for (size_t j=0; j<loc.Size(); ++j)
            {
                     if (!pre[i] && !pre[j]) { K11.PushEntry (loc[i], loc[j], K(i,j));  M11.PushEntry (loc[i], loc[j], M(i,j)); }
                else if (!pre[i] &&  pre[j]) { K12.PushEntry (loc[i], loc[j], K(i,j));  M12.PushEntry (loc[i], loc[j], M(i,j)); }
                else if ( pre[i] && !pre[j]) { K21.PushEntry (loc[i], loc[j], K(i,j));  M21.PushEntry (loc[i], loc[j], M(i,j)); }
                else if ( pre[i] &&  pre[j]) { K22.PushEntry (loc[i], loc[j], K(i,j));  M22.PushEntry (loc[i], loc[j], M(i,j)); }
            }
        }
    }
    // augment M11 and set m11
    for (size_t i=0; i<pDOFs.Size(); ++i) M11.PushEntry (pDOFs[i],pDOFs[i], 1.0);
    m11.Set (M11);
    //Mat_t D11;  m11.GetDense (D11);
    //std::cout << "M11 =\n" << PrintMatrix(D11);
}

inline void Solver::TgIncs (double t, double dt, Vec_t & dU, Vec_t & dF)
{
    // assemble global K matrix
    if (K11.Top()==0 || CteTg==false) AssembleK (); // not constant tangent matrices => non-linear problems

    // assemble dF and W (workspace) vectors
    for (size_t i=0; i<Dom.Nods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.Nods[i]->nDOF(); ++j)
        {
            long eq = Dom.Nods[i]->EQ[j];
            if (Dom.Nods[i]->pU[j]) // prescribed U
            {
                dF(eq) = 0.0;                   // clear dF2
                W (eq) = dt*Dom.Nods[i]->dU[j]; // set W2 equal to dU2
            }
            else
            {
                dF(eq) = dt*Dom.Nods[i]->dF[j]; // set dF1 equal to dF1
                W (eq) = dF(eq);                // set W1  equal to dF1
            }
        }
    }

    // W1 = dF1 - K12*dU2
    for (int k=0; k<K12.Top(); ++k) W(K12.Ai(k)) -= K12.Ax(k) * W(K12.Aj(k)); // W1 -= K12 * dU2

    // solve for dU1 (and dU2)
    UMFPACK::Solve (k11, W, dU); // inv(k11)*W = dU

    // calculate dF2
    for (int k=0; k<K21.Top(); ++k) dF(K21.Ai(k)) += K21.Ax(k) * dU(K21.Aj(k)); // dF2 += K21 * dU1
    for (int k=0; k<K22.Top(); ++k) dF(K22.Ai(k)) += K22.Ax(k) * dU(K22.Aj(k)); // dF2 += K22 * dU2
}

inline void Solver::TransTgIncs (double t, double dt, Vec_t & dU, Vec_t & dF)
{
    // assemble global K and M matrices
    if (K11.Top()==0 || CteTg==false) AssembleKandM(); // not constant tangent matrices => non-linear problems

    // assemble W (workspace) and U vectors
    for (size_t i=0; i<Dom.Nods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.Nods[i]->nDOF(); ++j)
        {
            long eq = Dom.Nods[i]->EQ[j];
            if (Dom.Nods[i]->pU[j]) // prescribed U
            {
                W(eq) = 0.0; // set W2 equal to V2==0
                U(eq) = Dom.Nods[i]->dU[j];
                //double alp = Dom.Nods[i]->dU[j]/tsw;
                //U(eq) = (t<tsw ? alp*t : Dom.Nods[i]->dU[j]);
                //W(eq) = Dom.Nods[i]->GetPrescV(j,t); // set W2 equal to V2
                //U(eq) = Dom.Nods[i]->GetPrescU(j,t); // set U2 equal to U2
            }
            else
            {
                W(eq) = Dom.Nods[i]->dF[j];
                //W(eq) = F(eq); // set W1 equal to F1
                //double alp = Dom.Nods[i]->dF[j]/tsw;
                //W(eq) = (t<tsw ? alp*t : Dom.Nods[i]->dF[j]);
                //W(eq) = Dom.Nods[i]->GetPrescF(j,t); // set W1 equal to F1
            }
        }
    }

    //cout << "W = \n" << PrintVector(W);

    // W1 = F1 - M12*V2 - K12*U2 - K11*U1
    for (int k=0; k<M12.Top(); ++k) W(M12.Ai(k)) -= M12.Ax(k) * W(M12.Aj(k)); // W1 -= M12 * V2
    for (int k=0; k<K12.Top(); ++k) W(K12.Ai(k)) -= K12.Ax(k) * U(K12.Aj(k)); // W1 -= K12 * U2
    //cout << "W = \n" << PrintVector(W);
    for (int k=0; k<K11.Top(); ++k) W(K11.Ai(k)) -= K11.Ax(k) * U(K11.Aj(k)); // W1 -= K11 * U1

    // solve for V1 (and V2)
    Vec_t V(NEQ);
    UMFPACK::Solve (m11, W, V); // inv(M11)*W = V

    // increments
    dU = V*dt;

    // dF2 = M21*V1 + M22*V2 + K21*Unew1 + K22*Unew2 - F2
    Vec_t Unew(U+dU);
    set_to_zero (dF);
    for (int k=0; k<M21.Top(); ++k) dF(M21.Ai(k)) += M21.Ax(k) * V   (M21.Aj(k)); // dF2 += M21 * V1
    for (int k=0; k<M22.Top(); ++k) dF(M22.Ai(k)) += M22.Ax(k) * V   (M22.Aj(k)); // dF2 += M22 * V2
    for (int k=0; k<K21.Top(); ++k) dF(K21.Ai(k)) += K21.Ax(k) * Unew(K21.Aj(k)); // dF2 += K21 * Unew1
    for (int k=0; k<K22.Top(); ++k) dF(K22.Ai(k)) += K22.Ax(k) * Unew(K22.Aj(k)); // dF2 += K22 * Unew2
    for (size_t i=0; i<pDOFs.Size(); ++i) dF(pDOFs[i]) -= F(pDOFs[i]);            // dF2 -= F2
}

inline void Solver::_initialize (bool Transient)
{
    // assign equation numbers
    NEQ = 0;
    pDOFs.Resize (0);
    for (size_t i=0; i<Dom.Nods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.Nods[i]->nDOF(); ++j)
        {
            Dom.Nods[i]->EQ[j] = NEQ;
            if (Dom.Nods[i]->pU[j]) pDOFs.Push (NEQ);
            NEQ++;
        }
    }

    // find total number of non-zero entries, including duplicates
    size_t K11_size = 0;
    size_t K12_size = 0;
    size_t K21_size = 0;
    size_t K22_size = 0;
    for (size_t k=0; k<Dom.Eles.Size(); ++k)
    {
        Array<size_t> loc; // location array
        Array<bool>   pre; // prescribed U ?
        Dom.Eles[k]->GetLoc (loc, pre);
        for (size_t i=0; i<loc.Size(); ++i)
        for (size_t j=0; j<loc.Size(); ++j)
        {
                 if (!pre[i] && !pre[j]) K11_size++;
            else if (!pre[i] &&  pre[j]) K12_size++;
            else if ( pre[i] && !pre[j]) K21_size++;
            else if ( pre[i] &&  pre[j]) K22_size++;
        }
    }

    // allocate triplets
    K12.AllocSpace (NEQ,NEQ,K12_size);
    K21.AllocSpace (NEQ,NEQ,K21_size);
    K22.AllocSpace (NEQ,NEQ,K22_size);
    if (Transient)
    {
        K11.AllocSpace (NEQ,NEQ,K11_size);
        M11.AllocSpace (NEQ,NEQ,K11_size+pDOFs.Size()); // augmented
        M12.AllocSpace (NEQ,NEQ,K12_size);
        M21.AllocSpace (NEQ,NEQ,K21_size);
        M22.AllocSpace (NEQ,NEQ,K22_size);
    }
    else
    {
        K11.AllocSpace (NEQ,NEQ,K11_size+pDOFs.Size()); // augmented
    }

    // resize workspace
    W.change_dim (NEQ);

    // build U and F
    U.change_dim (NEQ);
    F.change_dim (NEQ);
    for (size_t i=0; i<Dom.Nods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.Nods[i]->nDOF(); ++j)
        {
            long eq = Dom.Nods[i]->EQ[j];
            U(eq)   = Dom.Nods[i]->U [j];
            F(eq)   = Dom.Nods[i]->F [j];
            /*
            if (Transient)
            {
                if (Dom.Nods[i]->pU[j]) U(eq) = Dom.Nods[i]->dU[j];
                else                    F(eq) = Dom.Nods[i]->dF[j];
            }
            */
        }
    }

    // build F_int
    F_int.change_dim (NEQ);
    for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->InitFint (F_int);

    //std::cout << "F     = \n" << PrintVector(F);
    //std::cout << "F_int = \n" << PrintVector(F_int);
}

inline void Solver::_FE_update (p2TgIncs TgIncs, double tf)
{
    // auxiliar vectors
    Vec_t dU_fe(NEQ), dF_fe(NEQ);

    double dt = (tf-Time)/nSS;
    for (size_t i=0; i<nSS; ++i)
    {
        // calculate tangent increments
        (this->*TgIncs) (Time, dt, dU_fe, dF_fe);

        // update elements
        for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->UpdateState (dU_fe, &F_int);

        // update U, F, and Time
        U    += dU_fe;
        F    += dF_fe;
        Time += dt;
    }
}

inline void Solver::_ME_update (p2TgIncs TgIncs, double tf)
{
    // auxiliar vectors
    Vec_t dU_fe(NEQ), dU_tm(NEQ), dU_me(NEQ), U_me(NEQ), U_dif(NEQ);
    Vec_t dF_fe(NEQ), dF_tm(NEQ), dF_me(NEQ), F_me(NEQ), F_dif(NEQ);

    // for each pseudo time T
    double T  = 0.0;
    double dT = dTini;
    double Dt = tf-Time;
    size_t k  = 0;
    for (k=0; k<maxSS; ++k)
    {
        // exit point
        if (T>=1.0) break;

        // backup state of elements
        for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->BackupState ();

        // time increment
        double dt = Dt*dT;

        // FE state
        (this->*TgIncs) (Time, dt, dU_fe, dF_fe);
        for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->UpdateState (dU_fe);

        // ME state
        (this->*TgIncs) (Time+dt, dt, dU_tm, dF_tm);
        dU_me = 0.5*(dU_fe + dU_tm);
        dF_me = 0.5*(dF_fe + dF_tm);
        U_me  = U + dU_me;
        F_me  = F + dF_me;

        // local error
        U_dif = 0.5*(dU_tm - dU_fe);
        F_dif = 0.5*(dF_tm - dF_fe);
        double U_err = Norm(U_dif)/(1.0+Norm(U_me));
        double F_err = Norm(F_dif)/(1.0+Norm(F_me));
        double error = U_err + F_err;

        // step multiplier
        double m = (error>0.0 ? 0.9*sqrt(STOL/error) : mMax);

        // restore state of elements
        for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->RestoreState ();

        // update
        if (error<STOL)
        {
            for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->UpdateState (dU_me, &F_int);
            T    += dT;
            U     = U_me;
            F     = F_me;
            Time += dt;
            if (true && T<1.0)
            {
                Vec_t R(F-F_int);
                double norm_R = Norm(R);
                double tol_R  = 1.0e-9;
                std::cout << Util::_6_3 << Time << (norm_R>tol_R?"[1;31m":"[1;32m") << Util::_8s << norm_R << "[0m\n";
                Dom.OutResults (Time);
            }
            if (m>mMax) m = mMax;
        }
        else if (m<mMin) m = mMin;

        // change next step size
        dT = m * dT;

        // check for last increment
        if (dT>1.0-T) dT = 1.0-T;
    }
    if (k>=maxSS) throw new Fatal("Solver:_ME_update: Modified-Euler (global integration) did not converge for %d steps",k);
}

}; // namespace FEM

#endif // MECHSYS_FEM_SOLVER_H
