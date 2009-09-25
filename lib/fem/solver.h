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
    enum Scheme_t  { FE_t, ME_t };     ///< Steady time integration scheme
    enum TScheme_t { SS11_t, SS22_t }; ///< Transient time integration scheme

    // Constructor
    Solver (Domain const & Dom);

    // Methods
    void Solve         (size_t NDiv=1);
    void TransSolve    (double tf, double dt, double dtOut);
    void AssembleK     ();
    void AssembleKandM (double dt);
    void TgIncs        (double dT, Vec_t & dU, Vec_t & dF);

    // Data
    Domain const & Dom;   ///< Domain
    double         Time;  ///< Current time (t)
    size_t         NEQ;   ///< Total number of equations (DOFs)
    Array<size_t>  pDOFs; ///< prescribed DOFs (known equations)

    // Triplets and sparse matrices
    Sparse::Triplet<double,int> K11,K12,K21,K22; ///< Stiffness matrices
    Sparse::Triplet<double,int> M11,M12,M21,M22; ///< Mass matrices
    Sparse::Triplet<double,int> A11;             ///< M11 + Theta*dt*K11
    Sparse::Matrix <double,int> K11mat;          ///< Augmented K11 in compressed-column format
    Sparse::Matrix <double,int> A11mat;          ///< Augmented A11 in compressed-column format

    // Vectors
    Vec_t U, F, F_int, W; // U, F, F_int, and Workspace
    Vec_t U0, F0, V;      // Transient: U0, F0, and V

    // Constants for integration
    Scheme_t  Scheme;  ///< Scheme: FE_t (Forward-Euler), ME_t (Modified-Euler)
    size_t    nSS;     ///< FE: number of substeps
    double    STOL;    ///< ME:
    double    dTini;   ///< ME:
    double    mMin;    ///< ME:
    double    mMax;    ///< ME:
    size_t    maxSS;   ///< ME:
    bool      CteTg;   ///< Constant tangent matrices (linear problems) => K and M will be calculated once
    TScheme_t TScheme; ///< Transient scheme
    double    Theta;   ///< Transient scheme constant

private:
    void _initialize   (bool Transient=false); ///< Initialize global matrices and vectors
    void _FE_update    (double tf);            ///< (Forward-Euler) Update Time and elements to tf
    void _ME_update    (double tf);            ///< (Modified-Euler) Update Time and elements to tf
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Solver::Solver (Domain const & TheDom)
    : Dom     (TheDom),
      Time    (0.0),
      Scheme  (ME_t),
      nSS     (1),
      STOL    (1.0e-5),
      dTini   (1.0),
      mMin    (0.1),
      mMax    (10.0),
      maxSS   (2000),
      CteTg   (false),
      TScheme (SS11_t),
      Theta   (2./3.)
{
}

inline void Solver::Solve (size_t NDiv)
{
    // initialize global matrices and vectors
    _initialize ();

    // residual
    Vec_t R(F-F_int);
    double norm_R = Norm(R);
    double tol_R  = 1.0e-9;

    // output initial state
    std::cout << "\n[1;37m--- Stage solution --------- (" << (Scheme==FE_t?"FE":"ME") << ") --------------------------------------------\n";
    std::cout << Util::_6_3 << "Time" <<                                         Util::_8s <<"Norm(R)" << "[0m\n";
    std::cout << Util::_6_3 <<  Time  << (norm_R>tol_R?"[1;31m":"[1;32m") << Util::_8s << norm_R   << "[0m\n";
    Dom.OutResults (Time);

    // solve
    double t0   = Time;     // current time
    double tf   = t0 + 1.0; // final time
    double Dt   = tf - t0;  // total time increment
    double dt   = Dt/NDiv;  // global timestep
    double tout = t0 + dt;  // time for output
    for (size_t inc=0; inc<NDiv; ++inc)
    {
        // update Time and elements to tout
             if (Scheme==FE_t) _FE_update (tout);
        else if (Scheme==ME_t) _ME_update (tout);
        else throw new Fatal("Solver::Solve: Time integration scheme invalid");

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

inline void Solver::TransSolve (double tf, double dt, double dtOut)
{
    // initialize global matrices and vectors
    _initialize (/*Transient*/true);

    // residual
    Vec_t R(F-F_int);
    double norm_R = Norm(R);
    double tol_R  = 1.0e-9;

    // output initial state
    std::cout << "\n[1;37m--- Stage solution --------- (" << (Scheme==FE_t?"FE":"ME") << ") --------------------------------------------\n";
    std::cout << Util::_6_3 << "Time" <<                                         Util::_8s <<"Norm(R)" << "[0m\n";
    std::cout << Util::_6_3 <<  Time  << (norm_R>tol_R?"[1;31m":"[1;32m") << Util::_8s << norm_R   << "[0m\n";
    Dom.OutResults (Time);

    // auxiliar variables
    Vec_t F1(NEQ), U2(NEQ), G10(NEQ), G1(NEQ), dU(NEQ);
    set_to_zero (F1);
    set_to_zero (U2);

    // build U2 and F1
    for (size_t i=0; i<Dom.Nods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.Nods[i]->nDOF(); ++j)
        {
            long eq = Dom.Nods[i]->EQ[j];
            if (Dom.Nods[i]->pU[j]) U2(eq) = U0(eq) + Dom.Nods[i]->DU[j];
            else                    F1(eq) = F0(eq) + Dom.Nods[i]->DF[j];
        }
    }

    // initialize U
    //for (size_t i=0; i<pDOFs.Size(); ++i) U(pDOFs[i]) = U2(pDOFs[i]);

    // initialize G10
    set_to_zero (G10);

    // solve
    double t0   = Time;       // current time
    double tout = t0 + dtOut; // time for output
    while (Time<tf)
    {
        // update elements
        if (TScheme==SS11_t)
        {
            // assemble K, M and A
            AssembleKandM (dt);

            // calc G1 = F1 - K12*U2
            G1 = F1;
            for (int k=0; k<K12.Top(); ++k) G1(K12.Ai(k)) -= K12.Ax(k) * U2(K12.Aj(k));

            // calc W = G1bar - K11*U
            W = G10 + Theta*(G1 - G10);
            for (int k=0; k<K11.Top(); ++k) W(K11.Ai(k)) -= K11.Ax(k) * U(K11.Aj(k));

            //Vec_t tmp(G1-G10);
            //std::cout << "G1-G10  = \n" << PrintVector(tmp);

            // calc V and U
            UMFPACK::Solve (A11mat, W, V); // inv(A11mat)*W = V
            for (size_t i=0; i<pDOFs.Size(); ++i) V(pDOFs[i]) = 0.0;
            dU = dt*V;
            U += dU;
            for (size_t i=0; i<pDOFs.Size(); ++i) U(pDOFs[i]) = U2(pDOFs[i]);
            for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->UpdateState (dU, &F_int);

            // calc F
            F = F1;
            for (int k=0; k<M21.Top(); ++k) F(M21.Ai(k)) += M21.Ax(k) * V(M21.Aj(k)); // F2 += M21 * V1
            for (int k=0; k<K21.Top(); ++k) F(K21.Ai(k)) += K21.Ax(k) * U(K21.Aj(k)); // F2 += K21 * U1
            for (int k=0; k<K22.Top(); ++k) F(K22.Ai(k)) += K22.Ax(k) * U(K22.Aj(k)); // F2 += K22 * U2
            //std::cout << "F     = \n" << PrintVector(F);
            
            // update G10
            G10 = G1;
        }
        else throw new Fatal("Solver::TransSolve: Time integration scheme invalid");

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

        // output
        if (Time>=tout)
        {
            // residual
            R      = F - F_int;
            norm_R = Norm(R);

            // output
            std::cout << Util::_6_3 << Time << (norm_R>tol_R?"[1;31m":"[1;32m") << Util::_8s << norm_R << "[0m\n";
            Dom.OutResults (Time);
            tout += dtOut;
        }

        // next time
        Time += dt;
    }

    // last output
    std::cout << Util::_6_3 << Time << (norm_R>tol_R?"[1;31m":"[1;32m") << Util::_8s << norm_R << "[0m\n";
    Dom.OutResults (Time);
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
    // augment K11 and set K11mat
    for (size_t i=0; i<pDOFs.Size(); ++i) K11.PushEntry (pDOFs[i],pDOFs[i], 1.0);
    K11mat.Set (K11);
    //Mat_t D11;  K11mat.GetDense (D11);
    //std::cout << "K11 =\n" << PrintMatrix(D11);
}

inline void Solver::AssembleKandM (double dt)
{
    K11.ResetTop(); // reset top (position to insert new values) => clear triplet
    K12.ResetTop();
    K21.ResetTop();
    K22.ResetTop();
    M11.ResetTop(); // reset top (position to insert new values) => clear triplet
    M12.ResetTop();
    M21.ResetTop();
    M22.ResetTop();
    A11.ResetTop(); // reset top (position to insert new values) => clear triplet
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
                     if (!pre[i] && !pre[j]) { K11.PushEntry (loc[i], loc[j], K(i,j));  M11.PushEntry (loc[i], loc[j], M(i,j));  A11.PushEntry (loc[i], loc[j], M(i,j)+Theta*dt*K(i,j)); }
                else if (!pre[i] &&  pre[j]) { K12.PushEntry (loc[i], loc[j], K(i,j));  M12.PushEntry (loc[i], loc[j], M(i,j)); }
                else if ( pre[i] && !pre[j]) { K21.PushEntry (loc[i], loc[j], K(i,j));  M21.PushEntry (loc[i], loc[j], M(i,j)); }
                else if ( pre[i] &&  pre[j]) { K22.PushEntry (loc[i], loc[j], K(i,j));  M22.PushEntry (loc[i], loc[j], M(i,j)); }
            }
        }
    }
    // augment A11 and set A11mat
    for (size_t i=0; i<pDOFs.Size(); ++i) A11.PushEntry (pDOFs[i],pDOFs[i], 1.0);
    A11mat.Set (A11);
}

inline void Solver::TgIncs (double dT, Vec_t & dU, Vec_t & dF)
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
                W (eq) = dT*Dom.Nods[i]->DU[j]; // set W2 equal to dU2
            }
            else
            {
                dF(eq) = dT*Dom.Nods[i]->DF[j]; // set dF1 equal to dF1
                W (eq) = dF(eq);                // set W1  equal to dF1
            }
        }
    }

    // W1 = dF1 - K12*dU2
    for (int k=0; k<K12.Top(); ++k) W(K12.Ai(k)) -= K12.Ax(k) * W(K12.Aj(k)); // W1 -= K12 * dU2

    // solve for dU1 (and dU2)
    UMFPACK::Solve (K11mat, W, dU); // inv(K11mat)*W = dU

    // calculate dF2
    for (int k=0; k<K21.Top(); ++k) dF(K21.Ai(k)) += K21.Ax(k) * dU(K21.Aj(k)); // dF2 += K21 * dU1
    for (int k=0; k<K22.Top(); ++k) dF(K22.Ai(k)) += K22.Ax(k) * dU(K22.Aj(k)); // dF2 += K22 * dU2
}

inline void Solver::_initialize (bool Transient)
{
    // assign equation numbers and set pDOFs
    NEQ = 0;
    pDOFs.Resize (0); // prescribed DOFs
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
        M11.AllocSpace (NEQ,NEQ,K11_size);
        M12.AllocSpace (NEQ,NEQ,K12_size);
        M21.AllocSpace (NEQ,NEQ,K21_size);
        M22.AllocSpace (NEQ,NEQ,K22_size);
        A11.AllocSpace (NEQ,NEQ,K11_size+pDOFs.Size()); // augmented
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
        }
    }

    // build V, U0 and F0
    if (Transient)
    {
        V .change_dim (NEQ);
        U0.change_dim (NEQ);
        F0.change_dim (NEQ);
        set_to_zero   (V);
        U0 = U;
        F0 = F;
    }

    // build F_int
    F_int.change_dim (NEQ);
    for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->InitFint (F_int);

    //std::cout << "F     = \n" << PrintVector(F);
    //std::cout << "F_int = \n" << PrintVector(F_int);
}

inline void Solver::_FE_update (double tf)
{
    // auxiliar vectors
    Vec_t dU_fe(NEQ), dF_fe(NEQ);

    double dt = (tf-Time)/nSS;
    for (size_t i=0; i<nSS; ++i)
    {
        // calculate tangent increments
        TgIncs (dt, dU_fe, dF_fe);

        // update elements
        for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->UpdateState (dU_fe, &F_int);

        // update U, F, and Time
        U    += dU_fe;
        F    += dF_fe;
        Time += dt;
    }
}

inline void Solver::_ME_update (double tf)
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
        TgIncs (dt, dU_fe, dF_fe);
        for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->UpdateState (dU_fe);

        // ME state
        TgIncs (dt, dU_tm, dF_tm);
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
