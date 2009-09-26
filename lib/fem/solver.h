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
    enum Scheme_t  { FE_t, ME_t };         ///< Steady time integration scheme
    enum TScheme_t { SS11_t };             ///< Transient time integration scheme
    enum DScheme_t { SS22_t };             ///< Dynamic time integration scheme
    enum Damping_t { None_t, Rayleigh_t }; ///< Damping type

    // Constructor
    Solver (Domain const & Dom);

    // Methods
    void Solve        (size_t NDiv=1);
    void TransSolve   (double tf, double dt, double dtOut);
    void DynSolve     (double tf, double dt, double dtOut);
    void AssembleKA   ();                                         ///< A = K11
    void AssembleKMA  (double Coef1, double Coef2);               ///< A = Coef1*M + Coef2*K
    void AssembleKCMA (double Coef1, double Coef2, double Coef3); ///< A = Coef1*M + Coef2*C + Coef3*K
    void TgIncs       (double dT, Vec_t & dU, Vec_t & dF);

    // Data
    Domain const & Dom;   ///< Domain
    double         Time;  ///< Current time (t)
    size_t         NEQ;   ///< Total number of equations (DOFs)
    Array<size_t>  pDOFs; ///< prescribed DOFs (known equations)

    // Triplets and sparse matrices
    Sparse::Triplet<double,int> K11,K12,K21,K22; ///< Stiffness matrices
    Sparse::Triplet<double,int> C11,C12,C21,C22; ///< Damping matrices
    Sparse::Triplet<double,int> M11,M12,M21,M22; ///< Mass matrices
    Sparse::Triplet<double,int> A11;             ///< A=K  or  A=C1*M+C2*K  or  A=C1*M+C2*C+C3*K
    Sparse::Matrix <double,int> A11mat;          ///< Augmented A11 in compressed-column format

    // Vectors
    Vec_t U, F, F_int, W;  // U, F, F_int, and Workspace
    Vec_t V, U2, F1;       // Transient/dynamic variables

    // Constants for integration
    Scheme_t  Scheme;  ///< Scheme: FE_t (Forward-Euler), ME_t (Modified-Euler)
    size_t    nSS;     ///< FE: number of substeps
    double    STOL;    ///< ME:
    double    dTini;   ///< ME:
    double    mMin;    ///< ME:
    double    mMax;    ///< ME:
    size_t    maxSS;   ///< ME:
    bool      CteTg;   ///< Constant tangent matrices (linear problems) => K will be calculated once
    TScheme_t TScheme; ///< Transient scheme
    double    Theta;   ///< Transient scheme constant
    DScheme_t DScheme; ///< Dynamic scheme
    Damping_t DampTy;  ///< Damping type
    double    DampAlp; ///< Rayleigh damping alpha coefficient
    double    DampBet; ///< Rayleigh damping beta coefficient
    double    DynTh1;  ///< Dynamic coefficient Theta 1
    double    DynTh2;  ///< Dynamic coefficient Theta 2

private:
    void _initialize (bool Transient=false); ///< Initialize global matrices and vectors
    void _FE_update  (double tf);            ///< (Forward-Euler) Update Time and elements to tf
    void _ME_update  (double tf);            ///< (Modified-Euler) Update Time and elements to tf
    void _calc_Fbar  (double t, Vec_t & Fb); ///< Calc Fbar = F1 - K12*U2
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
      Theta   (2./3.),
      DScheme (SS22_t),
      DampTy  (None_t),
      DampAlp (0.5),
      DampBet (0.5),
      DynTh1  (0.5),
      DynTh2  (0.5)
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
    Vec_t dU(NEQ);

    // solve
    double t0   = Time;       // current time
    double tout = t0 + dtOut; // time for output
    while (Time<tf)
    {
        // update elements
        if (TScheme==SS11_t)
        {
            AssembleKMA    (1.0, Theta*dt); // A = M + Theta*dt*K
            _calc_Fbar     (Time, W);       // W = Fbar = F1 - K12*U2
            SubMult        (K11, U, W);     // W -= K11*U
            UMFPACK::Solve (A11mat, W, V);  // V = inv(A11mat)*W
            dU = dt*V;
            U += dU;
            for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->UpdateState (dU, &F_int);
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

inline void Solver::DynSolve (double tf, double dt, double dtOut)
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
    Vec_t Ubar(NEQ), A(NEQ), dU(NEQ);

    // solve
    double t0   = Time;       // current time
    double tout = t0 + dtOut; // time for output
    while (Time<tf)
    {
        // update elements
        if (DScheme==SS22_t)
        {
            Ubar = U + DynTh1*dt*V;
            if (DampTy==None_t) // no damping
            {
                AssembleKMA     (1.0, 0.5*DynTh2*dt*dt); // A = M + 0.5*Th2*(dt*dt)*K
                _calc_Fbar      (Time, W);               // W = Fbar = F1 - K12*U2
                Sparse::SubMult (K11, Ubar, W);          // W -= K11*Ubar
                UMFPACK::Solve  (A11mat, W, A);          // A = inv(A11mat)*W
                dU = dt*V + (0.5*dt*dt)*A;
                U += dU;
                V += dt*A;
                for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->UpdateState (dU, &F_int);
            }
            else throw new Fatal("Solver::DynSolve: Damping not implemented yet");
        }
        else throw new Fatal("Solver::DynSolve: Time integration scheme invalid");

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

inline void Solver::AssembleKA ()
{
    A11.ResetTop(); // reset top (position to insert new values) => clear triplet
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
                     if (!pre[i] && !pre[j]) A11.PushEntry (loc[i], loc[j], K(i,j));
                else if (!pre[i] &&  pre[j]) K12.PushEntry (loc[i], loc[j], K(i,j));
                else if ( pre[i] && !pre[j]) K21.PushEntry (loc[i], loc[j], K(i,j));
                else if ( pre[i] &&  pre[j]) K22.PushEntry (loc[i], loc[j], K(i,j));
            }
        }
    }
    // augment A11 and set A11mat
    for (size_t i=0; i<pDOFs.Size(); ++i) A11.PushEntry (pDOFs[i],pDOFs[i], 1.0);
    A11mat.Set (A11);
}

inline void Solver::AssembleKMA (double C1, double C2)
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
                     if (!pre[i] && !pre[j]) { A11.PushEntry (loc[i], loc[j], C1*M(i,j) + C2*K(i,j));
                                               K11.PushEntry (loc[i], loc[j], K(i,j));  M11.PushEntry (loc[i], loc[j], M(i,j)); }
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

inline void Solver::AssembleKCMA (double C1, double C2, double C3)
{
    if (DampTy!=Rayleigh_t) throw new Fatal("Solver::AssembleKCMA: Only Rayleigh damping implemented yet");
    K11.ResetTop(); // reset top (position to insert new values) => clear triplet
    K12.ResetTop();
    K21.ResetTop();
    K22.ResetTop();
    C11.ResetTop(); // reset top (position to insert new values) => clear triplet
    C12.ResetTop();
    C21.ResetTop();
    C22.ResetTop();
    M11.ResetTop(); // reset top (position to insert new values) => clear triplet
    M12.ResetTop();
    M21.ResetTop();
    M22.ResetTop();
    A11.ResetTop(); // reset top (position to insert new values) => clear triplet
    for (size_t k=0; k<Dom.Eles.Size(); ++k)
    {
        // calc K and M
        Mat_t         K, M; // matrices
        Array<size_t> loc;  // location array
        Array<bool>   pre;  // prescribed U ?
        Dom.Eles[k]->CalcK  (K);
        Dom.Eles[k]->CalcM  (M);
        Dom.Eles[k]->GetLoc (loc, pre);
        // calc C
        Mat_t C(K.num_rows(),K.num_cols());
        C = DampAlp*M + DampBet*K;
        // set K, C, M, and A matrices
        for (size_t i=0; i<loc.Size(); ++i)
        {
            for (size_t j=0; j<loc.Size(); ++j)
            {
                     if (!pre[i] && !pre[j]) { A11.PushEntry (loc[i], loc[j], C1*M(i,j) + C2*C(i,j) + C3*K(i,j));
                                               K11.PushEntry (loc[i], loc[j], K(i,j));  C11.PushEntry (loc[i], loc[j], C(i,j));  M11.PushEntry (loc[i], loc[j], M(i,j)); }
                else if (!pre[i] &&  pre[j]) { K12.PushEntry (loc[i], loc[j], K(i,j));  C12.PushEntry (loc[i], loc[j], C(i,j));  M12.PushEntry (loc[i], loc[j], M(i,j)); }
                else if ( pre[i] && !pre[j]) { K21.PushEntry (loc[i], loc[j], K(i,j));  C21.PushEntry (loc[i], loc[j], C(i,j));  M21.PushEntry (loc[i], loc[j], M(i,j)); }
                else if ( pre[i] &&  pre[j]) { K22.PushEntry (loc[i], loc[j], K(i,j));  C22.PushEntry (loc[i], loc[j], C(i,j));  M22.PushEntry (loc[i], loc[j], M(i,j)); }
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
    if (K11.Top()==0 || CteTg==false) AssembleKA (); // not constant tangent matrices => non-linear problems

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

    // calc dU and dF
    Sparse::SubMult (K12, W, W);     // W1  -= K12*dU2
    UMFPACK::Solve  (A11mat, W, dU); // dU   = inv(A11mat)*W
    Sparse::AddMult (K21, dU, dF);   // dF2 += K21*dU1
    Sparse::AddMult (K22, dU, dF);   // dF2 += K22*dU2
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
    A11.AllocSpace (NEQ,NEQ,K11_size+pDOFs.Size()); // augmented
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
        if (DampTy!=None_t)
        {
            C11.AllocSpace (NEQ,NEQ,K11_size);
            C12.AllocSpace (NEQ,NEQ,K12_size);
            C21.AllocSpace (NEQ,NEQ,K21_size);
            C22.AllocSpace (NEQ,NEQ,K22_size);
        }
    }

    // resize workspace
    W.change_dim (NEQ);

    // build U and F (and V, U2, F1, Fbar)
    U.change_dim (NEQ);
    F.change_dim (NEQ);
    if (Transient)
    {
        V .change_dim (NEQ);
        U2.change_dim (NEQ);
        F1.change_dim (NEQ);
        set_to_zero   (V);
        set_to_zero   (U2);
        set_to_zero   (F1);
    }
    for (size_t i=0; i<Dom.Nods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.Nods[i]->nDOF(); ++j)
        {
            long eq = Dom.Nods[i]->EQ[j];
            U(eq)   = Dom.Nods[i]->U [j];
            F(eq)   = Dom.Nods[i]->F [j];
            if (Transient)
            {
                if (Dom.Nods[i]->pU[j]) U2(eq) = U(eq) + Dom.Nods[i]->DU[j];
                else                    F1(eq) = F(eq) + Dom.Nods[i]->DF[j];
            }
        }
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

inline void Solver::_calc_Fbar (double t, Vec_t & Fb)
{
    // calc Fbar = F1 - K12*U2
    Fb = F1;
    for (int k=0; k<K12.Top(); ++k)
        Fb(K12.Ai(k)) -= K12.Ax(k) * U2(K12.Aj(k));
    if (Dom.NDim==3)
    {
        for (size_t i=0; i<Dom.NodsF.Size(); ++i)
        {
            (*Dom.CalcF[i]) (t, Fb(Dom.NodsF[i]->EQ[Dom.NodsF[i]->FMap("fx")]),
                                Fb(Dom.NodsF[i]->EQ[Dom.NodsF[i]->FMap("fy")]),
                                Fb(Dom.NodsF[i]->EQ[Dom.NodsF[i]->FMap("fz")]));
        }
    }
    else
    {
        double fz = 0.0;
        for (size_t i=0; i<Dom.NodsF.Size(); ++i)
        {
            (*Dom.CalcF[i]) (t, Fb(Dom.NodsF[i]->EQ[Dom.NodsF[i]->FMap("fx")]),
                                Fb(Dom.NodsF[i]->EQ[Dom.NodsF[i]->FMap("fy")]),
                                fz);
        }
    }
}

}; // namespace FEM

#endif // MECHSYS_FEM_SOLVER_H
