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

using std::cout;
using std::endl;

namespace FEM
{

class Solver
{
public:
    // enum
    enum Scheme_t  { FE_t, ME_t, NR_t };   ///< Steady time integration scheme: Forward-Euler, Modified-Euler, Newton-Rhapson
    enum TScheme_t { SS11_t };             ///< Transient time integration scheme: (Single step/O1/1st order)
    enum DScheme_t { SS22_t, GN22_t };     ///< Dynamic time integration scheme: (Single step/O2/2nd order), (Generalized Newmark/O2/2nd order)
    enum Damping_t { None_t, Rayleigh_t }; ///< Damping type: none, Rayleigh type (C=alp*M+bet*K)

    // typedefs
    typedef void (*pDbgFun) (Solver const & Sol, void * DbgDat); ///< Pointer to Debug function

    // Constructor
    Solver (Domain const & Dom, pDbgFun DbgFun=NULL, void * DbgDat=NULL);

    // Methods
    void Solve        (size_t NInc=1, Array<double> * Weights=NULL); ///< Solve steady/equilibrium equation
    void TransSolve   (double tf, double dt, double dtOut);          ///< Solve transient equation
    void DynSolve     (double tf, double dt, double dtOut);          ///< Solve dynamic equation
    void AssembleKA   ();                                            ///< A = K11
    void AssembleKMA  (double Coef1, double Coef2);                  ///< A = Coef1*M + Coef2*K
    void AssembleKCMA (double Coef1, double Coef2, double Coef3);    ///< A = Coef1*M + Coef2*C + Coef3*K
    void TgIncs       (double dT, Vec_t & dU, Vec_t & dF);           ///< Tangent increments: dU = inv(K)*dF
    void Initialize   (bool Transient=false);                        ///< Initialize global matrices and vectors

    // Data
    Domain const & Dom;    ///< Domain
    pDbgFun        DbgFun; ///< Debug function
    void         * DbgDat; ///< Debug data
    double         Time;   ///< Current time (t)
    size_t         Inc;    ///< Current increment
    size_t         Stp;    ///< Current (sub) step
    size_t         It;     ///< Current iteration
    size_t         NEQ;    ///< Total number of equations (DOFs)
    Array<long>    uDOFs;  ///< unknown DOFs
    Array<long>    pDOFs;  ///< prescribed DOFs (known equations)
    double         NormR;  ///< Euclidian norm of residual (R)
    double         TolR;   ///< Tolerance for the norm of residual

    // Triplets and sparse matrices
    Sparse::Triplet<double,int> K11,K12,K21,K22; ///< Stiffness matrices
    Sparse::Triplet<double,int> C11,C12,C21,C22; ///< Damping matrices
    Sparse::Triplet<double,int> M11,M12,M21,M22; ///< Mass matrices
    Sparse::Triplet<double,int> A11;             ///< A=K  or  A=C1*M+C2*K  or  A=C1*M+C2*C+C3*K

    // Vectors
    Vec_t R;        // Residual
    Vec_t F, F_int; // External and internal forces
    Vec_t W, U;     // Workspace, displacement
    Vec_t DU2, DF1; // Prescribed DU and DF
    Vec_t V, A;     // (Transient/Dynamic) velocity and acceleration

    // Constants for integration
    Scheme_t  Scheme;  ///< Scheme: FE_t (Forward-Euler), ME_t (Modified-Euler)
    size_t    nSS;     ///< FE and NR: number of substeps
    double    STOL;    ///< ME:
    double    dTini;   ///< ME:
    double    mMin;    ///< ME:
    double    mMax;    ///< ME:
    size_t    MaxSS;   ///< ME:
    bool      CteTg;   ///< Constant tangent matrices (linear problems) => K will be calculated once
    bool      ModNR;   ///< Modified Newton-Rhapson ?
    size_t    MaxIt;   ///< Max iterations (for Newton-Rhapson)
    TScheme_t TScheme; ///< Transient scheme
    double    Theta;   ///< Transient scheme constant
    DScheme_t DScheme; ///< Dynamic scheme
    Damping_t DampTy;  ///< Damping type
    double    DampAlp; ///< Rayleigh damping alpha coefficient
    double    DampBet; ///< Rayleigh damping beta coefficient
    double    DynTh1;  ///< Dynamic coefficient Theta 1
    double    DynTh2;  ///< Dynamic coefficient Theta 2

private:
    void _FE_update  (double tf); ///< (Forward-Euler)  Update Time and elements to tf
    void _ME_update  (double tf); ///< (Modified-Euler) Update Time and elements to tf
    void _NR_update  (double tf); ///< (Newton-Rhapson) Update Time and elements to tf
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Solver::Solver (Domain const & TheDom, pDbgFun TheDbgFun, void * TheDbgDat)
    : Dom     (TheDom),
      DbgFun  (TheDbgFun),
      DbgDat  (TheDbgDat),
      Time    (0.0),
      Inc     (0),
      Stp     (0),
      It      (0),
      TolR    (1.0e-9),
      Scheme  (ME_t),
      nSS     (1),
      STOL    (1.0e-5),
      dTini   (1.0),
      mMin    (0.1),
      mMax    (10.0),
      MaxSS   (2000),
      CteTg   (false),
      ModNR   (false),
      MaxIt   (20),
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

inline void Solver::Solve (size_t NInc, Array<double> * Weights)
{
    // info
    double start = std::clock();

    // initialize global matrices and vectors
    Initialize ();

    // output initial state
    std::cout << "\n[1;37m--- Stage solution -----------------------------------------------------------\n";
    std::cout << Util::_6_3 << "Time" <<                                       Util::_8s <<"Norm(R)" << "[0m\n";
    std::cout << Util::_6_3 <<  Time  << (NormR>TolR?"[1;31m":"[1;32m") << Util::_8s << NormR    << "[0m\n";
    Dom.OutResults (Time);

    // weights
    bool   del_weights = false;
    double sum_weights = 0.0;
    if (Weights==NULL)
    {
        Weights = new Array<double>;
        Weights->Resize (NInc);
        for (size_t i=0; i<NInc; ++i) (*Weights)[i] = 1.0/NInc;
        del_weights = true;
    }
    else if (NInc!=Weights->Size()) throw new Fatal("Solver::Solve: Array with weights must have size equal to NDiv (%d)",NInc);
    for (size_t i=0; i<NInc; ++i) sum_weights += (*Weights)[i];
    if (fabs(sum_weights-1.0)>1.0e-9) throw new Fatal("Solver::Solver: Sum of weights must be equal to 1.0 (Sum(W)=%g is invalid)",sum_weights);

    // solve
    double t0 = Time;     // current time
    double tf = t0 + 1.0; // final time
    double Dt = tf - t0;  // total time increment
    double dt, tout;      // timestep and time for output
    for (Inc=0; Inc<NInc; ++Inc)
    {
        // timestep
        dt   = (*Weights)[Inc]*Dt; // timestep
        tout = Time + dt;          // time for output

        // update Time and elements to tout
        String str;
        if      (Scheme==FE_t) { _FE_update (tout);  str.Printf("Forward-Euler (FE): nss = %d",Stp); }
        else if (Scheme==ME_t) { _ME_update (tout);  str.Printf("Modified-Euler (ME): nss = %d",Stp); }
        else if (Scheme==NR_t) { _NR_update (tout);  str.Printf("Newton-Rhapson (NR): nit = %d",It); }
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
        R     = F - F_int;
        NormR = Norm(R);

        // output
        std::cout << Util::_6_3 << Time << (NormR>TolR?"[1;31m":"[1;32m") << Util::_8s << NormR << "[0m    " << str << "\n";
        Dom.OutResults (Time);

        // next tout
        tout = Time + dt;
    }

    // info
    double total = std::clock() - start;
    std::cout << Util::_reset << "[1;36m Time elapsed = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";

    // clean up
    if (del_weights) delete Weights;
}

inline void Solver::TransSolve (double tf, double dt, double dtOut)
{
}

inline void Solver::DynSolve (double tf, double dt, double dtOut)
{
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
    // augment A11
    for (size_t i=0; i<pDOFs.Size(); ++i) A11.PushEntry (pDOFs[i],pDOFs[i], 1.0);
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
    // augment M11 and A11
    for (size_t i=0; i<pDOFs.Size(); ++i)
    {
        M11.PushEntry (pDOFs[i],pDOFs[i], 1.0);
        A11.PushEntry (pDOFs[i],pDOFs[i], 1.0);
    }
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
    // augment A11
    for (size_t i=0; i<pDOFs.Size(); ++i) A11.PushEntry (pDOFs[i],pDOFs[i], 1.0);
}

inline void Solver::TgIncs (double dT, Vec_t & dU, Vec_t & dF)
{
    // assemble global K matrix
    if (K11.Top()==0 || CteTg==false) AssembleKA (); // not constant tangent matrices => non-linear problems

    // assemble dF and W (workspace) vectors
    for (size_t i=0; i<uDOFs.Size(); ++i)
    {
        dF(uDOFs[i]) = dT*DF1(uDOFs[i]); // set dF1 equal to dT*DF1
        W (uDOFs[i]) = dF(uDOFs[i]);     // set W1  equal to dF1
    }
    for (size_t i=0; i<pDOFs.Size(); ++i)
    {
        dF(pDOFs[i]) = 0.0;              // clear dF2
        W (pDOFs[i]) = dT*DU2(pDOFs[i]); // set W2 equal to dT*DU2
    }

    // calc dU and dF
    Sparse::SubMult (K12,  W,  W); // W1  -= K12*dU2
    UMFPACK::Solve  (A11,  W, dU); // dU   = inv(A11)*W
    Sparse::AddMult (K21, dU, dF); // dF2 += K21*dU1
    Sparse::AddMult (K22, dU, dF); // dF2 += K22*dU2
}

inline void Solver::Initialize (bool Transient)
{
    // assign equation numbers
    NEQ = 0;
    uDOFs.Resize (0); // unknown DOFs
    pDOFs.Resize (0); // prescribed DOFs
    for (size_t i=0; i<Dom.Nods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.Nods[i]->nDOF(); ++j)
        {
            Dom.Nods[i]->EQ[j] = NEQ;
            if (Dom.Nods[i]->pU[j]) pDOFs.Push (NEQ);
            else                    uDOFs.Push (NEQ);
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
        M11.AllocSpace (NEQ,NEQ,K11_size+pDOFs.Size());
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

    // initialize variables
    R    .change_dim (NEQ);  set_to_zero (R);
    F    .change_dim (NEQ);  set_to_zero (F);
    F_int.change_dim (NEQ);  set_to_zero (F_int);
    W    .change_dim (NEQ);  set_to_zero (W);
    U    .change_dim (NEQ);  set_to_zero (U);
    DU2  .change_dim (NEQ);  set_to_zero (DU2);
    DF1  .change_dim (NEQ);  set_to_zero (DF1);
    if (Transient)
    {
        V.change_dim (NEQ);  set_to_zero (V);
        A.change_dim (NEQ);  set_to_zero (A);
    }

    // initialize F_int (must come before setting up of F(ext))
    for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->InitFint (F_int);

    // set variables
    for (size_t i=0; i<Dom.Nods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.Nods[i]->nDOF(); ++j)
        {
            long eq = Dom.Nods[i]->EQ[j];
            U (eq)  = Dom.Nods[i]->U [j];
            F (eq)  = Dom.Nods[i]->F [j];
            if (Dom.Nods[i]->pU[j]) DU2(eq) = Dom.Nods[i]->DU[j];
            else                    DF1(eq) = Dom.Nods[i]->DF[j];
        }
    }

    // calc residual
    R     = F - F_int;
    NormR = Norm(R);
    //std::cout << "F     = \n" << PrintVector(F);
    //std::cout << "F_int = \n" << PrintVector(F_int);
}

inline void Solver::_FE_update (double tf)
{
    // auxiliar vectors
    Vec_t dU(NEQ), dF(NEQ);

    double dt = (tf-Time)/nSS;
    for (Stp=0; Stp<nSS; ++Stp)
    {
        // calculate tangent increments
        TgIncs (dt, dU, dF);

        // update elements
        for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->UpdateState (dU, &F_int);

        // update U, F, and Time
        U    += dU;
        F    += dF;
        Time += dt;
    }
}

inline void Solver::_ME_update (double tf)
{
    // auxiliar vectors
    Vec_t dU_fe(NEQ), dU_tm(NEQ), dU_me(NEQ), U_me(NEQ), U_dif(NEQ);
    Vec_t dF_fe(NEQ), dF_tm(NEQ), dF_me(NEQ), F_me(NEQ), F_dif(NEQ);

    // for each pseudo time T
    double T   = 0.0;
    double dT  = dTini;
    double Dt  = tf-Time;
    for (Stp=0; Stp<MaxSS; ++Stp)
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
            if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
            for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->UpdateState (dU_me, &F_int);
            T    += dT;
            U     = U_me;
            F     = F_me;
            Time += dt;
            if (m>mMax) m = mMax;
        }
        else if (m<mMin) m = mMin;

        // change next step size
        dT = m * dT;

        // check for last increment
        if (dT>1.0-T) dT = 1.0-T;
    }
    if (Stp>=MaxSS) throw new Fatal("Solver:_ME_update: Modified-Euler (global integration) did not converge for %d steps",Stp);
}

inline void Solver::_NR_update (double tf)
{
    // auxiliar vectors
    Vec_t dU(NEQ), dF(NEQ);

    size_t max_it = It;
    double dt     = (tf-Time)/nSS;
    for (Stp=0; Stp<nSS; ++Stp)
    {
        // calculate tangent increments
        TgIncs (dt, dU, dF);

        // update elements
        for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->UpdateState (dU, &F_int);

        // update U, F, and Time
        U    += dU;
        F    += dF;
        Time += dt;

        // residual
        R     = F - F_int;
        NormR = Norm(R);

        // iterations
        It = 0;
        while (NormR>TolR && It<MaxIt)
        {
            // debug
            if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);

            // assemble global K matrix
            if (!ModNR) AssembleKA ();

            // clear unbalanced forces related to supports
            for (size_t i=0; i<pDOFs.Size(); ++i) R(pDOFs[i]) = 0.0;

            // calc corrector dU
            UMFPACK::Solve (A11, R, dU); // dU = inv(A11)*R

            // update elements and displacements
            for (size_t i=0; i<Dom.Eles.Size(); ++i) Dom.Eles[i]->UpdateState (dU, &F_int);
            U += dU;

            // residual
            R     = F - F_int;
            NormR = Norm(R);

            // next iteration
            It++;
        }
        if (It>=MaxIt) throw new Fatal("Solver::_NR_update: Newton-Rhapson did not converge after %d iterations",It);
        if (It>max_it) max_it = It;
    }

    // return max number of iterations
    It = max_it;
}

}; // namespace FEM

#endif // MECHSYS_FEM_SOLVER_H
