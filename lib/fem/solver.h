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
#include <iostream> // for cout

// Blitz++
#include <blitz/tinyvec-et.h>

// MPI
#ifdef HAS_MPI
  #include <mpi.h>
#endif

// MechSys
#include <mechsys/fem/node.h>
#include <mechsys/fem/element.h>
#include <mechsys/fem/domain.h>
#include <mechsys/linalg/sparse_triplet.h>
#include <mechsys/linalg/sparse_matrix.h>
#include <mechsys/linalg/umfpack.h>
#include <mechsys/linalg/mumps.h>
#include <mechsys/util/stopwatch.h>

namespace FEM
{

class Solver
{
public:
    // enum
    enum Scheme_t  { FE_t, ME_t, NR_t };             ///< Steady time integration scheme: Forward-Euler, Modified-Euler, Newton-Rhapson
    enum TScheme_t { SS11_t };                       ///< Transient time integration scheme: (Single step/O1/1st order)
    enum DScheme_t { GN22_t };                       ///< Dynamic time integration scheme: (Single step/O2/2nd order), (Generalized Newmark/O2/2nd order)
    enum Damping_t { None_t, Rayleigh_t, HMCoup_t }; ///< Damping type: none, Rayleigh type (C=alp*M+bet*K), HydroMechCoupling

    // typedefs
    typedef void (*pOutFun) (Solver const & Sol, void * OutDat); ///< Pointer to output function

    // Constructor
    Solver (Domain const & Dom, pOutFun OutFun=NULL, void * OutDat=NULL,
                                pOutFun DbgFun=NULL, void * DbgDat=NULL); ///< Allocate solver object

    // Methods
    void Solve          (size_t NInc=1, char const * FileKey=NULL);                      ///< Solve quasi-static problem
    void TransSolve     (double tf, double dt, double dtOut, char const * FileKey=NULL); ///< Solve transient problem
    void DynSolve       (double tf, double dt, double dtOut, char const * FileKey=NULL); ///< Solve dynamic problem
    void AssembleKA     ();                                                              ///< A = K11
    void AssembleKMA    (double Coef1, double Coef2);                                    ///< A = Coef1*M + Coef2*K
    void AssembleKCMA   (double Coef1, double Coef2, double Coef3);                      ///< A = Coef1*M + Coef2*C + Coef3*K
    void TgIncs         (double dT, Vec_t & dU, Vec_t & dF);                             ///< Tangent increments: dU = inv(K)*dF
    void UpdateElements (Vec_t const & dU, bool CalcFint);                               ///< Update elements
    void Initialize     (bool Transient=false);                                          ///< Initialize global matrices and vectors
    void SetScheme      (char const * StrScheme);                                        ///< Set solution scheme: 'FE', 'ME', 'NR'
    void SetIncsW       (size_t NInc, bool NonLinWei=false);                             ///< Set weights for quasi-static problem (If Weights.Size()==0: generate weights)
    bool ResidOK        () const;                                                        ///< Check if the residual is OK

    // Data (read-only)
    Domain  const & Dom;      ///< Domain
    pOutFun         OutFun;   ///< Output function (called during output)
    void          * OutDat;   ///< Debug data (to be used with either OutFun or DbgFun)
    pOutFun         DbgFun;   ///< Debug function (called everytime for some internal (update) methods)
    void          * DbgDat;   ///< Debug data (to be used with either OutFun or DbgFun)
    bool            Root;     ///< Root processor ?
    double          Time;     ///< Current time (t)
    size_t          Inc;      ///< Current increment
    size_t          IdxOut;   ///< Counter for generating VTU files in DynSolve
    size_t          Stp;      ///< Current (sub) step
    size_t          It;       ///< Current iteration
    size_t          NEq;      ///< Total number of equations (DOFs)
    size_t          NLag;     ///< Number of Lagrange multipliers
    Array<long>     pEQ;      ///< prescribed equations
    Array<bool>     pU;       ///< prescribed U
    double          NormR;    ///< Euclidian norm of residual (R)
    double          MaxNormF; ///< Max(Norm(F), Norm(Fint))
    Array<Node*>    ActNods;  ///< Active nodes
    Array<Element*> ActEles;  ///< Active elements

    // Data (read-write)
    Scheme_t      Scheme;   ///< Scheme: FE_t (Forward-Euler), ME_t (Modified-Euler)
    bool          CalcWork; ///< Calc work done == twice the stored (elastic) strain energy ?
    size_t        nSS;      ///< FE and NR: number of substeps
    double        STOL;     ///< ME:
    double        dTini;    ///< ME:
    double        dTlast;   ///< ME:
    double        mMin;     ///< ME:
    double        mMax;     ///< ME:
    size_t        MaxSS;    ///< ME:
    bool          SSOut;    ///< SubSteps ouput in ME ?
    bool          CteTg;    ///< Constant tangent matrices (linear problems) => K will be calculated once
    bool          ModNR;    ///< Modified Newton-Rhapson ?
    double        TolR;     ///< Tolerance for the norm of residual
    bool          CorR;     ///< Correct residual ?
    size_t        MaxIt;    ///< Max iterations (for Newton-Rhapson)
    TScheme_t     TScheme;  ///< Transient scheme
    double        Theta;    ///< Transient scheme constant
    DScheme_t     DScheme;  ///< Dynamic scheme
    Damping_t     DampTy;   ///< Damping type
    double        DampAm;   ///< Rayleigh damping Am coefficient (C = Am*M + Ak*K)
    double        DampAk;   ///< Rayleigh damping Ak coefficient (C = Am*M + Ak*K)
    double        DynTh1;   ///< Dynamic coefficient Theta 1
    double        DynTh2;   ///< Dynamic coefficient Theta 2
    Array<double> IncsW;    ///< Increments weights used in Solve

    // Triplets and sparse matrices
    Sparse::Triplet<double,int> K11,K12,K21,K22; ///< Stiffness matrices
    Sparse::Triplet<double,int> C11,C12,C21,C22; ///< Damping matrices
    Sparse::Triplet<double,int> M11,M12,M21,M22; ///< Mass matrices
    Sparse::Triplet<double,int> A11;             ///< A=K  or  A=C1*M+C2*K  or  A=C1*M+C2*C+C3*K

    // Vectors
    Vec_t R;            ///< Residual
    Vec_t F, F_int;     ///< External and internal forces
    Vec_t W, U;         ///< Workspace, displacement
    Vec_t V, A;         ///< (Transient/Dynamic) velocity and acceleration
    Vec_t dU, dF;       ///< Increments of U and F
    Vec_t Us, Vs;       ///< starred variables (for GN22)
    Vec_t TmpVec;       ///< Temporary vector (for parallel Allreduce)

private:
    void _set_A_Lag   ();                          ///< Set A matrix due to Lagrange multipliers
    void _cal_resid   (bool WithAccel=false);      ///< Calculate residual
    void _cor_resid   (Vec_t & dU, Vec_t & dF);    ///< Correct residual
    void _FE_update   (double tf);                 ///< (Forward-Euler)  Update Time and elements to tf
    void _ME_update   (double tf);                 ///< (Modified-Euler) Update Time and elements to tf
    void _NR_update   (double tf);                 ///< (Newton-Rhapson) Update Time and elements to tf
    void _presc_F     (double t);                  ///< Calculate prescribed F(t=Time)
    void _GN22_update (double tf, double dt);      ///< (Generalized-Newmark) Update Time and elements to tf
    void _time_print  (char const * Comment=NULL); ///< Print timestep data
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Solver::Solver (Domain const & TheDom, pOutFun TheOutFun, void * TheOutDat, pOutFun TheDbgFun, void * TheDbgDat)
    : Dom     (TheDom),
      OutFun  (TheOutFun),
      OutDat  (TheOutDat),
      DbgFun  (TheDbgFun),
      DbgDat  (TheDbgDat),
      Root    (true),
      Time    (0.0),
      Inc     (0),
      IdxOut  (0),
      Stp     (0),
      It      (0),
      Scheme  (ME_t),
      CalcWork(false),
      nSS     (1),
      STOL    (1.0e-5),
      dTini   (1.0),
      dTlast  (-1.0),
      mMin    (0.1),
      mMax    (10.0),
      MaxSS   (2000),
      SSOut   (false),
      CteTg   (false),
      ModNR   (false),
      TolR    (1.0e-7),
      CorR    (true),
      MaxIt   (20),
      TScheme (SS11_t),
      Theta   (2./3.),
      DScheme (GN22_t),
      DampTy  (None_t),
      DampAm  (0.005),
      DampAk  (0.5),
      DynTh1  (0.5),
      DynTh2  (0.5)
{
#if HAS_MPI
    if (FEM::Domain::PARA && MPI::COMM_WORLD.Get_rank()!=0) Root = false;
#endif
}

inline void Solver::Solve (size_t NInc, char const * FileKey)
{
    // info
    Util::Stopwatch stopwatch(/*only_root*/FEM::Domain::PARA);

    // initialize global matrices and vectors
    Initialize ();

    // output initial state
    if (Root)
    {
        if      (Scheme==FE_t) _time_print ("Quasi-static --- FE");
        else if (Scheme==ME_t) _time_print ("Quasi-static --- ME");
        else if (Scheme==NR_t) _time_print ("Quasi-static --- NR");
    }
    if (IdxOut==0)
    {
        Dom.OutResults (Time, F_int);
        if (OutFun!=NULL) (*OutFun) ((*this), OutDat);
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
    }

    // weights
    if (IncsW.Size()!=NInc) SetIncsW (NInc);

    // solve
    double t0 = Time;     // current time
    double tf = t0 + 1.0; // final time
    double Dt = tf - t0;  // total time increment
    double dt, tout;      // timestep and time for output
    for (Inc=0; Inc<NInc; ++Inc)
    {
        // timestep
        dt   = IncsW[Inc]*Dt; // timestep
        tout = Time + dt;     // time for output

        // update U, F, Time and elements to tout
        if      (Scheme==FE_t) _FE_update (tout);
        else if (Scheme==ME_t) _ME_update (tout);
        else if (Scheme==NR_t) _NR_update (tout);

        // update nodes to tout
        for (size_t i=0; i<ActNods.Size(); ++i)
        {
            for (size_t j=0; j<ActNods[i]->nDOF(); ++j)
            {
                long eq = ActNods[i]->EQ[j];
                ActNods[i]->U[j] = U(eq);
                ActNods[i]->F[j] = F(eq);
            }
        }

        // output
        IdxOut++;
        if (Root) _time_print ();
        if (Scheme!=ME_t) Dom.OutResults (Time, F_int);
        else if (!SSOut)  Dom.OutResults (Time, F_int);
        if (OutFun!=NULL) (*OutFun) ((*this), OutDat);

        // write VTU
        if (FileKey!=NULL)
        {
            String fkey;
            fkey.Printf  ("%s_%08d", FileKey, IdxOut);
            Dom.WriteVTU (fkey.CStr());
        }

        // next tout
        tout = Time + dt;
    }
}

inline void Solver::TransSolve (double tf, double dt, double dtOut, char const * FileKey)
{
}

inline void Solver::DynSolve (double tf, double dt, double dtOut, char const * FileKey)
{
    // info
    Util::Stopwatch stopwatch(/*only_root*/FEM::Domain::PARA);

    // initialize global matrices and vectors
    Initialize (/*Transient*/true);

    // output initial state
    if (Root)
    {
        if (DScheme==GN22_t) _time_print ("Dynamic ------ GN22");
    }
    if (IdxOut==0)
    {
        Dom.OutResults (Time, F_int);
        if (OutFun!=NULL) (*OutFun) ((*this), OutDat);
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
    }

    // solve
    double tout = Time + dtOut; // time for output
    while (Time<tf)
    {
        // update U, F, Time and elements to tout
        if (DScheme==GN22_t) _GN22_update (tout,dt);

        // update nodes to tout
        for (size_t i=0; i<ActNods.Size(); ++i)
        {
            for (size_t j=0; j<ActNods[i]->nDOF(); ++j)
            {
                long eq = ActNods[i]->EQ[j];
                ActNods[i]->U[j] = U(eq);
                ActNods[i]->F[j] = F(eq);
            }
        }

        // output
        IdxOut++;
        if (Root) _time_print ();
        Dom.OutResults (Time, F_int);
        if (OutFun!=NULL) (*OutFun) ((*this), OutDat);

        // write VTU
        if (FileKey!=NULL)
        {
            String fkey;
            fkey.Printf  ("%s_%08d", FileKey, IdxOut);
            Dom.WriteVTU (fkey.CStr());
        }

        // next tout
        tout = Time + dtOut;
    }
}

inline void Solver::AssembleKA ()
{
    if (CteTg && A11.Top()>0) return; // constant tangent matrices => linear problems
    A11.ResetTop(); // reset top (position to insert new values) => clear triplet
    K12.ResetTop();
    K21.ResetTop();
    K22.ResetTop();
    for (size_t k=0; k<ActEles.Size(); ++k)
    {
        Mat_t         K;   // K matrix
        Array<size_t> loc; // location array
        ActEles[k]->CalcK  (K);
        ActEles[k]->GetLoc (loc);
        for (size_t i=0; i<loc.Size(); ++i)
        {
            for (size_t j=0; j<loc.Size(); ++j)
            {
                     if (!pU[loc[i]] && !pU[loc[j]]) A11.PushEntry (loc[i], loc[j], K(i,j));
                else if (!pU[loc[i]] &&  pU[loc[j]]) K12.PushEntry (loc[i], loc[j], K(i,j));
                else if ( pU[loc[i]] && !pU[loc[j]]) K21.PushEntry (loc[i], loc[j], K(i,j));
                else if ( pU[loc[i]] &&  pU[loc[j]]) K22.PushEntry (loc[i], loc[j], K(i,j));
            }
        }
    }

    /*
    Sparse::Matrix<double,int> k11(A11), k12(K12), k21(K21), k22(K22);
    Mat_t kk11, kk12, kk21, kk22, kk;
    k11.GetDense(kk11), k12.GetDense(kk12), k21.GetDense(kk21), k22.GetDense(kk22);
    kk = kk11 + kk12 + kk21 + kk22;
    std::cout << "K   = \n" << PrintMatrix(kk);
    std::cout << "K11 = \n" << PrintMatrix(kk11);
    */

    // augment A11
    for (size_t i=0; i<pEQ.Size(); ++i) A11.PushEntry (pEQ[i],pEQ[i], 1.0);
    _set_A_Lag ();

    /*
    Sparse::Matrix<double,int> Asp(A11);
    Mat_t Adense;
    Asp.GetDense(Adense);
    std::cout << "\nA =\n" << PrintMatrix(Adense,"%12.2f");
    std::cout << "\ndet(A) = " << UMFPACK::Det(Asp) << "\n\n";
    */
}

inline void Solver::AssembleKMA (double C1, double C2)
{
    if (CteTg && A11.Top()>0) return; // constant tangent matrices => linear problems
    K11.ResetTop(); // reset top (position to insert new values) => clear triplet
    K12.ResetTop();
    K21.ResetTop();
    K22.ResetTop();
    M11.ResetTop(); // reset top (position to insert new values) => clear triplet
    M12.ResetTop();
    M21.ResetTop();
    M22.ResetTop();
    A11.ResetTop(); // reset top (position to insert new values) => clear triplet
    for (size_t k=0; k<ActEles.Size(); ++k)
    {
        Mat_t         K, M; // matrices
        Array<size_t> loc;  // location array
        ActEles[k]->CalcK  (K);
        ActEles[k]->CalcM  (M);
        ActEles[k]->GetLoc (loc);
        for (size_t i=0; i<loc.Size(); ++i)
        {
            for (size_t j=0; j<loc.Size(); ++j)
            {
                     if (!pU[loc[i]] && !pU[loc[j]]) { A11.PushEntry (loc[i], loc[j], C1*M(i,j) + C2*K(i,j));
                                                       K11.PushEntry (loc[i], loc[j], K(i,j));  M11.PushEntry (loc[i], loc[j], M(i,j)); }
                else if (!pU[loc[i]] &&  pU[loc[j]]) { K12.PushEntry (loc[i], loc[j], K(i,j));  M12.PushEntry (loc[i], loc[j], M(i,j)); }
                else if ( pU[loc[i]] && !pU[loc[j]]) { K21.PushEntry (loc[i], loc[j], K(i,j));  M21.PushEntry (loc[i], loc[j], M(i,j)); }
                else if ( pU[loc[i]] &&  pU[loc[j]]) { K22.PushEntry (loc[i], loc[j], K(i,j));  M22.PushEntry (loc[i], loc[j], M(i,j)); }
            }
        }
    }
    // augment A11
    for (size_t i=0; i<pEQ.Size(); ++i) A11.PushEntry (pEQ[i],pEQ[i], 1.0);
}

inline void Solver::AssembleKCMA (double C1, double C2, double C3)
{
    if (CteTg && A11.Top()>0) return; // constant tangent matrices => linear problems
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
    for (size_t k=0; k<ActEles.Size(); ++k)
    {
        // matrices
        Mat_t M, C, K;
        if      (DampTy==HMCoup_t) ActEles[k]->CalcKCM (K, C, M);
        else if (DampTy==Rayleigh_t)
        {
            ActEles[k]->CalcK (K);
            ActEles[k]->CalcM (M);
            C = DampAm*M + DampAk*K;
        }
        // set K, C, M, and A matrices
        Array<size_t> loc;
        ActEles[k]->GetLoc (loc);
        for (size_t i=0; i<loc.Size(); ++i)
        {
            for (size_t j=0; j<loc.Size(); ++j)
            {
                     if (!pU[loc[i]] && !pU[loc[j]]) { A11.PushEntry (loc[i], loc[j], C1*M(i,j) + C2*C(i,j) + C3*K(i,j));
                                                       K11.PushEntry (loc[i], loc[j], K(i,j));  C11.PushEntry (loc[i], loc[j], C(i,j));  M11.PushEntry (loc[i], loc[j], M(i,j)); }
                else if (!pU[loc[i]] &&  pU[loc[j]]) { K12.PushEntry (loc[i], loc[j], K(i,j));  C12.PushEntry (loc[i], loc[j], C(i,j));  M12.PushEntry (loc[i], loc[j], M(i,j)); }
                else if ( pU[loc[i]] && !pU[loc[j]]) { K21.PushEntry (loc[i], loc[j], K(i,j));  C21.PushEntry (loc[i], loc[j], C(i,j));  M21.PushEntry (loc[i], loc[j], M(i,j)); }
                else if ( pU[loc[i]] &&  pU[loc[j]]) { K22.PushEntry (loc[i], loc[j], K(i,j));  C22.PushEntry (loc[i], loc[j], C(i,j));  M22.PushEntry (loc[i], loc[j], M(i,j)); }
            }
        }
    }
    // augment A11
    for (size_t i=0; i<pEQ.Size(); ++i) A11.PushEntry (pEQ[i],pEQ[i], 1.0);
}

inline void Solver::TgIncs (double dT, Vec_t & dU, Vec_t & dF)
{
    // assemble global K matrix
    AssembleKA ();

    // set prescribed dF
    set_to_zero (dF);
    set_to_zero (W);
    for (NodBCs_t::const_iterator p=Dom.pF.begin(); p!=Dom.pF.end(); ++p)
    {
        for (IntDbl_t::const_iterator q=p->second.first.begin(); q!=p->second.first.end(); ++q)
        {
            size_t idof = q->first;
            long   eq   = p->first->EQ[idof];
            if (!pU[eq]) // set dF for unknown variables only
            {
                dF(eq) = dT*q->second;
                W (eq) = dT*q->second; // set W1 equal to dF1
            }
        }
    }

    // set prescribed dU
    set_to_zero (dU);
    for (NodBCs_t::const_iterator p=Dom.pU.begin(); p!=Dom.pU.end(); ++p)
    {
        for (IntDbl_t::const_iterator q=p->second.first.begin(); q!=p->second.first.end(); ++q)
        {
            size_t idof = q->first;
            long   eq   = p->first->EQ[idof];
            W(eq) = dT*q->second; // set W2 equal to dU2
        }
    }

    // correct W
    Sparse::SubMult (K12, W, W); // W1 -= K12*dU2

    // calc dU and dF
    if (FEM::Domain::PARA) MUMPS  ::Solve (A11, W, dU); // dU = inv(A11)*W
    else                   UMFPACK::Solve (A11, W, dU); // dU = inv(A11)*W

    // calc dF2
    Sparse::AddMult (K21, dU, dF); // dF2 += K21*dU1
    Sparse::AddMult (K22, dU, dF); // dF2 += K22*dU2

#ifdef HAS_MPI
    if (FEM::Domain::PARA)
    {
        // join dF
        //std::cout << "Proc # " << MPI::COMM_WORLD.Get_rank() << " BEFORE: dF = " << PrintVector(dF);
        MPI::COMM_WORLD.Allreduce (dF.data, TmpVec.data, NEq, MPI::DOUBLE, MPI::SUM);
        dF = TmpVec;
        //std::cout << "Proc # " << MPI::COMM_WORLD.Get_rank() << " AFTER:  dF = " << PrintVector(dF);
    }
    //else std::cout << "dF = " << PrintVector(dF);
#endif
}

inline void Solver::UpdateElements (Vec_t const & dU, bool CalcFint)
{
    if (CalcFint)
    {
        if (FEM::Domain::PARA)
        {
#ifdef HAS_MPI
            Vec_t dFint(NEq), dFint_tmp(NEq);
            set_to_zero (dFint);
            for (size_t i=0; i<ActEles.Size(); ++i) ActEles[i]->UpdateState (dU, &dFint);
            MPI::COMM_WORLD.Allreduce (dFint.data, dFint_tmp.data, NEq, MPI::DOUBLE, MPI::SUM);
            F_int += dFint_tmp;
#endif
        }
        else
        {
            for (size_t i=0; i<ActEles.Size(); ++i) ActEles[i]->UpdateState (dU, &F_int);
        }
    }
    else
    {
        for (size_t i=0; i<ActEles.Size(); ++i) ActEles[i]->UpdateState (dU);
    }
}

inline void Solver::Initialize (bool Transient)
{
    // info
    Util::Stopwatch stopwatch(/*only_root*/FEM::Domain::PARA);
    if (Root) printf("\n%s--- Solver --- initializing --------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    if (FEM::Domain::PARA)
    {
#ifdef HAS_MPI
        int my_id  = MPI::COMM_WORLD.Get_rank();
        int nprocs = MPI::COMM_WORLD.Get_size();

        // compute equation numbers corresponding to local DOFs of elements
        NEq = 0;
        for (size_t i=0; i<Dom.Eles.Size(); ++i)
        {
            if (Dom.Eles[i]->Active) Dom.Eles[i]->IncNLocDOF (NEq);
        }

        // compute equation numbers
        for (size_t i=0; i<Dom.Nods.Size(); ++i)
        {
            if (Dom.Nods[i]->NShares>0)
            {
                for (size_t j=0; j<Dom.Nods[i]->nDOF(); ++j)
                {
                    // only the domain with smallest ID will set EQ number
                    int min_part_id = Dom.Nods[i]->Vert.PartIDs[Dom.Nods[i]->Vert.PartIDs.Min()];
                    if (min_part_id==my_id) NEq++;
                }
            }
        }

        // communicate the number of DOFs numbered
        Array<int> send_alloc_dofs(nprocs); // number of dofs just allocated in this processor
        Array<int> recv_alloc_dofs(nprocs); // number of dofs just allocated in this processor
        send_alloc_dofs.SetValues (0);
        send_alloc_dofs[my_id] = NEq;
        MPI::COMM_WORLD.Scan (send_alloc_dofs.GetPtr(), recv_alloc_dofs.GetPtr(), nprocs, MPI::INT, MPI::SUM);

        // correct equation numbers
        NEq = 0;
        if (my_id>0)
        {
            for (int i=0; i<my_id; ++i) NEq += recv_alloc_dofs[i];
        }

        // assign equation numbers corresponding to local DOFs of elements
        for (size_t i=0; i<Dom.Eles.Size(); ++i)
        {
            if (Dom.Eles[i]->Active) Dom.Eles[i]->IncNLocDOF (NEq);
        }

        // assign equation numbers and set active nodes
        ActNods.Resize(0);
        for (size_t i=0; i<Dom.Nods.Size(); ++i)
        {
            if (Dom.Nods[i]->NShares>0)
            {
                for (size_t j=0; j<Dom.Nods[i]->nDOF(); ++j)
                {
                    int min_part_id = Dom.Nods[i]->Vert.PartIDs[Dom.Nods[i]->Vert.PartIDs.Min()];
                    if (min_part_id==my_id) // only the domain with smallest ID will set EQ number
                    {
                        Dom.Nods[i]->EQ[j] = NEq;
                        NEq++;
                    }
                    else Dom.Nods[i]->EQ[j] = -1;
                }
                ActNods.Push (Dom.Nods[i]);
            }
        }

        // post messages
        const int TAG_SENT_EQ = 1000;
        for (int i=my_id+1; i<nprocs; ++i)
        {
            Array<int> inter_eq; // equation of interface DOFs
            for (size_t j=0; j<Dom.InterNodes.Size(); ++j)
            {
                int  min_part_id       =  Dom.InterNodes[j]->Vert.PartIDs[Dom.InterNodes[j]->Vert.PartIDs.Min()]; // the smallest proc is the one supposed to send always
                bool do_send_to_proc_i = (Dom.InterNodes[j]->Vert.PartIDs.Find(i)>=0); // found processor on the interface and has higher id than me
                if (my_id==min_part_id && do_send_to_proc_i)
                {
                    for (size_t k=0; k<Dom.InterNodes[j]->nDOF(); ++k)
                        inter_eq.Push (Dom.InterNodes[j]->EQ[k]);
                }
            }
            MPI::Request req_send = MPI::COMM_WORLD.Isend (inter_eq.GetPtr(), inter_eq.Size(), MPI::INT, i, TAG_SENT_EQ);
        }

        // receive messages
        MPI::Status status;
        for (int i=0; i<my_id; ++i)
        {
            MPI::COMM_WORLD.Probe (MPI::ANY_SOURCE, TAG_SENT_EQ, status);
            int source = status.Get_source();
            int count  = status.Get_count(MPI::INT);
            Array<int> inter_eq(count); // equation of interface DOFs
            MPI::COMM_WORLD.Recv (inter_eq.GetPtr(), count, MPI::INT, source, TAG_SENT_EQ);

            int m = 0;
            for (size_t j=0; j<Dom.InterNodes.Size(); ++j)
            {
                int min_part_id = Dom.InterNodes[j]->Vert.PartIDs[Dom.InterNodes[j]->Vert.PartIDs.Min()]; // the smallest proc is the one supposed to send always
                if (source==min_part_id)
                {
                    for (size_t k=0; k<Dom.InterNodes[j]->nDOF(); ++k)
                    {
                        Dom.InterNodes[j]->EQ[k] = inter_eq[m];
                        m++;
                    }
                }
            }
        }
        
        // set NEq in all procs  
        if (my_id==nprocs-1) // I'm the last processor and the only one who knows the total num equations
        {
            NEq = 0;
            for (int i=0; i<nprocs; ++i) NEq += recv_alloc_dofs[i];
        }
        MPI::COMM_WORLD.Bcast (&NEq, 1, MPI::INT, /*from*/nprocs-1); // last processor will broadcast to everyone else
#else
        throw new Fatal("Solver::Initialize: Parallel code is not available (not compiled)");
#endif
    }
    else
    {
        // assign equation numbers corresponding to local DOFs of elements
        NEq = 0;
        for (size_t i=0; i<Dom.Eles.Size(); ++i)
        {
            if (Dom.Eles[i]->Active) Dom.Eles[i]->IncNLocDOF (NEq);
        }

        // assign equation numbers and set active nodes
        ActNods.Resize(0);
        for (size_t i=0; i<Dom.Nods.Size(); ++i)
        {
            if (Dom.Nods[i]->NShares>0)
            {
                for (size_t j=0; j<Dom.Nods[i]->nDOF(); ++j)
                {
                    Dom.Nods[i]->EQ[j] = NEq;
                    NEq++;
                }
                ActNods.Push (Dom.Nods[i]);
            }
        }
    }

    // prescribed equations and prescribed U
    pEQ.Resize    (0);
    pU .Resize    (NEq);
    pU .SetValues (false);
    for (NodBCs_t::const_iterator p=Dom.pU.begin(); p!=Dom.pU.end(); ++p)
    {
        for (IntDbl_t::const_iterator q=p->second.first.begin(); q!=p->second.first.end(); ++q)
        {
            size_t idof = q->first;
            long   eq   = p->first->EQ[idof];
            pU[eq] = true;
            pEQ.Push (eq);
        }
    }

    // number of Lagrange multipliers
    size_t nlag_pins = 0;//Dom.Msh.Pins.size()*Dom.NDim; // number of Lag mult due to pins
    size_t nlag_insu = Dom.InclSupport.size();       // number of Lag mult due to inclined supports
    NLag = nlag_pins + nlag_insu;                    // total number of Lagrange multipliers
    NEq += NLag;

    // number of extra non-zero values due to Lagrange multipliers
    size_t nzlag = nlag_pins*2*Dom.NDim + nlag_insu*2*Dom.NDim;

    // find total number of non-zero entries, including duplicates, and assign active elements
    size_t K11_size = 0;
    size_t K12_size = 0;
    size_t K21_size = 0;
    size_t K22_size = 0;
    ActEles.Resize(0);
    for (size_t k=0; k<Dom.Eles.Size(); ++k)
    {
        if (Dom.Eles[k]->Active)
        {
            Array<size_t> loc; // location array
            Dom.Eles[k]->GetLoc (loc);
            for (size_t i=0; i<loc.Size(); ++i)
            for (size_t j=0; j<loc.Size(); ++j)
            {
                     if (!pU[loc[i]] && !pU[loc[j]]) K11_size++;
                else if (!pU[loc[i]] &&  pU[loc[j]]) K12_size++;
                else if ( pU[loc[i]] && !pU[loc[j]]) K21_size++;
                else if ( pU[loc[i]] &&  pU[loc[j]]) K22_size++;
            }
            ActEles.Push (Dom.Eles[k]);
        }
    }

    // allocate triplets
    A11.AllocSpace (NEq,NEq,K11_size+pEQ.Size()+nzlag); // augmented
    K12.AllocSpace (NEq,NEq,K12_size);
    K21.AllocSpace (NEq,NEq,K21_size);
    K22.AllocSpace (NEq,NEq,K22_size);
    if (Transient)
    {
        K11.AllocSpace (NEq,NEq,K11_size);
        M11.AllocSpace (NEq,NEq,K11_size);
        M12.AllocSpace (NEq,NEq,K12_size);
        M21.AllocSpace (NEq,NEq,K21_size);
        M22.AllocSpace (NEq,NEq,K22_size);
        if (DampTy!=None_t)
        {
            C11.AllocSpace (NEq,NEq,K11_size);
            C12.AllocSpace (NEq,NEq,K12_size);
            C21.AllocSpace (NEq,NEq,K21_size);
            C22.AllocSpace (NEq,NEq,K22_size);
        }
    }

    // initialize variables
    R     .change_dim (NEq);  set_to_zero (R);
    F     .change_dim (NEq);  set_to_zero (F);
    F_int .change_dim (NEq);  set_to_zero (F_int);
    W     .change_dim (NEq);  set_to_zero (W);
    U     .change_dim (NEq);  set_to_zero (U);
    dU    .change_dim (NEq);  set_to_zero (dU);
    dF    .change_dim (NEq);  set_to_zero (dF);
    TmpVec.change_dim (NEq);  set_to_zero (TmpVec);
    if (Transient)
    {
        V .change_dim (NEq);  set_to_zero (V);
        A .change_dim (NEq);  set_to_zero (A);
        Us.change_dim (NEq);  set_to_zero (Us);
        Vs.change_dim (NEq);  set_to_zero (Vs);
    }

    // set variables
    for (size_t i=0; i<ActNods.Size(); ++i)
    {
        for (size_t j=0; j<ActNods[i]->nDOF(); ++j)
        {
            long  eq  = ActNods[i]->EQ[j];
            U    (eq) = ActNods[i]->U [j];
            F    (eq) = ActNods[i]->F [j];
            F_int(eq) = ActNods[i]->F [j];
        }
    }

    // calc residual
    _cal_resid ();

    // info
    if (Root)
    {
        printf("%s  Num of DOFs (NEq)  = %d%s\n", TERM_CLR2, NEq, TERM_RST);
        printf("%s  Num of non-zeros   = %d%s\n", TERM_CLR2, K11_size+pEQ.Size()+nzlag, TERM_RST);
    }
}

inline void Solver::SetScheme (char const * StrScheme)
{
    if      (strcmp(StrScheme,"FE")==0) Scheme = FE_t;
    else if (strcmp(StrScheme,"ME")==0) Scheme = ME_t;
    else if (strcmp(StrScheme,"NR")==0) Scheme = NR_t;
    else throw new Fatal("Solver::SetScheme: Key '%s' is invalid. The following keys are availabe: 'FE', 'ME', 'NR'",StrScheme);
}

inline void Solver::SetIncsW (size_t NInc, bool NonLinWei)
{
    IncsW.Resize (NInc);
    if (NonLinWei)
    {
        double delx = 1.0/NInc;
        for (size_t i=0; i<NInc; ++i)
        {
            double xi = 1.0-delx*i;
            double xj = 1.0-delx*(i+1);
            double yi = pow(xi,3.0);
            double yj = pow(xj,3.0);
            IncsW[i] = yi-yj;
        }
    }
    else
    {
        for (size_t i=0; i<NInc; ++i) IncsW[i] = 1.0/NInc;
    }
}

inline bool Solver::ResidOK () const
{
    if (MaxNormF<1.0e-8) // particular case when MaxNormF is very small
    {
        return (NormR<TolR);
    }
    else return (NormR<TolR*MaxNormF);
}

inline void Solver::_set_A_Lag ()
{
    // set equations corresponding to Lagrange multipliers
    long eqlag = NEq - NLag;

    // pins
    /*
    if (Dom.Msh.Pins.size()>0)
    {
        for (Mesh::Pin_t::const_iterator p=Dom.Msh.Pins.begin(); p!=Dom.Msh.Pins.end(); ++p)
        {
            Node const & nod0 = (*ActNods[p->first->ID]);
            for (size_t i=0; i<p->second.Size(); ++i)
            {
                Node const & nod1 = (*ActNods[p->second[i]->ID]);
                for (int j=0; j<Dom.NDim; ++j)
                {
                    long eq0 = nod0.EQ[nod0.UMap(Dom.DisplKeys[j])];
                    long eq1 = nod1.EQ[nod1.UMap(Dom.DisplKeys[j])];
                    A11.PushEntry (eq0,eqlag,1.0);   A11.PushEntry (eq1,eqlag,-1.0);
                    A11.PushEntry (eqlag,eq0,1.0);   A11.PushEntry (eqlag,eq1,-1.0);
                    eqlag++;
                }
            }
        }
    }
    */

    // inclined supports
    if (Dom.InclSupport.size()>0)
    {
        for (InclSupport_t::const_iterator p=Dom.InclSupport.begin(); p!=Dom.InclSupport.end(); ++p)
        {
            Node const & nod = (*p->first);
            double s = sin(p->second); // sin(alpha)
            double c = cos(p->second); // cos(alpha)
            long eq0 = nod.EQ[nod.UMap(Dom.DisplKeys[0])]; // ~ ux
            long eq1 = nod.EQ[nod.UMap(Dom.DisplKeys[1])]; // ~ uy
            A11.PushEntry (eqlag, eq0,  s);
            A11.PushEntry (eqlag, eq1, -c);
            A11.PushEntry (eq0, eqlag,  s);
            A11.PushEntry (eq1, eqlag, -c);
            eqlag++;
        }
    }
}

inline void Solver::_cal_resid (bool WithAccel)
{
    // calculate residual
    R = F - F_int;

    //std::cout << "\n######################################   Before\n";
    //std::cout << "F     = " << PrintVector(F,     "%10.3f");
    //std::cout << "F_int = " << PrintVector(F_int, "%10.3f");
    //std::cout << "R     = " << PrintVector(R,     "%10.3f");
    
    // number of the first equation corresponding to Lagrange multipliers
    long eqlag = NEq - NLag;

    // clear forces due to pins
    /*
    if (Dom.Msh.Pins.size()>0)
    {
        for (Mesh::Pin_t::const_iterator p=Dom.Msh.Pins.begin(); p!=Dom.Msh.Pins.end(); ++p)
        {
            Node const & nod0 = (*ActNods[p->first->ID]);
            for (size_t i=0; i<p->second.Size(); ++i)
            {
                Node const & nod1 = (*ActNods[p->second[i]->ID]);
                for (int j=0; j<Dom.NDim; ++j)
                {
                    long eq0 = nod0.EQ[nod0.UMap(Dom.DisplKeys[j])];
                    long eq1 = nod1.EQ[nod1.UMap(Dom.DisplKeys[j])];
                    //F(eq0) += -U(eqlag);
                    //F(eq1) +=  U(eqlag);
                    R(eq0) = 0.0;
                    R(eq1) = 0.0;
                    eqlag++;
                }
            }
        }
    }
    */

    // clear forces due to inclined supports
    if (Dom.InclSupport.size()>0)
    {
        for (InclSupport_t::const_iterator p=Dom.InclSupport.begin(); p!=Dom.InclSupport.end(); ++p)
        {
            Node const & nod = (*p->first);
            //double s = sin(p->second); // sin(alpha)
            //double c = cos(p->second); // cos(alpha)
            long eq0 = nod.EQ[nod.UMap(Dom.DisplKeys[0])]; // ~ ux
            long eq1 = nod.EQ[nod.UMap(Dom.DisplKeys[1])]; // ~ uy
            //F(eq0) += -U(eqlag)*s;
            //F(eq1) +=  U(eqlag)*c;
            R(eq0) = 0.0;
            R(eq1) = 0.0;

            //std::cout << Util::_12_4<< -U(eqlag)*s << "  " << Util::_12_4<< U(eqlag)*c << std::endl;

            eqlag++;
        }
    }

    // clear forces due to supports
    //for (size_t i=0; i<pEQ.Size(); ++i) R(pEQ[i]) = 0.0;

    /*
    std::cout << "\n######################################   After\n";
    std::cout << "F     = " << PrintVector(F,     "%10.3f");
    std::cout << "F_int = " << PrintVector(F_int, "%10.3f");
    std::cout << std::endl;
    */

    if (WithAccel)
    {
        Sparse::SubMult (M11, A, R); // R -= M11*A
        if (DampTy!=None_t) Sparse::SubMult (C11, V, R); // R -= C11*V
        //std::cout << "V  = " << PrintVector(V);
        //std::cout << "A  = " << PrintVector(A);
        //std::cout << "R  = " << PrintVector(R);
    }
    NormR    = Norm(R);
    MaxNormF = Util::Max (Norm(F), Norm(F_int));
}

inline void Solver::_cor_resid (Vec_t & dU, Vec_t & dF)
{
    if (!CorR) return;

    // iterations
    size_t it = 0;
    while (!ResidOK() && it<MaxIt)
    {
        // debug
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);

        // assemble global K matrix
        if (!ModNR) AssembleKA ();

        if (FEM::Domain::PARA)
        {
#ifdef HAS_MPI
            // set workspace: R
            set_to_zero (dF);
            for (size_t i=0; i<pEQ.Size(); ++i)
            {
                dF(pEQ[i]) = -R(pEQ[i]); // dF2 = -R2
            }

            if (FEM::Domain::PARA && MPI::COMM_WORLD.Get_rank()>0) R = 1.0;

            // set workspace: R
            for (size_t i=0; i<pEQ.Size(); ++i)
            {
                R (pEQ[i]) = 0.0;        // R2  = 0
            }

            for (size_t i=0; i<Dom.InterNodes.Size(); ++i)
            {
                for (size_t j=0; j<Dom.InterNodes[i]->nDOF(); ++j)
                {
                    long   eq   = Dom.InterNodes[i]->EQ[j];
                    size_t nint = Dom.InterNodes[i]->Vert.PartIDs.Size();
                    dF(eq) /= static_cast<double>(nint);
                }
            }

            // solve
            MUMPS::Solve (A11, R,  dU, /*Prod*/true); // dU1 = inv(A11)*R1

            // calc dF
            Sparse::AddMult (K21, dU, dF); // dF2 += K21*dU1  =>  dF2 = K21*dU1 - R2

            // join dF
            MPI::COMM_WORLD.Allreduce (dF.data, TmpVec.data, NEq, MPI::DOUBLE, MPI::SUM);
            dF = TmpVec;
#endif
        }
        else
        {
            // calc corrector dU
            set_to_zero (dF);
            for (size_t i=0; i<pEQ.Size(); ++i)
            {
                dF(pEQ[i]) = -R(pEQ[i]); // dF2 = -R2
                R (pEQ[i]) = 0.0;        // R2  = 0
            }

            // solve
            UMFPACK::Solve (A11, R,  dU); // dU1 = inv(A11)*R1

            // calc dF
            Sparse::AddMult (K21, dU, dF); // dF2 += K21*dU1  =>  dF2 = K21*dU1 - R2
        }

        //std::cout << "Norm(dU) = " << Util::_8s << Norm(dU)    << std::endl;
        //std::cout << "Norm(R)  = " << Util::_8s << Norm(R)     << std::endl;
        //std::cout << "Norm(F)  = " << Util::_8s << Norm(F)     << std::endl;
        //std::cout << "Norm(Fi) = " << Util::_8s << Norm(F_int) << std::endl;

        // update
        UpdateElements (dU, /*CalcFint*/true);
        U += dU;
        F += dF;

        // residual
        _cal_resid ();

        // next iteration
        it++;
    }
    if (it>=MaxIt) throw new Fatal("Solver::_cor_resid: Residual correction did not converge after %d iterations.\n\t(%e>=%e)  (NormR>=TolR*MaxNormF).\n\tNormR         = %g\n\tMaxNormF      = %g\n\tTolR          = %e\n\tTolR*MaxNormF = %e",it,NormR,TolR*MaxNormF,NormR,MaxNormF,TolR,TolR*MaxNormF);
    if (it>It) It = it;
}

inline void Solver::_FE_update (double tf)
{
    double dt = (tf-Time)/nSS;
    for (Stp=0; Stp<nSS; ++Stp)
    {
        // calculate tangent increments
        TgIncs (dt, dU, dF);

        // update elements
        UpdateElements (dU, /*CalcFint*/true);

        // update U, F, and Time
        U    += dU;
        F    += dF;
        Time += dt;

        // debug
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
    }

    // residual
    _cal_resid ();
}

inline void Solver::_ME_update (double tf)
{
    // auxiliar vectors
    Vec_t dU_fe(NEq), dU_tm(NEq), dU_me(NEq), U_me(NEq), U_dif(NEq);
    Vec_t dF_fe(NEq), dF_tm(NEq), dF_me(NEq), F_me(NEq), F_dif(NEq);

    // for each pseudo time T
    double T   = 0.0;
    double dT  = (dTlast>0.0 ? dTlast : dTini);
    double Dt  = tf-Time;
    for (Stp=0; Stp<MaxSS; ++Stp)
    {
        // exit point
        if (T>=1.0) break;

        // backup state of elements
        for (size_t i=0; i<ActEles.Size(); ++i) ActEles[i]->BackupState ();

        // time increment
        double dt = Dt*dT;

        // FE state
        TgIncs (dt, dU_fe, dF_fe);
        UpdateElements (dU_fe, /*CalcFint*/false);

        // ME state
        TgIncs (dt, dU_tm, dF_tm);
        dU_me = 0.5*(dU_fe + dU_tm);
        dF_me = 0.5*(dF_fe + dF_tm);
        U_me  = U + dU_me;
        F_me  = F + dF_me;

        // local error
        U_dif = 0.5*(dU_tm - dU_fe);
        F_dif = 0.5*(dF_tm - dF_fe);
        for (size_t i=NEq-NLag; i<NEq; ++i) { U_dif(i)=0.0; F_dif(i)=0.0; } // ignore equations corresponding to Lagrange multipliers
        double U_err = Norm(U_dif)/(1.0+Norm(U_me));
        double F_err = Norm(F_dif)/(1.0+Norm(F_me));
        double error = U_err + F_err;

        // step multiplier
        double m = (error>0.0 ? 0.9*sqrt(STOL/error) : mMax);

        // restore state of elements
        for (size_t i=0; i<ActEles.Size(); ++i) ActEles[i]->RestoreState ();

        // update
        if (error<STOL)
        {
            UpdateElements (dU_me, /*CalcFint*/true);
            T    += dT;
            U     = U_me;
            F     = F_me;
            Time += dt;
            if (m>mMax) m = mMax;
            if (SSOut || (DbgFun!=NULL))
            {
                // update nodes
                for (size_t i=0; i<ActNods.Size(); ++i)
                {
                    for (size_t j=0; j<ActNods[i]->nDOF(); ++j)
                    {
                        long eq = ActNods[i]->EQ[j];
                        ActNods[i]->U[j] = U(eq);
                        ActNods[i]->F[j] = F(eq);
                    }
                }
            }
            if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
            if (SSOut) Dom.OutResults (Time, F_int);
        }
        else if (m<mMin) m = mMin;

        // change next step size
        dT = m * dT;

        // check for last increment
        if (dT>1.0-T) dT = 1.0-T;
        else dTlast = dT;
    }
    if (Stp>=MaxSS) throw new Fatal("Solver:_ME_update: Modified-Euler (global integration) did not converge for %d steps",Stp);

    // correct residual
    _cal_resid ();
    _cor_resid (dU_me, dF_me);
}

inline void Solver::_NR_update (double tf)
{
    It = 0;
    double dt = (tf-Time)/nSS;
    for (Stp=0; Stp<nSS; ++Stp)
    {
        // calculate tangent increments
        TgIncs (dt, dU, dF);

        // update elements
        UpdateElements (dU, /*CalcFint*/true);

        // update U, F, and Time
        U    += dU;
        F    += dF;
        Time += dt;

        // residual
        _cal_resid ();
        _cor_resid (dU, dF);

        // debug
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
    }
}

inline void Solver::_GN22_update (double tf, double dt)
{
    // constants
    const double c1 = dt*dt*(1.0-DynTh2)/2.0;
    const double c2 = dt*(1.0-DynTh1);
    const double c3 = 2.0/(dt*dt*DynTh2);
    const double c4 = 2.0*DynTh1/(dt*DynTh2);

    while (Time<tf)
    {
        // predictor
        Us = U + dt*V + c1*A;
        Vs = V + c2*A;
        A  = c3*(U - Us);
        V  = Vs + (DynTh1*dt)*A;

        // new F
        for (NodBCs_t::const_iterator p=Dom.pF.begin(); p!=Dom.pF.end(); ++p)
        {
            for (IntDbl_t::const_iterator q=p->second.first.begin(); q!=p->second.first.end(); ++q)
            {
                size_t idof = q->first;
                long   eq   = p->first->EQ[idof];
                if (!pU[eq]) // set dF for unknown variables only
                {
                    double dFdt = (*p->second.second)(Time+DynTh1*dt);
                    F(eq) += dFdt*dt;
                    //F(eq) = (*p->second.second)(Time+DynTh1*dt);
                }
            }
        }
        double normF = Norm(F);

        // iterations
        for (It=0; It<MaxIt; ++It)
        {
            // new F and residual
            R = F - F_int;
            //std::cout << "F    = " << PrintVector(F,     "%8.2f");
            //std::cout << "Fint = " << PrintVector(F_int, "%8.2f");
            //std::cout << "R    = " << PrintVector(R,     "%8.2f");
            //std::cout << std::endl;

            // assemble Amat
            if (DampTy==None_t) AssembleKMA  (c3,     1.0);  // A = c3*M        + K
            else                AssembleKCMA (c3, c4, 1.0);  // A = c3*M + c4*C + K

            // solve for dU
            Sparse::SubMult (M11, A, R);  if (DampTy!=None_t)  // R -= M11*A
            Sparse::SubMult (C11, V, R);                       // R -= C11*V
            UMFPACK::Solve  (A11, R, dU);                      // dU = inv(A11)*R

            // update elements
            UpdateElements (dU, /*CalcFint*/true);
            for (size_t i=0; i<pEQ.Size(); ++i) F_int(pEQ[i]) = 0.0; // clear internal forces related to supports

            // update state
            U += dU;
            A  = c3*(U - Us);
            V  = Vs + (DynTh1*dt)*A;

            // check convergence
            NormR    = Norm(R);
            MaxNormF = Util::Max (normF, Norm(F_int));
            if (ResidOK()) break;
        }
        if (It>=MaxIt) throw new Fatal("Solver::_GN22_update: Generalized-Newmark (GN22) did not converge after %d iterations",It);

        // next time step
        Time += dt;

        // debug
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
    }
}

inline void Solver::_time_print (char const * Comment)
{
    if (Comment!=NULL)
    {
        printf("\n%s--- Stage solution --- %s -----------------------------------------%s\n",TERM_CLR1,Comment,TERM_RST);
        printf("%s%10s  %12s  %4s  %4s%s\n",TERM_CLR2,"Time","Norm(R)","NSS","NIT",TERM_RST);
    }
    printf("%10.6f  %s%8e%s  %4zd  %4zd\n",Time,(ResidOK()?TERM_GREEN:TERM_RED),NormR,TERM_RST,Stp,It);
}

}; // namespace FEM

#endif // MECHSYS_FEM_SOLVER_H
