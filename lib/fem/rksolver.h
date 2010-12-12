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

#ifndef MECHSYS_FEM_RKSOLVER_H
#define MECHSYS_FEM_RKSOLVER_H

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
#include <mechsys/numerical/odesolver.h>

namespace FEM
{

class RKSolver
{
public:
    // enum
    enum EqType_t  { Steady_t, Trans_t, Dyn_t };     ///< Equation type (internal)
    enum Damping_t { None_t, Rayleigh_t, HMCoup_t }; ///< Damping type: none, Rayleigh type (C=alp*M+bet*K), HydroMechCoupling

    // Constructor
    RKSolver (Domain & dom);

    // Methods
    void SteadySolve (int NInc,                           char const * FKey=NULL) { EqType=Steady_t;  Solve(Dom.Time+1.0, 1.0/NInc, 1.0/NInc, FKey); }
    void TransSolve  (double tf, double dt, double dtOut, char const * FKey=NULL) { EqType=Trans_t;   Solve(tf, dt, dtOut, FKey); }
    void DynSolve    (double tf, double dt, double dtOut, char const * FKey=NULL) { EqType=Dyn_t;     Solve(tf, dt, dtOut, FKey); }

    // Data (read/write)
    bool          WithInfo; ///< Print information ?
    bool          LinProb;  ///< Linear problems => constant matrices
    Damping_t     DampTy;   ///< Damping type
    double        DampAm;   ///< Rayleigh damping Am coefficient (C = Am*M + Ak*K)
    double        DampAk;   ///< Rayleigh damping Ak coefficient (C = Am*M + Ak*K)
    String        RKScheme; ///< Runge-Kutta scheme
    double        STOL;     ///< Error tolerance for RK scheme
    SDPair      * Steps;    ///< Non-linear steps

    // Data (read-only)
    Domain      & Dom;      ///< Domain
    size_t        NEq;      ///< Total number of equations (DOFs)
    size_t        NIv;      ///< Total number of internal variables of elements
    size_t        NLag;     ///< Number of Lagrange multipliers
    Array<int>    pEQ;      ///< prescribed equations;  size = NEq2
    Array<bool>   pU;       ///< prescribed U;  size = NEq (all equations)
    EqType_t      EqType;   ///< Equation type

    // Auxiliary data
    Array<int>    Eq1_to_V1;
    Array<int>    Eq1_to_U1;
    Array<int>    Eq2_to_F2;
    Vec_t         A,V,U,F;

    // Triplets
    Sparse::Triplet<double,int> K11,K12,K21,K22; ///< Stiffness matrices
    Sparse::Triplet<double,int> C11,C12,C21,C22; ///< Damping matrices
    Sparse::Triplet<double,int> M11,M12,M21,M22; ///< Mass matrices

    // Internal methods
    void Solve      (double tf, double dt, double dtOut, char const * FKey);
    int  SteadyFunc (double t, double const Y[], double dYdt[]);
    int  TransFunc  (double t, double const Y[], double dYdt[]);
    int  DynFunc    (double t, double const Y[], double dYdt[]);
    void Initialize ();
    void Assemble   ();

    // Auxiliary methods
    static double Timestep (int i, int a=7, double L=100.0, int sch=0, double m=2.0)
    { return (i<a ? pow(2.0,i)/L : (sch==0 ? pow(2.0,i-a) : pow(2.0,a)/L+m*(i-a))); }
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation ///////


inline RKSolver::RKSolver (Domain & dom)
    : WithInfo (true),
      LinProb  (false),
      DampTy   (None_t),
      DampAm   (0.01),
      DampAk   (0.01),
      RKScheme ("RK4I"),
      STOL     (1.0e-5),
      Dom      (dom),
      EqType   (Steady_t)
{
}

inline void RKSolver::Solve (double tf, double dt, double dtOut, char const * FKey)
{
    // initialize 
    Initialize ();
    int neq2 = pEQ.Size(); // number of prescribed equations (2)
    int neq1 = NEq - neq2; // number of unknown equations (1)

    // pointer to function and number of equations
    typedef int (RKSolver::*pFun) (double t, double const Y[], double dYdt[]);
    pFun pfun;
    int  neq;
    if (EqType==Steady_t)
    {
        pfun = &RKSolver::SteadyFunc;
        neq  = NEq;
    }
    else if (EqType==Trans_t)
    {
        pfun = &RKSolver::TransFunc;
        neq  = neq1;
    }
    else if (EqType==Dyn_t)
    {
        pfun = &RKSolver::DynFunc;
        neq  = 2*neq1;
    }

    // allocate ode solver
    Numerical::ODESolver<RKSolver> ode(this, pfun, neq, RKScheme.CStr(), STOL, dt);

    // set initial values
    ode.t    = Dom.Time;
    int rkeq = 0;
    if (EqType==Steady_t)
    {
        for (size_t i=0; i<Dom.ActNods.Size();     ++i)
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            int eq = Dom.ActNods[i]->Eq(j);
            if (!pU[eq])
            {
                Eq1_to_U1[eq]        = rkeq++;
                ode.Y[Eq1_to_U1[eq]] = 0.0;//Dom.ActNods[i]->U(j); // U1
            }
            else
            {
                Eq2_to_F2[eq]        = rkeq++;
                ode.Y[Eq2_to_F2[eq]] = 0.0;//Dom.ActNods[i]->F(j); // F2
            }
        }
    }
    else if (EqType==Trans_t)
    {
        for (size_t i=0; i<Dom.ActNods.Size();     ++i)
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            int eq = Dom.ActNods[i]->Eq(j);
            if (!pU[eq])
            {
                Eq1_to_U1[eq]        = rkeq++;
                ode.Y[Eq1_to_U1[eq]] = 0.0;//Dom.ActNods[i]->U(j); // U1
            }
        }
    }
    else if (EqType==Dyn_t)
    {
        for (size_t i=0; i<Dom.ActNods.Size();     ++i)
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            int eq = Dom.ActNods[i]->Eq(j);
            if (!pU[eq])
            {
                Eq1_to_V1[eq]        = rkeq++;
                Eq1_to_U1[eq]        = rkeq++;
                ode.Y[Eq1_to_V1[eq]] = 0.0;//Dom.ActNods[i]->V(j); // V1
                ode.Y[Eq1_to_U1[eq]] = 0.0;//Dom.ActNods[i]->U(j); // U1
            }
        }
    }

    // first output
    printf("\n%s--- Stage solution --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    printf("%s%10s%s\n",TERM_CLR2,"Time",TERM_RST);

    // solve
    double tout = Dom.Time + dtOut;
    while (Dom.Time<tf)
    {
        // evolve from Time to tout
        ode.Evolve (tout);
        Dom.Time = ode.t;

        // output
        //IdxOut++;
        if (WithInfo) printf("%10.6f\n",Dom.Time);
        //Dom.OutResults (IdxOut, Dom.Time, FileKey);
        tout += dtOut;
    }
}

inline int RKSolver::SteadyFunc (double t, double const Y[], double dYdt[])
{
    for (size_t i=0; i<Dom.ActNods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.ActNods[i]->NPF(); ++j)
        {
            int eq1 = Dom.ActNods[i]->EqPF (j);    // eq1
            F(eq1)  = Dom.ActNods[i]->PF   (j, t); // F1
        }
        for (size_t j=0; j<Dom.ActNods[i]->NPU(); ++j)
        {
            int eq2 = Dom.ActNods[i]->EqPU (j);    // eq2
            U(eq2)  = Dom.ActNods[i]->PU   (j, t); // U2
        }
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            int eq = Dom.ActNods[i]->Eq(j);
            if (!pU[eq]) U(eq) = Y[Eq1_to_U1[eq]]; // U1
            else
            {
                F(eq) = Y[Eq2_to_F2[eq]]; // F2
                V(eq) = 0.0;              // dF2dT
            }
        }
    }

    Assemble        ();
    Sparse::SubMult (K12, U, F);  // F1 -= K12*U2
    UMFPACK::Solve  (K11, F, V);  // V1  = inv(K11)*F1 == dU1dT
    Sparse::AddMult (K21, V, V);  // V2 += K21*V1
    Sparse::AddMult (K22, U, V);  // V2 += K22*U1 == dF2dT

    for (size_t i=0; i<Dom.ActNods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            int eq = Dom.ActNods[i]->Eq(j);
            if (!pU[eq]) dYdt[Eq1_to_U1[eq]] = V(eq); // dU1dT
            else         dYdt[Eq2_to_F2[eq]] = V(eq); // dF2dT
        }
    }

    return GSL_SUCCESS;
}

inline int RKSolver::TransFunc (double t, double const Y[], double dYdt[])
{
    for (size_t i=0; i<Dom.ActNods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.ActNods[i]->NPF(); ++j)
        {
            int eq1 = Dom.ActNods[i]->EqPF (j);    // eq1
            F(eq1)  = Dom.ActNods[i]->PF   (j, t); // F1
        }
        for (size_t j=0; j<Dom.ActNods[i]->NPU(); ++j)
        {
            int eq2 = Dom.ActNods[i]->EqPU (j);    // eq2
            V(eq2)  = Dom.ActNods[i]->PU   (j, t); // V2
            U(eq2)  = Dom.ActNods[i]->PU   (j, t); // U2
        }
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            int eq = Dom.ActNods[i]->Eq(j);
            if (!pU[eq]) U(eq) = Y[Eq1_to_U1[eq]]; // U1
        }
    }

    Assemble        ();
    Sparse::SubMult (M12, V, F);  // F1 -= M12*V2
    Sparse::SubMult (K11, U, F);  // F1 -= K11*U1
    Sparse::SubMult (K12, U, F);  // F1 -= K12*U2
    UMFPACK::Solve  (M11, F, V);  // V1  = inv(M11)*F1

    for (size_t i=0; i<Dom.ActNods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            int eq = Dom.ActNods[i]->Eq(j);
            if (!pU[eq]) dYdt[Eq1_to_U1[eq]] = V(eq); // V1
        }
    }

    return GSL_SUCCESS;
}

inline int RKSolver::DynFunc (double t, double const Y[], double dYdt[])
{
    for (size_t i=0; i<Dom.ActNods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.ActNods[i]->NPF(); ++j)
        {
            int eq1 = Dom.ActNods[i]->EqPF (j);    // eq1
            F(eq1)  = Dom.ActNods[i]->PF   (j, t); // F1
        }
        for (size_t j=0; j<Dom.ActNods[i]->NPU(); ++j)
        {
            int eq2 = Dom.ActNods[i]->EqPU (j);    // eq2
            A(eq2)  = Dom.ActNods[i]->PU   (j, t); // A2
            V(eq2)  = Dom.ActNods[i]->PU   (j, t); // V2
            U(eq2)  = Dom.ActNods[i]->PU   (j, t); // U2
        }
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            int eq = Dom.ActNods[i]->Eq(j);
            if (!pU[eq])
            {
                V(eq)  = Y[Eq1_to_V1[eq]]; // V1
                U(eq)  = Y[Eq1_to_U1[eq]]; // U1
            }
        }
    }

    Assemble        ();
    Sparse::SubMult (M12, A, F); if (DampTy!=None_t) { // F1 -= M12*A2
    Sparse::SubMult (C11, V, F);                       // F1 -= C11*V1
    Sparse::SubMult (C12, V, F); }                     // F1 -= C12*V2
    Sparse::SubMult (K11, U, F);                       // F1 -= K11*U1
    Sparse::SubMult (K12, U, F);                       // F1 -= K12*U2
    UMFPACK::Solve  (M11, F, A);                       // A1  = inv(M11)*F1

    for (size_t i=0; i<Dom.ActNods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            int eq = Dom.ActNods[i]->Eq(j);
            if (!pU[eq])
            {
                dYdt[Eq1_to_V1[eq]] = A(eq); // A1
                dYdt[Eq1_to_U1[eq]] = V(eq); // V1
            }
        }
    }

    return GSL_SUCCESS;
}

inline void RKSolver::Initialize ()
{
    // assign equation numbers corresponding to local DOFs of elements
    NEq = 0;
    //for (size_t i=0; i<Dom.ActEles.Size(); ++i) Dom.ActEles[i]->IncNLocDOF (NEq);

    // assign equation numbers
    for (size_t i=0; i<Dom.ActNods.Size(); ++i)
    {
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            Dom.ActNods[i]->Eq(j) = NEq;
            NEq++;
        }
    }

    // prescribed equations and prescribed U
    pEQ.Resize    (0);
    pU .Resize    (NEq);
    pU .SetValues (false);
    for (size_t i=0; i<Dom.NodsWithPU.Size(); ++i)
    {
        Node const * nod = Dom.NodsWithPU[i];
        for (size_t j=0; j<nod->NPU(); ++j)
        {
            int eq = nod->EqPU(j);
            pU[eq] = true;
            pEQ.Push (eq);
        }
    }

    // number of Lagrange multipliers
    size_t nlag_pins = Dom.NodsPins.Size()*Dom.NDim; // number of Lag mult due to pins
    size_t nlag_insu = Dom.NodsIncSup.Size();        // number of Lag mult due to inclined supports
    NLag = nlag_pins + nlag_insu;                    // total number of Lagrange multipliers
    NEq += NLag;

    // number of extra non-zero values due to Lagrange multipliers
    size_t nzlag = nlag_pins*2*Dom.NDim + nlag_insu*2*Dom.NDim;

    // find total number of non-zero entries, including duplicates
    NIv = 0; // total number of internal variables of elements
    size_t K11_size = 0;
    size_t K12_size = 0;
    size_t K21_size = 0;
    size_t K22_size = 0;
    for (size_t k=0; k<Dom.ActEles.Size(); ++k)
    {
        Array<size_t> loc; // location array
        Dom.ActEles[k]->GetLoc (loc);
        for (size_t i=0; i<loc.Size(); ++i)
        for (size_t j=0; j<loc.Size(); ++j)
        {
                 if (!pU[loc[i]] && !pU[loc[j]]) K11_size++;
            else if (!pU[loc[i]] &&  pU[loc[j]]) K12_size++;
            else if ( pU[loc[i]] && !pU[loc[j]]) K21_size++;
            else if ( pU[loc[i]] &&  pU[loc[j]]) K22_size++;
        }
        NIv += Dom.ActEles[k]->NIVs();
    }

    // allocate triplets and vectors
    U  .change_dim (NEq);
    F  .change_dim (NEq);
    V  .change_dim (NEq);
    K12.AllocSpace (NEq,NEq,K12_size);
    K21.AllocSpace (NEq,NEq,K21_size);
    K22.AllocSpace (NEq,NEq,K22_size);
    if (EqType==Steady_t)
    {
        K11.AllocSpace (NEq,NEq,K11_size+pEQ.Size()+nzlag); // augmented
    }
    else
    {
        K11.AllocSpace (NEq,NEq,K11_size);
        M11.AllocSpace (NEq,NEq,K11_size+pEQ.Size()+nzlag); // augmented
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
        if (EqType==Dyn_t) A.change_dim (NEq);
    }

    // RK maps
    Eq1_to_V1.Resize (NEq);   Eq1_to_V1.SetValues (-1);
    Eq1_to_U1.Resize (NEq);   Eq1_to_U1.SetValues (-1);
    Eq2_to_F2.Resize (NEq);   Eq2_to_F2.SetValues (-1);

    // info
    if (WithInfo)
    {
        printf("%s  Num of DOFs (NEq)  = %zd%s\n", TERM_CLR2, NEq, TERM_RST);
        printf("%s  Num of non-zeros   = %zd%s\n", TERM_CLR2, K11_size+pEQ.Size()+nzlag, TERM_RST);
    }
}

inline void RKSolver::Assemble ()
{
    Sparse::Triplet<double,int> * A11;
    if (EqType==Steady_t)
    {
        if (LinProb && K11.Top()>0) return;
        K11.ResetTop();
        K12.ResetTop();
        K21.ResetTop();
        K22.ResetTop();
        for (size_t k=0; k<Dom.ActEles.Size(); ++k)
        {
            Mat_t         K;
            Array<size_t> loc;
            Dom.ActEles[k]->CalcK  (K);
            Dom.ActEles[k]->GetLoc (loc);
            for (size_t i=0; i<loc.Size(); ++i)
            {
                for (size_t j=0; j<loc.Size(); ++j)
                {
                    if      (!pU[loc[i]] && !pU[loc[j]]) { K11.PushEntry(loc[i], loc[j], K(i,j)); }
                    else if (!pU[loc[i]] &&  pU[loc[j]]) { K12.PushEntry(loc[i], loc[j], K(i,j)); }
                    else if ( pU[loc[i]] && !pU[loc[j]]) { K21.PushEntry(loc[i], loc[j], K(i,j)); }
                    else if ( pU[loc[i]] &&  pU[loc[j]]) { K22.PushEntry(loc[i], loc[j], K(i,j)); }
                }
            }
        }
        A11 = &K11;
    }
    else if (EqType==Trans_t)
    {
        if (LinProb && M11.Top()>0) return;
        K11.ResetTop();   M11.ResetTop();
        K12.ResetTop();   M12.ResetTop();
        K21.ResetTop();   M21.ResetTop();
        K22.ResetTop();   M22.ResetTop();
        for (size_t k=0; k<Dom.ActEles.Size(); ++k)
        {
            Mat_t         K, M;
            Array<size_t> loc;
            Dom.ActEles[k]->CalcK  (K);
            Dom.ActEles[k]->CalcM  (M);
            Dom.ActEles[k]->GetLoc (loc);
            for (size_t i=0; i<loc.Size(); ++i)
            {
                for (size_t j=0; j<loc.Size(); ++j)
                {
                    if      (!pU[loc[i]] && !pU[loc[j]]) { K11.PushEntry(loc[i], loc[j], K(i,j));  M11.PushEntry(loc[i], loc[j], M(i,j)); }
                    else if (!pU[loc[i]] &&  pU[loc[j]]) { K12.PushEntry(loc[i], loc[j], K(i,j));  M12.PushEntry(loc[i], loc[j], M(i,j)); }
                    else if ( pU[loc[i]] && !pU[loc[j]]) { K21.PushEntry(loc[i], loc[j], K(i,j));  M21.PushEntry(loc[i], loc[j], M(i,j)); }
                    else if ( pU[loc[i]] &&  pU[loc[j]]) { K22.PushEntry(loc[i], loc[j], K(i,j));  M22.PushEntry(loc[i], loc[j], M(i,j)); }
                }
            }
        }
        A11 = &M11;
    }
    else if (EqType==Dyn_t)
    {
        if (LinProb && M11.Top()>0) return;
        K11.ResetTop();   M11.ResetTop();
        K12.ResetTop();   M12.ResetTop();
        K21.ResetTop();   M21.ResetTop();
        K22.ResetTop();   M22.ResetTop();
        if (DampTy!=None_t)
        {
            C11.ResetTop();
            C12.ResetTop();
            C21.ResetTop();
            C22.ResetTop();
        }
        for (size_t k=0; k<Dom.ActEles.Size(); ++k)
        {
            Mat_t         K, C, M;
            Array<size_t> loc;
            Dom.ActEles[k]->CalcK  (K);
            Dom.ActEles[k]->CalcM  (M);
            Dom.ActEles[k]->GetLoc (loc);
            if (DampTy==Rayleigh_t) C = DampAm*M + DampAk*K;
            for (size_t i=0; i<loc.Size(); ++i)
            {
                for (size_t j=0; j<loc.Size(); ++j)
                {
                    if      (!pU[loc[i]] && !pU[loc[j]]) { K11.PushEntry(loc[i], loc[j], K(i,j));  M11.PushEntry(loc[i], loc[j], M(i,j));  if (DampTy!=None_t) { C11.PushEntry(loc[i], loc[j], C(i,j)); } }
                    else if (!pU[loc[i]] &&  pU[loc[j]]) { K12.PushEntry(loc[i], loc[j], K(i,j));  M12.PushEntry(loc[i], loc[j], M(i,j));  if (DampTy!=None_t) { C12.PushEntry(loc[i], loc[j], C(i,j)); } }
                    else if ( pU[loc[i]] && !pU[loc[j]]) { K21.PushEntry(loc[i], loc[j], K(i,j));  M21.PushEntry(loc[i], loc[j], M(i,j));  if (DampTy!=None_t) { C21.PushEntry(loc[i], loc[j], C(i,j)); } }
                    else if ( pU[loc[i]] &&  pU[loc[j]]) { K22.PushEntry(loc[i], loc[j], K(i,j));  M22.PushEntry(loc[i], loc[j], M(i,j));  if (DampTy!=None_t) { C22.PushEntry(loc[i], loc[j], C(i,j)); } }
                }
            }
        }
        A11 = &M11;
    }

    // augment A11 matrix and set Lagrange multipliers (if any)
    for (size_t i=0; i<pEQ.Size(); ++i) A11->PushEntry (pEQ[i],pEQ[i], 1.0);

    // set equations corresponding to Lagrange multipliers
    int eqlag = NEq - NLag;

    // pins and inclined supports
    for (size_t i=0; i<Dom.NodsPins  .Size(); ++i) Dom.NodsPins  [i]->SetLagPin    (eqlag, (*A11));
    for (size_t i=0; i<Dom.NodsIncSup.Size(); ++i) Dom.NodsIncSup[i]->SetLagIncSup (eqlag, (*A11));
}


}; // namespace FEM

#endif // MECHSYS_FEM_RKSOLVER_H
