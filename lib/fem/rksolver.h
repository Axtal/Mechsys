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
    enum Damping_t { None_t, Rayleigh_t, HMCoup_t }; ///< Damping type: none, Rayleigh type (C=alp*M+bet*K), HydroMechCoupling

    // Constructor
    RKSolver (Domain & dom);

    // Methods
    void SteadySolve (int NInc=1,                         char const * FKey=NULL);
    void TransSolve  (double tf, double dt, double dtOut, char const * FKey=NULL);
    void DynSolve    (double tf, double dt, double dtOut, char const * FKey=NULL);

    // Data (read/write)
    bool          WithInfo; ///< Print information ?
    bool          LinProb;  ///< Linear problems => constant matrices
    Damping_t     DampTy;   ///< Damping type
    double        DampAm;   ///< Rayleigh damping Am coefficient (C = Am*M + Ak*K)
    double        DampAk;   ///< Rayleigh damping Ak coefficient (C = Am*M + Ak*K)
    String        RKScheme; ///< Runge-Kutta scheme
    double        STOL;     ///< Error tolerance for RK scheme

    // Data (read-only)
    Domain      & Dom;      ///< Domain
    size_t        NEq;      ///< Total number of equations (DOFs)
    size_t        NIv;      ///< Total number of internal variables of elements
    size_t        NLag;     ///< Number of Lagrange multipliers
    size_t        NnzLag;   ///< number of extra non-zero values due to Lagrange multipliers
    Array<int>    pEQ;      ///< prescribed equations;  size = NEq2
    Array<bool>   pU;       ///< prescribed U;  size = NEq (all equations)

    // Auxiliary data
    Array<int>    Eq1_to_V1;
    Array<int>    Eq1_to_U1;
    Vec_t         A,V,U,F;

    // Triplets
    size_t K11_size, K12_size, K21_size, K22_size; ///< Size of triplets
    Sparse::Triplet<double,int> K11,K12,K21,K22;   ///< Stiffness matrices
    Sparse::Triplet<double,int> C11,C12,C21,C22;   ///< Damping matrices
    Sparse::Triplet<double,int> M11,M12,M21,M22;   ///< Mass matrices

    // Internal methods
    int  SteadyFunc     (double t, double const Y[], double dYdt[]); ///< Steady problems: callback for RK solver
    int  TransFunc      (double t, double const Y[], double dYdt[]); ///< Transient problems: callback for RK solver
    int  DynFunc        (double t, double const Y[], double dYdt[]); ///< Dynamic problems: callback for RK solver
    void Initialize     ();                                          ///< Allocate memory
    void AssembleK      ();                                          ///< Assemble K
    void AssembleKM     ();                                          ///< Assemble K and M
    void AssembleKCM    ();                                          ///< Assemble K, C, and M
    void AugMatrix      (Sparse::Triplet<double,int> & A11);         ///< Augment A matrix with BCs and Lag multipliers
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation ///////


inline RKSolver::RKSolver (Domain & dom)
    : WithInfo (true),
      LinProb  (false),
      DampTy   (None_t),
      DampAm   (0.01),
      DampAk   (0.01),
      RKScheme ("RK4I"),
      STOL     (1.0e-3),
      Dom      (dom)
{
}

inline void RKSolver::SteadySolve (int NInc, char const * FKey)
{
    // initialize 
    Initialize ();

    // allocate triplets and vectors
    V  .change_dim (NEq);
    U  .change_dim (NEq);
    F  .change_dim (NEq);
    K11.AllocSpace (NEq,NEq,K11_size+pEQ.Size()+NnzLag); // augmented
    K12.AllocSpace (NEq,NEq,K12_size);
    K21.AllocSpace (NEq,NEq,K21_size);
    K22.AllocSpace (NEq,NEq,K22_size);

    // allocate ode solver
    Numerical::ODESolver<RKSolver> ode(this, &RKSolver::SteadyFunc, NEq-pEQ.Size(), RKScheme.CStr(), STOL, 1.0/NInc);

    // set initial values
    ode.t    = 0.0; // pseudo-time (0 <= ode.t <= 1)
    int rkeq = 0;
    for (size_t i=0; i<Dom.ActNods.Size();     ++i)
    for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
    {
        int eq = Dom.ActNods[i]->Eq(j);
        if (!Dom.ActNods[i]->pU(j))
        {
            Eq1_to_U1[eq]        = rkeq++;
            ode.Y[Eq1_to_U1[eq]] = Dom.ActNods[i]->U(j); // U1
        }
        // fill vectors with initial values
        U(eq) = Dom.ActNods[i]->U(j);
        F(eq) = Dom.ActNods[i]->F(j);
    }

    // first output
    printf("\n%s--- Stage solution --- Steady ---------------------------------%s\n",TERM_CLR1,TERM_RST);
    printf("%s%10s%s\n",TERM_CLR2,"Time",TERM_RST);

    // solve
    double Time0 = Dom.Time;
    double tout  = 1.0/NInc;
    while (ode.t<1.0)
    {
        // evolve
        ode.Evolve (tout);
        Dom.Time = Time0 + ode.t;

        // collect results into U,F and update nodes
        for (size_t i=0; i<Dom.ActNods.Size();     ++i)
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            int eq = Dom.ActNods[i]->Eq(j);
            if (!Dom.ActNods[i]->pU(j))
            {
                // fill vectors
                U(eq)  = ode.Y[Eq1_to_U1[eq]];                // U1
                F(eq) += ode.t*Dom.ActNods[i]->PFOrZero(j,0); // F1
                // set node
                Dom.ActNods[i]->U(j) = U(eq); // U1
                Dom.ActNods[i]->F(j) = F(eq); // F1
            }
            else
            {
                // fill vectors
                U(eq) += ode.t*Dom.ActNods[i]->PUIdxDOF(j,0); // U2
                F(eq)  = 0.0;                                 // F2
                // set node
                Dom.ActNods[i]->U(j) = U(eq); // U2
            }
        }

        // calculate F2
        AssembleK       ();
        Sparse::AddMult (K21, U, F); // F2 += K21*U1
        Sparse::AddMult (K22, U, F); // F2 += K22*U2

        // update F2 in nodes
        for (size_t i=0; i<Dom.NodsWithPU.Size();     ++i)
        for (size_t j=0; j<Dom.NodsWithPU[i]->NDOF(); ++j)
        {
            if (Dom.NodsWithPU[i]->pU(j))
            {
                int eq = Dom.NodsWithPU[i]->Eq(j);
                Dom.NodsWithPU[i]->F(j) = F(eq); // F2
            }
        }

        // output
        if (WithInfo) printf("%10.6f\n",Dom.Time);
        Dom.OutResults (FKey);
        tout += 1.0/NInc;
    }
}

inline int RKSolver::SteadyFunc (double t, double const Y[], double dYdt[])
{
    for (size_t i=0; i<Dom.ActNods.Size();     ++i)
    for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
    {
        int eq = Dom.ActNods[i]->Eq(j);
        if (!Dom.ActNods[i]->pU(j))
        {
            U(eq) = Y[Eq1_to_U1[eq]];              // U1
            F(eq) = Dom.ActNods[i]->PFOrZero(j,0); // DF1
            V(eq) = 0.0;                           // dU1dT
        }
        else U(eq) = Dom.ActNods[i]->PUIdxDOF(j,0); // DU2
    }

    AssembleK       ();
    Sparse::SubMult (K12, U, F);  // DF1 -= K12*DU2
    UMFPACK::Solve  (K11, F, V);  // V1   = inv(K11)*DF1 == dU1dT

    for (size_t i=0; i<Dom.ActNods.Size();     ++i)
    for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
    {
        if (!Dom.ActNods[i]->pU(j))
        {
            int eq = Dom.ActNods[i]->Eq(j);
            dYdt[Eq1_to_U1[eq]] = V(eq); // dU1dT
        }
    }

    return GSL_SUCCESS;
}

inline void RKSolver::TransSolve (double tf, double dt, double dtOut, char const * FKey)
{
    // initialize 
    Initialize ();

    // allocate triplets and vectors
    V  .change_dim (NEq);
    U  .change_dim (NEq);
    F  .change_dim (NEq);
    K11.AllocSpace (NEq,NEq,K11_size);
    K12.AllocSpace (NEq,NEq,K12_size);
    K21.AllocSpace (NEq,NEq,K21_size);
    K22.AllocSpace (NEq,NEq,K22_size);
    M11.AllocSpace (NEq,NEq,K11_size+pEQ.Size()+NnzLag); // augmented
    M12.AllocSpace (NEq,NEq,K12_size);
    M21.AllocSpace (NEq,NEq,K21_size);
    M22.AllocSpace (NEq,NEq,K22_size);

    // allocate ode solver
    Numerical::ODESolver<RKSolver> ode(this, &RKSolver::TransFunc, NEq-pEQ.Size(), RKScheme.CStr(), STOL, dt);

    // set initial values
    ode.t    = Dom.Time;
    int rkeq = 0;
    for (size_t i=0; i<Dom.ActNods.Size();     ++i)
    for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
    {
        if (!Dom.ActNods[i]->pU(j))
        {
            int eq = Dom.ActNods[i]->Eq(j);
            Eq1_to_U1[eq]        = rkeq++;
            ode.Y[Eq1_to_U1[eq]] = Dom.ActNods[i]->U(j); // U1
        }
    }

    // first output
    printf("\n%s--- Stage solution --- Transient ------------------------------%s\n",TERM_CLR1,TERM_RST);
    printf("%s%10s%s\n",TERM_CLR2,"Time",TERM_RST);

    // solve
    double tout = Dom.Time + dtOut;
    while (ode.t<tf)
    {
        // evolve
        ode.Evolve (tout);
        Dom.Time = ode.t;

        // collect results into V,U,F and update nodes
        for (size_t i=0; i<Dom.ActNods.Size();     ++i)
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            int eq = Dom.ActNods[i]->Eq(j);
            if (!Dom.ActNods[i]->pU(j))
            {
                // fill vectors
                U(eq) = ode.Y[Eq1_to_U1[eq]];              // U1
                F(eq) = Dom.ActNods[i]->PFOrZero(j,ode.t); // F1
                // set node
                Dom.ActNods[i]->U(j) = U(eq); // U1
                Dom.ActNods[i]->F(j) = F(eq); // F1
            }
            else
            {
                // fill vectors
                V(eq) = Dom.ActNods[i]->PVIdxDOF(j,ode.t); // V2
                U(eq) = Dom.ActNods[i]->PUIdxDOF(j,ode.t); // U2
                F(eq) = 0.0;                               // F2
                // set node
                Dom.ActNods[i]->U(j) = U(eq); // U2
            }
        }

        // calculate F2 TODO: should calculate V1 here
        AssembleKM      ();
        Sparse::AddMult (M21, V, F);  // F2 += M21*V1
        Sparse::AddMult (M22, V, F);  // F2 += M22*V2
        Sparse::AddMult (K21, U, F);  // F2 += K21*U1
        Sparse::AddMult (K22, U, F);  // F2 += K22*U2

        // update F2 in nodes
        for (size_t i=0; i<Dom.NodsWithPU.Size();     ++i)
        for (size_t j=0; j<Dom.NodsWithPU[i]->NDOF(); ++j)
        {
            if (Dom.NodsWithPU[i]->pU(j))
            {
                int eq = Dom.NodsWithPU[i]->Eq(j);
                Dom.NodsWithPU[i]->F(j) = F(eq); // F2
            }
        }

        // output
        if (WithInfo) printf("%10.6f\n",Dom.Time);
        Dom.OutResults (FKey);
        tout += dtOut;
    }
}

inline int RKSolver::TransFunc (double t, double const Y[], double dYdt[])
{
    for (size_t i=0; i<Dom.ActNods.Size();     ++i)
    for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
    {
        int eq = Dom.ActNods[i]->Eq(j);
        if (!Dom.ActNods[i]->pU(j))
        {
            U(eq) = Y[Eq1_to_U1[eq]];              // U1
            F(eq) = Dom.ActNods[i]->PFOrZero(j,t); // F1
        }
        else
        {
            V(eq) = Dom.ActNods[i]->PVIdxDOF(j,t); // V2
            U(eq) = Dom.ActNods[i]->PUIdxDOF(j,t); // U2
        }
    }

    AssembleKM      ();
    Sparse::SubMult (M12, V, F);  // F1 -= M12*V2
    Sparse::SubMult (K11, U, F);  // F1 -= K11*U1
    Sparse::SubMult (K12, U, F);  // F1 -= K12*U2
    UMFPACK::Solve  (M11, F, V);  // V1  = inv(M11)*F1

    for (size_t i=0; i<Dom.ActNods.Size();     ++i)
    for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
    {
        if (!Dom.ActNods[i]->pU(j))
        {
            int eq = Dom.ActNods[i]->Eq(j);
            dYdt[Eq1_to_U1[eq]] = V(eq); // V1
        }
    }

    return GSL_SUCCESS;
}

inline void RKSolver::DynSolve (double tf, double dt, double dtOut, char const * FKey)
{
    // initialize 
    Initialize ();

    // allocate triplets and vectors
    A  .change_dim (NEq);
    V  .change_dim (NEq);
    U  .change_dim (NEq);
    F  .change_dim (NEq);
    K11.AllocSpace (NEq,NEq,K11_size);
    K12.AllocSpace (NEq,NEq,K12_size);
    K21.AllocSpace (NEq,NEq,K21_size);
    K22.AllocSpace (NEq,NEq,K22_size);
    M11.AllocSpace (NEq,NEq,K11_size+pEQ.Size()+NnzLag); // augmented
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

    // allocate ode solver
    Numerical::ODESolver<RKSolver> ode(this, &RKSolver::DynFunc, 2*(NEq-pEQ.Size()), RKScheme.CStr(), STOL, dt);

    // set initial values
    ode.t    = Dom.Time;
    int rkeq = 0;
    for (size_t i=0; i<Dom.ActNods.Size();     ++i)
    for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
    {
        if (!Dom.ActNods[i]->pU(j))
        {
            int eq = Dom.ActNods[i]->Eq(j);
            Eq1_to_V1[eq]        = rkeq++;
            Eq1_to_U1[eq]        = rkeq++;
            ode.Y[Eq1_to_V1[eq]] = Dom.ActNods[i]->V(j); // V1
            ode.Y[Eq1_to_U1[eq]] = Dom.ActNods[i]->U(j); // U1
        }
    }

    // first output
    printf("\n%s--- Stage solution --- Dynamic --------------------------------%s\n",TERM_CLR1,TERM_RST);
    printf("%s%10s%s\n",TERM_CLR2,"Time",TERM_RST);

    // solve
    double tout = Dom.Time + dtOut;
    while (Dom.Time<tf)
    {
        // evolve
        ode.Evolve (tout);
        Dom.Time = ode.t;

        // collect results into A,V,U,F and update nodes
        for (size_t i=0; i<Dom.ActNods.Size();     ++i)
        for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
        {
            int eq = Dom.ActNods[i]->Eq(j);
            if (!Dom.ActNods[i]->pU(j))
            {
                // fill vectors
                V(eq) = ode.Y[Eq1_to_V1[eq]];              // V1
                U(eq) = ode.Y[Eq1_to_U1[eq]];              // U1
                F(eq) = Dom.ActNods[i]->PFOrZero(j,ode.t); // F1
                // set node
                Dom.ActNods[i]->V(j) = V(eq); // V1
                Dom.ActNods[i]->U(j) = U(eq); // U1
                Dom.ActNods[i]->F(j) = F(eq); // F1
            }
            else
            {
                // fill vectors
                A(eq) = Dom.ActNods[i]->PAIdxDOF(j,ode.t); // A2
                V(eq) = Dom.ActNods[i]->PVIdxDOF(j,ode.t); // V2
                U(eq) = Dom.ActNods[i]->PUIdxDOF(j,ode.t); // U2
                F(eq) = 0.0;                               // F2
                // set node
                Dom.ActNods[i]->V(j) = V(eq); // V2
                Dom.ActNods[i]->U(j) = U(eq); // U2
            }
        }

        // calculate F2 TODO: should calculate A1 here
        AssembleKCM     ();
        Sparse::AddMult (M21, A, F);                       // F2 += M21*A1
        Sparse::AddMult (M22, A, F); if (DampTy!=None_t) { // F2 += M22*A2
        Sparse::AddMult (C21, V, F);                       // F2 += C21*V1
        Sparse::AddMult (C22, V, F); }                     // F2 += C22*V2
        Sparse::AddMult (K21, U, F);                       // F2 += K21*U1
        Sparse::AddMult (K22, U, F);                       // F2 += K22*U2

        // update F2 in nodes
        for (size_t i=0; i<Dom.NodsWithPU.Size();     ++i)
        for (size_t j=0; j<Dom.NodsWithPU[i]->NDOF(); ++j)
        {
            if (Dom.NodsWithPU[i]->pU(j))
            {
                int eq = Dom.NodsWithPU[i]->Eq(j);
                Dom.NodsWithPU[i]->F(j) = F(eq); // F2
            }
        }

        // output
        if (WithInfo) printf("%10.6f\n",Dom.Time);
        Dom.OutResults (FKey);
        tout += dtOut;
    }
}

inline int RKSolver::DynFunc (double t, double const Y[], double dYdt[])
{
    for (size_t i=0; i<Dom.ActNods.Size();     ++i)
    for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
    {
        int eq = Dom.ActNods[i]->Eq(j);
        if (!Dom.ActNods[i]->pU(j))
        {
            V(eq) = Y[Eq1_to_V1[eq]];              // V1
            U(eq) = Y[Eq1_to_U1[eq]];              // U1
            F(eq) = Dom.ActNods[i]->PFOrZero(j,t); // F1
        }
        else
        {
            A(eq) = Dom.ActNods[i]->PAIdxDOF(j,t); // A2
            V(eq) = Dom.ActNods[i]->PVIdxDOF(j,t); // V2
            U(eq) = Dom.ActNods[i]->PUIdxDOF(j,t); // U2
        }
    }

    AssembleKCM     ();
    Sparse::SubMult (M12, A, F); if (DampTy!=None_t) { // F1 -= M12*A2
    Sparse::SubMult (C11, V, F);                       // F1 -= C11*V1
    Sparse::SubMult (C12, V, F); }                     // F1 -= C12*V2
    Sparse::SubMult (K11, U, F);                       // F1 -= K11*U1
    Sparse::SubMult (K12, U, F);                       // F1 -= K12*U2
    UMFPACK::Solve  (M11, F, A);                       // A1  = inv(M11)*F1

    for (size_t i=0; i<Dom.ActNods.Size();     ++i)
    for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
    {
        if (!Dom.ActNods[i]->pU(j))
        {
            int eq = Dom.ActNods[i]->Eq(j);
            dYdt[Eq1_to_V1[eq]] = A(eq); // A1
            dYdt[Eq1_to_U1[eq]] = V(eq); // V1
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
    for (size_t i=0; i<Dom.ActNods.Size();     ++i)
    for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
    {
        Dom.ActNods[i]->Eq(j) = NEq;
        NEq++;
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
    NnzLag = nlag_pins*2*Dom.NDim + nlag_insu*2*Dom.NDim;

    // find total number of non-zero entries, including duplicates
    NIv      = 0; // total number of internal variables of elements
    K11_size = 0;
    K12_size = 0;
    K21_size = 0;
    K22_size = 0;
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

    // RK maps
    Eq1_to_V1.Resize (NEq);   Eq1_to_V1.SetValues (-1);
    Eq1_to_U1.Resize (NEq);   Eq1_to_U1.SetValues (-1);

    // info
    if (WithInfo)
    {
        printf("%s  Num of DOFs (NEq)  = %zd%s\n", TERM_CLR2, NEq, TERM_RST);
        printf("%s  Num of non-zeros   = %zd%s\n", TERM_CLR2, K11_size+pEQ.Size()+NnzLag, TERM_RST);
    }
}

inline void RKSolver::AssembleK ()
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
        for (size_t j=0; j<loc.Size(); ++j)
        {
            if      (!pU[loc[i]] && !pU[loc[j]]) { K11.PushEntry(loc[i], loc[j], K(i,j)); }
            else if (!pU[loc[i]] &&  pU[loc[j]]) { K12.PushEntry(loc[i], loc[j], K(i,j)); }
            else if ( pU[loc[i]] && !pU[loc[j]]) { K21.PushEntry(loc[i], loc[j], K(i,j)); }
            else if ( pU[loc[i]] &&  pU[loc[j]]) { K22.PushEntry(loc[i], loc[j], K(i,j)); }
        }
    }
    AugMatrix (K11);
}

inline void RKSolver::AssembleKM ()
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
        for (size_t j=0; j<loc.Size(); ++j)
        {
            if      (!pU[loc[i]] && !pU[loc[j]]) { K11.PushEntry(loc[i], loc[j], K(i,j));  M11.PushEntry(loc[i], loc[j], M(i,j)); }
            else if (!pU[loc[i]] &&  pU[loc[j]]) { K12.PushEntry(loc[i], loc[j], K(i,j));  M12.PushEntry(loc[i], loc[j], M(i,j)); }
            else if ( pU[loc[i]] && !pU[loc[j]]) { K21.PushEntry(loc[i], loc[j], K(i,j));  M21.PushEntry(loc[i], loc[j], M(i,j)); }
            else if ( pU[loc[i]] &&  pU[loc[j]]) { K22.PushEntry(loc[i], loc[j], K(i,j));  M22.PushEntry(loc[i], loc[j], M(i,j)); }
        }
    }
    AugMatrix (M11);
}

inline void RKSolver::AssembleKCM ()
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
        for (size_t j=0; j<loc.Size(); ++j)
        {
            if      (!pU[loc[i]] && !pU[loc[j]]) { K11.PushEntry(loc[i], loc[j], K(i,j));  M11.PushEntry(loc[i], loc[j], M(i,j));  if (DampTy!=None_t) { C11.PushEntry(loc[i], loc[j], C(i,j)); } }
            else if (!pU[loc[i]] &&  pU[loc[j]]) { K12.PushEntry(loc[i], loc[j], K(i,j));  M12.PushEntry(loc[i], loc[j], M(i,j));  if (DampTy!=None_t) { C12.PushEntry(loc[i], loc[j], C(i,j)); } }
            else if ( pU[loc[i]] && !pU[loc[j]]) { K21.PushEntry(loc[i], loc[j], K(i,j));  M21.PushEntry(loc[i], loc[j], M(i,j));  if (DampTy!=None_t) { C21.PushEntry(loc[i], loc[j], C(i,j)); } }
            else if ( pU[loc[i]] &&  pU[loc[j]]) { K22.PushEntry(loc[i], loc[j], K(i,j));  M22.PushEntry(loc[i], loc[j], M(i,j));  if (DampTy!=None_t) { C22.PushEntry(loc[i], loc[j], C(i,j)); } }
        }
    }
    AugMatrix (M11);
}

inline void RKSolver::AugMatrix (Sparse::Triplet<double,int> & A11)
{
    // augment A11 matrix and set Lagrange multipliers (if any)
    for (size_t i=0; i<pEQ.Size(); ++i) A11.PushEntry (pEQ[i],pEQ[i], 1.0);

    // set equations corresponding to Lagrange multipliers
    int eqlag = NEq - NLag;

    // pins and inclined supports
    for (size_t i=0; i<Dom.NodsPins  .Size(); ++i) Dom.NodsPins  [i]->SetLagPin    (eqlag, A11);
    for (size_t i=0; i<Dom.NodsIncSup.Size(); ++i) Dom.NodsIncSup[i]->SetLagIncSup (eqlag, A11);
}


}; // namespace FEM

#endif // MECHSYS_FEM_RKSOLVER_H
