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
#include <mechsys/numerical/odesolver.h>

namespace FEM
{

class Solver
{
public:
    // enum
    enum Scheme_t  { FE_t, ME_t, NR_t };             ///< Steady time integration scheme: Forward-Euler, Modified-Euler, Newton-Rhapson
    enum TScheme_t { TH_t };                         ///< Transient time integration scheme: Theta
    enum DScheme_t { GN22_t, RK_t };                 ///< Dynamic time integration scheme: (Single step/O2/2nd order), (Generalized Newmark/O2/2nd order)
    enum Damping_t { None_t, Rayleigh_t, HMCoup_t }; ///< Damping type: none, Rayleigh type (C=alp*M+bet*K), HydroMechCoupling

    // typedefs
    typedef void (*pOutFun) (Solver const & Sol, void * OutDat); ///< Pointer to output function

    // Constructor
    Solver (Domain & Dom, pOutFun OutFun=NULL, void * OutDat=NULL,
                          pOutFun DbgFun=NULL, void * DbgDat=NULL); ///< Allocate solver object

    // Methods
    void Solve          (size_t NInc=1, char const * FileKey=NULL);                      ///< Solve quasi-static problem
    void TransSolve     (double tf, double dt, double dtOut, char const * FileKey=NULL, SDPair * Steps=NULL); ///< Solve transient problem
    void DynSolve       (double tf, double dt, double dtOut, char const * FileKey=NULL, SDPair * Steps=NULL); ///< Solve dynamic problem
    void AssembleKA     ();                                                              ///< A = K11
    void AssembleKMA    (double Coef1, double Coef2);                                    ///< A = Coef1*M + Coef2*K
    void AssembleKCMA   (double Coef1, double Coef2, double Coef3);                      ///< A = Coef1*M + Coef2*C + Coef3*K
    void TgIncs         (double dT, Vec_t & dU, Vec_t & dF);                             ///< Tangent increments: dU = inv(K)*dF
    void UpdateNodes    ();                                                              ///< Copy U and F values into Nodes' U and F structures
    void UpdateElements (Vec_t const & dU, bool CalcFint);                               ///< Update elements
    void Initialize     (bool Transient=false);                                          ///< Initialize global matrices and vectors
    void SetScheme      (char const * StrScheme);                                        ///< Set solution scheme: 'FE', 'ME', 'NR'
    void SetTransScheme (char const * StrScheme);                                        ///< Set transient scheme: 'TH'
    void SetIncsW       (size_t NInc, bool NonLinWei=false);                             ///< Set weights for quasi-static problem (If Weights.Size()==0: generate weights)
    bool ResidOK        () const;                                                        ///< Check if the residual is OK

    // Data (read-only)
    Domain      & Dom;      ///< Domain
    pOutFun       OutFun;   ///< Output function (called during output)
    void        * OutDat;   ///< Debug data (to be used with either OutFun or DbgFun)
    pOutFun       DbgFun;   ///< Debug function (called everytime for some internal (update) methods)
    void        * DbgDat;   ///< Debug data (to be used with either OutFun or DbgFun)
    size_t        Inc;      ///< Current increment
    size_t        Stp;      ///< Current (sub) step
    size_t        It;       ///< Current iteration
    size_t        NEq;      ///< Total number of equations (DOFs)
    size_t        NIv;      ///< Total number of internal variables of elements
    size_t        NLag;     ///< Number of Lagrange multipliers
    Array<long>   pEQ;      ///< prescribed equations
    Array<long>   pEQproc;  ///< prescribed equations: processor which will handle the equation (in case of shared nodes)
    Array<bool>   pU;       ///< prescribed U
    double        NormR;    ///< Euclidian norm of residual (R)
    double        MaxNormF; ///< Max(Norm(F), Norm(Fint))

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
    double        TransTh;  ///< Transient scheme constant (theta)
    DScheme_t     DScheme;  ///< Dynamic scheme
    Damping_t     DampTy;   ///< Damping type
    double        DampAm;   ///< Rayleigh damping Am coefficient (C = Am*M + Ak*K)
    double        DampAk;   ///< Rayleigh damping Ak coefficient (C = Am*M + Ak*K)
    double        DynTh1;   ///< Dynamic coefficient Theta 1
    double        DynTh2;   ///< Dynamic coefficient Theta 2
    Array<double> IncsW;    ///< Increments weights used in Solve
    bool          WithInfo; ///< Print information ?
    bool          WarnRes;  ///< Warning if residual exceeds limit ?
    double        WrnTol;   ///< Warning tolerance
    String        RKScheme; ///< Runge-Kutta scheme
    double        RKSTOL;   ///< Runge-Kutta tolerance

    // Triplets and sparse matrices
    Sparse::Triplet<double,int> K11,K12,K21,K22; ///< Stiffness matrices
    Sparse::Triplet<double,int> C11,C12,C21,C22; ///< Damping matrices
    Sparse::Triplet<double,int> M11,M12,M21,M22; ///< Mass matrices
    Sparse::Triplet<double,int> A11;             ///< A=K  or  A=C1*M+C2*K  or  A=C1*M+C2*C+C3*K

    // Vectors
    Vec_t R;            ///< Residual
    Vec_t F0;           ///< External force at the beginning of the stage
    Vec_t F, F_int;     ///< External and internal forces
    Vec_t W, U;         ///< Workspace, displacement
    Vec_t V, A;         ///< (Transient/Dynamic) velocity and acceleration
    Vec_t dU, dF;       ///< Increments of U and F
    Vec_t Us, Vs;       ///< starred variables (for GN22)
    Vec_t TmpVec;       ///< Temporary vector (for parallel Allreduce)

    // Auxiliary methods
    static double Timestep (int i, int a=7, double L=100.0, int sch=0, double m=2.0)
    { return (i<a ? pow(2.0,i)/L : (sch==0 ? pow(2.0,i-a) : pow(2.0,a)/L+m*(i-a))); }

private:
    void   _aug_and_set_A ();                          ///< Augment A matrix and set Lagrange multipliers if any
    void   _wrn_resid     ();                          ///< Warning message for large residuals (debugging)
    void   _cal_resid     ();                          ///< Calculate residual
    void   _cor_resid     (Vec_t & dU, Vec_t & dF);    ///< Correct residual
    void   _FE_update     (double tf);                 ///< (Forward-Euler)  Update Time and elements to tf
    void   _ME_update     (double tf, char const*FNK); ///< (Modified-Euler) Update Time and elements to tf (FNK=FilenameKey)
    void   _NR_update     (double tf);                 ///< (Newton-Rhapson) Update Time and elements to tf
    void   _TH_update     (double tf, double dt);      ///< (theta) Transient: update Time and elements to tf
    void   _GN22_update   (double tf, double dt);      ///< (Generalized-Newmark) Update Time and elements to tf
    void   _time_print    (char const * Comment=NULL); ///< Print timestep data
    void   _VUIV_to_Y     (double Y[]);
    void   _Y_to_VUIV     (double const Y[]);
    int    _RK_func       (double t, double const Y[], double dYdt[]);
    void   _RK_update     (double tf, double dt);      ///< Runge-Kutta update
    void   _debug_print_matrices (bool Stop=true);     ///< Debug method
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Solver::Solver (Domain & TheDom, pOutFun TheOutFun, void * TheOutDat, pOutFun TheDbgFun, void * TheDbgDat)
    : Dom      (TheDom),
      OutFun   (TheOutFun),
      OutDat   (TheOutDat),
      DbgFun   (TheDbgFun),
      DbgDat   (TheDbgDat),
      Inc      (0),
      Stp      (0),
      It       (0),
      Scheme   (ME_t),
      CalcWork (false),
      nSS      (1),
      STOL     (1.0e-5),
      dTini    (1.0),
      dTlast   (-1.0),
      mMin     (0.1),
      mMax     (10.0),
      MaxSS    (2000),
      SSOut    (false),
      CteTg    (false),
      ModNR    (false),
      TolR     (1.0e-7),
      CorR     (true),
      MaxIt    (20),
      TScheme  (TH_t),
      TransTh  (1./2.),
      DScheme  (GN22_t),
      DampTy   (None_t),
      DampAm   (0.005),
      DampAk   (0.5),
      DynTh1   (0.5),
      DynTh2   (0.5),
      WithInfo (true),
      WarnRes  (false),
      WrnTol   (1.0e-7),
      RKScheme ("RK4I"),
      RKSTOL   (1.0e-2)
{
#if HAS_MPI
    if (FEM::Domain::PARA && MPI::COMM_WORLD.Get_rank()!=0) WithInfo = false;
#endif
}

inline void Solver::Solve (size_t NInc, char const * FileKey)
{
    // info
    Util::Stopwatch stopwatch(/*activated*/WithInfo);

    // initialize global matrices and vectors
    Initialize ();

    // output initial state
    if (WithInfo)
    {
        if      (Scheme==FE_t) _time_print ("Quasi-static --- FE");
        else if (Scheme==ME_t) _time_print ("Quasi-static --- ME");
        else if (Scheme==NR_t) _time_print ("Quasi-static --- NR");
    }
    if (Dom.IdxOut==0)
    {
        Dom.OutResults (FileKey);
        if (OutFun!=NULL) (*OutFun) ((*this), OutDat);
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
    }

    // weights
    if (IncsW.Size()!=NInc) SetIncsW (NInc);

    // solve
    double t0 = Dom.Time; // current time
    double tf = t0 + 1.0; // final time
    double Dt = tf - t0;  // total time increment
    double dt, tout;      // timestep and time for output
    for (Inc=0; Inc<NInc; ++Inc)
    {
        // timestep
        dt   = IncsW[Inc]*Dt; // timestep
        tout = Dom.Time + dt; // time for output

        // update U, F, Time and elements to tout
        if      (Scheme==FE_t) _FE_update (tout);
        else if (Scheme==ME_t) _ME_update (tout, FileKey);
        else if (Scheme==NR_t) _NR_update (tout);

        // update nodes to tout
        UpdateNodes ();

        // output
        if (WithInfo) _time_print ();
        if (Scheme!=ME_t) Dom.OutResults (FileKey);
        else if (!SSOut)  Dom.OutResults (FileKey);
        if (OutFun!=NULL) (*OutFun) ((*this), OutDat);

        // next tout
        tout = Dom.Time + dt;
    }
}

inline void Solver::TransSolve (double tf, double dt, double DtOut, char const * FileKey, SDPair * Steps)
{
    // info
    Util::Stopwatch stopwatch(/*activated*/WithInfo);

    // initialize global matrices and vectors
    Initialize (/*Transient*/true);

    // output initial state
    if (WithInfo)
    {
        if (TScheme==TH_t) _time_print ("Transient ------ TH(theta)");
    }
    if (Dom.IdxOut==0)
    {
        Dom.OutResults (FileKey);
        if (OutFun!=NULL) (*OutFun) ((*this), OutDat);
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
    }

    // time for output
    double dtOut = DtOut;
    double tout  = Dom.Time + dtOut;

    // nonlinear timesteps
    int    nl_nsml = 7;     // number fo small timestep __sets__
    int    nl_sch  = 0;     // scheme for larger timesteps
    double nl_ll   = 100.0; // denominator
    double nl_m    = 2.0;   // multiplier for larger timesteps in sch==1
    int    nl_n    = 10;    // number of timesteps per __set__
    int    nl_i    = 0;     // timestep set index
    int    nl_k    = 0;     // current accumulated timesteps
    int    nl_K    = 0;     // current total number of timesteps
    bool   nl_stp  = false; // use nonlinear timesteps ?
    if (Steps!=NULL)
    {
        nl_nsml = static_cast<int>((*Steps)("nsml"));
        nl_n    = static_cast<int>((*Steps)("n"));
        if (Steps->HasKey("sch")) nl_sch = static_cast<int>((*Steps)("sch"));
        if (Steps->HasKey("ll" )) nl_ll  = (*Steps)("ll");
        if (Steps->HasKey("m"  )) nl_m   = (*Steps)("m");
        nl_stp  = true;
        dt      = Timestep (nl_i, nl_nsml, nl_ll, nl_sch, nl_m);
        dtOut   = (dtOut<dt ? dt : dtOut);
    }

    // solve
    if (TScheme==TH_t)
    {
        while (Dom.Time<tf)
        {
            // update U, F, Time and elements to tout
            _TH_update (tout,dt);

            // update nodes to tout
            UpdateNodes ();

            // output
            if (WithInfo) _time_print ();
            Dom.OutResults (FileKey);
            if (OutFun!=NULL) (*OutFun) ((*this), OutDat);

            // next tout
            tout = Dom.Time + dtOut;

            // next timestep set
            if (nl_stp)
            {
                nl_K++;
                nl_k++;
                if (nl_k==nl_n)
                {
                    nl_k = 0;
                    nl_i++;
                    dt    = Timestep (nl_i, nl_nsml, nl_ll, nl_sch, nl_m);
                    dtOut = (dtOut<dt ? dt : dtOut);
                }
            }
        }
    }
}

inline void Solver::DynSolve (double tf, double dt, double DtOut, char const * FileKey, SDPair * Steps)
{
    // info
    Util::Stopwatch stopwatch(/*activated*/WithInfo);

    // initialize global matrices and vectors
    Initialize (/*Transient*/true);

    // output initial state
    if (WithInfo)
    {
        if      (DScheme==GN22_t) _time_print ("Dynamic ------ GN22");
        else if (DScheme==RK_t)
        {
            String buf("Dynamic ------ "); buf.append(RKScheme);
            _time_print (buf.CStr());
        }
    }
    if (Dom.IdxOut==0)
    {
        Dom.OutResults (FileKey);
        if (OutFun!=NULL) (*OutFun) ((*this), OutDat);
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
    }

    // time for output
    double dtOut = DtOut;
    double tout  = Dom.Time + dtOut;

    // nonlinear timesteps
    int    nl_nsml = 7;     // number fo small timestep __sets__
    int    nl_sch  = 0;     // scheme for larger timesteps
    double nl_ll   = 100.0; // denominator
    double nl_m    = 2.0;   // multiplier for larger timesteps in sch==1
    int    nl_n    = 10;    // number of timesteps per __set__
    int    nl_i    = 0;     // timestep set index
    int    nl_k    = 0;     // current accumulated timesteps
    int    nl_K    = 0;     // current total number of timesteps
    bool   nl_stp  = false; // use nonlinear timesteps ?
    if (Steps!=NULL)
    {
        nl_nsml = static_cast<int>((*Steps)("nsml"));
        nl_n    = static_cast<int>((*Steps)("n"));
        if (Steps->HasKey("sch")) nl_sch = static_cast<int>((*Steps)("sch"));
        if (Steps->HasKey("ll" )) nl_ll  = (*Steps)("ll");
        if (Steps->HasKey("m"  )) nl_m   = (*Steps)("m");
        nl_stp  = true;
        dt      = Timestep (nl_i, nl_nsml, nl_ll, nl_sch, nl_m);
        dtOut   = (dtOut<dt ? dt : dtOut);
    }

    // solve
    if (DScheme==GN22_t)
    {
        while (Dom.Time<tf)
        {
            // update U, F, Time and elements to tout
            _GN22_update (tout,dt);

            // update nodes to tout
            UpdateNodes ();

            // output
            if (WithInfo) _time_print ();
            Dom.OutResults (FileKey);
            if (OutFun!=NULL) (*OutFun) ((*this), OutDat);

            // next tout
            tout = Dom.Time + dtOut;

            // next timestep set
            if (nl_stp)
            {
                nl_K++;
                nl_k++;
                if (nl_k==nl_n)
                {
                    nl_k = 0;
                    nl_i++;
                    dt    = Timestep (nl_i, nl_nsml, nl_ll, nl_sch, nl_m);
                    dtOut = (dtOut<dt ? dt : dtOut);
                }
            }
        }
    }
    else if (DScheme==RK_t)
    {
        // ODE
        if (CteTg) NIv = 0;
        size_t nvars = 2*NEq + NIv;
        Numerical::ODESolver<Solver> ode(this, &Solver::_RK_func, nvars, RKScheme.CStr(), RKSTOL, dt);

        // set initial state
        ode.t = Dom.Time;
        _VUIV_to_Y (ode.Y);

        // solve
        while (Dom.Time<tf)
        {
            // evolve from Time to tout
            ode.Evolve (tout);

            // update V, U and elements (if not CteTg)
            _Y_to_VUIV (ode.Y);

            // update elements if not already updated by _Y_to_VUIV
            if (CteTg)
            {
                double Dt = ode.t - Dom.Time;
                for (size_t i=0; i<Dom.ActEles.Size(); ++i)
                {
                    size_t niv = Dom.ActEles[i]->NIVs();
                    if (niv>0)
                    {
                        Vec_t rate;
                        Dom.ActEles[i]->CalcIVRate (Dom.Time, U, V, rate);
                        for (size_t j=0; j<niv; ++j) Dom.ActEles[i]->SetIV (j, rate(j)*Dt);
                    }
                }
            }

            // update nodes to tout
            UpdateNodes ();

            // new time
            Dom.Time = ode.t;

            // output
            if (WithInfo) _time_print ();
            Dom.OutResults (FileKey);
            if (OutFun!=NULL) (*OutFun) ((*this), OutDat);

            // next tout
            tout = Dom.Time + dtOut;
        }
    }
}

inline void Solver::AssembleKA ()
{
    if (CteTg && A11.Top()>0) return; // constant tangent matrices => linear problems
    A11.ResetTop(); // reset top (position to insert new values) => clear triplet
    K12.ResetTop();
    K21.ResetTop();
    K22.ResetTop();
    for (size_t k=0; k<Dom.ActEles.Size(); ++k)
    {
        Mat_t         K;   // K matrix
        Array<size_t> loc; // location array
        Dom.ActEles[k]->CalcK  (K);
        Dom.ActEles[k]->GetLoc (loc);
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

    // augment A matrix and set Lagrange multipliers (if any)
    _aug_and_set_A();
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
    for (size_t k=0; k<Dom.ActEles.Size(); ++k)
    {
        Mat_t         K, M; // matrices
        Array<size_t> loc;  // location array
        Dom.ActEles[k]->CalcK  (K);
        Dom.ActEles[k]->CalcM  (M);
        Dom.ActEles[k]->GetLoc (loc);
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

    // augment A matrix and set Lagrange multipliers (if any)
    _aug_and_set_A();
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
    for (size_t k=0; k<Dom.ActEles.Size(); ++k)
    {
        // matrices
        Mat_t M, C, K;
        if      (DampTy==HMCoup_t) Dom.ActEles[k]->CalcKCM (K, C, M);
        else if (DampTy==Rayleigh_t)
        {
            Dom.ActEles[k]->CalcK (K);
            Dom.ActEles[k]->CalcM (M);
            C = DampAm*M + DampAk*K;
        }
        // set K, C, M, and A matrices
        Array<size_t> loc;
        Dom.ActEles[k]->GetLoc (loc);
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

    // augment A matrix and set Lagrange multipliers (if any)
    _aug_and_set_A();
}

inline void Solver::TgIncs (double dT, Vec_t & dU, Vec_t & dF)
{
    // assemble global K matrix
    AssembleKA ();

    // set prescribed dF
    set_to_zero (dF);
    set_to_zero (W);
    for (size_t i=0; i<Dom.NodsWithPF.Size(); ++i)
    {
        Node * const nod = Dom.NodsWithPF[i];
        for (size_t j=0; j<nod->NPF(); ++j)
        {
            int eq = nod->EqPF(j);
            if (!pU[eq]) // set dF for unknown variables only
            {
                dF(eq) = dT*nod->PF(j, /*time*/0);
                W (eq) = dF(eq); // set W1 equal to dF1
            }
        }
    }

    // set prescribed dU
    set_to_zero (dU);
    for (size_t i=0; i<Dom.NodsWithPU.Size(); ++i)
    {
        Node * const nod = Dom.NodsWithPU[i];
        for (size_t j=0; j<nod->NPU(); ++j)
        {
            int eq = nod->EqPU(j);
            if (FEM::Domain::PARA)
            {
#ifdef HAS_MPI
                if (nod->Vert.PartIDs.TheMin()==MPI::COMM_WORLD.Get_rank())
                    W(eq) = dT*nod->PU(j, /*time*/0); // set W2 equal to dU2
#endif
            }
            else W(eq) = dT*nod->PU(j, /*time*/0);
            dU(eq) = dT*nod->PU(j, /*time*/0);
        }
    }

    // correct W
    Sparse::SubMult (K12, dU, W); // W1 -= K12*dU2

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
        MPI::COMM_WORLD.Allreduce (dF.data, TmpVec.data, NEq, MPI::DOUBLE, MPI::SUM);
        dF = TmpVec;
    }
#endif
}

inline void Solver::UpdateNodes ()
{
    for (size_t i=0; i<Dom.ActNods.Size(); ++i) Dom.ActNods[i]->SetUF (U,F);
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
            for (size_t i=0; i<Dom.ActEles.Size(); ++i) Dom.ActEles[i]->UpdateState (dU, &dFint);
            MPI::COMM_WORLD.Allreduce (dFint.data, dFint_tmp.data, NEq, MPI::DOUBLE, MPI::SUM);
            F_int += dFint_tmp;
#endif
        }
        else
        {
            for (size_t i=0; i<Dom.ActEles.Size(); ++i) Dom.ActEles[i]->UpdateState (dU, &F_int);
        }
    }
    else
    {
        for (size_t i=0; i<Dom.ActEles.Size(); ++i) Dom.ActEles[i]->UpdateState (dU);
    }
}

inline void Solver::Initialize (bool Transient)
{
    // info
    Util::Stopwatch stopwatch(/*activated*/WithInfo);
    if (WithInfo) printf("\n%s--- Solver --- initializing --------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    if (FEM::Domain::PARA)
    {
#ifdef HAS_MPI
        int my_id  = MPI::COMM_WORLD.Get_rank();
        int nprocs = MPI::COMM_WORLD.Get_size();

        // compute equation numbers corresponding to local DOFs of elements
        NEq = 0;
        for (size_t i=0; i<Dom.ActEles.Size(); ++i) Dom.ActEles[i]->IncNLocDOF (NEq);

        // compute equation numbers
        for (size_t i=0; i<Dom.ActNods.Size(); ++i)
        {
            for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
            {
                // only the domain with smallest ID will set EQ number
                int min_part_id = Dom.ActNods[i]->Vert.PartIDs.TheMin();
                if (min_part_id==my_id) NEq++;
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
        for (size_t i=0; i<Dom.ActEles.Size(); ++i) Dom.ActEles[i]->IncNLocDOF (NEq);

        // assign equation numbers
        for (size_t i=0; i<Dom.ActNods.Size(); ++i)
        {
            for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
            {
                int min_part_id = Dom.ActNods[i]->Vert.PartIDs.TheMin();
                if (min_part_id==my_id) // only the domain with smallest ID will set EQ number
                {
                    Dom.ActNods[i]->Eq(j) = NEq;
                    NEq++;
                }
            }
        }

        // post messages
        const int TAG_SENT_EQ = 1000;
        for (int i=my_id+1; i<nprocs; ++i)
        {
            Array<int> inter_eq; // equation of interface DOFs
            for (size_t j=0; j<Dom.InterNodes.Size(); ++j)
            {
                int  min_part_id       =  Dom.InterNodes[j]->Vert.PartIDs.TheMin(); // the smallest proc is the one supposed to send always
                bool do_send_to_proc_i = (Dom.InterNodes[j]->Vert.PartIDs.Find(i)>=0); // found processor on the interface and has higher id than me
                if (my_id==min_part_id && do_send_to_proc_i)
                {
                    for (size_t k=0; k<Dom.InterNodes[j]->NDOF(); ++k)
                        inter_eq.Push (Dom.InterNodes[j]->Eq(k));
                }
            }
            MPI::Request req_send = MPI::COMM_WORLD.Isend (inter_eq.GetPtr(), inter_eq.Size(), MPI::INT, i, TAG_SENT_EQ);
            req_send.Wait ();
        }

        // receive messages
        MPI::Status status;
#ifdef PARALLEL_DEBUG
        Array<int> assigned_equations;
#endif
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
                int min_part_id = Dom.InterNodes[j]->Vert.PartIDs.TheMin(); // the smallest proc is the one supposed to send always
                if (source==min_part_id)
                {
                    for (size_t k=0; k<Dom.InterNodes[j]->NDOF(); ++k)
                    {
#ifdef PARALLEL_DEBUG
                        if (assigned_equations.Find(inter_eq[m])>=0) throw new Fatal("problem during assignment: Node # %d, iDOF=%zd, eq=%d. Proc # %d got message from proc # %d",Dom.InterNodes[j]->Vert.ID,k,inter_eq[m],my_id,source);
                        assigned_equations.Push(inter_eq[m]);
#endif
                        Dom.InterNodes[j]->Eq(k) = inter_eq[m];
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
        for (size_t i=0; i<Dom.ActEles.Size(); ++i) Dom.ActEles[i]->IncNLocDOF (NEq);

        // assign equation numbers
        for (size_t i=0; i<Dom.ActNods.Size(); ++i)
        {
            for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
            {
                Dom.ActNods[i]->Eq(j) = NEq;
                NEq++;
            }
        }
    }

    // prescribed equations and prescribed U
    pEQ    .Resize    (0);
    pU     .Resize    (NEq);
    pU     .SetValues (false);
    pEQproc.Resize    (0);
    for (size_t i=0; i<Dom.NodsWithPU.Size(); ++i)
    {
        Node const * nod = Dom.NodsWithPU[i];
        for (size_t j=0; j<nod->NPU(); ++j)
        {
            int eq = nod->EqPU(j);
            pU[eq] = true;
            pEQ.Push (eq);
#ifdef HAS_MPI
            if (FEM::Domain::PARA) pEQproc.Push (nod->Vert.PartIDs.TheMin());
#endif
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
    F0    .change_dim (NEq);  set_to_zero (F0);
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
    for (size_t i=0; i<Dom.ActNods.Size(); ++i) Dom.ActNods[i]->GetUF (U,F);
    F_int = F;
    F0    = F;

    // calc residual
    _cal_resid ();

    // info
    if (WithInfo)
    {
        printf("%s  Num of DOFs (NEq)  = %zd%s\n", TERM_CLR2, NEq, TERM_RST);
        printf("%s  Num of non-zeros   = %zd%s\n", TERM_CLR2, K11_size+pEQ.Size()+nzlag, TERM_RST);
    }
}

inline void Solver::SetScheme (char const * StrScheme)
{
    if      (strcmp(StrScheme,"FE")==0) Scheme = FE_t;
    else if (strcmp(StrScheme,"ME")==0) Scheme = ME_t;
    else if (strcmp(StrScheme,"NR")==0) Scheme = NR_t;
    else throw new Fatal("Solver::SetScheme: Key '%s' is invalid. The following keys are availabe: 'FE', 'ME', 'NR'",StrScheme);
}

inline void Solver::SetTransScheme (char const * StrScheme)
{
    if (strcmp(StrScheme,"TH")==0) TScheme = TH_t;
    else throw new Fatal("Solver::SetTransScheme: Key '%s' is invalid. The following keys are availabe: 'TH'",StrScheme);
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

inline void Solver::_aug_and_set_A ()
{
    // augment A11
    if (FEM::Domain::PARA)
    {
#ifdef HAS_MPI
        for (size_t i=0; i<pEQ.Size(); ++i)
        {
            if (pEQproc[i]==MPI::COMM_WORLD.Get_rank()) A11.PushEntry (pEQ[i],pEQ[i], 1.0);
        }
#endif
    }
    else
    {
        for (size_t i=0; i<pEQ.Size(); ++i) A11.PushEntry (pEQ[i],pEQ[i], 1.0);
    }

    // set equations corresponding to Lagrange multipliers
    int eqlag = NEq - NLag;

    // pins and inclined supports
    for (size_t i=0; i<Dom.NodsPins  .Size(); ++i) Dom.NodsPins  [i]->SetLagPin    (eqlag, A11);
    for (size_t i=0; i<Dom.NodsIncSup.Size(); ++i) Dom.NodsIncSup[i]->SetLagIncSup (eqlag, A11);
}

inline void Solver::_wrn_resid ()
{
    for (size_t eq=0; eq<NEq; ++eq)
    {
        if (fabs(R(eq))>WrnTol)
        {
            for (size_t i=0; i<Dom.ActNods.Size();     ++i)
            for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
                if (Dom.ActNods[i]->Eq(j)==static_cast<int>(eq))
                    printf("%sWarning: Problem with residual: Vert # %4zd, eq=%6zd, F=%15.6e, F_int=%15.6e, R=%15.6e\n%s", TERM_RED, Dom.ActNods[i]->Vert.ID, eq, F(eq), F_int(eq), R(eq),TERM_RST);
                    //throw new Fatal("%sWarning: Problem with residual: Vert # %4ld, eq=%6zd, F=%15.6e, F_int=%15.6e, R=%15.6e\n%s", TERM_RED, Dom.ActNods[i]->Vert.ID, eq, F(eq), F_int(eq), R(eq),TERM_RST);
        }
    }
}

inline void Solver::_cal_resid ()
{
    // calculate residual
    R = F - F_int;

    //printf("\n######################################   Before   ####################################\n");
#ifdef DO_DEBUG
    _wrn_resid ();
#endif
    
    // number of the first equation corresponding to Lagrange multipliers
    int eqlag = NEq - NLag;

    // clear forces due to pins and inclined supports
    for (size_t i=0; i<Dom.NodsPins  .Size(); ++i) Dom.NodsPins  [i]->ClrRPin    (eqlag, R);
    for (size_t i=0; i<Dom.NodsIncSup.Size(); ++i) Dom.NodsIncSup[i]->ClrRIncSup (eqlag, R);

    // clear forces due to supports
    //for (size_t i=0; i<pEQ.Size(); ++i) R(pEQ[i]) = 0.0;

    //printf("\n######################################   After   #####################################\n");
    if (WarnRes) _wrn_resid ();

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
                for (size_t j=0; j<Dom.InterNodes[i]->NDOF(); ++j)
                {
                    long   eq   = Dom.InterNodes[i]->Eq(j);
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
    double dt = (tf-Dom.Time)/nSS;
    for (Stp=0; Stp<nSS; ++Stp)
    {
        // calculate tangent increments
        TgIncs (dt, dU, dF);

        // update elements
        UpdateElements (dU, /*CalcFint*/true);

        // update U, F, and Time
        U        += dU;
        F        += dF;
        Dom.Time += dt;

        // debug
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
    }

    // residual
    _cal_resid ();
}

inline void Solver::_ME_update (double tf, char const * FNK)
{
    // auxiliar vectors
    Vec_t dU_fe(NEq), dU_tm(NEq), dU_me(NEq), U_me(NEq), U_dif(NEq);
    Vec_t dF_fe(NEq), dF_tm(NEq), dF_me(NEq), F_me(NEq), F_dif(NEq);

    // for each pseudo time T
    double T   = 0.0;
    double dT  = (dTlast>0.0 ? dTlast : dTini);
    double Dt  = tf-Dom.Time;
    for (Stp=0; Stp<MaxSS; ++Stp)
    {
        // exit point
        if (T>=1.0) break;

        // backup state of elements
        for (size_t i=0; i<Dom.ActEles.Size(); ++i) Dom.ActEles[i]->BackupState ();

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
        for (size_t i=0; i<Dom.ActEles.Size(); ++i) Dom.ActEles[i]->RestoreState ();

        // update
        if (error<STOL)
        {
            UpdateElements (dU_me, /*CalcFint*/true);
            T        += dT;
            U         = U_me;
            F         = F_me;
            Dom.Time += dt;
            //_cal_resid ();
            //_cor_resid (dU_me, dF_me);
            if (m>mMax) m = mMax;
            if (SSOut || (DbgFun!=NULL)) UpdateNodes ();
            if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
            if (SSOut) Dom.OutResults (FNK);
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
    double dt = (tf-Dom.Time)/nSS;
    for (Stp=0; Stp<nSS; ++Stp)
    {
        // calculate tangent increments
        TgIncs (dt, dU, dF);

        // update elements
        UpdateElements (dU, /*CalcFint*/true);

        // update U, F, and Time
        U        += dU;
        F        += dF;
        Dom.Time += dt;

        // residual
        _cal_resid ();
        _cor_resid (dU, dF);

        // debug
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
    }
}

inline void Solver::_TH_update (double tf, double Dt)
{
    // timestep
    double dt = (Dom.Time+Dt>tf ? tf-Dom.Time : Dt);

    while (Dom.Time<tf)
    {
        // calculate W1 and F1
        set_to_zero (W);
        double t0 = Dom.Time;
        double t1 = Dom.Time + dt;
        for (size_t i=0; i<Dom.NodsWithPF.Size(); ++i)
        {
            Node * const nod = Dom.NodsWithPF[i];
            for (size_t j=0; j<nod->NPF(); ++j)
            {
                int eq = nod->EqPF(j);
                if (!pU[eq])
                {
                    W(eq) = (F0(eq) + nod->PF(j,t1)*TransTh + nod->PF(j,t0)*(1.0-TransTh)) * dt;
                }
            }
        }

        // set W2 with prescribed DeltaU = U(n+1) - U(n)
        for (size_t i=0; i<Dom.NodsWithPU.Size(); ++i)
        {
            Node * const nod = Dom.NodsWithPU[i];
            for (size_t j=0; j<nod->NPU(); ++j)
            {
                int eq = nod->EqPU(j);
                W(eq) = nod->PU(j, t1) - U(eq);
            }
        }

        // solve
        double c1 = dt*TransTh;
        AssembleKMA     (1.0, c1);         // A11 = M + dt*th*K
        Sparse::SubMult (dt, K11, U, W);   // W1 -=    dt*K11*U1(n)
        Sparse::SubMult (dt, K12, U, W);   // W1 -=    dt*K12*U2(n)
        Sparse::SubMult (    M12, W, W);   // W1 -=       M12*W2
        Sparse::SubMult (c1, K12, W, W);   // W1 -= dt*th*K12*W2
        UMFPACK::Solve  (A11, W, dU);      // dU  = inv(A11)*W
        UpdateElements  (dU, /*CalcFint*/false);

        // update
        U += dU;

        // next time step
        dt = (Dom.Time+dt>tf ? tf-Dom.Time : dt);
        Dom.Time += dt;

        // debug
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
    }
}

inline void Solver::_GN22_update (double tf, double Dt)
{
    // timestep
    double dt = (Dom.Time+Dt>tf ? tf-Dom.Time : Dt);

    // constants
    const double c1 = dt*dt*(1.0-DynTh2)/2.0;
    const double c2 = dt*(1.0-DynTh1);
    const double c3 = 2.0/(dt*dt*DynTh2);
    const double c4 = 2.0*DynTh1/(dt*DynTh2);

    while (Dom.Time<tf)
    {
        // set prescribed F
        F = F0;
        double tb = Dom.Time+DynTh1*dt;
        for (size_t i=0; i<Dom.NodsWithPF.Size(); ++i)
        {
            Node * const nod = Dom.NodsWithPF[i];
            for (size_t j=0; j<nod->NPF(); ++j)
            {
                int eq = nod->EqPF(j);
                if (!pU[eq]) F(eq) = F0(eq) + nod->PF(j, tb);
            }
        }
        for (size_t i=0; i<Dom.ActEles.Size(); ++i) Dom.ActEles[i]->AddToF (tb, F);
        double normF = Norm(F);

        // predictor
        Us = U + dt*V + c1*A;
        Vs = V + c2*A;
        A  = c3*(U - Us);
        V  = Vs + (DynTh1*dt)*A;

        // iterations
        for (It=0; It<MaxIt; ++It)
        {
            // residual
            R = F - F_int;
            for (size_t i=0; i<pEQ.Size(); ++i)
            {
                F(pEQ[i]) = 0.0; // F2 = 0
                R(pEQ[i]) = 0.0; // R2 = 0   clear residual corresponding to supports
            }

            // assemble Amat
            if (DampTy==None_t) AssembleKMA  (c3,     1.0);  // A = c3*M        + K
            else                AssembleKCMA (c3, c4, 1.0);  // A = c3*M + c4*C + K

            // solve for dU
            Sparse::SubMult (M11, A, R);  if (DampTy!=None_t)  // R -= M11*A
            Sparse::SubMult (C11, V, R);                       // R -= C11*V
            UMFPACK::Solve  (A11, R, dU);                      // dU = inv(A11)*R

            // set prescribed dU
            for (size_t i=0; i<Dom.NodsWithPU.Size(); ++i)
            {
                Node * const nod = Dom.NodsWithPU[i];
                for (size_t j=0; j<nod->NPU(); ++j)
                {
                    int eq = nod->EqPU(j);
                    dU(eq) = nod->PU(j,tb) - U(eq);
                }
            }

#ifdef DO_DEBUG
            double normdU = Norm(dU);
            if (Util::IsNan(normdU)) throw new Fatal("Solver::_GN22_update: normdU is NaN");
#endif

            // update elements
            UpdateElements (dU, /*CalcFint*/true);

#ifdef DO_DEBUG
            double normFint = Norm(F_int);
            if (Util::IsNan(normFint)) throw new Fatal("Solver::_GN22_update: normFint is NaN");
#endif

            // update state
            U += dU;
            A  = c3*(U - Us);
            V  = Vs + (DynTh1*dt)*A;

            // set prescribed U, V and A
            for (size_t i=0; i<Dom.NodsWithPU.Size(); ++i)
            {
                Node * const nod = Dom.NodsWithPU[i];
                for (size_t j=0; j<nod->NPU(); ++j)
                {
                    int eq = nod->EqPU(j);
                    U(eq) = nod->PU(j, tb);
                    V(eq) = nod->PV(j, tb);
                    A(eq) = nod->PA(j, tb);
                }
            }

            // calculate F2
            Sparse::AddMult (M21, A, F);                        // F2 += M21*A1
            Sparse::AddMult (M22, A, F);  if (DampTy!=None_t) { // F2 += M22*A2
            Sparse::AddMult (C21, V, F);                        // F2 += C21*V1
            Sparse::AddMult (C22, V, F); }                      // F2 += C22*V2
            Sparse::AddMult (K21, U, F);                        // F2 += K21*V1
            Sparse::AddMult (K22, U, F);                        // F2 += K22*V2

            // check convergence
            NormR    = Norm(R);
            MaxNormF = Util::Max (normF, Norm(F_int));
#ifdef DO_DEBUG
            if (Util::IsNan(NormR)) throw new Fatal("Solver::_GN22_update: NormR is NaN");
            printf("NormR = %g\n",NormR);
#endif
            if (ResidOK()) break;
        }
        if (It>=MaxIt) throw new Fatal("Solver::_GN22_update: Generalized-Newmark (GN22) did not converge after %d iterations (TolR=%g). NormR = %g",It,TolR,NormR);

        // next time step
        dt = (Dom.Time+dt>tf ? tf-Dom.Time : dt);
        Dom.Time += dt;

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
    printf("%10.6f  %s%8e%s  %4zd  %4zd\n",Dom.Time,(ResidOK()?TERM_GREEN:TERM_RED),NormR,TERM_RST,Stp,It);
}

inline void Solver::_VUIV_to_Y (double Y[])
{
    for (size_t i=0; i<NEq; ++i)
    {
        Y[i]     = V(i);
        Y[NEq+i] = U(i);
    }
    if (!CteTg)
    {
        for (size_t i=0; i<Dom.ActEles.Size();     ++i)
        for (size_t j=0; j<Dom.ActEles[i]->NIVs(); ++j)
            Y[2*NEq+j] = Dom.ActEles[i]->GetIV (j);
    }
}

inline void Solver::_Y_to_VUIV (double const Y[])
{
    for (size_t i=0; i<NEq; ++i)
    {
        V(i) = Y[i];
        U(i) = Y[NEq+i];
    }
    if (!CteTg)
    {
        for (size_t i=0; i<Dom.ActEles.Size(); ++i)
        for (size_t j=0; j<Dom.ActEles[i]->NIVs(); ++j)
            Dom.ActEles[i]->SetIV (j, Y[2*NEq+j]);
    }
}

inline int Solver::_RK_func (double t, double const Y[], double dYdt[])
{
    // get current U and V and set elements with current IVs (needs to be before the assembly)
    _Y_to_VUIV (Y);

    // prescribed U, V, and A
    set_to_zero (A);
    for (size_t i=0; i<Dom.NodsWithPU.Size(); ++i)
    {
        Node * const nod = Dom.NodsWithPU[i];
        for (size_t j=0; j<nod->NPU(); ++j)
        {
            int eq = nod->EqPU(j);
            U(eq) = nod->PU(j, t);
            V(eq) = nod->PV(j, t);
            A(eq) = nod->PA(j, t);
        }
    }

    // prescribed F
    F = F0;
    for (size_t i=0; i<Dom.NodsWithPF.Size(); ++i)
    {
        Node * const nod = Dom.NodsWithPF[i];
        for (size_t j=0; j<nod->NPF(); ++j)
        {
            int eq = nod->EqPF(j);
            if (!pU[eq]) F(eq) = F0(eq) + nod->PF(j, t);
        }
    }
    for (size_t i=0; i<Dom.ActEles.Size(); ++i) Dom.ActEles[i]->AddToF (t, F);

    // assembly
    if (DampTy==None_t) AssembleKMA  (1.0, 0.0);      // A11 = M11
    else                AssembleKCMA (1.0, 0.0, 0.0); // A11 = M11

    // set A(t) and V(t)
    // F = F1 - p1 - c1 - M12*A2
    //   = F1 - K11*U1 - K12*U2 - C11*V1 - C12*V2 - M12*A2
    Sparse::SubMult (K11, U, F);                               // F -= K11*U1
    Sparse::SubMult (K12, U, F);  if (DampTy!=None_t) {        // F -= K12*U2
    Sparse::SubMult (C11, V, F);                               // F -= C11*V1
    Sparse::SubMult (C12, V, F); }                             // F -= C12*V2
    Sparse::SubMult (M12, A, F);                               // F -= M12*A2
    for (size_t i=0; i<pEQ.Size(); ++i) F(pEQ[i]) = A(pEQ[i]); // F2 = A2
    UMFPACK::Solve (A11, F, A);                                // A = inv(M11)*F

    // set dYdt
    for (size_t i=0; i<NEq; ++i)
    {
        dYdt[i]     = A(i); // dVdt
        dYdt[NEq+i] = V(i); // dUdt
    }
    if (!CteTg)
    {
        for (size_t i=0; i<Dom.ActEles.Size(); ++i)
        {
            size_t niv = Dom.ActEles[i]->NIVs();
            if (niv>0)
            {
                Vec_t rate;
                Dom.ActEles[i]->CalcIVRate (t, U, V, rate);
                for (size_t j=0; j<niv; ++j) dYdt[2*NEq+j] = rate(j);
            }
        }
    }

    return GSL_SUCCESS;
}

inline void Solver::_debug_print_matrices (bool Stop)
{
    Sparse::Matrix<double,int> MM11(M11), MM12(M12), MM21(M21), MM22(M22);
    Sparse::Matrix<double,int> KK11(K11), KK12(K12), KK21(K21), KK22(K22);
    Mat_t mm11, mm12, mm21, mm22, mm;
    Mat_t kk11, kk12, kk21, kk22, kk;
    MM11.GetDense(mm11); MM12.GetDense(mm12); MM21.GetDense(mm21); MM22.GetDense(mm22);
    KK11.GetDense(kk11); KK12.GetDense(kk12); KK21.GetDense(kk21); KK22.GetDense(kk22);
    mm = mm11 + mm12 + mm21 + mm22;
    kk = kk11 + kk12 + kk21 + kk22;
    A11.WriteSMAT("A11");
    WriteSMAT    (mm,"M");
    WriteSMAT    (kk,"K");
    std::cout << std::endl << std::endl;
    printf("det(A11) = %g\n", UMFPACK::Det(A11));
    printf("det(M)   = %g\n", Det(mm));
    printf("det(K)   = %g\n", Det(kk));
    if (DampTy!=None_t)
    {
        Sparse::Matrix<double,int> CC11(C11), CC12(C12), CC21(C21), CC22(C22);
        Mat_t cc11, cc12, cc21, cc22, cc;
        CC11.GetDense(cc11); CC12.GetDense(cc12); CC21.GetDense(cc21); CC22.GetDense(cc22);
        cc = cc11 + cc12 + cc21 + cc22;
        printf("det(C)   = %g\n", Det(cc));
        WriteSMAT(cc,"C");  printf("Matrix <%sC.smat%s> written\n",TERM_CLR_BLUE_H,TERM_RST);
    }
    printf("Matrix <%sA11.smat%s> written\n",TERM_CLR_BLUE_H,TERM_RST);
    printf("Matrix <%sM.smat%s> written\n",TERM_CLR_BLUE_H,TERM_RST);
    printf("Matrix <%sK.smat%s> written\n",TERM_CLR_BLUE_H,TERM_RST);
    std::cout << std::endl;
    if (Stop) throw new Fatal("Solver::_debug_print_matrices   STOP");
}


}; // namespace FEM

#endif // MECHSYS_FEM_SOLVER_H
