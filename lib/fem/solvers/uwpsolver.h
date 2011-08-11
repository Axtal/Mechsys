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

#ifndef MECHSYS_FEM_UWPSOLVER_H
#define MECHSYS_FEM_UWPSOLVER_H

// MechSys
#include <mechsys/fem/solver.h>
#include <mechsys/fem/uwpelem.h>

using std::cout;
using std::endl;

namespace FEM
{

class UWPSolver : public Solver
{
public:
    // Constructor
    UWPSolver (Domain & Dom, SDPair const & Flags, pOutFun OutFun=NULL, void * OutDat=NULL,
                                                   pOutFun DbgFun=NULL, void * DbgDat=NULL); ///< Allocate solver object

    // Methods
    String Name       () const { return "UWPSolver"; }
    void   Initialize ();
    void   Assemble   (double Time, Vec_t const & V, Vec_t const & Ww, Vec_t const & Pw);
    void   RKFunc     (double Time, Vec_t const & V, Vec_t const & Ww, Vec_t const & Pw);
    void   Solve      (double tf, double dt, double dtOut, char const * FileKey=NULL); ///< Solve dynamic problem

    // Data
    Array<long>                   pEqU, pEqW, pEqP;   ///< prescribed equations
    Array<bool>                   pU,   pW,   pP;     ///< prescribed ?
    Vec_t                         V, Ww, Pw;          ///< Vectors
    Vec_t                         V_fe, Ww_fe, Pw_fe; ///< Vectors
    Vec_t                         V_me, Ww_me, Pw_me; ///< Vectors
    Vec_t                         V_er, Ww_er, Pw_er; ///< Vectors
    Vec_t                         DV,   DWw,   DPw;   ///< Vectors
    Vec_t                         Rv, Rww;            ///< Vectors
    Vec_t                         dVdt, dWwdt, dPwdt; ///< Vectors
    Vec_t                         f, fw, hw, cw, ew;  ///< Vectors
    Sparse::Triplet<double,int>   M, Mb, Mw, Hw;      ///< Matrices
    bool                          WithInfo;           ///< Print information ?
    double                        STOL, mMin, mMax, dTini, dTlast;
    size_t                        MaxSS;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline UWPSolver::UWPSolver (Domain & Dom, SDPair const & Flags, pOutFun OutFun, void * OutDat, pOutFun DbgFun, void * DbgDat)
    : Solver (Dom, Flags, OutFun, OutDat, DbgFun, DbgDat), WithInfo(true),
      STOL(1.0e-2), mMin(0.2), mMax(2.0), dTini(0.01), dTlast(-1.0), MaxSS(100)
{
}

inline void UWPSolver::Initialize ()
{
    // assign equation numbers
    size_t NEqU = 0;
    size_t NEqW = 0;
    size_t NEqP = 0;
    for (size_t i=0; i<Dom.ActNods.Size();     ++i)
    for (size_t j=0; j<Dom.ActNods[i]->NDOF(); ++j)
    {
        if (Dom.ActNods[i]->UKey(j)=="wwx" || Dom.ActNods[i]->UKey(j)=="wwy" || Dom.ActNods[i]->UKey(j)=="wwz")
        {
            Dom.ActNods[i]->Eq(j) = NEqW;
            NEqW++;
        }
        else if (Dom.ActNods[i]->UKey(j)=="pw")
        {
            Dom.ActNods[i]->Eq(j) = NEqP;
            NEqP++;
        }
        else
        {
            Dom.ActNods[i]->Eq(j) = NEqU;
            NEqU++;
        }
    }

    // prescribed equations
    pEqU.Resize(0);  pU.Resize(NEqU);  pU.SetValues(false);
    pEqW.Resize(0);  pW.Resize(NEqW);  pW.SetValues(false);
    pEqP.Resize(0);  pP.Resize(NEqP);  pP.SetValues(false);
    for (size_t i=0; i<Dom.NodsWithPU.Size(); ++i)
    {
        Node const * nod = Dom.NodsWithPU[i];
        for (size_t j=0; j<nod->NPU(); ++j)
        {
            int eq = nod->EqPU(j);
            if (nod->PUKey(j)=="wwx" || nod->PUKey(j)=="wwy" || nod->PUKey(j)=="wwz") { pEqW.Push(eq);  pW[eq]=true; }
            else if (nod->PUKey(j)=="pw")                                             { pEqP.Push(eq);  pP[eq]=true; }
            else                                                                      { pEqU.Push(eq);  pU[eq]=true; }
        }
    }
    //cout << "pEqU = " << pEqU << endl;
    //cout << "pEqW = " << pEqW << endl;
    //cout << "pEqP = " << pEqP << endl;

    // find total number of non-zero entries, including duplicates
    size_t Mb_size = 0;
    size_t Mw_size = 0;
    size_t Hw_size = 0;
    for (size_t k=0; k<Dom.ActEles.Size(); ++k)
    {
        if (Dom.ActEles[k]->IsUWP)
        {
            Dom.ActEles[k]->GetLoc ();
            for (size_t i=0; i<UWPElem::LocU .Size(); ++i)
            for (size_t j=0; j<UWPElem::LocU .Size(); ++j) Mb_size++;
            for (size_t i=0; i<UWPElem::LocWw.Size(); ++i)
            for (size_t j=0; j<UWPElem::LocWw.Size(); ++j) Mw_size++;
            for (size_t i=0; i<UWPElem::LocPw.Size(); ++i)
            for (size_t j=0; j<UWPElem::LocPw.Size(); ++j) Hw_size++;
        }
        else throw new Fatal("UWPElem::Initialize: this method does not work with elements that are not UWP yet");
    }

    // resize vectors and matrices
    V     . change_dim (NEqU);
    Ww    . change_dim (NEqW);
    Pw    . change_dim (NEqP);
    V_fe  . change_dim (NEqU);
    Ww_fe . change_dim (NEqW);
    Pw_fe . change_dim (NEqP);
    V_me  . change_dim (NEqU);
    Ww_me . change_dim (NEqW);
    Pw_me . change_dim (NEqP);
    V_er  . change_dim (NEqU);
    Ww_er . change_dim (NEqW);
    Pw_er . change_dim (NEqP);
    DV    . change_dim (NEqU);
    DWw   . change_dim (NEqW);
    DPw   . change_dim (NEqP);
    Rv    . change_dim (NEqU);
    Rww   . change_dim (NEqW);
    dVdt  . change_dim (NEqU);
    dWwdt . change_dim (NEqW);
    dPwdt . change_dim (NEqP);
    f     . change_dim (NEqU);
    fw    . change_dim (NEqW);
    hw    . change_dim (NEqW);
    cw    . change_dim (NEqU);
    ew    . change_dim (NEqP);
    M     . AllocSpace (NEqU,NEqU,Mb_size);
    Mb    . AllocSpace (NEqU,NEqU,Mb_size);
    Mw    . AllocSpace (NEqW,NEqW,Mw_size);
    Hw    . AllocSpace (NEqP,NEqP,Hw_size);
}

inline void UWPSolver::Assemble (double Time, Vec_t const & TheV, Vec_t const & TheWw, Vec_t const & ThePw)
{
    M .ResetTop ();
    Mb.ResetTop ();
    Mw.ResetTop ();
    Hw.ResetTop ();
    set_to_zero (f);
    set_to_zero (fw);
    set_to_zero (hw);
    set_to_zero (cw);
    set_to_zero (ew);
    for (size_t k=0; k<Dom.ActEles.Size(); ++k)
    {
        if (Dom.ActEles[k]->IsUWP)
        {
            Dom.ActEles[k]->GetLoc  ();
            Dom.ActEles[k]->ElemEqs (Time, TheV, TheWw, ThePw);

            for (size_t i=0; i<UWPElem::LocU .Size(); ++i)
            {
                f (UWPElem::LocU[i]) += UWPElem::f (i);
                cw(UWPElem::LocU[i]) += UWPElem::cw(i);
                for (size_t j=0; j<UWPElem::LocU.Size(); ++j)
                {
                    M.PushEntry (UWPElem::LocU[i], UWPElem::LocU[j], UWPElem::M(i,j));
                    if (!pU[UWPElem::LocU[i]] && !pU[UWPElem::LocU[j]]) Mb.PushEntry (UWPElem::LocU[i], UWPElem::LocU[j], UWPElem::M(i,j) - UWPElem::Mw(i,j));
                }
            }

            for (size_t i=0; i<UWPElem::LocWw.Size(); ++i)
            {
                fw(UWPElem::LocWw[i]) += UWPElem::fw(i);
                hw(UWPElem::LocWw[i]) += UWPElem::hw(i);
                for (size_t j=0; j<UWPElem::LocWw.Size(); ++j)
                {
                    if (!pW[UWPElem::LocWw[i]] && !pW[UWPElem::LocWw[j]]) Mw.PushEntry (UWPElem::LocWw[i], UWPElem::LocWw[j], UWPElem::Mw(i,j));
                }
            }

            for (size_t i=0; i<UWPElem::LocPw.Size(); ++i)
            {
                ew(UWPElem::LocPw[i]) += UWPElem::ew(i);
                for (size_t j=0; j<UWPElem::LocPw.Size(); ++j)
                {
                    if (!pP[UWPElem::LocPw[i]] && !pP[UWPElem::LocPw[j]]) Hw.PushEntry (UWPElem::LocPw[i], UWPElem::LocPw[j], UWPElem::Hw(i,j));
                }
            }

            //cout << "M  = \n" << PrintMatrix(UWPElem::M);
            //cout << "Mw = \n" << PrintMatrix(UWPElem::Mw);
            //cout << "Hw = \n" << PrintMatrix(UWPElem::Hw);
            //cout << "f  = \n" << PrintVector(UWPElem::f);
            //cout << "fw = \n" << PrintVector(UWPElem::fw);
            //cout << "hw = \n" << PrintVector(UWPElem::hw);
            //cout << "cw = \n" << PrintVector(UWPElem::cw);
            //cout << "ew = \n" << PrintVector(UWPElem::ew);
        }
        else
        {
            throw new Fatal("UWPElem::Assembly: this method does not work with elements that are not u-w-p yet");
            Array<size_t> loc;
            Mat_t         M;
            Dom.ActEles[k]->GetLoc (loc);
            Dom.ActEles[k]->CalcM  (M);
        }
    }

    for (size_t i=0; i<pEqU.Size(); ++i) Mb.PushEntry (pEqU[i],pEqU[i], 1.0);
    for (size_t i=0; i<pEqW.Size(); ++i) Mw.PushEntry (pEqW[i],pEqW[i], 1.0);
    for (size_t i=0; i<pEqP.Size(); ++i) Hw.PushEntry (pEqP[i],pEqP[i], 1.0);

    /*
    M .WriteSMAT("M");
    Mb.WriteSMAT("Mb");
    Mw.WriteSMAT("Mw");
    Hw.WriteSMAT("Hw");
    Sparse::Matrix<double,int> sM(M), sMb(Mb), sMw(Mw), sHw(Hw);
    Mat_t Mmat, Mbmat, Mwmat, Hwmat, Mi, Mbi, Mwi, Hwi;
    sM .GetDense (Mmat);
    sMb.GetDense (Mbmat);
    sMw.GetDense (Mwmat);
    sHw.GetDense (Hwmat);
    cout << "M  = \n"       << PrintMatrix(Mmat,  "%10.5f");
    cout << "Mb = \n"       << PrintMatrix(Mbmat, "%10.5f");
    cout << "Mw = \n"       << PrintMatrix(Mwmat, "%10.5f");
    cout << "Hw = \n"       << PrintMatrix(Hwmat, "%10.5f");
    Inv (Mmat,  Mi);
    Inv (Mbmat, Mbi);
    Inv (Mwmat, Mwi);
    Inv (Hwmat, Hwi);
    cout << "Inv(M)   = \n" << PrintMatrix(Mi,  "%10.5f");
    cout << "Inv(Mbi) = \n" << PrintMatrix(Mbi, "%10.5f");
    cout << "Inv(Mw)  = \n" << PrintMatrix(Mwi, "%10.5f");
    cout << "Inv(Hw)  = \n" << PrintMatrix(Hwi, "%10.1f");
    cout << "det(M)  = "    << UMFPACK::Det (sM)  << endl;
    cout << "det(Mb) = "    << UMFPACK::Det (sMb) << endl;
    cout << "det(Mw) = "    << UMFPACK::Det (sMw) << endl;
    cout << "det(Hw) = "    << UMFPACK::Det (sHw) << endl;
    */
}

inline void UWPSolver::RKFunc (double Time, Vec_t const & TheV, Vec_t const & TheWw, Vec_t const & ThePw)
{
    Assemble (Time, TheV, TheWw, ThePw);

    Rv  = f + hw - fw;
    Rww = f - cw;

    /*
    cout << "f     = " << PrintVector(f,      "%10.5f");
    cout << "fw    = " << PrintVector(fw,     "%10.5f");
    cout << "hw    = " << PrintVector(hw,     "%10.5f");
    cout << "cw    = " << PrintVector(cw,     "%10.5f");
    cout << "ew    = " << PrintVector(ew,     "%10.5f");
    cout << "Rv    = " << PrintVector(Rv,     "%10.5f");
    cout << "Rww   = " << PrintVector(Rww,    "%10.5f");
    */

    UMFPACK::Solve  (Mb, Rv,    dVdt);  // dVdt  = inv(Mb)*Rv

    for (size_t i=0; i<pEqU.Size(); ++i) dVdt(pEqU[i]) = 0.0;

    Sparse::SubMult (M,  dVdt,  Rww);   // Rww  -= M*dVdt
    UMFPACK::Solve  (Mw, Rww,   dWwdt); // dWwdt = inv(Mw)*Rww

    for (size_t i=0; i<pEqW.Size(); ++i) dWwdt(pEqW[i]) = 0.0;

    UMFPACK::Solve  (Hw, ew,    dPwdt); // dPwdt = inv(Hw)*ew

    // set prescribed U, V and A
    for (size_t i=0; i<Dom.NodsWithPU.Size(); ++i)
    {
        Node * const nod = Dom.NodsWithPU[i];
        for (size_t j=0; j<nod->NPU(); ++j)
        {
            if (nod->PUKey(j)=="pw")
            {
                int eq = nod->EqPU(j);
                dPwdt(eq) = nod->PV(j, Time);
            }
        }
    }

    /*
    cout << "Rv   (after) = " << PrintVector(Rv,     "%10.5f");
    cout << "Rww  (after) = " << PrintVector(Rww,    "%10.5f");
    cout << "dVdt  = "        << PrintVector(dVdt,   "%10.5f");
    cout << "dWwdt = "        << PrintVector(dWwdt,  "%10.5f");
    cout << "dPwdt = "        << PrintVector(dPwdt,  "%10.5f");
    cout << "V     = "        << PrintVector(V,      "%10.5f");
    cout << "Ww    = "        << PrintVector(Ww,     "%10.5f");
    cout << "Pw    = "        << PrintVector(Pw,     "%10.5f") << "\n";
    */
}

inline void UWPSolver::Solve (double tf, double, double dtOut, char const * FileKey)
{
    // info
    Util::Stopwatch stopwatch(WithInfo);

    // initialize global matrices and vectors
    Initialize ();

    // output initial state
    printf ("\n%s--- Stage solution --- u-w-p Solver -- Runge-Kutta --------------------------%s\n",TERM_CLR1,TERM_RST);
    printf ("%s%10s  %4s%s\n",TERM_CLR2,"Time","NSS",TERM_RST);
    printf ("%10.6f  %4zd\n",Dom.Time,0);
    if (Dom.IdxOut==0)
    {
        Dom.OutResults (FileKey);
        if (OutFun!=NULL) (*OutFun) (this, OutDat);
        if (DbgFun!=NULL) (*DbgFun) (this, DbgDat);
    }

    // time for output
    double tout = Dom.Time + dtOut;
    while (Dom.Time<tf)
    {
        // update to tout
        size_t stp = 0;
        double T   = 0.0;
        double dT  = (dTlast>0.0 ? dTlast : dTini);
        double Dt  = tout - Dom.Time;
        for (stp=0; stp<MaxSS; ++stp)
        {
            // exit
            if (T>=1.0) break;

            // time increment
            double dt = dT*Dt;

            // FE state
            RKFunc (Dom.Time, V, Ww, Pw);
            DV    = dVdt  * dt;
            DWw   = dWwdt * dt;
            DPw   = dPwdt * dt;
            V_fe  = V  + DV;
            Ww_fe = Ww + DWw;
            Pw_fe = Pw + DPw;
            for (size_t k=0; k<Dom.ActEles.Size(); ++k) Dom.ActEles[k]->Update (0, V, dPwdt, dt);

            // ME State
            RKFunc (Dom.Time+dt, V_fe, Ww_fe, Pw_fe);
            V_me  = V  + 0.5*DV  + (0.5*dt)*dVdt;
            Ww_me = Ww + 0.5*DWw + (0.5*dt)*dWwdt;
            Pw_me = Pw + 0.5*DPw + (0.5*dt)*dPwdt;
            for (size_t k=0; k<Dom.ActEles.Size(); ++k) Dom.ActEles[k]->Update (1, V_fe, dPwdt, dt);

            // error estimate
            V_er  = V_me  - V_fe;
            Ww_er = Ww_me - Ww_fe;
            Pw_er = Pw_me - Pw_fe;
            double V_error  = Norm(V_er)  / (1.0+Norm(V_me));
            double Ww_error = Norm(Ww_er) / (1.0+Norm(Ww_me));
            double Pw_error = Norm(Pw_er) / (1.0+Norm(Pw_me));
            double error    = Util::Max (V_error, Ww_error, Pw_error);

            // step multiplier
            double m = (error>0.0 ? 0.9*sqrt(STOL/error) : mMax);

            //cout << "error = " << error << ", elems_error = " << elems_error << endl;
            //cout << "error = " << error << endl;
            //break;

            // update
            //if (error<STOL)
            {
                T        += dT;
                Dom.Time += dt;
                V         = V_me;
                Ww        = Ww_me;
                Pw        = Pw_me;
                if (m>mMax) m = mMax;
            }
            //else
            //{
                //for (size_t k=0; k<Dom.ActEles.Size(); ++k) Dom.ActEles[k]->Restore ();
                //if (m<mMin) m = mMin;
            //}
            //dT = m * dT;
            if (dT>1.0-T) dT = 1.0-T;
            else dTlast = dT;
        }
        //if (stp==MaxSS) throw new Fatal("UWPSolver::Solve: Runge-Kutta (2nd order / ME) did not converge after %zd substeps",stp);

        // update nodes to tout
        
        // output
        printf ("%10.6f  %4zd\n",Dom.Time,stp);
        Dom.OutResults (FileKey);
        if (OutFun!=NULL) (*OutFun) (this, OutDat);

        // next tout
        tout = Dom.Time + dtOut;
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


Solver * UWPSolverMaker(Domain & Dom, SDPair const & Flags, Solver::pOutFun OutFun, void * OutDat, Solver::pOutFun DbgFun, void * DbgDat)
{
    return new UWPSolver (Dom, Flags, OutFun, OutDat, DbgFun, DbgDat);
}

int UWPSolverRegister()
{
    SolverFactory["UWPSolver"] = UWPSolverMaker;
    return 0;
}

int __UWPSolver_dummy_int  = UWPSolverRegister();

};

#endif // MECHSYS_FEM_UWPSOLVER_H
