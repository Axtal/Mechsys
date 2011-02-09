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

#ifdef HMSTRESSUPDATE_DECLARE

class HMStressUpdate
{
public:
    // static
    static Vec_t Iv; ///< Identity vector (NCo)

    // callbacks
    typedef void (*pDbgFun) (HMStressUpdate const & SU, void * UserData); ///< Pointer to debug function

    // Constructor & Destructor
     HMStressUpdate ();
    ~HMStressUpdate () { if (ODE!=NULL) delete ODE; }

    // Methods
    void SetModel (Model const * TheMdl, UnsatFlow const * TheFMdl);
    void Update   (Vec_t const & DEps, double Dpw, State * Sta, UnsatFlowState * FSta, Vec_t & DSig);
    void GetInfo  (std::ostream & os, bool Header=false) const;

    // Data
    Model     const * Mdl;  ///< Equilib model
    UnsatFlow const * FMdl; ///< Flow model
    pDbgFun           DbgFun;
    void            * DbgDat;

    // ODE solver
    Numerical::ODESolver<HMStressUpdate> * ODE;

    // Constants for integration
    double  STOL;
    double  dTini;
    double  mMin;
    double  mMax;
    size_t  MaxSS;
    bool    CDrift; ///< correct drift ?
    String  RKScheme;
    size_t  SS;     ///< num of substeps
    size_t  SSs;    ///< num of successful substeps
    size_t  DCit;   ///< drift correction iterations
    size_t  DCitEl; ///< drift correction (elastic)
    size_t  ncp;    ///< num of stress components
    size_t  niv;    ///< num of internal variables
    size_t  neq;    ///< num of equations for ODE
    Vec_t   dsig;
    Vec_t   deps;
    Vec_t   divs;
    Vec_t   eps0;
    double  pc0;
    double  n0;
    double  dpc;
    double  dev;

    // State
    EquilibState   * sta;
    UnsatFlowState * fsta;

    // auxiliary methods
    size_t GetMax (size_t PrevMax, size_t NewMax)
    {
        if (NewMax>PrevMax) return NewMax;
        else                return PrevMax;
    }

private:
    // Auxiliary methods
    int  _RK_fun    (double t, double const Y[], double dYdt[]);
    void _RK_up_fun (double t, double Y[]);
};

#endif // HMSTRESSUPDATE_DECLARE


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


#ifdef HMSTRESSUPDATE_IMPLEMENT

Vec_t UnsatFlow::HMStressUpdate::Iv;

#include <mechsys/numerical/odesolver.h>

inline UnsatFlow::HMStressUpdate::HMStressUpdate ()
    : Mdl      (NULL),
      FMdl     (NULL),
      DbgFun   (NULL),
      DbgDat   (NULL),
      ODE      (NULL),
      STOL     (1.0e-5),
      dTini    (1.0),
      mMin     (0.1),
      mMax     (10.0),
      MaxSS    (2000),
      CDrift   (true),
      RKScheme ("ME"),
      SS       (0),
      SSs      (0),
      DCit     (0),
      DCitEl   (0),
      ncp      (0),
      niv      (0),
      neq      (0),
      sta      (NULL),
      fsta     (NULL)
{
}

inline void UnsatFlow::HMStressUpdate::SetModel (Model const * TheMdl, UnsatFlow const * TheFMdl)
{
    // model
    Mdl  = TheMdl;
    FMdl = TheFMdl;

    // constants
    if (Mdl!=NULL)
    {
        ncp = Mdl->NCps;
        niv = Mdl->NIvs;
        neq = ncp + niv + 1; // +1 => Sw

        dsig.change_dim (ncp);
        deps.change_dim (ncp);
        divs.change_dim (niv);

        if (size(Iv)==0)
        {
            Iv.change_dim(ncp);
            Iv(0) = 1.;
            Iv(1) = 1.;
            Iv(2) = 1.;
        }
    }
    else neq = 1; // Sw

    // ODE
    ODE = new Numerical::ODESolver<HMStressUpdate> (this, &HMStressUpdate::_RK_fun, neq, "ME", STOL, dTini);
    ODE->UpFun = &HMStressUpdate::_RK_up_fun;
    ODE->mMin  = mMin;
    ODE->mMax  = mMax;
    ODE->MaxSS = MaxSS;
}

inline void UnsatFlow::HMStressUpdate::Update (Vec_t const & DEps, double Dpw, State * Sta, UnsatFlowState * FSta, Vec_t & DSig)
{
    // check
    if (Mdl ==NULL) throw new Fatal("UnsatFlow::HMStressUpdate::Update: SetModel must be called before calling this method (equilib mdl is null)");
    if (FMdl==NULL) throw new Fatal("UnsatFlow::HMStressUpdate::Update: SetModel must be called before calling this method (unsatflow mdl is null)");

    // current state
    sta  = static_cast<EquilibState*>(Sta);
    fsta = FSta;
    DSig = sta->Sig; // temporary copy to calculate increment later

    // driving increments
    eps0 = sta->Eps;
    deps = DEps;
    pc0  = FSta->pc;
    n0   = FSta->n;
    dpc  = (-Dpw);
    dev  = Calc_ev(deps);

    // loading condition
    double aint;
    sta->Ldg = Mdl->LoadCond (sta, deps, aint);

    // initial state
    ODE->t = 0.0;
    for (size_t i=0; i<ncp; ++i) ODE->Y[    i] = sta->Sig(i);
    for (size_t i=0; i<niv; ++i) ODE->Y[ncp+i] = sta->Ivs(i);
    ODE->Y[ncp+niv] = FSta->Sw;

    // evolve
    ODE->Evolve (1.0);

    // return total stress increment
    DSig = sta->Sig - DSig;

    // debug
    if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
}

inline void UnsatFlow::HMStressUpdate::GetInfo (std::ostream & os, bool Header) const
{
    String buf;
    if (Header)
    {
        if (ODE==NULL)
        {
            os << "\n" << TERM_BLACK_WHITE << "----------------------------- HMStressUpdate/Scheme: RK ----------------------------" << TERM_RST << "\n\n";
            return;
        }
        os << "\n" << TERM_BLACK_WHITE << "----------------------------- HMStressUpdate/Scheme: RK(" << ODE->Scheme << ") ----------------------------" << TERM_RST << "\n\n";
        os << "STOL   = " << ODE->STOL   << std::endl;
        os << "dTini  = " << ODE->dTini  << std::endl;
        os << "mMin   = " << ODE->mMin   << std::endl;
        os << "mMax   = " << ODE->mMax   << std::endl;
        os << "MaxSS  = " << ODE->MaxSS  << std::endl;
        os << "CDrift = " << CDrift      << std::endl;
        buf.Printf("\n%6s %6s %6s %6s %12s %12s\n", "Scheme", "SS", "SSs", "DCit", "T", "dT");
        os << buf;
    }
    if (ODE==NULL) return;
    buf.Printf("%6s %6zd %6zd %6zd %12.8f %12.8f\n", "ME", ODE->SS, ODE->SSs, DCit, ODE->T, ODE->dT);
    os << buf;
}

inline int UnsatFlow::HMStressUpdate::_RK_fun (double t, double const Y[], double dYdt[])
{
    // current strain, stress, and internal values
    sta->Eps = eps0 + t*deps;
    fsta->pc = pc0  + t*dpc;
    fsta->n  = n0   + t*dev;
    for (size_t i=0; i<ncp; ++i) sta->Sig(i) = Y[    i];
    for (size_t i=0; i<niv; ++i) sta->Ivs(i) = Y[ncp+i];
    fsta->Sw = Y[ncp+niv];

    // dSw and dchi
    double dSw, dchi;
    FMdl->TgIncs (fsta, -dpc, dev, dSw, dchi);

    // constitutive stresses
    sta->Sig += (FMdl->chi*(-fsta->pc))*Iv;

    // tangent increments
    Mdl->TgIncs (sta, deps, dsig, divs); // the model expects effective stresses and returns effective stress increments

    // total stress increments
    dsig -= (FMdl->chi*(-dpc) + (-fsta->pc)*dchi)*Iv;

    for (size_t i=0; i<ncp; ++i) dYdt[    i] = dsig(i);
    for (size_t i=0; i<niv; ++i) dYdt[ncp+i] = divs(i);
    dYdt[ncp+niv] = dSw;

    // return success
    return GSL_SUCCESS;
}

inline void UnsatFlow::HMStressUpdate::_RK_up_fun (double t, double Y[])
{
    // current strain, stress, and internal values
    sta->Eps = eps0 + t*deps;
    fsta->pc = pc0  + t*dpc;
    fsta->n  = n0   + t*dev;
    for (size_t i=0; i<ncp; ++i) sta->Sig(i) = Y[    i];
    for (size_t i=0; i<niv; ++i) sta->Ivs(i) = Y[ncp+i];
    fsta->Sw = Y[ncp+niv];

    // correct drift
    if (CDrift)
    {
        DCit = GetMax (DCit, Mdl->CorrectDrift (sta));
        for (size_t i=0; i<ncp; ++i) Y[    i] = sta->Sig(i);
        for (size_t i=0; i<niv; ++i) Y[ncp+i] = sta->Ivs(i);
    }

    // debug function
    if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
}

#endif // HMSTRESSUPDATE_IMPLEMENT
