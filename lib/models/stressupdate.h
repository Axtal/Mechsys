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

#ifdef STRESSUPDATE_DECLARE

class StressUpdate
{
public:
    // callbacks
    typedef void (*pDbgFun) (StressUpdate const & SU, void * UserData); ///< Pointer to debug function

    // enum
    enum Scheme_t { ME_t, SingleFE_t, RK_t }; ///< Integration scheme

    // Constructor & Destructor
     StressUpdate ();
    ~StressUpdate () { if (ODE!=NULL) delete ODE;  if (sta_1!=NULL) delete sta_1;  if (sta_ME!=NULL) delete sta_ME; }

    // Methods
    void SetModel  (Model const * TheMdl);
    void SetScheme (String const & Name);
    void Update    (Vec_t const & DEps, State * Sta, Vec_t & DSig);
    void GetInfo   (std::ostream & os, bool Header=false) const;

    // Data
    Model const * Mdl;
    pDbgFun       DbgFun;
    void        * DbgDat;

    // ODE solver (if sheme==RK)
    Numerical::ODESolver<StressUpdate> * ODE;

    // Constants for integration
    Scheme_t       Scheme; ///< Scheme: ME_t (Modified-Euler)
    double         STOL;
    double         dTini;
    double         mMin;
    double         mMax;
    size_t         MaxSS;
    bool           CDrift; ///< correct drift ?
    double         Error;
    String         RKScheme;
    double         T;
    double         dT;
    size_t         SS;     ///< num of substeps
    size_t         SSs;    ///< num of successful substeps
    size_t         DCit;   ///< drift correction iterations
    size_t         DCitEl; ///< drift correction (elastic)
    size_t         ncp;    ///< num of stress components
    size_t         niv;    ///< num of internal variables
    size_t         neq;    ///< ncp + niv
    Vec_t          dsig;
    Vec_t          deps;
    Vec_t          divs;
    Vec_t          eps0;
    EquilibState * sta;

    // auxiliary variables
    EquilibState * sta_1;         // intermediate state
    EquilibState * sta_ME;        // Modified-Euler state
    Vec_t deps_1, dsig_1, divs_1; // intermediate increments
    Vec_t dsig_2, divs_2;         // ME increments
    Vec_t sig_dif;                // ME - FE stress difference

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

#endif // STRESSUPDATE_DECLARE


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


#ifdef STRESSUPDATE_IMPLEMENT

#include <mechsys/numerical/odesolver.h>

inline Model::StressUpdate::StressUpdate ()
    : Mdl      (NULL),
      DbgFun   (NULL),
      DbgDat   (NULL),
      ODE      (NULL),
      Scheme   (ME_t),
      STOL     (1.0e-5),
      dTini    (1.0),
      mMin     (0.1),
      mMax     (10.0),
      MaxSS    (2000),
      CDrift   (true),
      Error    (0.0),
      RKScheme ("ME"),
      T        (0.0),
      dT       (dTini),
      SS       (0),
      SSs      (0),
      DCit     (0),
      DCitEl   (0),
      sta      (NULL)
{
}

inline void Model::StressUpdate::SetModel (Model const * TheMdl)
{
    Mdl = TheMdl;
    ncp = Mdl->NCps;
    niv = Mdl->NIvs;
    neq = ncp + niv;
    ODE = new Numerical::ODESolver<StressUpdate> (this, &StressUpdate::_RK_fun, neq, "ME", STOL, dTini);
    ODE->UpFun = &StressUpdate::_RK_up_fun;
    ODE->mMin  = mMin;
    ODE->mMax  = mMax;

    // auxiliary variables
    dsig.change_dim (ncp);
    deps.change_dim (ncp);
    divs.change_dim (niv);

    // ME variables
    sta_1  = new EquilibState (Mdl->NDim);
    sta_ME = new EquilibState (Mdl->NDim);
    deps_1.change_dim(ncp);  dsig_1.change_dim(ncp);  divs_1.change_dim(niv);
    dsig_2.change_dim(ncp);  divs_2.change_dim(niv);
    sig_dif.change_dim(ncp);
}

inline void Model::StressUpdate::SetScheme (String const & Name)
{
    if      (Name=="ME")       Scheme = ME_t;
    else if (Name=="SingleFE") Scheme = SingleFE_t;
    else if (Name=="RK")       Scheme = RK_t;
    else throw new Fatal("StressUpdate::SetScheme: Scheme named %s is invalid",Name.CStr());
}

inline void Model::StressUpdate::Update (Vec_t const & DEps, State * Sta, Vec_t & DSig)
{
    // check
    if (Mdl==NULL) throw new Fatal("Model::StressUpdate::Update: SetModel must be called before calling this method");

    // scheme
    if (Mdl->Prms.HasKey("newsu")) { if ((int)Mdl->Prms("newsu")) Scheme = RK_t; }

    // current state
    sta  = static_cast<EquilibState*>(Sta);
    DSig = sta->Sig; // temporary copy to calculate increment later

    if (Scheme==SingleFE_t) // without intersection detection (should be used for linear elasticity only)
    {
        deps = DEps;
        Mdl->TgIncs (sta, deps, dsig, divs);
        sta->Eps += deps;
        sta->Sig += dsig;
        sta->Ivs += divs;
        if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
    }
    else if (Scheme==ME_t)
    {
        // loading-unloading ?
        double aint = -1.0; // no intersection
        bool   ldg  = Mdl->LoadCond (sta, DEps, aint);

        // set loading flag
        sta    -> Ldg = ldg;
        sta_1  -> Ldg = ldg;
        sta_ME -> Ldg = ldg;

        // with intersection ?
        if (aint>0.0 && aint<1.0)
        {
            // update to intersection
            deps = aint*DEps;
            Mdl->TgIncs (sta, deps, dsig, divs);
            sta->Eps += deps;
            sta->Sig += dsig;
            sta->Ivs += divs;
            deps = fabs(1.0-aint)*DEps; // remaining of DEps to be applied

            // change loading flag
            ldg           = true;
            sta    -> Ldg = ldg;
            sta_1  -> Ldg = ldg;
            sta_ME -> Ldg = ldg;

            // drift correction
            //if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
            if (CDrift) DCitEl = GetMax (DCitEl, Mdl->CorrectDrift(sta));

            // debug
            if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
        }
        else deps = DEps; // update with full DEps

        // for each pseudo time T
        T      = 0.0;
        dT     = dTini;
        SS     = 0;
        SSs    = 0;
        DCitEl = 0;
        DCit   = 0;
        for (SS=0; SS<MaxSS; ++SS)
        {
            // exit point
            if (T>=1.0) break;

            // FE and ME increments
            deps_1 = dT*deps;
            Mdl->TgIncs (sta, deps_1, dsig_1, divs_1);
            sta_1->Eps  = sta->Eps + deps_1;
            sta_1->Sig  = sta->Sig + dsig_1;
            sta_1->Ivs  = sta->Ivs + divs_1;
            Mdl->TgIncs (sta_1, deps_1, dsig_2, divs_2);
            sta_ME->Sig = sta->Sig + 0.5*(dsig_1+dsig_2);
            sta_ME->Ivs = sta->Ivs + 0.5*(divs_1+divs_2);

            // local error estimate
            sig_dif = sta_ME->Sig - sta_1->Sig;
            double sig_err = Norm(sig_dif)/(1.0+Norm(sta_ME->Sig));
            double ivs_err = 0.0;
            for (size_t i=0; i<niv; ++i) ivs_err += fabs(sta_ME->Ivs(i)-sta_1->Ivs(i))/(1.0+fabs(sta_ME->Ivs(i)));
            Error = sig_err + ivs_err;

            // step multiplier
            double m = (Error>0.0 ? 0.9*sqrt(STOL/Error) : mMax);

            // update
            if (Error<STOL)
            {
                // update state
                SSs++;
                T += dT;
                sta->Eps = sta_1  -> Eps;
                sta->Sig = sta_ME -> Sig;
                sta->Ivs = sta_ME -> Ivs;

                // drift correction
                //if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
                if (CDrift) DCit = GetMax (DCit, Mdl->CorrectDrift (sta));

                // update stress path in model
                Mdl->UpdatePath (sta, deps_1, Vec_t(0.5*(dsig_1+dsig_2)));

                // limit change on stepsize
                if (m>mMax) m = mMax;

                // debug
                if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
            }
            else if (m<mMin) m = mMin;

            // change next step size
            dT = m * dT;

            // check for last increment
            if (dT>1.0-T) dT = 1.0-T;
        }
        if (SS>=MaxSS) throw new Fatal("StressUpdate::Update: Modified-Euler (local) did not converge after %d substeps",SS);
    }
    else if (Scheme==RK_t)
    {
        // initial strain and increment
        eps0 = sta->Eps;
        deps = DEps;

        // loading condition
        double aint;
        sta->Ldg = Mdl->LoadCond (sta, deps, aint);

        // initial state
        ODE->t = 0.0;
        for (size_t i=0; i<ncp; ++i) ODE->Y[    i] = sta->Sig(i);
        for (size_t i=0; i<niv; ++i) ODE->Y[ncp+i] = sta->Ivs(i);

        // evolve
        ODE->Evolve (1.0);

        // final state
        for (size_t i=0; i<ncp; ++i) sta->Sig(i) = ODE->Y[    i];
        for (size_t i=0; i<niv; ++i) sta->Ivs(i) = ODE->Y[ncp+i];
        sta->Eps = eps0 + deps;
    }
    else throw new Fatal("StressUpdate::Update: Scheme is not available yet");

    // return total stress increment
    DSig = sta->Sig - DSig;

    // debug
    if (DbgFun!=NULL) (*DbgFun) ((*this), DbgDat);
}

inline void Model::StressUpdate::GetInfo (std::ostream & os, bool Header) const
{
    String buf;
    if (Scheme==ME_t)
    {
        if (Header)
        {
            os << "\n" << TERM_BLACK_WHITE << "----------------------------- StressUpdate/Scheme: ME ------------------------------" << TERM_RST << "\n\n";
            os << "STOL   = " << STOL   << std::endl;
            os << "dTini  = " << dTini  << std::endl;
            os << "mMin   = " << mMin   << std::endl;
            os << "mMax   = " << mMax   << std::endl;
            os << "MaxSS  = " << MaxSS  << std::endl;
            os << "CDrift = " << CDrift << std::endl;
            buf.Printf("\n%6s %6s %6s %6s %6s %12s %12s %16s\n", "Scheme", "SS", "SSs", "DCitEl", "DCit", "T", "dT", "Error");
            os << buf;
        }
        buf.Printf("%6s %6zd %6zd %6zd %6zd %12.8f %12.8f %16.8e\n", "ME", SS, SSs, DCitEl, DCit, T, dT, Error);
        os << buf;
    }
    else if (Scheme==SingleFE_t)
    {
        if (Header) os << "\n" << TERM_BLACK_WHITE << "----------------------------- StressUpdate/Scheme: SingleFE ------------------------" << TERM_RST << "\n\n";
    }
    else if (Scheme==RK_t)
    {
        if (Header)
        {
            if (ODE==NULL)
            {
                os << "\n" << TERM_BLACK_WHITE << "----------------------------- StressUpdate/Scheme: RK ----------------------------" << TERM_RST << "\n\n";
                return;
            }
            os << "\n" << TERM_BLACK_WHITE << "----------------------------- StressUpdate/Scheme: RK(" << ODE->Scheme << ") ----------------------------" << TERM_RST << "\n\n";
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
}

inline int Model::StressUpdate::_RK_fun (double t, double const Y[], double dYdt[])
{
    // current strain, stress, and internal values
    sta->Eps = eps0 + t*deps;
    for (size_t i=0; i<ncp; ++i) sta->Sig(i) = Y[    i];
    for (size_t i=0; i<niv; ++i) sta->Ivs(i) = Y[ncp+i];

    // tangent increments
    Mdl->TgIncs (sta, deps, dsig, divs);
    for (size_t i=0; i<ncp; ++i) dYdt[    i] = dsig(i);
    for (size_t i=0; i<niv; ++i) dYdt[ncp+i] = divs(i);

    // return success
    return GSL_SUCCESS;
}

inline void Model::StressUpdate::_RK_up_fun (double t, double Y[])
{
    // current strain, stress, and internal values
    sta->Eps = eps0 + t*deps;
    for (size_t i=0; i<ncp; ++i) sta->Sig(i) = Y[    i];
    for (size_t i=0; i<niv; ++i) sta->Ivs(i) = Y[ncp+i];

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

#endif // STRESSUPDATE_IMPLEMENT
