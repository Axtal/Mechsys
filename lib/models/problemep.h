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

#ifndef MECHSYS_PROBLEMEP_H
#define MECHSYS_PROBLEMEP_H

// MechSys
#include <mechsys/models/model.h>
#include <mechsys/models/equilibstate.h>
#include <mechsys/numerical/root.h>

class ProblemEP : public Model
{
public:
    // enums
    enum FCrit_t { VM_t, DP_t, MC_t, MN_t }; ///< Failure criterion type

    // Constructor & Destructor
    ProblemEP (int NDim, SDPair const & Prms, bool DerivedModel=false);
    virtual ~ProblemEP () {}

    // Derived methods
    virtual void   TgIncs       (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const;
    virtual void   Stiffness    (State const * Sta, Mat_t & D)                                const;
    virtual size_t CorrectDrift (State       * Sta)                                           const;
    virtual bool   LoadCond     (State const * Sta, Vec_t const & DEps, double & alpInt)      const;

    // Internal methods to be overloaded by derived classes
    virtual void   InitIvs   (SDPair const & Ini, State * Sta)      const;
    virtual void   Gradients (Vec_t const & Sig, Vec_t const & Ivs) const;
    virtual void   FlowRule  (Vec_t const & Sig, Vec_t const & Ivs) const;
    virtual void   Hardening (Vec_t const & Sig, Vec_t const & Ivs) const;
    virtual double YieldFunc (Vec_t const & Sig, Vec_t const & Ivs) const;
    virtual double CalcE     (Vec_t const & Sig, Vec_t const & Ivs) const { return E; }
    virtual void   ELStiff   (Vec_t const & Sig, Vec_t const & Ivs) const;

    // Constants
    bool    Derived;    ///< Derived model (such as CamClay)
    double  E;          ///< Young
    double  nu;         ///< Poisson
    FCrit_t FC;         ///< Failure criterion: VM:Von-Mises
    double  kVM;        ///< von Mises coefficient
    double  kDP;        ///< Drucker-Prager coefficient
    double  kMN;        ///< Matsuoka-Nakai coefficient
    double  Hb;         ///< Hardening coefficient H_bar
    bool    NonAssoc;   ///< Non-associated flow rule ?
    double  sphi;       ///< Sin(phi) friction angle
    double  spsi;       ///< Sin(psi) dilatancy angle
    double  FTol;       ///< Tolerance to be used when finding the intersection
    double  DCFTol;     ///< Drift correction ys function tolerance
    size_t  DCMaxIt;    ///< Drift correction max iterations
    double  pTol;       ///< Tolerance for minimum poct
    double  qTol;       ///< Tolerance for minimum qoct
    bool    NewSU;      ///< New stress update ?
    double  BetSU;      ///< Beta coefficient for new stress update
    double  AlpSU;      ///< Alpha coefficient for new stress update
    Vec_t   I;          ///< Idendity tensor

    // State data (mutable/scratch-pad)
    mutable Vec_t Sig0, Ivs0, SigA, DSigTr;    ///< Variables for yield crossing detection
    mutable Vec_t DEpsEl;                      ///< Elastic strain increment
    mutable Vec_t DEpsPl;                      ///< Plastic strain increment == Lam*W
    mutable Vec_t V;                           ///< NCps: Gradient of the yield surface
    mutable Vec_t W, devW;                     ///< NCps: Plastic flow rule direction
    mutable Vec_t Y;                           ///< NIvs: Derivative of the yield surface w.r.t internal variables
    mutable Vec_t H;                           ///< NIvs: Hardening coefficients, one for each internal variable
    mutable Mat_t De;                          ///< Elastic stiffness
    mutable Mat_t Dep;                         ///< Elastoplastic stiffness
    mutable Vec_t VDe;                         ///< V*De
    mutable Vec_t DeW;                         ///< De*W
    mutable Vec_t s;                           ///< Deviator of sigma
    mutable Vec_t dthdsig;                     ///< Derivative of theta w.r.t sigma
    mutable Vec_t dgdsig;                      ///< derivative of g w.r.t sigma
    mutable Vec_t dI1dsig,dI2dsig,dI3dsig;     ///< derivative of characteristic invariants
    mutable double p, q, t, th, g, I1, I2, I3; ///< Invariants

    // Methods for yield surface crossing
    double Fa    (double Alp, void*) { SigA=Sig0+Alp*DSigTr;  return YieldFunc(SigA,Ivs0); }
    double dFada (double Alp, void*) { SigA=Sig0+Alp*DSigTr;  Gradients(SigA,Ivs0);  return dot(V,DSigTr); }

    // Auxiliary methods
    void Calc_pq     (Vec_t const & Sig) const { p=Calc_poct(Sig);  q=Calc_qoct(Sig); }
    void Calc_pqg    (Vec_t const & Sig) const { OctInvs(Sig,p,q,t);  th=asin(t)/3.0;  g=Util::SQ2*sphi/(Util::SQ3*cos(th)-sphi*sin(th)); }
    void Calc_dgdsig (Vec_t const & Sig, bool Potential=false) const
    {
        double sinp = (Potential ? spsi : sphi);
        OctInvs (Sig, p,q,t,th,s, qTol, &dthdsig);
        g = Util::SQ2*sinp/(Util::SQ3*cos(th)-sinp*sin(th));
        double dgdth = g*(Util::SQ3*sin(th)+sinp*cos(th))/(Util::SQ3*cos(th)-sinp*sin(th));
        dgdsig = dgdth * dthdsig;
    }
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline ProblemEP::ProblemEP (int NDim, SDPair const & Prms, bool Deriv)
    : Model (NDim,Prms,/*niv*/3,"ProblemEP"), Derived(Deriv),
      E(0.0), nu(0.0), FC(VM_t), kVM(0.0), Hb(0.0), NonAssoc(false), FTol(1.0e-7), DCFTol(1.0e-8), DCMaxIt(10),
      pTol(1.0e-6), qTol(1.0e-7), NewSU(false)
{
    // resize scratchpad arrays
    Sig0   .change_dim (NCps);
    Ivs0   .change_dim (NIvs);
    SigA   .change_dim (NCps);
    DSigTr .change_dim (NCps);
    DEpsEl .change_dim (NCps);
    V      .change_dim (NCps);
    W      .change_dim (NCps);
    devW   .change_dim (NCps);
    De     .change_dim (NCps,NCps);
    Dep    .change_dim (NCps,NCps);
    VDe    .change_dim (NCps);
    DeW    .change_dim (NCps);
    s      .change_dim (NCps);
    dthdsig.change_dim (NCps);
    dgdsig .change_dim (NCps);
    dI1dsig.change_dim (NCps);
    dI2dsig.change_dim (NCps);
    dI3dsig.change_dim (NCps);
    I      .change_dim (NCps);
    set_to_zero(I);
    I(0) = 1.0;
    I(1) = 1.0;
    I(2) = 1.0;

    if (!Derived) // for instance, CamClay
    {
        // parameters
        if (!Prms.HasKey("E"))  throw new Fatal("ProblemEP::ProblemEP: Young modulus (E) must be provided");
        if (!Prms.HasKey("nu")) throw new Fatal("ProblemEP::ProblemEP: Poisson coefficient (nu) must be provided");
        E  = Prms("E");
        nu = Prms("nu");
        if (Prms.HasKey("DP")) FC = DP_t;
        if (Prms.HasKey("MC")) FC = MC_t;
        if (Prms.HasKey("MN")) FC = MN_t;
        if (FC==VM_t)
        {
            if      (Prms.HasKey("sY")) kVM = sqrt(2.0/3.0)*Prms("sY");
            else if (Prms.HasKey("c"))  kVM = (GTy==psa_t ? sqrt(2.0)*Prms("c") : 2.0*sqrt(2.0/3.0)*Prms("c"));
            else throw new Fatal("ProblemEP::ProblemEP: With VM (von Mises), either sY (uniaxial yield stress) or c (undrained cohesion) must be provided");
        }
        else
        {
            if (!Prms.HasKey("phi")) throw new Fatal("ProblemEP::ProblemEP: friction angle phi (degrees) must be provided");
            double c       = (Prms.HasKey("c") ? Prms("c") : 0.0);
            double phi_deg = Prms("phi");
            double phi_rad = phi_deg*Util::PI/180.0;
            double psi_rad = 0.0;
            if (phi_deg<1.0e-3) throw new Fatal("ProblemEP::ProblemEP: Friction angle (phi [deg]) must be greater than zero (1.0e-3). phi=%g is invalid",phi_deg);
            if (c<0.0)          throw new Fatal("ProblemEP::ProblemEP: 'cohesion' must be greater than zero. c=%g is invalid",c);
            if (Prms.HasKey("psi"))
            {
                NonAssoc = true;
                psi_rad  = Prms("psi")*Util::PI/180.0;
            }
            sphi = sin(phi_rad);
            spsi = sin(psi_rad);
            if (FC==DP_t) kDP = 2.0*sqrt(2.0)*sphi/(3.0-sphi);
            if (FC==MN_t) kMN = 9.0+8.0*pow(tan(phi_rad),2.0);
        }

        // hardening
        if (Prms.HasKey("Hp")) Hb = (2.0/3.0)*Prms("Hp"); // Hp=H_prime

        // internal values
        NIvs = 3;
        Y.change_dim (NIvs);
        H.change_dim (NIvs);
        IvNames.Push ("z0");
        IvNames.Push ("evp");
        IvNames.Push ("edp");

        // name
        switch (FC)
        {
            case VM_t: { Name = "ProblemEP(VM)"; break; }
            case DP_t: { Name = "ProblemEP(DP)"; break; }
            case MC_t: { Name = "ProblemEP(MC)"; break; }
            case MN_t: { Name = "ProblemEP(MN)"; break; }
        }

        // set model in stress update
        SUp.SetModel (this);
    }

    // new stress update parmeters
    NewSU = (Prms.HasKey("newsu") ? (int)Prms("newsu") : false);
    if (NewSU)
    {
        BetSU = Prms("betsu");
        AlpSU = Prms("alpsu");
    }
}

inline void ProblemEP::TgIncs (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // zero internal values (if any)
    DIvs.change_dim (NIvs);
    set_to_zero     (DIvs);

    // De: elastic stiffness
    ELStiff (sta->Sig, sta->Ivs);

    // increments
    if (sta->Ldg)
    {
        // limit tensile strength
        /*
        if (FC!=VM_t)
        {
            p = Calc_poct (sta->Sig);
            if (p<pTol)
            {
                printf("ProblemEP::TgIncs: __tensile_loading__\n");
                if (!Derived) // plastic strains
                {
                    DIvs(1) = Calc_ev (DEps); // devp
                    DIvs(2) = Calc_ed (DEps); // dedp
                }
                set_to_zero (DSig);
                return;
            }
        }
        */

        // gradients, flow rule, hardening, and hp
        Gradients (sta->Sig, sta->Ivs);
        FlowRule  (sta->Sig, sta->Ivs);
        Hardening (sta->Sig, sta->Ivs);
        double hp = (NIvs>0 ? Y(0)*H(0) : 0.0);

        // plastic multiplier
        Mult (V, De, VDe);
        double phi = dot(VDe,W) - hp;
        double Lam = dot(VDe,DEps)/phi;

        // increments
        DEpsPl = Lam*W;
        DEpsEl = DEps - DEpsPl;
        DSig   = De*DEpsEl;

        // increment of internal values
        DIvs = Lam*H;

        // plastic strains
        if (!Derived)
        {
            DIvs(1) = Calc_ev (DEpsPl); // devp
            DIvs(2) = Calc_ed (DEpsPl); // dedp
        }
    }
    else
    {
        // stress increment
        DSig = De*DEps;

        // new stress update
        if (NewSU) DIvs(0) = -dot(V, DSig) / Y(0);
    }

    // correct strain increment for plane stress
    if (GTy==pse_t) DEps(2) = -nu*(DSig(0)+DSig(1))/CalcE(sta->Sig, sta->Ivs);
}

inline void ProblemEP::Stiffness (State const * Sta, Mat_t & D) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // De: elastic stiffness
    ELStiff (sta->Sig, sta->Ivs);

    // stiffness
    if (sta->Ldg)
    {
        // limit tensile strength
        /*
        if (FC!=VM_t)
        {
            p = Calc_poct (sta->Sig);
            if (p<pTol)
            {
                D.change_dim (NCps,NCps);
                set_to_zero  (D);
                return;
            }
        }
        */

        // gradients, flow rule, hardening, and hp
        Gradients (sta->Sig, sta->Ivs);
        FlowRule  (sta->Sig, sta->Ivs);
        Hardening (sta->Sig, sta->Ivs);
        double hp = (NIvs>0 ? Y(0)*H(0) : 0.0);

        // auxiliary vectors
        Mult (V, De, VDe);
        double phi = dot(VDe,W) - hp;
        DeW = De*W;

        // elastoplastic stiffness
        D.change_dim (NCps, NCps);
        for (size_t i=0; i<NCps; ++i)
        for (size_t j=0; j<NCps; ++j)
            D(i,j) = De(i,j) - DeW(i)*VDe(j)/phi;
    }
    else D = De;

    // plane stress
    if (GTy==pse_t)
    {
        for (size_t i=0; i<NCps; ++i)
        {
            D(2,i) = 0.0;
            D(i,2) = 0.0;
        }
    }
}

inline size_t ProblemEP::CorrectDrift (State * Sta) const
{
    // state
    EquilibState * sta = static_cast<EquilibState *>(Sta);

    // limit tensile strength
    /*
    if (FC!=VM_t)
    {
        p = Calc_poct (sta->Sig);
        if (p<pTol) { p = pTol; return 666; }
    }
    */

    // iterations
    double fnew = YieldFunc (sta->Sig, sta->Ivs);
    size_t it   = 0;
    while (fnew>DCFTol && it<DCMaxIt)
    {
        // gradients, flow rule, hardening, and hp
        Gradients (sta->Sig, sta->Ivs);
        FlowRule  (sta->Sig, sta->Ivs);
        Hardening (sta->Sig, sta->Ivs);
        double hp = (NIvs>0 ? Y(0)*H(0) : 0.0);

        // elastic stiffness
        if (it==0) ELStiff (sta->Sig, sta->Ivs);

        // auxiliary vectors
        Mult (V, De, VDe);
        double dgam = fnew/(dot(VDe,W)-hp);
        DeW = De*W;

        // update stress and ivs (only)
        sta->Sig -= dgam*DeW;
        sta->Ivs += dgam*H;
        fnew = YieldFunc (sta->Sig, sta->Ivs);

        // check convergence
        if (fabs(fnew)<DCFTol) break;
        it++;
    }

    // check number of iterations
    if (it>=DCMaxIt) throw new Fatal("ProblemEP::CorrectDrift: Yield surface drift correction did not converge after %d iterations (fnew=%g, DCFTol=%g)",it,fnew,DCFTol);
    return it;
}

inline bool ProblemEP::LoadCond (State const * Sta, Vec_t const & DEps, double & alpInt) const
{
    // default return values
    alpInt   = -1.0;  // no intersection
    bool ldg = false; // => unloading (or crossing)

    // current state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // trial state
    ELStiff (sta->Sig, sta->Ivs);
    DSigTr = De * DEps;
    SigA   = sta->Sig + DSigTr;

    // yield function values
    double f    = YieldFunc (sta->Sig, sta->Ivs);
    double f_tr = YieldFunc (SigA,     sta->Ivs);

    // numerator of Lagrange multiplier
    Gradients (sta->Sig, sta->Ivs);
    double numL = dot(V, DSigTr);
    printf("f=%g,  ftr=%g,  numL=%g\n",f,f_tr,numL);

    // new stress update
    if (NewSU)
    {
        q = Calc_qoct (sta->Sig);
        if (q>qTol)
        {
            if (numL>0.0) ldg = true;
        }
        else
        {
            double qf = Calc_qoct (SigA);
            double Dq = qf - q;
            if (Dq>0.0) ldg = true;
        }
        printf(">>>>>>> (new su) >>>>>>>>>>> f=%g,  ftr=%g,  ldg=%d\n",f,f_tr,ldg);
        return ldg;
    }

    // going outside
    if (f_tr>0.0)
    {
        bool crossing = false;
        if (f<-FTol && f_tr>FTol) // (f<0.0) does not work, (f<-FTol) does not work
        {
            Sig0     = sta->Sig;
            Ivs0     = sta->Ivs;
            crossing = true;
        }
        else if (numL<-qTol) // (numL<0.0) does not work, since it may be equal to -0.0 => neutral loading
        {
            // moving stress state slighlty to the inside of the yield surface in order to have f=-FTol
            double df   = -FTol - f;
            double dalp = df / numL; // numL = dfdalpha
            Sig0        = sta->Sig + dalp*DSigTr;
            Ivs0        = sta->Ivs;
            crossing    = true;
            double fnew = YieldFunc (Sig0, Ivs0);
            if (fnew>0.0) throw new Fatal("ProblemEP::LoadCond: __internal_error__: correction for numL<0 failed (dalp=%g, fnew=%g)",dalp,fnew);
            //throw new Fatal("ProblemEP::LoadCond: Strain increment is too large (f=%g, f_tr=%g, numL=%g). Crossing and going all the way through the yield surface to the other side.",f,f_tr,numL);
        }
        if (crossing)
        {
            Numerical::Root<ProblemEP> root(const_cast<ProblemEP*>(this), &ProblemEP::Fa, &ProblemEP::dFada);
            //root.Scheme = "Newton";
            //root.Verbose = true;
            alpInt = root.Solve (0.0, 1.0);
            if (alpInt<0) throw new Fatal("ProblemEP::LoadCond: alpInt=%g must be positive",alpInt);
        }
        else ldg = true;
    }

    // return true if there is loading. returns false (unloading) with intersection
    return ldg;
}

inline void ProblemEP::InitIvs (SDPair const & Ini, State * Sta) const
{
    // initialize state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini, NIvs);

    // internal variables
    sta->Ivs(0) = 1.0; // z0
    sta->Ivs(1) = 0.0; // evp
    sta->Ivs(2) = 0.0; // edp

    // new stress update
    if (NewSU)
    {
        switch (FC)
        {
            case VM_t:
            {
                q = Calc_qoct (sta->Sig);
                sta->Ivs(0) = q/kVM;
                break;
            }
            case DP_t:
            {
                Calc_pq (sta->Sig);
                sta->Ivs(0) = q/(p*kDP);
                break;
            }
            case MC_t:
            {
                Calc_pqg (sta->Sig);
                sta->Ivs(0) = q/(p*g);
                break;
            }
            case MN_t:
            {
                CharInvs (sta->Sig, I1,I2,I3);
                sta->Ivs(0) = I1*I2/(I3*kMN);
                break;
            }
        }
    }

    // check initial yield function
    double f = YieldFunc (sta->Sig, sta->Ivs);
    if (f>FTol)           throw new Fatal("ProblemEP:InitIvs: stress point (sig=(%g,%g,%g,%g]) is outside yield surface (f=%g) with z0=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,f,sta->Ivs(0));
    if (NewSU && f<-FTol) throw new Fatal("ProblemEP:InitIvs: stress point (sig=(%g,%g,%g,%g]) is outside yield surface (f=%g) with z0=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,f,sta->Ivs(0));
}

inline void ProblemEP::Gradients (Vec_t const & Sig, Vec_t const & Ivs) const
{
    // derivative of internal values
    Y(0) = 0.0; // dfdz0
    Y(1) = 0.0; // dfdz1
    Y(2) = 0.0; // dfdz2

    // new stress update
    if (NewSU) Y(0) = -1.0;

    // limit tensile strength
    /*
    if (FC!=VM_t)
    {
        p = Calc_poct (Sig);
        if (p<pTol) { V = I; return; }
    }
    */

    switch (FC)
    {
        case VM_t:
        {
            OctInvs (Sig, p,q,s, qTol);
            V = s/(q*kVM);
            break;
        }
        case DP_t:
        {
            OctInvs (Sig, p,q,s, qTol);
            V = (1.0/(p*q*kDP))*s + (q/(p*p*kDP*Util::SQ3))*I;
            break;
        }
        case MC_t:
        {
            Calc_dgdsig (Sig);
            V = (q/(g*p*p*Util::SQ3))*I + (1.0/(p*q*g))*s - (q/(p*g*g))*dgdsig;
            break;
        }
        case MN_t:
        {
            CharInvs (Sig, I1,I2,I3, dI1dsig,dI2dsig,dI3dsig);
            V = (I2/(I3*kMN))*dI1dsig + (I1/(I3*kMN))*dI2dsig - (I1*I2/(I3*I3*kMN))*dI3dsig;
            break;
        }
    }
}

inline void ProblemEP::FlowRule (Vec_t const & Sig, Vec_t const & Ivs) const
{
    switch (FC)
    {
        case VM_t: { W = V; break; }
        case DP_t: { W = V; break; }
        case MC_t:
        {
            if (NonAssoc)
            {
                Calc_dgdsig (Sig, true);
                //W = (q/(g*p*p*Util::SQ3))*I + (1.0/(p*q*g))*s - (q/(p*g*g))*dgdsig;
                W = (g/Util::SQ3)*I + (1.0/q)*s - p*dgdsig;
            }
            else W = V;
            break;
        }
        case MN_t: { W = V; break; }
    }
}

inline void ProblemEP::Hardening (Vec_t const & Sig, Vec_t const & Ivs) const
{
    // internal values
    H(0) = 0.0;
    H(1) = 0.0;
    H(2) = 0.0;

    // new stress update
    if (NewSU)
    {
        double F=0;
        switch (FC)
        {
            case VM_t:
            {
                q = Calc_qoct (Sig);
                F = q/kVM - 1.0;
                break;
            }
            case DP_t:
            {
                Calc_pq (Sig);
                F = q/(p*kDP) - 1.0;
                break;
            }
            case MC_t:
            {
                Calc_pqg (Sig);
                F = q/(p*g) - 1.0;
                break;
            }
            case MN_t:
            {
                CharInvs (Sig, I1,I2,I3);
                F = I1*I2/(I3*kMN) - 1.0;
                break;
            }
        }
        double H1 = 0.0;
        if (F>0.0) F = 0.0;
        H(0) = AlpSU + (H1-AlpSU)*exp(BetSU*F);
    }
}

inline double ProblemEP::YieldFunc (Vec_t const & Sig, Vec_t const & Ivs) const
{
    // limit tensile strength
    /*
    if (FC!=VM_t)
    {
        p = Calc_poct (Sig);
        if (p<pTol) { return FTol; }
    }
    */

    double f=0;
    switch (FC)
    {
        case VM_t:
        {
            q = Calc_qoct (Sig);
            f = q/kVM - Ivs(0);
            break;
        }
        case DP_t:
        {
            Calc_pq (Sig);
            f = q/(p*kDP) - Ivs(0);
            break;
        }
        case MC_t:
        {
            Calc_pqg (Sig);
            f = q/(p*g) - Ivs(0);
            break;
        }
        case MN_t:
        {
            CharInvs (Sig, I1,I2,I3);
            f = I1*I2/(I3*kMN) - Ivs(0);
            break;
        }
    }
    return f;
}

inline void ProblemEP::ELStiff (Vec_t const & Sig, Vec_t const & Ivs) const
{
    if (NDim==2)
    {
        if (GTy==pse_t)
        {
            double c = CalcE(Sig,Ivs)/(1.0-nu*nu);
            De = c,    c*nu, 0.0,        0.0,
                 c*nu, c,    0.0,        0.0,
                 0.0,  0.0,  0.0,        0.0,
                 0.0,  0.0,  0.0, c*(1.0-nu);
        }
        else if (GTy==psa_t || GTy==axs_t)
        {
            double c = CalcE(Sig,Ivs)/((1.0+nu)*(1.0-2.0*nu));
            De = c*(1.0-nu),       c*nu ,      c*nu ,            0.0,
                      c*nu ,  c*(1.0-nu),      c*nu ,            0.0,
                      c*nu ,       c*nu , c*(1.0-nu),            0.0,
                       0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu);
        }
        else throw new Fatal("ProblemEP::Stiffness: 2D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
    }
    else
    {
        if (GTy==d3d_t)
        {
            double c = CalcE(Sig,Ivs)/((1.0+nu)*(1.0-2.0*nu));
            De = c*(1.0-nu),       c*nu ,      c*nu ,            0.0,            0.0,            0.0,
                      c*nu ,  c*(1.0-nu),      c*nu ,            0.0,            0.0,            0.0,
                      c*nu ,       c*nu , c*(1.0-nu),            0.0,            0.0,            0.0,
                       0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu),            0.0,            0.0,
                       0.0 ,        0.0 ,       0.0 ,            0.0, c*(1.0-2.0*nu),            0.0,
                       0.0 ,        0.0 ,       0.0 ,            0.0,            0.0, c*(1.0-2.0*nu);
        }
        else throw new Fatal("ProblemEP::Stiffness: 3D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
    }
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * ProblemEPMaker(int NDim, SDPair const & Prms, Model const * O) { return new ProblemEP(NDim,Prms); }

int ProblemEPRegister()
{
    ModelFactory   ["ProblemEP"] = ProblemEPMaker;
    MODEL.Set      ("ProblemEP", (double)MODEL.Keys.Size());
    MODEL_PRM_NAMES["ProblemEP"].Resize (14);
    MODEL_PRM_NAMES["ProblemEP"] = "E", "nu", "sY", "c", "phi", "Hp", "psi", "VM", "DP", "MC", "MN", "newsu", "betsu", "alpsu";
    MODEL_IVS_NAMES["ProblemEP"].Resize (0);
    return 0;
}

int __ProblemEP_dummy_int = ProblemEPRegister();


#endif // MECHSYS_PROBLEMEP_H
