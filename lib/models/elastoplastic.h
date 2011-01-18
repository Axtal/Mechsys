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

#ifndef MECHSYS_ELASTOPLASTIC_H
#define MECHSYS_ELASTOPLASTIC_H

// MechSys
#include <mechsys/models/model.h>
#include <mechsys/models/equilibstate.h>
#include <mechsys/numerical/root.h>

class ElastoPlastic : public Model
{
public:
    // struct
    struct AlphaData
    {
        AlphaData (int NDim, Mat_t const & De, EquilibState const & Sta, Vec_t const & DEps)
        {
            DSig_tr     = De * DEps;
            StaAlp      = new EquilibState (NDim);
            StaAlp->Ivs = Sta.Ivs;
            Sig         = Sta.Sig;
        }
        ~AlphaData () { delete StaAlp; }
        Vec_t          Sig;
        Vec_t          DSig_tr;
        EquilibState * StaAlp;
    };

    // enums
    enum FCrit_t { VM_t, DP_t, MC_t, MN_t }; ///< Failure criterion type

    // Constructor & Destructor
    ElastoPlastic (int NDim, SDPair const & Prms, bool DerivedModel=false);
    virtual ~ElastoPlastic () {}

    // Derived methods
    virtual void TgIncs       (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const;
    void         Stiffness    (State const * Sta, Mat_t & D)                                const;
    virtual bool LoadCond     (State const * Sta, Vec_t const & DEps, double & alpInt)      const;
    void         CorrectDrift (State       * Sta)                                           const;

    // Internal methods
    virtual void ELStiff (EquilibState const * Sta) const;

    // Internal methods to be overloaded by derived classes
    virtual void   InitIvs   (SDPair const & Ini, State * Sta) const;
    virtual void   Gradients (EquilibState const * Sta)        const;
    virtual void   FlowRule  (EquilibState const * Sta)        const;
    virtual void   Hardening (EquilibState const * Sta)        const;
    virtual double YieldFunc (EquilibState const * Sta)        const;
    virtual double CalcE     (EquilibState const * Sta)        const { return E; }

    // Constants
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
    double  CDFtol;     ///< Tolerance to be used in correct drift
    double  qTol;       ///< Tolerance for minium qoct
    bool    NewSU;      ///< New stress update ?
    double  BetSU;      ///< Beta coefficient for new stress update
    double  AlpSU;      ///< Alpha coefficient for new stress update
    Vec_t   I;          ///< Idendity tensor

    // State data (mutable/scratch-pad)
    mutable Vec_t V;       ///< NCps: Gradient of the yield surface
    mutable Vec_t W, devW; ///< NCps: Plastic flow rule direction
    mutable Vec_t Y;       ///< NIvs: Derivative of the yield surface w.r.t internal variables
    mutable Vec_t H;       ///< NIvs: Hardening coefficients, one for each internal variable
    mutable Mat_t De;      ///< Elastic stiffness
    mutable Mat_t Dep;     ///< Elastoplastic stiffness
    mutable Vec_t VDe;     ///< V*De
    mutable Vec_t DeW;     ///< De*W
    mutable Vec_t s;       ///< Deviator of sigma
    mutable Vec_t dthdsig; ///< Derivative of theta w.r.t sigma
    mutable Vec_t dgdsig;  ///< derivative of g w.r.t sigma
    mutable Vec_t dI1dsig,dI2dsig,dI3dsig; ///< derivative of characteristic invariants
    mutable double p, q, t, th, g, I1, I2, I3; ///< Invariants


    // Methods for yield surface crossing
    double Falpha  (double Alp, void * UserData);
    double dFalpha (double Alp, void * UserData);

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


inline ElastoPlastic::ElastoPlastic (int NDim, SDPair const & Prms, bool Derived)
    : Model (NDim,Prms,"ElastoPlastic"),
      E(0.0), nu(0.0), FC(VM_t), kVM(0.0), Hb(0.0), NonAssoc(false), FTol(1.0e-8), CDFtol(1.0e-9), qTol(1.0e-8), NewSU(false)
{
    // resize scratchpad arrays
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
        if (!Prms.HasKey("E"))  throw new Fatal("ElastoPlastic::ElastoPlastic: Young modulus (E) must be provided");
        if (!Prms.HasKey("nu")) throw new Fatal("ElastoPlastic::ElastoPlastic: Poisson coefficient (nu) must be provided");
        E  = Prms("E");
        nu = Prms("nu");
        if (Prms.HasKey("DP")) FC = DP_t;
        if (Prms.HasKey("MC")) FC = MC_t;
        if (Prms.HasKey("MN")) FC = MN_t;
        if (FC==VM_t)
        {
            if      (Prms.HasKey("sY")) kVM = sqrt(2.0/3.0)*Prms("sY");
            else if (Prms.HasKey("c"))  kVM = (GTy==psa_t ? sqrt(2.0)*Prms("c") : 2.0*sqrt(2.0/3.0)*Prms("c"));
            else throw new Fatal("ElastoPlastic::ElastoPlastic: With VM (von Mises), either sY (uniaxial yield stress) or c (undrained cohesion) must be provided");
        }
        else
        {
            if (!Prms.HasKey("phi")) throw new Fatal("ElastoPlastic::ElastoPlastic: friction angle phi (degrees) must be provided");
            double c       = (Prms.HasKey("c") ? Prms("c") : 0.0);
            double phi_deg = Prms("phi");
            double phi_rad = phi_deg*Util::PI/180.0;
            double psi_rad = 0.0;
            if (phi_deg<1.0e-3) throw new Fatal("ElastoPlastic::ElastoPlastic: Friction angle (phi [deg]) must be greater than zero (1.0e-3). phi=%g is invalid",phi_deg);
            if (c<0.0)          throw new Fatal("ElastoPlastic::ElastoPlastic: 'cohesion' must be greater than zero. c=%g is invalid",c);
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
            case VM_t: { Name = "ElastoPlastic(VM)"; break; }
            case DP_t: { Name = "ElastoPlastic(DP)"; break; }
            case MC_t: { Name = "ElastoPlastic(MC)"; break; }
            case MN_t: { Name = "ElastoPlastic(MN)"; break; }
        }
    }

    // new stress update parmeters
    NewSU = (Prms.HasKey("newsu") ? (int)Prms("newsu") : false);
    if (NewSU)
    {
        BetSU = Prms("betsu");
        AlpSU = Prms("alpsu");
    }
}

inline void ElastoPlastic::InitIvs (SDPair const & Ini, State * Sta) const
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
    double f = YieldFunc (sta);
    if (f>FTol)           throw new Fatal("ElastoPlastic:InitIvs: stress point (sig=(%g,%g,%g,%g]) is outside yield surface (f=%g) with z0=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,f,sta->Ivs(0));
    if (NewSU && f<-FTol) throw new Fatal("ElastoPlastic:InitIvs: stress point (sig=(%g,%g,%g,%g]) is outside yield surface (f=%g) with z0=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,f,sta->Ivs(0));
}

inline void ElastoPlastic::TgIncs (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // De: elastic stiffness
    ELStiff (sta);

    // increments
    DIvs.change_dim (NIvs);
    if (sta->Ldg)
    {
        // gradients, flow rule, hardening, and hp
        Gradients (sta);
        FlowRule  (sta);
        Hardening (sta);
        double hp = (NIvs>0 ? Y(0)*H(0) : 0.0);

        // plastic multiplier
        Mult (V, De, VDe);
        double phi = dot(VDe,W) - hp;
        double gam = dot(VDe,DEps)/phi;

        // stress increment
        Vec_t deps_elastic(DEps-gam*W);
        DSig = De*deps_elastic;

        // increment of internal values
        for (size_t i=0; i<NIvs; ++i) DIvs(i) = gam*H(i);
    }
    else
    {
        DSig = De*DEps;
        for (size_t i=0; i<NIvs; ++i) DIvs(i) = 0.0;

        // new stress update
        if (NewSU) DIvs(0) = -dot(V, DSig) / Y(0);
    }

    // correct strain increment for plane stress
    if (GTy==pse_t) DEps(2) = -nu*(DSig(0)+DSig(1))/CalcE(sta);
}

inline void ElastoPlastic::Stiffness (State const * Sta, Mat_t & D) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // De: elastic stiffness
    ELStiff (sta);

    // stiffness
    if (sta->Ldg)
    {
        // gradients, flow rule, hardening, and hp
        Gradients (sta);
        FlowRule  (sta);
        Hardening (sta);
        double hp = (NIvs>0 ? Y(0)*H(0) : 0.0);

        // auxiliar vectors
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

inline bool ElastoPlastic::LoadCond (State const * Sta, Vec_t const & DEps, double & alpInt) const
{
    // default return values
    alpInt = -1.0;    // no intersection
    bool ldg = false; // => unloading

    // current state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // elastic stiffness
    ELStiff (sta);

    // trial state
    Vec_t dsig_tr(De * DEps);
    EquilibState sta_tr(NDim);
    sta_tr = (*sta);
    sta_tr.Sig += dsig_tr;

    // yield function values
    double f    = YieldFunc (sta);
    double f_tr = YieldFunc (&sta_tr);

    // numerator of Lagrange multiplier
    Gradients (sta);
    double numL = dot(V, dsig_tr);

    // new stress update
    if (NewSU)
    {
        q = Calc_qoct (sta->Sig);
        if (q>qTol)
        {
            if (numL>0.0) ldg = true;
            //if (f_tr>0.0 && numL<0.0) throw new Fatal("ElastoPlastic::LoadCond (new update): Strain increment is too large (f=%g, f_tr=%g, numL=%g). Crossing and going all the way through the yield surface to the other side.",f,f_tr,numL);
            //if (f_tr>0.0 && numL<0.0) printf("ElastoPlastic::LoadCond (new update): Strain increment is too large (f=%g, f_tr=%g, numL=%g). Crossing and going all the way through the yield surface to the other side.\n",f,f_tr,numL);
        }
        else
        {
            double qf = Calc_qoct (sta_tr.Sig);
            double Dq = qf - q;
            if (Dq>0.0) ldg = true;
        }
        return ldg;
    }

    // going outside
    if (f_tr>0.0)
    {
        ldg = true;
        bool crossing = false;
        if (f<-FTol) crossing = true; // works
        //if (f<0.0) crossing = true; // does not work
        //else if (numL<0.0) // crossing to the other side
        //{
            //f = -1.0e-10;
            //crossing = true;
        //}
        if (crossing)
        {
            ldg = false;
            AlphaData dat(NDim, De, (*sta), DEps);
            Numerical::Root<ElastoPlastic> root(const_cast<ElastoPlastic*>(this), &ElastoPlastic::Falpha, &ElastoPlastic::dFalpha);
            //root.Scheme = "Newton";
            //root.Verbose = true;
            alpInt = root.Solve (0.0, 1.0, NULL, &dat);
            if (alpInt<0) throw new Fatal("ElastoPlastic::LoadCond: alpInt=%g must be positive",alpInt);
        }
        else if (numL<0.0) throw new Fatal("ElastoPlastic::LoadCond: Strain increment is too large (f=%g, f_tr=%g, numL=%g). Crossing and going all the way through the yield surface to the other side.",f,f_tr,numL);
    }

    // return true if there is loading
    // with intersection, return false (unloading)
    return ldg;
}

inline void ElastoPlastic::ELStiff (EquilibState const * Sta) const
{
    if (NDim==2)
    {
        if (GTy==pse_t)
        {
            double c = CalcE(Sta)/(1.0-nu*nu);
            De = c,    c*nu, 0.0,        0.0,
                 c*nu, c,    0.0,        0.0,
                 0.0,  0.0,  0.0,        0.0,
                 0.0,  0.0,  0.0, c*(1.0-nu);
        }
        else if (GTy==psa_t || GTy==axs_t)
        {
            double c = CalcE(Sta)/((1.0+nu)*(1.0-2.0*nu));
            De = c*(1.0-nu),       c*nu ,      c*nu ,            0.0,
                      c*nu ,  c*(1.0-nu),      c*nu ,            0.0,
                      c*nu ,       c*nu , c*(1.0-nu),            0.0,
                       0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu);
        }
        else throw new Fatal("ElastoPlastic::Stiffness: 2D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
    }
    else
    {
        if (GTy==d3d_t)
        {
            double c = CalcE(Sta)/((1.0+nu)*(1.0-2.0*nu));
            De = c*(1.0-nu),       c*nu ,      c*nu ,            0.0,            0.0,            0.0,
                      c*nu ,  c*(1.0-nu),      c*nu ,            0.0,            0.0,            0.0,
                      c*nu ,       c*nu , c*(1.0-nu),            0.0,            0.0,            0.0,
                       0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu),            0.0,            0.0,
                       0.0 ,        0.0 ,       0.0 ,            0.0, c*(1.0-2.0*nu),            0.0,
                       0.0 ,        0.0 ,       0.0 ,            0.0,            0.0, c*(1.0-2.0*nu);
        }
        else throw new Fatal("ElastoPlastic::Stiffness: 3D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
    }
}

inline void ElastoPlastic::Gradients (EquilibState const * Sta) const
{
    Y(0) = 0.0; // dfdz0
    Y(1) = 0.0; // dfdz1
    Y(2) = 0.0; // dfdz2

    switch (FC)
    {
        case VM_t:
        {
            OctInvs (Sta->Sig, p,q,s, qTol);
            V = s/(q*kVM);
            break;
        }
        case DP_t:
        {
            OctInvs (Sta->Sig, p,q,s, qTol);
            V = (1.0/(p*q*kDP))*s + (q/(p*p*kDP*Util::SQ3))*I;
            break;
        }
        case MC_t:
        {
            Calc_dgdsig (Sta->Sig);
            V = (q/(g*p*p*Util::SQ3))*I + (1.0/(p*q*g))*s - (q/(p*g*g))*dgdsig;
            break;
        }
        case MN_t:
        {
            CharInvs (Sta->Sig, I1,I2,I3, dI1dsig,dI2dsig,dI3dsig);
            V = (I2/(I3*kMN))*dI1dsig + (I1/(I3*kMN))*dI2dsig - (I1*I2/(I3*I3*kMN))*dI3dsig;
            break;
        }
    }

    // new stress update
    if (NewSU) Y(0) = -1.0;
}

inline void ElastoPlastic::FlowRule (EquilibState const * Sta) const
{
    switch (FC)
    {
        case VM_t: { W = V; break; }
        case DP_t: { W = V; break; }
        case MC_t:
        {
            if (NonAssoc)
            {
                Calc_dgdsig (Sta->Sig, true);
                //W = (q/(g*p*p*Util::SQ3))*I + (1.0/(p*q*g))*s - (q/(p*g*g))*dgdsig;
                W = (g/Util::SQ3)*I + (1.0/q)*s - p*dgdsig;
            }
            else W = V;
            break;
        }
        case MN_t: { W = V; break; }
    }
}

inline void ElastoPlastic::Hardening (EquilibState const * Sta) const
{
    Dev (W, devW);
    H(0) = 0.0;
    H(1) = Tra  (W);
    H(2) = Norm (devW);

    // new stress update
    if (NewSU)
    {
        double F;
        switch (FC)
        {
            case VM_t:
            {
                q = Calc_qoct (Sta->Sig);
                F = q/kVM - 1.0;
                break;
            }
            case DP_t:
            {
                Calc_pq (Sta->Sig);
                F = q/(p*kDP) - 1.0;
                break;
            }
            case MC_t:
            {
                Calc_pqg (Sta->Sig);
                F = q/(p*g) - 1.0;
                break;
            }
            case MN_t:
            {
                CharInvs (Sta->Sig, I1,I2,I3);
                F = I1*I2/(I3*kMN) - 1.0;
                break;
            }
        }
        double H1 = 0.0;
        if (F>0.0) F = 0.0;
        H(0) = AlpSU + (H1-AlpSU)*exp(BetSU*F);
        //printf("f=%g, F=%g, exp(BetSU*f)=%g, H(0)=%g\n",YieldFunc(Sta),F,exp(BetSU*F),H(0));
    }
}

inline double ElastoPlastic::YieldFunc (EquilibState const * Sta) const
{
    double f;
    switch (FC)
    {
        case VM_t:
        {
            q = Calc_qoct (Sta->Sig);
            f = q/kVM - Sta->Ivs(0);
            break;
        }
        case DP_t:
        {
            Calc_pq (Sta->Sig);
            f = q/(p*kDP) - Sta->Ivs(0);
            break;
        }
        case MC_t:
        {
            Calc_pqg (Sta->Sig);
            f = q/(p*g) - Sta->Ivs(0);
            break;
        }
        case MN_t:
        {
            CharInvs (Sta->Sig, I1,I2,I3);
            f = I1*I2/(I3*kMN) - Sta->Ivs(0);
            break;
        }
    }
    return f;
}

inline void ElastoPlastic::CorrectDrift (State * Sta) const
{
    EquilibState * sta = static_cast<EquilibState *>(Sta);
    double fnew  = YieldFunc (sta);
    size_t it    = 0;
    size_t maxIt = 10;
    Vec_t  VDe(NCps), DeW(NCps);
    while (fnew>CDFtol && it<maxIt)
    {
        Gradients (sta);
        FlowRule  (sta);
        Hardening (sta);
        double hp = (NIvs>0 ? Y(0)*H(0) : 0.0);
        if (it==0) ELStiff (sta);
        Mult (V, De, VDe);
        double dgam = fnew/(dot(VDe,W)-hp);
        DeW = De*W;
        sta->Sig -= dgam*DeW;
        sta->Ivs += dgam*H;
        fnew = YieldFunc (sta);
        if (fabs(fnew)<CDFtol) break;
        it++;
    }
    if (it>=maxIt) throw new Fatal("ElastoPlastic::CorrectDrift: Yield surface drift correction did not converge after %d iterations (CDFtol=%g)",it,CDFtol);
}

inline double ElastoPlastic::Falpha (double Alp, void * UserData)
{
    AlphaData const & dat = (*static_cast<AlphaData const *>(UserData));
    dat.StaAlp->Sig = dat.Sig + Alp * dat.DSig_tr;
    return YieldFunc (dat.StaAlp);
}

inline double ElastoPlastic::dFalpha (double Alp, void * UserData)
{
    AlphaData const & dat = (*static_cast<AlphaData const *>(UserData));
    dat.StaAlp->Sig = dat.Sig + Alp * dat.DSig_tr;
    Gradients (dat.StaAlp);
    return dot (V, dat.DSig_tr);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * ElastoPlasticMaker(int NDim, SDPair const & Prms) { return new ElastoPlastic(NDim,Prms); }

int ElastoPlasticRegister()
{
    ModelFactory   ["ElastoPlastic"] = ElastoPlasticMaker;
    MODEL.Set      ("ElastoPlastic", (double)MODEL.Keys.Size());
    MODEL_PRM_NAMES["ElastoPlastic"].Resize (14);
    MODEL_PRM_NAMES["ElastoPlastic"] = "E", "nu", "sY", "c", "phi", "Hp", "psi", "VM", "DP", "MC", "MN", "newsu", "betsu", "alpsu";
    MODEL_IVS_NAMES["ElastoPlastic"].Resize (0);
    return 0;
}

int __ElastoPlastic_dummy_int = ElastoPlasticRegister();


#endif // MECHSYS_ELASTOPLASTIC_H
