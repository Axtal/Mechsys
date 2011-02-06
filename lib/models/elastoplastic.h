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
#include <mechsys/models/smpinvs.h>
#include <mechsys/numerical/root.h>

class ElastoPlastic : public Model
{
public:
    // enums
    enum FCrit_t { VM_t, MC_t, GE_t }; ///< Failure criterion type

    // Constructor & Destructor
    ElastoPlastic (int NDim, SDPair const & Prms, size_t NIv=3, char const * Name="ElastoPlastic(VM)", bool DerivedModel=false);
    virtual ~ElastoPlastic () {}

    // Derived methods
    virtual void   TgIncs       (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const;
    virtual void   Stiffness    (State const * Sta, Mat_t & D)                                const;
    virtual size_t CorrectDrift (State       * Sta)                                           const;
    virtual bool   LoadCond     (State const * Sta, Vec_t const & DEps, double & alpInt)      const;

    // Internal methods to be overloaded by derived classes
    virtual void   InitIvs   (SDPair const & Ini, State * Sta)                            const;
    virtual void   Gradients (Vec_t const & Sig, Vec_t const & Ivs, bool Potential=false) const;
    virtual void   Hardening (Vec_t const & Sig, Vec_t const & Ivs)                       const;
    virtual double YieldFunc (Vec_t const & Sig, Vec_t const & Ivs)                       const;
    virtual void   ExtraDIvs (Vec_t & DIvs)                                               const;
    virtual double CalcE     (Vec_t const & Sig, Vec_t const & Ivs)                       const { return E; }
    virtual void   ELStiff   (Vec_t const & Sig, Vec_t const & Ivs)                       const;

    // Constants
    bool    Derived;    ///< Derived model (such as CamClay)
    double  E;          ///< Young
    double  nu;         ///< Poisson
    FCrit_t FC;         ///< Failure criterion: VM:Von-Mises
    double  kVM;        ///< von Mises coefficient
	double  kGE;        ///< General FC coefficient
	double  kGEpot;     ///< General FC coefficient (for potential/flow rule)
    double  pTol;       ///< Tolerance for minimum poct
    double  Hb;         ///< Hardening coefficient H_bar
    double  sphi;       ///< Sin(phi) friction angle
    double  spsi;       ///< Sin(psi) dilatancy angle
    double  FTol;       ///< Tolerance to be used when finding the intersection
    double  DCFTol;     ///< Drift correction ys function tolerance
    size_t  DCMaxIt;    ///< Drift correction max iterations
    double  qTol;       ///< Tolerance for minimum qoct
    bool    NewSU;      ///< New stress update ?
    double  BetSU;      ///< Beta coefficient for new stress update
    double  AlpSU;      ///< Alpha coefficient for new stress update
    Vec_t   I;          ///< Idendity tensor
    Mat_t   Psd, IdyI;  ///< Sym-dev and I dy I 4th order tensors

    // State data (mutable/scratch-pad)
    mutable double  pc, pm, r2;                ///< Variables for smoothing with arc
    mutable SMPInvs SMP;                       ///< SMP invariants
    mutable Vec_t SigB;                        ///< SigmaBar: shifted sigma for cohesion
    mutable Vec_t Sig0, Ivs0, SigA, DSigTr;    ///< Variables for yield crossing detection
    mutable Vec_t DEpsEl;                      ///< Elastic strain increment
    mutable Vec_t DEpsPl;                      ///< Plastic strain increment == Lam*W
    mutable Vec_t V;                           ///< NCps: Gradient of the yield surface
    mutable Vec_t W;                           ///< NCps: Plastic flow rule direction
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
    void Calc_pq     (Vec_t const & Sig) const              { p=Calc_poct(Sig);  q=Calc_qoct(Sig); }
    void Calc_pqt    (Vec_t const & Sig) const              { OctInvs(Sig,p,q,t);  th=asin(t)/3.0; }
    void Calc_pqg    (Vec_t const & Sig, double sinp) const { OctInvs(Sig,p,q,t);  th=asin(t)/3.0;  g=Util::SQ2*sinp/(Util::SQ3*cos(th)-sinp*sin(th)); }
    void Calc_dgdsig (Vec_t const & Sig, double sinp) const
    {
        OctInvs (Sig, p,q,t,th,s, qTol, &dthdsig);
        g = Util::SQ2*sinp/(Util::SQ3*cos(th)-sinp*sin(th));
        double dgdth = g*(Util::SQ3*sin(th)+sinp*cos(th))/(Util::SQ3*cos(th)-sinp*sin(th));
        dgdsig = dgdth * dthdsig;
    }

    void Calc_arc (double M) const
    {
        double den = 1.0 + M*M;
        pc = pTol/(1.0-M/sqrt(den));
        pm = pc/den;
        r2 = pow(pc*M,2.0)/den;
    }
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline ElastoPlastic::ElastoPlastic (int NDim, SDPair const & Prms, size_t NIv, char const * Name, bool Deriv)
    : Model    (NDim, Prms, NIv, Name),
      Derived  (Deriv),
      E        (0.0),
      nu       (0.0),
      FC       (VM_t),
      kVM      (0.0),
      pTol     (1.0e-5),
      Hb       (0.0),
      FTol     (1.0e-7),
      DCFTol   (1.0e-8),
      DCMaxIt  (10),
      qTol     (1.0e-7),
      NewSU    (false)
{
    // resize scratchpad arrays
    SigB   .change_dim (NCps);
    Sig0   .change_dim (NCps);
    Ivs0   .change_dim (NIvs);
    SigA   .change_dim (NCps);
    DSigTr .change_dim (NCps);
    DEpsEl .change_dim (NCps);
    V      .change_dim (NCps);
    W      .change_dim (NCps);
    Y      .change_dim (NIvs);
    H      .change_dim (NIvs);
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
    Calc_Psd  (NCps, Psd);
    Calc_IdyI (NCps, IdyI);

    // reference pressure
    if (Prms.HasKey("ptol")) pTol = Prms("ptol");

    if (!Derived) // for instance, CamClay
    {
        // parameters
        if (!Prms.HasKey("E"))  throw new Fatal("ElastoPlastic::ElastoPlastic: Young modulus (E) must be provided");
        if (!Prms.HasKey("nu")) throw new Fatal("ElastoPlastic::ElastoPlastic: Poisson coefficient (nu) must be provided");
        E  = Prms("E");
        nu = Prms("nu");
        if (Prms.HasKey("MC")) { Name="ElastoPlastic(MC)";  FC = MC_t;                     }
        if (Prms.HasKey("DP")) { Name="ElastoPlastic(DP)";  FC = GE_t;  SMP.b = 0.0;       }
        if (Prms.HasKey("MN")) { Name="ElastoPlastic(MN)";  FC = GE_t;  SMP.b = 0.5;       }
        if (Prms.HasKey("GE")) { Name="ElastoPlastic(GE)";  FC = GE_t;  SMP.b = Prms("b"); }
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
            double psi_rad = phi_rad;
            if (phi_deg<1.0e-3) throw new Fatal("ElastoPlastic::ElastoPlastic: Friction angle (phi [deg]) must be greater than zero (1.0e-3). phi=%g is invalid",phi_deg);
            if (c<0.0)          throw new Fatal("ElastoPlastic::ElastoPlastic: 'cohesion' must be greater than zero. c=%g is invalid",c);
            if (Prms.HasKey("psi")) psi_rad = Prms("psi")*Util::PI/180.0;
            sphi = sin(phi_rad);
            spsi = sin(psi_rad);
            if (FC==GE_t)
            {
                // yield surface
                double A = (1.0+sphi)/(1.0-sphi);
                SigA = -A, -1., -1., 0., 0., 0.;
                SMP.Calc (SigA, false);
                kGE = SMP.sq/SMP.sp;
                // potential
                A      = (1.0+spsi)/(1.0-spsi);
                SigA   = -A, -1., -1., 0., 0., 0.;
                SMP.Calc (SigA, false);
                kGEpot = SMP.sq/SMP.sp;
            }
        }

        // hardening
        if (Prms.HasKey("Hp")) Hb = (2.0/3.0)*Prms("Hp"); // Hp=H_prime

        // internal values
        IvNames = "z0", "evp", "edp";
    }

    // new stress update parmeters
    NewSU = (Prms.HasKey("newsu") ? (int)Prms("newsu") : false);
    if (NewSU)
    {
        BetSU = Prms("betsu");
        AlpSU = Prms("alpsu");
    }
}

inline void ElastoPlastic::TgIncs (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const
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
        // gradients, flow rule, hardening, and hp
        Gradients (sta->Sig, sta->Ivs);
        Gradients (sta->Sig, sta->Ivs, true);
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

        // plastic strains => extra DIvs
        ExtraDIvs (DIvs);
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

inline void ElastoPlastic::Stiffness (State const * Sta, Mat_t & D) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // De: elastic stiffness
    ELStiff (sta->Sig, sta->Ivs);

    // stiffness
    if (sta->Ldg)
    {
        // gradients, flow rule, hardening, and hp
        Gradients (sta->Sig, sta->Ivs);
        Gradients (sta->Sig, sta->Ivs, true);
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

inline size_t ElastoPlastic::CorrectDrift (State * Sta) const
{
    // state
    EquilibState * sta = static_cast<EquilibState *>(Sta);

    // iterations
    double fnew = YieldFunc (sta->Sig, sta->Ivs);
    size_t it   = 0;
    //if (NewSU) fnew = fabs(fnew); // TODO: should do this, but it does not converge
    while (fnew>DCFTol && it<DCMaxIt)
    {
        // gradients, flow rule, hardening, and hp
        Gradients (sta->Sig, sta->Ivs);
        Gradients (sta->Sig, sta->Ivs, true);
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
        //if (NewSU) fnew = fabs(fnew); // TODO: should do this, but it does not converge

        // check convergence
        if (fabs(fnew)<DCFTol) break;
        it++;
    }

    // check number of iterations
    if (it>=DCMaxIt) throw new Fatal("ElastoPlastic::CorrectDrift: Yield surface drift correction did not converge after %d iterations (fnew=%g, DCFTol=%g)",it,fnew,DCFTol);
    return it;
}

inline bool ElastoPlastic::LoadCond (State const * Sta, Vec_t const & DEps, double & alpInt) const
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

    // new stress update
    if (NewSU)
    {
        if (numL>0.0) ldg = true;
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
            if (fnew>0.0) throw new Fatal("ElastoPlastic::LoadCond: __internal_error__: correction for numL<0 failed (dalp=%g, fnew=%g)",dalp,fnew);
            //throw new Fatal("ElastoPlastic::LoadCond: Strain increment is too large (f=%g, f_tr=%g, numL=%g). Crossing and going all the way through the yield surface to the other side.",f,f_tr,numL);
        }
        if (crossing)
        {
            Numerical::Root<ElastoPlastic> root(const_cast<ElastoPlastic*>(this), &ElastoPlastic::Fa, &ElastoPlastic::dFada);
            //root.Scheme = "Newton";
            //root.Verbose = true;
            alpInt = root.Solve (0.0, 1.0);
            if (alpInt<0) throw new Fatal("ElastoPlastic::LoadCond: alpInt=%g must be positive",alpInt);
        }
        else ldg = true;
    }

    // return true if there is loading. returns false (unloading) with intersection
    return ldg;
}

inline void ElastoPlastic::InitIvs (SDPair const & Ini, State * Sta) const
{
    // initialize state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini, NIvs);

    // internal variables
    sta->Ivs(0) = 0.0; // z0
    sta->Ivs(1) = 0.0; // evp
    sta->Ivs(2) = 0.0; // edp

    // new stress update
    if (NewSU)
    {
        q = Calc_qoct (sta->Sig);
        if (q>qTol) throw new Fatal("ElastoPlastic::InitIvs: qoct==%g must be smaller than qTol==%g for the time being",q,qTol);
    }
    else
    {
        switch (FC)
        {
            case VM_t: { sta->Ivs(0) = 1.0;  break; }
            case MC_t: { sta->Ivs(0) = sphi; break; }
            case GE_t: { sta->Ivs(0) = kGE;  break; }
        }
    }

    // check initial yield function
    double f = YieldFunc (sta->Sig, sta->Ivs);
    if (f>FTol)           throw new Fatal("ElastoPlastic:InitIvs: stress point (sig=(%g,%g,%g,%g]) is outside yield surface (f=%g) with z0=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,f,sta->Ivs(0));
    if (NewSU && f<-FTol) throw new Fatal("ElastoPlastic:InitIvs: stress point (sig=(%g,%g,%g,%g]) is outside yield surface (f=%g) with z0=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,f,sta->Ivs(0));
}

inline void ElastoPlastic::Gradients (Vec_t const & Sig, Vec_t const & Ivs, bool Potential) const
{
    Vec_t * VorW = &V;
    if (Potential) VorW = &W;
    else
    {
        Y(0) = 0.0; // dfdz0
        Y(1) = 0.0; // dfdz1
        Y(2) = 0.0; // dfdz2
    }

    // pertubation
    SigA = Sig;
    if (NewSU)
    {
        q = Calc_qoct (SigA);
        if (q<qTol)
        {
            if (SigA(3)<0.0) SigA(3) -= qTol*Util::SQ2;
            else             SigA(3) += qTol*Util::SQ2;
        }
    }

    switch (FC)
    {
        case VM_t:
        {
            OctInvs (SigA, p,q,s, qTol);
            (*VorW) = s/(q*kVM);
            if (NewSU && !Potential) Y(0) = -1.0;
            break;
        }
        case MC_t:
        {
            // gradients: V
            Calc_dgdsig (SigA, (Potential ? spsi : Ivs(0)));
            Calc_arc    (g);
            double dfdp = -g;
            double dfdq = 1.0;
            double dfdg = -p;
            if (p<pm)
            {
                dfdp = 2.0*(p-pc);
                dfdq = 2.0*q;
                dfdg = 2.0*g*(r2-pc*pc)/(1.0+g*g);
            }
            if (q>qTol) (*VorW) = (-dfdp/Util::SQ3)*I + (dfdq/q)*s + dfdg*dgdsig;
            else        (*VorW) = (-dfdp/Util::SQ3)*I;

            // gradients: y0
            if (NewSU && !Potential)
            {
                double dgdz0 = (Util::SQ2+g*sin(th))/(Util::SQ3*cos(th)-Ivs(0)*sin(th));
                Y(0) = dfdg * dgdz0;
            }
            break;
        }
        case GE_t:
        {
            // gradients: V
            SMP.Calc (SigA, true);
            double M = (Potential ? kGEpot : Ivs(0));
            Calc_arc (M);
            double dfdp = -M;
            double dfdq = 1.0;
            if (SMP.sp<pm)
            {
                dfdp = 2.0*(SMP.sp-pc);
                dfdq = 2.0*SMP.sq;
            }
            (*VorW) = dfdp*SMP.dspdSig + dfdq*SMP.dsqdSig;

            // gradients: y0
            if (NewSU && !Potential)
            {
                if (SMP.sp<pm) Y(0) = 2.0*M*(r2-pc*pc)/(1.0+M*M);
                else           Y(0) = -SMP.sp;
            }
            break;
        }
    }
}

inline void ElastoPlastic::Hardening (Vec_t const & Sig, Vec_t const & Ivs) const
{
    // internal values
    H(0) = 0.0;
    H(1) = 0.0;
    H(2) = 0.0;

    // new stress update
    if (NewSU)
    {
        double D = 0.0;
        switch (FC)
        {
            case VM_t:
            {
                q = Calc_qoct (Sig);
                D = q/kVM - 1.0;
                break;
            }
            case MC_t:
            {
                Calc_pqg (Sig, sphi);
                Calc_arc (g);
                if (p<pm) D = q*q + pow(p-pc,2.0) - r2;
                else      D = q - g*p;
                break;
            }
            case GE_t:
            {
                SMP.Calc (Sig, false);
                Calc_arc (kGE);
                if (SMP.sp<pm) D = SMP.sq*SMP.sq + pow(SMP.sp-pc,2.0) - r2;
                else           D = SMP.sq - kGE*SMP.sp;
                break;
            }
        }
        double H1 = 0.0;
        if (D>0.0) D = 0.0;
        H(0) = AlpSU + (H1-AlpSU)*exp(BetSU*D);
    }
}

inline double ElastoPlastic::YieldFunc (Vec_t const & Sig, Vec_t const & Ivs) const
{
    double f = 1e+8;
    switch (FC)
    {
        case VM_t:
        {
            q = Calc_qoct (Sig);
            f = q/kVM - Ivs(0);
            break;
        }
        case MC_t:
        {
            Calc_pqg (Sig, Ivs(0));
            Calc_arc (g);
            if (p<pm) f = q*q + pow(p-pc,2.0) - r2;
            else      f = q - g*p;
            break;
        }
        case GE_t:
        {
            SMP.Calc (Sig, false);
            Calc_arc (Ivs(0));
            if (SMP.sp<pm) f = SMP.sq*SMP.sq + pow(SMP.sp-pc,2.0) - r2;
            else           f = SMP.sq - Ivs(0)*SMP.sp;
        }
    }
    return f;
}

inline void ElastoPlastic::ExtraDIvs (Vec_t & DIvs) const
{
    if (Derived) return;
    DIvs(1) = Calc_ev (DEpsPl); // devp
    DIvs(2) = Calc_ed (DEpsPl); // dedp
}

inline void ElastoPlastic::ELStiff (Vec_t const & Sig, Vec_t const & Ivs) const
{
    /*
    if (NLElastic)
    {
        double K = Calc_K (E, nu);
        double G = Calc_G (E, nu);
        De = (2.0*G)*Psd + K*IdyI;
    }
    */

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
        else throw new Fatal("ElastoPlastic::Stiffness: 2D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
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
        else throw new Fatal("ElastoPlastic::Stiffness: 3D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
    }
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * ElastoPlasticMaker(int NDim, SDPair const & Prms) { return new ElastoPlastic(NDim,Prms); }

int ElastoPlasticRegister()
{
    ModelFactory   ["ElastoPlastic"] = ElastoPlasticMaker;
    MODEL.Set      ("ElastoPlastic", (double)MODEL.Keys.Size());
    MODEL_PRM_NAMES["ElastoPlastic"].Resize (17);
    MODEL_PRM_NAMES["ElastoPlastic"] = "E", "nu", "sY", "c", "phi", "Hp", "psi", "VM", "DP", "MC", "MN", "GE", "b", "ptol", "newsu", "betsu", "alpsu";
    MODEL_IVS_NAMES["ElastoPlastic"].Resize (0);
    return 0;
}

int __ElastoPlastic_dummy_int = ElastoPlasticRegister();


#endif // MECHSYS_ELASTOPLASTIC_H
