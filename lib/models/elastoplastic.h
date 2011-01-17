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
    enum FCrit_t { VM_t, DP_t, MC_t, MN_t, AN_t, MNnl_t }; ///< Failure criterion type

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
    double  kY;         ///< Coefficient for yielding
    double  Hb;         ///< Hardening coefficient H_bar
    bool    NonAssoc;   ///< Non-associated flow rule ? (for Mohr-Coulomb)
    double  sphi;       ///< Sin(phi) friction angle
    double  spsi;       ///< Sin(psi) dilatancy angle
    double  cbar;       ///< Cohesion_bar
    double  ftol;       ///< Tolerance to be used when finding the intersection
    double  kDP;        ///< Drucker-Prager coefficient
    double  kMN;        ///< Matsuoka-Nakai coefficient
    bool    NewSU;      ///< New stress update ?
    double  BetSU;      ///< Beta coefficient for new stress update
    double  Gbar;       ///< G for new stress update
    Vec_t   I;          ///< Idendity tensor

    // State data (mutable/scratch-pad)
    mutable Vec_t V;    ///< NCps: Gradient of the yield surface
    mutable Vec_t W;    ///< NCps: Plastic flow rule direction
    mutable Vec_t Y;    ///< NIvs: Derivative of the yield surface w.r.t internal variables
    mutable Vec_t H;    ///< NIvs: Hardening coefficients, one for each internal variable
    mutable Mat_t De;   ///< Elastic stiffness
    mutable Mat_t Dep;  ///< Elastoplastic stiffness
    mutable Vec_t VDe;  ///< V*De
    mutable Vec_t DeW;  ///< De*W

    // Methods for yield surface crossing
    double Falpha  (double Alp, void * UserData);
    double dFalpha (double Alp, void * UserData);

private:
    // Auxiliary methods
    void _MC_grads (EquilibState const * Sta, double SinPhiOrSinPsi, Vec_t & VorW) const; ///< Mohr-Coulomb gradients of YS or potential
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline ElastoPlastic::ElastoPlastic (int NDim, SDPair const & Prms, bool Derived)
    : Model (NDim,Prms,"ElastoPlastic"),
      E(0.0), nu(0.0), FC(VM_t), kY(0.0), Hb(0.0), NonAssoc(false), ftol(1.0e-8), NewSU(false)
{
    // resize scratchpad arrays
    V  .change_dim (NCps);
    W  .change_dim (NCps);
    De .change_dim (NCps,NCps);
    Dep.change_dim (NCps,NCps);
    VDe.change_dim (NCps);
    DeW.change_dim (NCps);
    I  .change_dim (NCps);
    set_to_zero(I);
    I(0) = 1.0;
    I(1) = 1.0;
    I(2) = 1.0;

    if (!Derived) // for instance, CamClay
    {
        // parameters
        if (!Prms.HasKey("E"))  throw new Fatal("ElastoPlastic::ElastoPlastic: Young modulus (E) must be provided");
        if (!Prms.HasKey("nu")) throw new Fatal("ElastoPlastic::ElastoPlastic: Poisson's coefficient (nu) must be provided");
        E  = Prms("E");
        nu = Prms("nu");
        if (Prms.HasKey("DP"))   FC = DP_t;
        if (Prms.HasKey("MC"))   FC = MC_t;
        if (Prms.HasKey("MN"))   FC = MN_t;
        if (Prms.HasKey("AN"))   FC = AN_t;
        if (Prms.HasKey("MNnl")) FC = MNnl_t;
        if (FC==VM_t)
        {
            if      (Prms.HasKey("sY")) kY = sqrt(2.0/3.0)*Prms("sY");
            else if (Prms.HasKey("c"))  kY = (GTy==psa_t ? sqrt(2.0)*Prms("c") : 2.0*sqrt(2.0/3.0)*Prms("c"));
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
            cbar = sqrt(3.0)*c/tan(phi_rad);
            //ftol = 1.0e-5;
            if (FC==DP_t) kDP = 2.0*sqrt(2.0)*sphi/(3.0-sphi);
            if (FC==MN_t || FC==MNnl_t) kMN = 9.0+8.0*pow(tan(phi_rad),2.0);
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
    }

    // TODO: new stress update parmeters
    NewSU  = (Prms.HasKey("newsu") ? (int)Prms("newsu") : false);
    BetSU  = (Prms.HasKey("betsu") ? Prms("betsu") : 10.0);
    Gbar   = (Prms.HasKey("Gbar")  ? Prms("Gbar")  : 0.5*E/(1.0+nu));
}

inline void ElastoPlastic::InitIvs (SDPair const & Ini, State * Sta) const
{
    // initialize state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini, NIvs);

    // internal variables
    sta->Ivs(0) = kY;  // size of YS (von Mises only)
    sta->Ivs(1) = 0.0; // evp
    sta->Ivs(2) = 0.0; // edp

    // TODO: new stress update
    if (NewSU)
    {
        double p = Calc_poct (sta->Sig);
        double q = Calc_qoct (sta->Sig);
        if (FC==MN_t || FC==MNnl_t)
        {
            double I1,I2,I3;
            CharInvs (sta->Sig, I1,I2,I3);
            sta->Ivs(0) = I1*I2/I3;
        }
        else if (FC==DP_t) sta->Ivs(0) = q/p;
        else sta->Ivs(0) = q;
    }

    // check initial yield function
    double f = YieldFunc (sta);
    if (f>1.0e-8) throw new Fatal("ElastoPlastic:InitIvs: stress point (sig=(%g,%g,%g,%g]) is outside yield surface (f=%g) with z0=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,f,sta->Ivs(0));
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

        // TODO: new stress update
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

    // TODO: new stress update
    if (NewSU)
    {
        double q = Calc_qoct (sta->Sig);
        if (q>1.0e-7)
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
        if (f<-1.0e-8 && f_tr>0.0) crossing = true; // works
        //if (f<0.0 && f_tr>0.0) crossing = true; // does not work
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

    double qoct = Calc_qoct (Sta->Sig);
    Vec_t s;
    Dev (Sta->Sig, s);
    if (qoct<1.0e-8)
    {
        Vec_t sig(Sta->Sig);
        sig(2) += 1.0e-5;
        Dev (sig, s);
        qoct = Calc_qoct (sig);
        if (qoct<1.0e-8)
        {
            std::ostringstream oss;
            oss << "Sig = "      << PrintVector(Sta->Sig);
            oss << "sig = "      << PrintVector(sig);
            oss << "dev(sig) = " << PrintVector(s);
            throw new Fatal("ElastoPlastic::Gradients:: __internal_error__ qoct=%g is too small\n%s",qoct,oss.str().c_str());
        }
    }

    if (FC==VM_t)
    {
        V    = s/qoct;
        Y(0) = -1.0; // dfdz0
    }
    else if (FC==DP_t) V = s/qoct + (kDP/Util::SQ3)*I;
    else if (FC==MC_t) _MC_grads (Sta, sphi, V);
    else if (FC==MN_t)
    {
        double I1,I2,I3;
        Vec_t dI1ds,dI2ds,dI3ds;
        CharInvs (Sta->Sig, I1,I2,I3, dI1ds,dI2ds,dI3ds);
        V = (I2/I3)*dI1ds + (I1/I3)*dI2ds - (I1*I2/pow(I3,2.0))*dI3ds;
    }

    // TODO: new stress update
    if (NewSU)
    {
        double poct = Calc_poct (Sta->Sig);
        if (FC==DP_t) Y(0) = -poct;
        else          Y(0) = -1.0; // dfdz0
    }
}

inline void ElastoPlastic::FlowRule (EquilibState const * Sta) const
{
    if (FC==MC_t)
    {
        if (NonAssoc) _MC_grads (Sta, spsi, W);
        else W = V;
    }
    else W = V; 
}

inline void ElastoPlastic::Hardening (EquilibState const * Sta) const
{
    Vec_t dev_W;
    Dev (W, dev_W);
    if (FC==VM_t) H(0) = Hb; 
    else          H(0) = 0.0;
    H(1) = Tra  (W);
    H(2) = Norm (dev_W);

    // TODO: new stress update
    if (NewSU)
    {
        double k = 0;
        if      (FC==VM_t) k = kY;
        else if (FC==DP_t) k = kDP;
        else if (FC==MN_t) k = kMN;
        else if (FC==MC_t)
        {
            double p, q, t;
            OctInvs (Sta->Sig, p, q, t);
            double th = asin(t)/3.0;
            double g  = sqrt(2.0)*sphi/(sqrt(3.0)*cos(th)-sphi*sin(th));
            k = (p + cbar)*g;
        }
        else throw new Fatal("ElastoPlastic::Hardening: NewSU FC not implemented");


        double qoct = Calc_qoct (Sta->Sig);
        Vec_t s;
        Dev (Sta->Sig, s);
        if (qoct<1.0e-8)
        {
            Vec_t sig(Sta->Sig);
            sig(2) += 1.0e-5;
            Dev (sig, s);
            qoct = Calc_qoct (sig);
            if (qoct<1.0e-8)
            {
                std::ostringstream oss;
                oss << "Sig = "      << PrintVector(Sta->Sig);
                oss << "sig = "      << PrintVector(sig);
                oss << "dev(sig) = " << PrintVector(s);
                throw new Fatal("ElastoPlastic::Hardening:: __internal_error__ qoct=%g is too small\n%s",qoct,oss.str().c_str());
            }
        }


        double D = 2.0*k/(k+Sta->Ivs(0))-1.0;
        if (D<0.0) D = 0.0;
        double m = 1.0-exp(-BetSU*D);
        //H(0) = Gbar*m*Norm(dev_W);
        H(0) = Gbar*m;//*Norm(s);
        printf("z0=%g, H0=%g, Norm(devW)=%g, D=%g, m=%g\n",Sta->Ivs(0),H(0),Norm(dev_W),D,m);
    }
}

inline double ElastoPlastic::YieldFunc (EquilibState const * Sta) const
{
    if (FC==VM_t)
    {
        double q = Calc_qoct (Sta->Sig);
        return q - Sta->Ivs(0);
    }
    else if (FC==DP_t)
    {
        double p = Calc_poct (Sta->Sig);
        double q = Calc_qoct (Sta->Sig);
        if (NewSU) return q - p*Sta->Ivs(0); // TODO: new stress update
        return q - p*kDP;
    }
    else if (FC==MC_t)
    {
        double p, q, t;
        OctInvs (Sta->Sig, p, q, t);
        double th = asin(t)/3.0;
        double g  = sqrt(2.0)*sphi/(sqrt(3.0)*cos(th)-sphi*sin(th));
        if (NewSU) return q - Sta->Ivs(0); // TODO: new stress update
        return q - (p + cbar)*g;
    }
    else if (FC==MN_t)
    {
        double I1,I2,I3;
        CharInvs (Sta->Sig, I1,I2,I3);
        if (NewSU) return I1*I2/I3 - Sta->Ivs(0); // TODO: new stress update
        return I1*I2/I3 - kMN;
    }
    return 0;
}

inline void ElastoPlastic::CorrectDrift (State * Sta) const
{
    EquilibState * sta = static_cast<EquilibState *>(Sta);
    double fnew  = YieldFunc (sta);
    size_t it    = 0;
    size_t maxIt = 10;
    double tol   = 1.0e-8;
    Vec_t  VDe(NCps), DeW(NCps);
    while (fnew>tol && it<maxIt)
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
        if (fabs(fnew)<tol) break;
        it++;
    }
    if (it>=maxIt) throw new Fatal("ElastoPlastic::CorrectDrift: Yield surface drift correction did not converge after %d iterations",it);
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

inline void ElastoPlastic::_MC_grads (EquilibState const * Sta, double sinp, Vec_t & VorW) const
{
    // eigenvalues and eigenprojectors
    Vec3_t L, v0, v1, v2;
    Vec_t  P0,P1,P2;
    EigenProj (Sta->Sig, L, v0, v1, v2, P0, P1, P2);

    // oct invariants and its derivatives w.r.t principal values (L)
    double p,q,t;
    Vec3_t dpdL,dqdL,dtdL;
    OctInvs (L, p, q, t, dpdL, dqdL, dtdL);

    // derivatives of f w.r.t. oct invariants
    double th   = asin(t)/3.0;
    double g    = sqrt(2.0)*sinp/(sqrt(3.0)*cos(th)-sinp*sin(th));
    double dfdp = -g;
    double dfdq = 1.0;
    double dfdt = 0.0;
    if (t>-0.999 && t<0.999)
    {
        double dgdth = g*(sqrt(3.0)*sin(th)+sinp*cos(th))/(sqrt(3.0)*cos(th)-sinp*sin(th));
        double dfdth = -(p+cbar)*dgdth;
        double dthdt = 1.0/(3.0*sqrt(1.0-t*t));
        dfdt  = dfdth*dthdt;
    }

    // gradient w.r.t principal values (L)
    Vec3_t dfdL(dfdp*dpdL + dfdq*dqdL + dfdt*dtdL);

    // gradient w.r.t sig
    VorW = dfdL(0)*P0 + dfdL(1)*P1 + dfdL(2)*P2;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * ElastoPlasticMaker(int NDim, SDPair const & Prms) { return new ElastoPlastic(NDim,Prms); }

int ElastoPlasticRegister()
{
    ModelFactory   ["ElastoPlastic"] = ElastoPlasticMaker;
    MODEL.Set      ("ElastoPlastic", (double)MODEL.Keys.Size());
    MODEL_PRM_NAMES["ElastoPlastic"].Resize (15);
    MODEL_PRM_NAMES["ElastoPlastic"] = "E", "nu", "sY", "c", "phi", "Hp", "psi", "VM", "DP", "MC", "MN", "AN", "newsu", "betsu", "Gbar";
    MODEL_IVS_NAMES["ElastoPlastic"].Resize (0);
    return 0;
}

int __ElastoPlastic_dummy_int = ElastoPlasticRegister();


#endif // MECHSYS_ELASTOPLASTIC_H
