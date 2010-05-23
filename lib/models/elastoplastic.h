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

using std::cout;
using std::endl;

/** Failure criteria names. */
SDPair FAILCRIT;

class ElastoPlastic : public Model
{
public:
    // enums
    enum FCrit_t { VM_t, DP_t, MC_t }; ///< Failure criterion type

    // Constructor
    ElastoPlastic (int NDim, SDPair const & Prms, bool DerivedModel=false);

    // Derived methods
    void   Stiffness    (State const * Sta, Mat_t & D, Vec_t * h=NULL, Vec_t * d=NULL) const;
    bool   LoadCond     (State const * Sta, Vec_t const & DEps, double & alpInt)       const;
    void   CorrectDrift (State       * Sta)                                            const;
    double CalcDEz      (State const * Sta, Vec_t const & DSig)                        const;

    // Internal methods
    void ELStiff (EquilibState const * Sta)                                 const;
    void EPStiff (EquilibState const * Sta, Vec_t * h=NULL, Vec_t * d=NULL) const;

    // Internal methods to be overloaded by derived classes
    virtual void   InitIvs   (SDPair const & Ini, State * Sta) const;
    virtual void   Gradients (EquilibState const * Sta)        const;
    virtual void   FlowRule  (EquilibState const * Sta)        const;
    virtual void   Hardening (EquilibState const * Sta)        const;
    virtual double YieldFunc (EquilibState const * Sta)        const;
    virtual double CalcE     (EquilibState const * Sta)        const { return E; }

    // Constants
    double  E;        ///< Young
    double  nu;       ///< Poisson
    FCrit_t FC;       ///< Failure criterion: VM:Von-Mises
    double  kY;       ///< Coefficient for yielding
    double  Hb;       ///< Hardening coefficient H_bar
    double  sphi;     ///< Sin(phi) friction angle
    double  spsi;     ///< Sin(psi) dilatancy angle
    double  cbar;     ///< Cohesion_bar
    double  ftol;     ///< Tolerance to be used when finding the intersection
    bool    NonAssoc; ///< Non-associated flow rule ? (for Mohr-Coulomb)

    // State data (mutable/scratch-pad)
    mutable Vec_t  V;    ///< Gradient of the yield surface
    mutable Vec_t  Y;    ///< Derivative of the yield surface w.r.t internal variables
    mutable Vec_t  W;    ///< Plastic flow rule direction
    mutable Vec_t  H;    ///< Hardening coefficients, one for each internal variable
    mutable Mat_t  De;   ///< Elastic stiffness
    mutable Mat_t  Dep;  ///< Elastoplastic stiffness

    // Auxiliar methods
private:
    void _MC_grads (EquilibState const * Sta, double SinPhiOrSinPsi, Vec_t & VorW) const; ///< Mohr-Coulomb gradients of YS or potential
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline ElastoPlastic::ElastoPlastic (int NDim, SDPair const & Prms, bool Derived)
    : Model (NDim,Prms,"ElastoPlastic"),
      E(0.0), nu(0.0), ftol(1.0e-9), NonAssoc(false)
{
    // number of stress/strain components
    NCps = 2*NDim;
    V  .change_dim (NCps);
    W  .change_dim (NCps);
    De .change_dim (NCps,NCps);
    Dep.change_dim (NCps,NCps);

    if (!Derived)
    {
        // parameters
        E  = Prms("E");
        nu = Prms("nu");
        FC = VM_t;
        kY = 0.0;
        Hb = 0.0;
        if (Prms.HasKey("fc"))
        {
            String fc;
            FAILCRIT.Val2Key (Prms("fc"), fc);
            if      (fc=="VM") FC = VM_t;
            //else if (fc=="DP") FC = DP_t;
            else if (fc=="MC") FC = MC_t;
            else throw new Fatal("ElastoPlastic::ElastoPlastic: Failure criterion fc=%s is not available",fc.CStr());
        }

        if (Prms.HasKey("NonAssoc")) NonAssoc = static_cast<bool>(Prms("NonAssoc"));

        if (FC==VM_t)
        {
            // constants
            if      (Prms.HasKey("sY")) kY = sqrt(2.0/3.0)*Prms("sY");
            else if (Prms.HasKey("cu")) kY = (GTy==psa_t ? sqrt(2.0)*Prms("cu") : 2.0*sqrt(2.0/3.0)*Prms("cu"));
            else throw new Fatal("ElastoPlastic::ElastoPlastic: With fc=VM (von Mises), either sY (uniaxial yield stress) or cu (undrained cohesion) must be provided");

            // internal values
            NIvs = 3;
            Y.change_dim (NIvs);
            H.change_dim (NIvs);
            IvNames.Push ("z0");
            IvNames.Push ("evp");
            IvNames.Push ("edp");
        }
        if (FC==MC_t)
        {
            // constants
            spsi = (Prms.HasKey("psi") ? sin(Prms("psi")*Util::PI/180.0) : 0.0);
            sphi = sin(Prms("phi")*Util::PI/180.0);
            cbar = sqrt(3.0)*Prms("c")/tan(Prms("phi")*Util::PI/180.0);
            ftol = 1.0e-5;

            // internal values
            NIvs = 2;
            Y.change_dim (NIvs);
            H.change_dim (NIvs);
            IvNames.Push ("evp");
            IvNames.Push ("edp");
        }
        if (Prms.HasKey("Hp")) Hb = (2.0/3.0)*Prms("Hp"); // H_prime
    }
}

inline void ElastoPlastic::InitIvs (SDPair const & Ini, State * Sta) const
{
    // initialize state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini, NIvs);

    // internal variables
    if (FC==VM_t)
    {
        sta->Ivs(0) = kY;  // size of YS
        sta->Ivs(1) = 0.0; // evp
        sta->Ivs(2) = 0.0; // edp
    }
    else
    {
        sta->Ivs(0) = 0.0; // evp
        sta->Ivs(1) = 0.0; // edp
    }

    // check initial yield function
    double f = YieldFunc (sta);
    if (f>1.0e-8) throw new Fatal("ElastoPlastic:InitIvs: stress point (sig=(%g,%g,%g,%g]) is outside yield surface (f=%g) with z0=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,f,sta->Ivs(0));
}

inline void ElastoPlastic::Stiffness (State const * Sta, Mat_t & D, Vec_t * h, Vec_t * d) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // elastic stiffness
    ELStiff (sta);

    // stiffness
    if (sta->Ldg)
    {
        EPStiff (sta, h, d);
        D = Dep;
        //cout << "Dep =\n" << PrintMatrix(D);
    }
    else
    {
    	D = De;
    	if (h!=NULL && d !=NULL)
    	{
    		h->change_dim (NIvs);
    		d->change_dim (NCps);
    		set_to_zero   ((*h));
    		set_to_zero   ((*d));
    	}
        //cout << "De =\n" << PrintMatrix(D);
    }

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

    //cout << "f = " << f << "   f_tr = " << f_tr << "\n";

    // going outside
    if (f_tr>0.0)
    {
        ldg = true;
        if (f*f_tr<0.0) // with crossing
        {
            size_t k     = 0;
            size_t maxIt = 10;
            double tol   = ftol;
            alpInt       = f/(f-f_tr);
            sta_tr.Sig   = sta->Sig + alpInt*dsig_tr;
            for (k=0; k<maxIt; ++k)
            {
                f_tr = YieldFunc (&sta_tr);
                if (fabs(f_tr)<tol) break;
                Gradients (&sta_tr);
                alpInt    += (-f_tr/dot(V,dsig_tr));
                sta_tr.Sig = sta->Sig + alpInt*dsig_tr;
            }
            if (k>=maxIt) throw new Fatal("ElastoPlastic::LoadCond: Newton-Rhapson (for calculating intersection) did not converge after %d iterations",k);
        }
        //else std::cout << "ElastoPlastic::LoadCond: f=" << f << "  f_tr=" << f_tr << "\n";
    }

    // return true if there is loading (also when there is intersection)
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

inline void ElastoPlastic::EPStiff (EquilibState const * Sta, Vec_t * h, Vec_t * d) const
{
    // gradients, flow rule, hardening, and hp
    Gradients (Sta);
    FlowRule  (Sta);
    Hardening (Sta);
    double hp = dot(Y,H);

    // auxiliar vectors
    Vec_t VDe;
    Mult (V, De, VDe);
    double phi = dot(VDe,W) - hp;
    Vec_t DeW(De*W);

    // elastoplastic stiffness
    for (size_t i=0; i<NCps; ++i)
    for (size_t j=0; j<NCps; ++j)
        Dep(i,j) = De(i,j) - DeW(i)*VDe(j)/phi;

    // internal values stiffness
    if (NIvs>0)
    {
        if (h!=NULL && d!=NULL)
        {
            (*d) = VDe/phi;
            (*h) = H;
        }
    }
}

inline void ElastoPlastic::Gradients (EquilibState const * Sta) const
{
    if (FC==VM_t)
    {
        double qoct = Calc_qoct (Sta->Sig);
        Vec_t s;
        Dev (Sta->Sig, s);
        V = s/qoct;
        Y(0) = -1.0; // dfdz0
        Y(1) = 0.0;
        Y(2) = 0.0;
    }
    else if (FC==DP_t)
    {
        //dfdp = -self.kdp
        //dfdq = 1.0
    }
    else if (FC==MC_t)
    {
        _MC_grads (Sta, sphi, V);
        Y(0) = 0.0;
        Y(1) = 0.0;
    }
}

inline void ElastoPlastic::FlowRule (EquilibState const * Sta) const
{ 
    if (FC==VM_t) W = V; 
    else if (FC==MC_t)
    {
        if (NonAssoc) _MC_grads (Sta, spsi, W);
        else W = V;
    }
}

inline void ElastoPlastic::Hardening (EquilibState const * Sta) const
{
    Vec_t dev_W;
    Dev (W, dev_W);
    if (FC==VM_t)
    {
        H(0) = Hb; 
        H(1) = Tra  (W);
        H(2) = Norm (dev_W);
    }
    else if (FC==MC_t)
    {
        H(0) = Tra  (W);
        H(1) = Norm (dev_W);
    }
}

inline double ElastoPlastic::YieldFunc (EquilibState const * Sta) const
{
    if (FC==VM_t)
    {
        double qoct = Calc_qoct (Sta->Sig);
        return qoct - Sta->Ivs(0);
    }
    else if (FC==DP_t)
    {
        // q - (p + self.cbar)*self.kdp
    }
    else if (FC==MC_t)
    {
        double p, q, t;
        OctInvs (Sta->Sig, p, q, t);
        double th = asin(t)/3.0;
        double g  = sqrt(2.0)*sphi/(sqrt(3.0)*cos(th)-sphi*sin(th));
        return q - (p + cbar)*g;
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
    while (fnew>tol)
    {
        Gradients (sta);
        FlowRule  (sta);
        Hardening (sta);
        double hp = dot(Y,H);
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

inline double ElastoPlastic::CalcDEz (State const * Sta, Vec_t const & DSig) const
{
    if (GTy!=pse_t) throw new Fatal("ElastoPlastic::CalcDEz: %dD: This method is not available for GeometryType = %s",NDim,GTypeToStr(GTy).CStr());
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);
    return -nu*(DSig(0)+DSig(1))/CalcE(sta);
}

inline void ElastoPlastic::_MC_grads (EquilibState const * Sta, double sinp, Vec_t & VorW) const
{
    // eigenvalues and eigenprojectors
    Vec3_t L;
    Vec_t  P0,P1,P2;
    EigenProj (Sta->Sig, L, P0, P1, P2);

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
    ModelFactory["ElastoPlastic"] = ElastoPlasticMaker;
    MODEL.Set ("ElastoPlastic", (double)MODEL.Keys.Size());
    double sz = FAILCRIT.Keys.Size();
    FAILCRIT.Set ("VM DP MC", sz, (sz+1.0), (sz+2.0));
    return 0;
}

int __ElastoPlastic_dummy_int = ElastoPlasticRegister();


#endif // MECHSYS_ELASTOPLASTIC_H
