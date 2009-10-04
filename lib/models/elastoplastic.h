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
#include "models/model.h"
#include "models/equilibstate.h"

using std::cout;
using std::endl;

class ElastoPlastic : public Model
{
public:
    // enums
    enum FCrit_t { VM_t, DP_t, MC_t }; ///< Failure criterion type

    // Constructor
    ElastoPlastic (int NDim, SDPair const & Prms);

    // Methods
    void InitIvs   (SDPair const & Ini, State * Sta) const;
    void Stiffness (State const * Sta, Mat_t & D, Array<Vec_t> * d=NULL) const;

    // Internal methods
    void   ELStiff   (EquilibState const * Sta, Mat_t & D)                         const;
    void   EPStiff   (EquilibState const * Sta, Mat_t const & De, Vec_t const & V,
                      Vec_t const & Y, Mat_t & Dep, Array<Vec_t> & dep)            const;
    void   Gradients (EquilibState const * Sta, Vec_t & V, Vec_t & Y)              const;
    void   FlowRule  (EquilibState const * Sta, Vec_t const & V, Vec_t & W)        const { W = V; }
    void   Hardening (EquilibState const * Sta, Vec_t const & W, Vec_t & H)        const {}
    double YieldFunc (EquilibState const * Sta)                                    const;
    double FailCrit  (EquilibState const * Sta) const { return YieldFunc (Sta); }

    // Data
    double  E;
    double  nu;
    double  sY; ///< Uniaxial yield stress
    FCrit_t FC; ///< Failure criterion: VM:Von-Mises
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline ElastoPlastic::ElastoPlastic (int NDim, SDPair const & Prms)
    : Model (NDim,Prms,"ElastoPlastic"), FC(VM_t)
{
    E  = Prms("E");
    nu = Prms("nu");
    sY = Prms("sY");
}

inline void ElastoPlastic::InitIvs (SDPair const & Ini, State * Sta) const
{
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini);
}

inline void ElastoPlastic::Stiffness (State const * Sta, Mat_t & D, Array<Vec_t> * d) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // elastic stiffness
    Mat_t De;
    ELStiff (sta, De);

    // stiffness
    if (sta->Loading)
    {
        Vec_t V, Y;
        Gradients (sta, V, Y);
        EPStiff   (sta, De, V, Y, D, (*d));
    }
    else
    {
        D = De;
    }
}

inline void ElastoPlastic::ELStiff (EquilibState const * Sta, Mat_t & D) const
{
    if (NDim==2)
    {
        D.change_dim (4,4);
        if (GTy==pse_t)
        {
            double c = E/(1.0-nu*nu);
            D = c,    c*nu, 0.0,        0.0,
                c*nu, c,    0.0,        0.0,
                0.0,  0.0,  0.0,        0.0,
                0.0,  0.0,  0.0, c*(1.0-nu);
        }
        else if (GTy==psa_t || GTy==axs_t)
        {
            double c = E/((1.0+nu)*(1.0-2.0*nu));
            D = c*(1.0-nu),       c*nu ,      c*nu ,            0.0,
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
            D.change_dim (6,6);
            double c = E/((1.0+nu)*(1.0-2.0*nu));
            D = c*(1.0-nu),       c*nu ,      c*nu ,            0.0,            0.0,            0.0,
                     c*nu ,  c*(1.0-nu),      c*nu ,            0.0,            0.0,            0.0,
                     c*nu ,       c*nu , c*(1.0-nu),            0.0,            0.0,            0.0,
                      0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu),            0.0,            0.0,
                      0.0 ,        0.0 ,       0.0 ,            0.0, c*(1.0-2.0*nu),            0.0,
                      0.0 ,        0.0 ,       0.0 ,            0.0,            0.0, c*(1.0-2.0*nu);
        }
        else throw new Fatal("ElastoPlastic::Stiffness: 3D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
    }
}

inline void ElastoPlastic::EPStiff (EquilibState const * Sta, Mat_t const & De, Vec_t const & V, Vec_t const & Y, Mat_t & Dep, Array<Vec_t> & dep) const
{
    Vec_t W, H;
    FlowRule  (Sta, V, W);
    Hardening (Sta, W, H);
    double hp  = 0.0;
    size_t niv = size(Sta->Ivs);
    for (size_t i=0; i<niv; ++i) hp += Y(i)*H(i);
    Vec_t VDe;
    Mult (V, De, VDe);
    double phi = dot(VDe,W) - hp;
    Vec_t DeW(De*W);
    Dyad (DeW, VDe, Dep);
    Dep = De - (1.0/phi)*Dep;
    if (niv>0)
    {
        dep.Resize (niv);
        for (size_t i=0; i<niv; ++i) dep[i] = (H(i)/phi)*VDe;
    }
}

inline void ElastoPlastic::Gradients (EquilibState const * Sta, Vec_t & V, Vec_t & Y) const
{
    // eigenvalues and eigenprojectors
    Vec3_t L;
    Vec_t  P0,P1,P2;
    EigenProj (Sta->Sig, L, P0, P1, P2);

    // derivatives w.r.t. (oct) invariants
    double dfdp=0., dfdq=0., dfdt=0.;
    if (FC==VM_t)
    {
        dfdp = 0.0;
        dfdq = 1.0;
    }
    else if (FC==DP_t)
    {
        /*
        dfdp = -self.kdp
        dfdq = 1.0
        */
    }
    else if (FC==MC_t)
    {
        /*
        th   = arcsin(t)/3.0
        g    = self.g(th)
        dfdp = -g
        dfdq = 1.0
        if t>-0.999 and t<0.999:
            dgdth = g*(sqrt(3.0)*sin(th)+self.sphi*cos(th))/(sqrt(3.0)*cos(th)-self.sphi*sin(th))
            dfdth = -(p+self.cbar)*dgdth
            dthdt = 1.0/(3.0*sqrt(1.0-t**2.0))
            dfdt  = dfdth*dthdt
        */
    }

    // derivatives of (oct) invariants w.r.t principal values (L)
    double p,q,t;
    Vec3_t dpdL,dqdL,dtdL;
    OctInvs (L, p, q, t, dpdL, dqdL, dtdL);

    // gradient w.r.t principal values (L)
    Vec3_t dfdL;
    dfdL = dfdp*dpdL + dfdq*dqdL + dfdt*dtdL;

    // gradient w.r.t sig
    V = dfdL(0)*P0 + dfdL(1)*P1 + dfdL(2)*P2;
}

inline double ElastoPlastic::YieldFunc (EquilibState const * Sta) const
{
    double p, q, t, f;
    OctInvs (Sta->Sig, p, q, t);
    if (FC==VM_t)
    {
        f = q - Util::SQ3*sY;
    }
    else if (FC==DP_t)
    {
        //f  = q - (p + self.cbar)*self.kdp
    }
    else if (FC==MC_t)
    {
        //th   = arcsin(t)/3.0
        //f    = q - (p + self.cbar)*self.g(th)
    }
    return f;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * ElastoPlasticMaker(int NDim, SDPair const & Prms) { return new ElastoPlastic(NDim,Prms); }

int ElastoPlasticRegister()
{
    ModelFactory["ElastoPlastic"] = ElastoPlasticMaker;
    MODEL.Set ("ElastoPlastic", (double)MODEL.Keys.Size());
    return 0;
}

int __ElastoPlastic_dummy_int = ElastoPlasticRegister();


#endif // MECHSYS_ELASTOPLASTIC_H
