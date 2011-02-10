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

#ifndef MECHSYS_LINELASTIC_H
#define MECHSYS_LINELASTIC_H

// MechSys
#include <mechsys/models/model.h>
#include <mechsys/models/equilibstate.h>

class LinElastic : public Model
{
public:
    // Constructor
    LinElastic (int NDim, SDPair const & Prms);

    // Methods
    void InitIvs   (SDPair const & Ini, State * Sta)                             const;
    void Rate      (State const * Sta, Vec_t const & DEpsDt, Vec_t & DSigDt, Vec_t & DIvsDt) const;
    void TgIncs    (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const;
    void InvTgIncs (State const * Sta, Vec_t & DSig, Vec_t & DEps, Vec_t & DIvs) const;
    void Stiffness (State const * Sta, Mat_t & D)                                const { D = De; }

    // Data
    double E;
    double nu;
    Mat_t  De;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline LinElastic::LinElastic (int NDim, SDPair const & Prms)
    : Model (NDim,Prms,/*niv*/0,"LinElastic")
{
    // parameters
    E  = Prms("E");
    nu = Prms("nu");

    // elastic stiffness
    if (NDim==2)
    {
        De.change_dim (4,4);
        if (GTy==pse_t)
        {
            double c = E/(1.0-nu*nu);
            De = c,    c*nu, 0.0,        0.0,
                 c*nu, c,    0.0,        0.0,
                 0.0,  0.0,  0.0,        0.0,
                 0.0,  0.0,  0.0, c*(1.0-nu);
        }
        else if (GTy==psa_t || GTy==axs_t)
        {
            double c = E/((1.0+nu)*(1.0-2.0*nu));
            De = c*(1.0-nu),       c*nu ,      c*nu ,            0.0,
                      c*nu ,  c*(1.0-nu),      c*nu ,            0.0,
                      c*nu ,       c*nu , c*(1.0-nu),            0.0,
                       0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu);
        }
        else throw new Fatal("LinElastic::Stiffness: 2D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
    }
    else
    {
        if (GTy==d3d_t)
        {
            De.change_dim (6,6);
            double c = E/((1.0+nu)*(1.0-2.0*nu));
            De = c*(1.0-nu),       c*nu ,      c*nu ,            0.0,            0.0,            0.0,
                      c*nu ,  c*(1.0-nu),      c*nu ,            0.0,            0.0,            0.0,
                      c*nu ,       c*nu , c*(1.0-nu),            0.0,            0.0,            0.0,
                       0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu),            0.0,            0.0,
                       0.0 ,        0.0 ,       0.0 ,            0.0, c*(1.0-2.0*nu),            0.0,
                       0.0 ,        0.0 ,       0.0 ,            0.0,            0.0, c*(1.0-2.0*nu);
        }
        else throw new Fatal("LinElastic::Stiffness: 3D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
    }

    // set model in stress update
    SUp.SetModel (this);
}

inline void LinElastic::InitIvs (SDPair const & Ini, State * Sta) const
{
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini);
}

inline void LinElastic::InvTgIncs (State const * Sta, Vec_t & DSig, Vec_t & DEps, Vec_t & DIvs) const
{
    Mat_t Ce(NCps,NCps);
    Inv (De, Ce);
    DEps = Ce*DSig;
}

inline void LinElastic::Rate (State const * Sta, Vec_t const & DEpsDt, Vec_t & DSigDt, Vec_t & DIvsDt) const
{
    DSigDt = De*DEpsDt;
    if (GTy==pse_t) throw new Fatal("LinElastic::Rate: this method does not work with plane-stress (pse)");
}

inline void LinElastic::TgIncs (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const
{
    DSig = De*DEps;
    if (GTy==pse_t) DEps(2) = -nu*(DSig(0)+DSig(1))/E;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * LinElasticMaker(int NDim, SDPair const & Prms, Model const * O) { return new LinElastic(NDim,Prms); }

int LinElasticRegister()
{
    ModelFactory   ["LinElastic"] = LinElasticMaker;
    MODEL.Set      ("LinElastic", (double)MODEL.Keys.Size());
    MODEL_PRM_NAMES["LinElastic"].Resize(2);
    MODEL_PRM_NAMES["LinElastic"] = "E", "nu";
    MODEL_IVS_NAMES["LinElastic"].Resize(0);
    return 0;
}

int __LinElastic_dummy_int = LinElasticRegister();


#endif // MECHSYS_LINELASTIC_H
