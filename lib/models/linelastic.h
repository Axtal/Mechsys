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
#include "models/model.h"
#include "models/equilibstate.h"

using std::cout;
using std::endl;

class LinElastic : public Model
{
public:
    // Constructor
    LinElastic (int NDim, SDPair const & Prms);

    // Methods
    void   InitIvs   (SDPair const & Ini, State * Sta)        const;
    void   Stiffness (State const * Sta, Mat_t & D, 
                      Array<double> * h=NULL, Vec_t * d=NULL) const;
    double CalcDEz   (Vec_t const & DSig) const;

    // Data
    double E;
    double nu;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline LinElastic::LinElastic (int NDim, SDPair const & Prms)
    : Model (NDim,Prms,"LinElastic")
{
    E  = Prms("E");
    nu = Prms("nu");
}

inline void LinElastic::InitIvs (SDPair const & Ini, State * Sta) const
{
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini);
}

inline void LinElastic::Stiffness (State const * Sta, Mat_t & D, Array<double> * h, Vec_t * d) const
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
        else throw new Fatal("LinElastic::Stiffness: 2D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
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
        else throw new Fatal("LinElastic::Stiffness: 3D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
    }
}

inline double LinElastic::CalcDEz (Vec_t const & DSig) const
{
    if (GTy!=pse_t) throw new Fatal("LinElastic::CalcDEz: %dD: This method is not available for GeometryType = %s",NDim,GTypeToStr(GTy).CStr());
    return -nu*(DSig(0)+DSig(1))/E;
}

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * LinElasticMaker(int NDim, SDPair const & Prms) { return new LinElastic(NDim,Prms); }

int LinElasticRegister()
{
    ModelFactory["LinElastic"] = LinElasticMaker;
    MODEL.Set ("LinElastic", (double)MODEL.Keys.Size());
    return 0;
}

int __LinElastic_dummy_int = LinElasticRegister();


#endif // MECHSYS_LINELASTIC_H
