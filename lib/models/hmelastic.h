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

#ifndef MECHSYS_HMELASTIC_H
#define MECHSYS_HMELASTIC_H

// MechSys
#include <mechsys/models/model.h>
#include <mechsys/models/hmstate.h>

using std::cout;
using std::endl;

class HMElastic : public Model
{
public:
    // Constructor
    HMElastic (int NDim, SDPair const & Prms);

    // Methods
    void InitIvs   (SDPair const & Ini, State * Sta)                                           const;
    void Stiffness (State const * Sta, Mat_t & D, Vec_t & Dw, Vec_t * h=NULL, Vec_t * d=NULL)  const;
    void Hydraulic (State const * Sta, Mat_t & Kw, double & ChiW, double & InvQs) const;

    // Data
    double E;
    double nu;
    Mat_t  Kw;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline HMElastic::HMElastic (int NDim, SDPair const & Prms)
    : Model (NDim,Prms,"HMElastic")
{
    NCps = 2*NDim;
    E    = Prms("E");
    nu   = Prms("nu");
    double gamW = Prms("gamW");
    double kw  = (Prms.HasKey("kw")  ? Prms("kw")  : 1.0);
    double kwx = (Prms.HasKey("kwx") ? Prms("kwx") : kw);
    double kwy = (Prms.HasKey("kwy") ? Prms("kwy") : kw);
    double kwz = (Prms.HasKey("kwz") ? Prms("kwz") : kw);
    Kw.change_dim (NDim,NDim);
    if (NDim==2)
    {
        Kw = kwx, 0.0,
             0.0, kwy;
    }
    else
    {
        Kw = kwx, 0.0, 0.0,
             0.0, kwy, 0.0,
             0.0, 0.0, kwz;
    }
    Kw /= gamW;
}

inline void HMElastic::InitIvs (SDPair const & Ini, State * Sta) const
{
    HMState * sta = static_cast<HMState*>(Sta);
    sta->Init (Ini);
}

inline void HMElastic::Stiffness (State const * Sta, Mat_t & D, Vec_t & Dw, Vec_t * h, Vec_t * d) const
{
    Dw.change_dim(NCps);
    set_to_zero(Dw);

    if (NDim==2)
    {
        D.change_dim (4,4);
        if (GTy==psa_t || GTy==axs_t)
        {
            double c = E/((1.0+nu)*(1.0-2.0*nu));
            D = c*(1.0-nu),       c*nu ,      c*nu ,            0.0,
                     c*nu ,  c*(1.0-nu),      c*nu ,            0.0,
                     c*nu ,       c*nu , c*(1.0-nu),            0.0,
                      0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu);
        }
        else throw new Fatal("HMElastic::Stiffness: 2D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
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
        else throw new Fatal("HMElastic::Stiffness: 3D: This model is not available for GeometryType = %s",GTypeToStr(GTy).CStr());
    }
}

inline void HMElastic::Hydraulic (State const * Sta, Mat_t & TheKw, double & ChiW, double & InvQs) const
{
    ChiW  = 1.0;
    InvQs = 0.0;
    TheKw = Kw;
}

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * HMElasticMaker(int NDim, SDPair const & Prms) { return new HMElastic(NDim,Prms); }

int HMElasticRegister()
{
    ModelFactory["HMElastic"] = HMElasticMaker;
    MODEL.Set ("HMElastic", (double)MODEL.Keys.Size());
    return 0;
}

int __HMElastic_dummy_int = HMElasticRegister();


#endif // MECHSYS_HMELASTIC_H
