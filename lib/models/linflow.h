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

#ifndef MECHSYS_LINFLOW_H
#define MECHSYS_LINFLOW_H

// MechSys
#include "models/model.h"
#include "models/flowstate.h"

class LinFlow : public Model
{
public:
    // Constructor
    LinFlow (int NDim, SDPair const & Prms);

    // Methods
    void InitIvs   (SDPair const & Ini, State * Sta) const;
    void Stiffness (State const * Sta, Mat_t & D, Array<Vec_t> * d=NULL) const;

    // Data
    double kx;
    double ky;
    double kz;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline LinFlow::LinFlow (int NDim, SDPair const & Prms)
    : Model (NDim,Prms,"LinFlow")
{
    if (Prms.HasKey("k"))
    {
        double k = Prms("k");
        kx = k;
        ky = k;
        kz = k;
    }
    else
    {
        kx = Prms("kx");
        ky = Prms("ky");  if (NDim==3)
        kz = Prms("kz");
    }
}

inline void LinFlow::InitIvs (SDPair const & Ini, State * Sta) const
{
    FlowState * sta = static_cast<FlowState*>(Sta);
    sta->Init (Ini);
}

inline void LinFlow::Stiffness (State const * Sta, Mat_t & D, Array<Vec_t> * d) const
{
    D.change_dim (NDim,NDim);
    if (NDim==2)
    {
        D =  kx, 0.0,
            0.0,  ky;
    }
    else if (NDim==3)
    {
        D =  kx, 0.0, 0.0,
            0.0,  ky, 0.0,
            0.0, 0.0,  kz;
    }
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * LinFlowMaker(int NDim, SDPair const & Prms) { return new LinFlow(NDim,Prms); }

int LinFlowRegister()
{
    ModelFactory["LinFlow"] = LinFlowMaker;
    MODEL.Set ("LinFlow", (double)MODEL.Keys.Size());
    return 0;
}

int __LinFlow_dummy_int = LinFlowRegister();


#endif // MECHSYS_LINFLOW_H
