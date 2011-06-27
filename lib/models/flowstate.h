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

#ifndef MECHSYS_FLOWSTATE_H
#define MECHSYS_FLOWSTATE_H

// Std Lib
#include <iostream>
#include <cmath>   // for sqrt

// MechSys
#include <mechsys/models/model.h>
#include <mechsys/linalg/matvec.h>

class FlowState : public State
{
public:
	// Constructor
	FlowState (int NDim);

	// Methods
	void   Init    (SDPair const & Ini, size_t NIvs=0);
    void   Backup  () { VelBkp=Vel; GraBkp=Gra; IvsBkp=Ivs; }
    void   Restore () { Vel=VelBkp; Gra=GraBkp; Ivs=IvsBkp; }
    size_t PckSize () const { return 2*size(Vel)+size(Ivs); }
    void   Pack    (Array<double>       & V) const;
    void   Unpack  (Array<double> const & V);
    void   Output  (SDPair & KeysVals) const;

	// Data
	Vec_t Vel, VelBkp; ///< Velocity
	Vec_t Gra, GraBkp; ///< Gradient
	Vec_t Ivs, IvsBkp; ///< Internal values
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline FlowState::FlowState (int NDim)
	: State(NDim)
{
    Vel   .change_dim(NDim);  set_to_zero(Vel   );
    Gra   .change_dim(NDim);  set_to_zero(Gra   );
    VelBkp.change_dim(NDim);  set_to_zero(VelBkp);
    GraBkp.change_dim(NDim);  set_to_zero(GraBkp);
}

inline void FlowState::Init (SDPair const & Ini, size_t NIvs)
{
	if (Ini.HasKey("vx")) Vel(0) = Ini("vx");
	if (Ini.HasKey("vy")) Vel(1) = Ini("vy");
	if (num_rows(Vel)>2)
	{
        if (Ini.HasKey("vz")) Vel(2) = Ini("vz");
	}
	else
	{
		bool error = false;
		String key;
		if (Ini.HasKey("vz")) { error=true; key="vz"; }
		if (error) throw new Fatal("FlowState::Init: For a 2D state, there are only 4 stress components. %s is not available",key.CStr());
	}
    if (NIvs>0)
    {
        Ivs.change_dim (NIvs);
        set_to_zero (Ivs);
        for (size_t i=0; i<NIvs; ++i)
        {
            String buf;
            buf.Printf ("z%d",i);
            if (Ini.HasKey(buf)) Ivs(i) = Ini(buf);
        }
        IvsBkp.change_dim (NIvs);
        IvsBkp = Ivs;
    }
}

inline void FlowState::Pack (Array<double> & V) const
{
    size_t nrw = size(Vel);
    size_t niv = size(Ivs);
    V.Resize (2*nrw + niv);
    for (size_t i=0; i<nrw; ++i)
    {
        V[    i] = Vel(i);
        V[nrw+i] = Gra(i);
    }
    for (size_t i=0; i<niv; ++i) V[2*nrw+i] = Ivs(i);
}

inline void FlowState::Unpack (Array<double> const & V)
{
    if (V.Size()!=PckSize()) throw new Fatal("FlowState::Unpack: Size of given vector (%zd) is different of correct size of Pack (%zd)",V.Size(),PckSize());
    size_t nrw = size(Vel);
    size_t niv = size(Ivs);
    for (size_t i=0; i<nrw; ++i)
    {
        Vel(i) = V[    i];
        Gra(i) = V[nrw+i];
    }
    for (size_t i=0; i<niv; ++i) Ivs(i) = V[2*nrw+i];
}

inline void FlowState::Output (SDPair & KeysVals) const
{
    size_t ndim = size(Vel);
    if (ndim==2)
    {
        KeysVals.Set("vx vy  gx gy",
                     Vel(0), Vel(1),
                     Gra(0), Gra(1));
    }
    else
    {
        KeysVals.Set("vx vy vz  gx gy gz",
                     Vel(0), Vel(1), Vel(2),
                     Gra(0), Gra(1), Gra(2));
    }
}


#endif // MECHSYS_FLOWSTATE_H
