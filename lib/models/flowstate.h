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
#include "models/model.h"
#include "linalg/matvec.h"

class FlowState : public State
{
public:
	// Constructor
	FlowState (int NDim);

	// Methods
	void Init    (SDPair const & Ini);
    void Backup  () { VelBkp=Vel; GraBkp=Gra; IvsBkp=Ivs; }
    void Restore () { Vel=VelBkp; Gra=GraBkp; Ivs=IvsBkp; }

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

inline void FlowState::Init (SDPair const & Ini)
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
}

#endif // MECHSYS_FLOWSTATE_H