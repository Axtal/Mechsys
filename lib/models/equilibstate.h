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

#ifndef MECHSYS_EQUILIBSTATE_H
#define MECHSYS_EQUILIBSTATE_H

// Std Lib
#include <iostream>
#include <cmath>   // for sqrt

// MechSys
#include <mechsys/models/model.h>
#include <mechsys/linalg/matvec.h>

class EquilibState : public State
{
public:
	// Constructor
	EquilibState (int NDim);

	// Methods
	void Init    (SDPair const & Ini, size_t NIvs=0);
    void Backup  () { SigBkp=Sig; EpsBkp=Eps; IvsBkp=Ivs; }
    void Restore () { Sig=SigBkp; Eps=EpsBkp; Ivs=IvsBkp; }

    // Operators
    void operator= (EquilibState const & Another);

	// Data
	Vec_t Sig, SigBkp; ///< Stress
	Vec_t Eps, EpsBkp; ///< Strain
	Vec_t Ivs, IvsBkp; ///< Internal values
    bool  Ldg;         ///< Loading ?
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline EquilibState::EquilibState (int NDim)
	: State(NDim), Ldg(false)
{
    int ncomp = NDim*2; // number of stress/strain components
    Sig   .change_dim(ncomp);  set_to_zero(Sig   );
    Eps   .change_dim(ncomp);  set_to_zero(Eps   );
    SigBkp.change_dim(ncomp);  set_to_zero(SigBkp);
    EpsBkp.change_dim(ncomp);  set_to_zero(EpsBkp);
}

inline void EquilibState::Init (SDPair const & Ini, size_t NIvs)
{
	if (Ini.HasKey("sx"))  Sig(0) = Ini("sx");
	if (Ini.HasKey("sy"))  Sig(1) = Ini("sy");
	if (Ini.HasKey("sz"))  Sig(2) = Ini("sz");
	if (Ini.HasKey("sxy")) Sig(3) = Ini("sxy")*sqrt(2.0);
	if (num_rows(Sig)>4)
	{
		if (Ini.HasKey("syz")) Sig(4) = Ini("syz")*sqrt(2.0);
		if (Ini.HasKey("sxz")) Sig(5) = Ini("sxz")*sqrt(2.0);
	}
	else
	{
		bool error = false;
		String key;
		if (Ini.HasKey("syz")) { error=true; key="syz"; }
		if (Ini.HasKey("sxz")) { error=true; key="sxz"; }
		if (error) throw new Fatal("EquilibState::Init: For a 2D state, there are only 4 stress components. %s is not available",key.CStr());
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

inline void EquilibState::operator= (EquilibState const & A)
{
    Sig = A.Sig;  SigBkp = A.SigBkp;
    Eps = A.Eps;  EpsBkp = A.EpsBkp;
    Ivs = A.Ivs;  IvsBkp = A.IvsBkp;
    Ldg = A.Ldg;
}

#endif // MECHSYS_EQUILIBSTATE_H
