/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
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

#ifndef MECHSYS_UNSATFLOWSTATE_H
#define MECHSYS_UNSATFLOWSTATE_H

// MechSys
#include <mechsys/models/model.h>

class UnsatFlowState : public State
{
public:
    // Constructor
    UnsatFlowState (int NDim) : State(NDim) {}

    // Methods
    void   Init    (SDPair const & Ini, size_t NIvs=0);
    void   Backup  () { throw new Fatal("UnsatFlowState::Backup: this method is not available yet"); }
    void   Restore () { throw new Fatal("UnsatFlowState::Restore: this method is not available yet"); }
    size_t PckSize () const { return 4; }
    void   Pack    (Array<double>       & V) const;
    void   Unpack  (Array<double> const & V);

    // Data
    double n;      ///< Porosity
    double Sw;     ///< Saturation
    double pc;     ///< Capillary pressure (pc = pa-pw = -pw)
    bool   Drying; ///< Drying ? (otherwise wetting) : for history-dependent models (hysteresis)
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void UnsatFlowState::Init (SDPair const & Ini, size_t NIvs)
{
    n      = Ini("n");
    pc     = Ini("pc");
    Sw     = Ini("Sw");
    Drying = false;
}

inline void UnsatFlowState::Pack (Array<double> & V) const
{
    V.Resize (PckSize());
    V[0] = n;
    V[1] = Sw;
    V[2] = pc;
    V[3] = Drying;
}

inline void UnsatFlowState::Unpack (Array<double> const & V)
{
    if (V.Size()!=PckSize()) throw new Fatal("UnsatFlowState::Unpack: Size of given vector (%zd) is different of correct size of Pack (%zd)",V.Size(),PckSize());
    n      = V[0];
    Sw     = V[1];
    pc     = V[2];
    Drying = V[3];
}


#endif
