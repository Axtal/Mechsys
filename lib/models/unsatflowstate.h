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
    // static
    static Array<String> Keys;

    // Constructor
    UnsatFlowState (int NDim) : State(NDim) {}

    // Methods
    void   Init    (SDPair const & Ini, size_t NIvs=0);
    void   Backup  () { throw new Fatal("UnsatFlowState::Backup: this method is not available yet"); }
    void   Restore () { throw new Fatal("UnsatFlowState::Restore: this method is not available yet"); }
    size_t PckSize () const { return 6; }
    void   Pack    (Array<double>       & V) const;
    void   Unpack  (Array<double> const & V);

    // Auxiliary methods
    void Output (SDPair & KeysVals) const;

    // Data
    double n;      ///< Porosity
    double Sw;     ///< Saturation
    double pc;     ///< Capillary pressure (pc = pa-pw = -pw)
    bool   Drying; ///< Drying ? (otherwise wetting) : for history-dependent models (hysteresis)
    double RhoW;   ///< Real water density
    double RhoS;   ///< Real solids density
};

Array<String> UnsatFlowState::Keys;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void UnsatFlowState::Init (SDPair const & Ini, size_t NIvs)
{
    n      = Ini("n");
    Sw     = Ini("Sw");
    pc     = Ini("pc");
    Drying = false;
    RhoW   = Ini("RhoW");
    RhoS   = Ini("RhoS");

    if (Keys.Size()==0)
    {
        Keys.Resize(6);
        Keys = "n", "Sw", "pc", "drying", "RhoW", "RhoS";
    }
}

inline void UnsatFlowState::Pack (Array<double> & V) const
{
    V.Resize (PckSize());
    V[0] = n;
    V[1] = Sw;
    V[2] = pc;
    V[3] = Drying;
    V[4] = RhoW;
    V[5] = RhoS;
}

inline void UnsatFlowState::Unpack (Array<double> const & V)
{
    if (V.Size()!=PckSize()) throw new Fatal("UnsatFlowState::Unpack: Size of given vector (%zd) is different of correct size of Pack (%zd)",V.Size(),PckSize());
    n      = V[0];
    Sw     = V[1];
    pc     = V[2];
    Drying = V[3];
    RhoW   = V[4];
    RhoS   = V[5];
}

inline void UnsatFlowState::Output (SDPair & KeysVals) const
{
    KeysVals.Set ("n Sw pc drying RhoW RhoS", n, Sw, pc, (double)Drying, RhoW, RhoS);
}


#endif
