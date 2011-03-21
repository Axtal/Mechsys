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

#ifndef MECHSYS_NEOHOOKEAN_H
#define MECHSYS_NEOHOOKEAN_H

// MechSys
#include <mechsys/models/model.h>
#include <mechsys/models/equilibstate.h>

class NeoHookean : public Model
{
public:
    // Constructor
    NeoHookean (int NDim, SDPair const & Prms);

    // Methods
    void InitIvs   (SDPair const & Ini, State * Sta) const;
    void UpdateSta (ATensor2 const & F, State * Sta) const;

    // Data
    double  lam;     ///< Lame constant
    double  G;       ///< Shear modulus
    mutable Vec_t b; ///< left Cauchy-Green tensor = F*F^T
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline NeoHookean::NeoHookean (int NDim, SDPair const & Prms)
    : Model (NDim,Prms,/*niv*/0,"NeoHookean")
{
    // parameters
    double E  = Prms("E");
    double nu = Prms("nu");
    lam = Calc_lam (E, nu);
    G   = Calc_G   (E, nu);
    b.change_dim (NCps);

    // set UseUpdateSta
    UseUpdateSta = true;
}

inline void NeoHookean::InitIvs (SDPair const & Ini, State * Sta) const
{
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini);
}

inline void NeoHookean::UpdateSta (ATensor2 const & F, State * Sta) const
{
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    CalcLCauchyGreen (F, NCps, b); // b = F*Ft
    b -= I; // b = F*Ft - I
    double J = Det (F);
    sta->Sig = (lam*log(J)/J)*I + (G/J)*b;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * NeoHookeanMaker(int NDim, SDPair const & Prms, Model const * O) { return new NeoHookean(NDim,Prms); }

int NeoHookeanRegister()
{
    ModelFactory   ["NeoHookean"] = NeoHookeanMaker;
    MODEL.Set      ("NeoHookean", (double)MODEL.Keys.Size());
    MODEL_PRM_NAMES["NeoHookean"].Resize(2);
    MODEL_PRM_NAMES["NeoHookean"] = "E", "nu";
    MODEL_IVS_NAMES["NeoHookean"].Resize(0);
    return 0;
}

int __NeoHookean_dummy_int = NeoHookeanRegister();


#endif // MECHSYS_NEOHOOKEAN_H
