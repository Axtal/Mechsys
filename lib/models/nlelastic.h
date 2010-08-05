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

#ifndef MECHSYS_NLELASTIC_H
#define MECHSYS_NLELASTIC_H

// MechSys
#include <mechsys/models/model.h>
#include <mechsys/models/equilibstate.h>

using std::cout;
using std::endl;

class NLElastic : public Model
{
public:
    // Constructor
    NLElastic (int NDim, SDPair const & Prms);

    // Methods
    void InitIvs   (SDPair const & Ini, State * Sta)                             const;
    void TgIncs    (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const;
    void InvTgIncs (State const * Sta, Vec_t & DSig, Vec_t & DEps, Vec_t & DIvs) const;
    void Stiffness (State const * Sta, Mat_t & D)                                const;

    // Parameters
    double K0;  ///< Initial bulk modulus 
    double G0;  ///< Initial shear modulus 
    double alp; ///< Nonlinear parameter
    double bet; ///< Nonlinear parameter

    // Data
    Mat_t Psd, IdyI; ///< Isotropic tensors

    // Scratchpad
    mutable Mat_t De;  ///< Stiffness
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline NLElastic::NLElastic (int NDim, SDPair const & Prms)
    : Model (NDim,Prms,"NLElastic")
{
    if (GTy==pse_t) throw new Fatal("NLElastic::NLElastic: This model is not available for plane-stress (pse)");

    // parameters
    K0  = Prms("K0");
    G0  = Prms("G0");
    alp = Prms("alp");
    bet = Prms("bet");

    // isotropic tensors
    Calc_Psd  (NCps, Psd);
    Calc_IdyI (NCps, IdyI);

    // elastic stiffness
    De.change_dim (NCps,NCps);
}

inline void NLElastic::InitIvs (SDPair const & Ini, State * Sta) const
{
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini);
}

inline void NLElastic::TgIncs (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const
{
    Stiffness (Sta, De);
    DSig = De*DEps;
}

inline void NLElastic::InvTgIncs (State const * Sta, Vec_t & DSig, Vec_t & DEps, Vec_t & DIvs) const
{
    /*
    Mat_t Ce(NCps,NCps);
    Stiffness (Sta, De);
    Inv (De, Ce);
    DEps = Ce*DSig;
    */

    EquilibState const * sta = static_cast<EquilibState const *>(Sta);
    double ev = Calc_ev (sta->Eps);
    double ed = Calc_ed (sta->Eps);
    double K  = K0*exp(-alp*ev*100.0);
    double G  = G0*exp(-bet*ed*100.0);
    double ck = (K<1.0e-14 ? 1.0e+14 : (1.0/(9.0*K)));
    double cg = (G<1.0e-14 ? 1.0e+14 : (1.0/(2.0*G)));
    Mat_t Ce(NCps,NCps);
    Ce   = cg*Psd + ck*IdyI;
    DEps = Ce*DSig;

    //cout << ev << "   " << K << endl;
}

inline void NLElastic::Stiffness (State const * Sta, Mat_t & D) const
{
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);
    double ev = Calc_ev (sta->Eps);
    double ed = Calc_ed (sta->Eps);
    double K  = K0*exp(-alp*ev*100.0);
    double G  = G0*exp(-bet*ed*100.0);
    D = (2.0*G)*Psd + K*IdyI;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * NLElasticMaker(int NDim, SDPair const & Prms) { return new NLElastic(NDim,Prms); }

int NLElasticRegister()
{
    ModelFactory["NLElastic"] = NLElasticMaker;
    MODEL.Set ("NLElastic", (double)MODEL.Keys.Size());
    return 0;
}

int __NLElastic_dummy_int = NLElasticRegister();


#endif // MECHSYS_NLELASTIC_H
