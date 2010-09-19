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
#include <mechsys/util/string.h>
#include <mechsys/models/model.h>
#include <mechsys/linalg/matvec.h>

class EquilibState : public State
{
public:
    // static
    static Array<String> Keys;

    // Constructor
    EquilibState (int NDim);

    // Methods
    void Init    (SDPair const & Ini, size_t NIvs=0);
    void Backup  () { SigBkp=Sig; EpsBkp=Eps; IvsBkp=Ivs; LdgBkp=Ldg; }
    void Restore () { Sig=SigBkp; Eps=EpsBkp; Ivs=IvsBkp; Ldg=LdgBkp; }

    // Auxiliar methods
    void Output (std::ostream & os, bool WithHeader=false, char const * NF="%13g") const;

    // Operators
    void operator= (EquilibState const & Another);

    // Data
    Vec_t Sig, SigBkp; ///< Stress
    Vec_t Eps, EpsBkp; ///< Strain
    Vec_t Ivs, IvsBkp; ///< Internal values
    bool  Ldg, LdgBkp; ///< Loading ?
};

Array<String> EquilibState::Keys;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline EquilibState::EquilibState (int NDim)
    : State(NDim), Ldg(false), LdgBkp(false)
{
    int ncomp = NDim*2; // number of stress/strain components
    Sig   .change_dim(ncomp);  set_to_zero(Sig   );
    Eps   .change_dim(ncomp);  set_to_zero(Eps   );
    SigBkp.change_dim(ncomp);  set_to_zero(SigBkp);
    EpsBkp.change_dim(ncomp);  set_to_zero(EpsBkp);

    if (Keys.Size()==0)
    {
        Keys.Resize (2*ncomp+4);
        if (NDim==2) 
        {
            Keys = "sx",   "sy",   "sz", "sxy",
                   "ex",   "ey",   "ez", "exy",
                   "pcam", "qcam", "ev", "ed";
        }
        else
        {
            Keys = "sx",   "sy",   "sz", "sxy", "syz", "szx",
                   "ex",   "ey",   "ez", "exy", "eyz", "ezx",
                   "pcam", "qcam", "ev", "ed";
        }
    }
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
    if (Ini.HasKey("zero"))
    {
        set_to_zero (Sig   );
        set_to_zero (Eps   );
        set_to_zero (SigBkp);
        set_to_zero (EpsBkp);
        if (NIvs>0)
        {
            set_to_zero (Ivs   );
            set_to_zero (IvsBkp);
        }
    }
}

inline void EquilibState::Output (std::ostream & os, bool WithHeader, char const * NF) const
{
    size_t ncp = size(Sig);
    size_t niv = size(Ivs);
    String buf;
    if (WithHeader)
    {
        String nf, str;
        nf.TextFmt(NF);
        if (ncp>4)
        {
            str.Printf(nf.CStr(),"sx");  buf.append(str);
            str.Printf(nf.CStr(),"sy");  buf.append(str);
            str.Printf(nf.CStr(),"sz");  buf.append(str);
            str.Printf(nf.CStr(),"sxy"); buf.append(str);
            str.Printf(nf.CStr(),"syz"); buf.append(str);
            str.Printf(nf.CStr(),"szx"); buf.append(str);
            str.Printf(nf.CStr(),"ex");  buf.append(str);
            str.Printf(nf.CStr(),"ey");  buf.append(str);
            str.Printf(nf.CStr(),"ez");  buf.append(str);
            str.Printf(nf.CStr(),"exy"); buf.append(str);
            str.Printf(nf.CStr(),"eyz"); buf.append(str);
            str.Printf(nf.CStr(),"ezx"); buf.append(str);
        }
        else
        {
            str.Printf(nf.CStr(),"sx");  buf.append(str);
            str.Printf(nf.CStr(),"sy");  buf.append(str);
            str.Printf(nf.CStr(),"sz");  buf.append(str);
            str.Printf(nf.CStr(),"sxy"); buf.append(str);
            str.Printf(nf.CStr(),"ex");  buf.append(str);
            str.Printf(nf.CStr(),"ey");  buf.append(str);
            str.Printf(nf.CStr(),"ez");  buf.append(str);
            str.Printf(nf.CStr(),"exy"); buf.append(str);
        }
        for (size_t i=0; i<niv; ++i)
        { 
            str.Printf("z%d",i);
            str.Printf(nf.CStr(), str.CStr()); buf.append(str);
        }
        os << buf << "\n";
    }
    for (size_t i=0; i<ncp; ++i) { buf.Printf(NF,Sig(i)); os<<buf; }
    for (size_t i=0; i<ncp; ++i) { buf.Printf(NF,Eps(i)); os<<buf; }
    for (size_t i=0; i<niv; ++i) { buf.Printf(NF,Ivs(i)); os<<buf; }
    os << "\n";
}

inline void EquilibState::operator= (EquilibState const & A)
{
    Sig = A.Sig;  SigBkp = A.SigBkp;
    Eps = A.Eps;  EpsBkp = A.EpsBkp;
    Ivs = A.Ivs;  IvsBkp = A.IvsBkp;
    Ldg = A.Ldg;  LdgBkp = A.LdgBkp;
}

#endif // MECHSYS_EQUILIBSTATE_H
