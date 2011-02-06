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

#ifndef MECHSYS_MODEL_H
#define MECHSYS_MODEL_H

// Std Lib
#include <iostream>
#include <fstream>

// MechSys
#include <mechsys/geomtype.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/numerical/odesolver.h>

class State
{
public:
    State (int NDim) {}
    virtual ~State () {}
    virtual void   Init    (SDPair const & Ini, size_t NIvs=0) =0;
    virtual void   Backup  () =0;
    virtual void   Restore () =0;
    virtual size_t PckSize ()                        const =0; ///< Size of pack
    virtual void   Pack    (Array<double>       & V) const =0; ///< Pack all values into V
    virtual void   Unpack  (Array<double> const & V)       =0; ///< Unpack all values from V
};

#include <mechsys/models/equilibstate.h>

class Model
{
public:
    // Constructor & Destructor
    Model (int NDim, SDPair const & Prms, size_t NIvs=0, char const * Name="__unnamed_model__"); ///< NDim:space dimension, Prms:parameters
    virtual ~Model () {}

    // Methods
    virtual void   InitIvs      (SDPair const & Ini, State * Sta)                              const =0;
    virtual void   TgIncs       (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs)  const { throw new Fatal("Model::TgIncs: This method is not available in this model (%s)",Name.CStr()); }
    virtual void   InvTgIncs    (State const * Sta, Vec_t & DSig, Vec_t & DEps, Vec_t & DIvs)  const { throw new Fatal("Model::InvTgIncs: This method is not available in this model (%s)",Name.CStr()); }
    virtual void   Stiffness    (State const * Sta, Mat_t & D)                                 const { throw new Fatal("Model::Stiffness: This method is not available in this model (%s)",Name.CStr()); }
    virtual void   Stiffness    (State const * Sta, Mat_t & D, Vec_t & Dw)                     const { throw new Fatal("Model::Stiffness(Dw): This method is not available in this model (%s)",Name.CStr()); }
    virtual void   Hydraulic    (State const * Sta, Mat_t & Kw, double & ChiW, double & InvQs) const { throw new Fatal("Model::Hydraulic: This method is not available in this model (%s)",Name.CStr()); }
    virtual bool   LoadCond     (State const * Sta, Vec_t const & DEps, double & alpInt)       const { alpInt=-1; return false; }
    virtual size_t CorrectDrift (State * Sta)                                                  const { return 0; }
    virtual void   UpdatePath   (State const * Sta, Vec_t const & DEps, Vec_t const & DSig)    const {}

    // Data
    int           NDim;    ///< Space dimension: 2 or 3
    SDPair        Prms;    ///< Parameters
    GeomType      GTy;     ///< Geometry type
    String        Name;    ///< Model name. Ex: LinElastic
    size_t        NCps;    ///< Number of stress/strain components
    size_t        NIvs;    ///< Number of internal values
    Array<String> IvNames; ///< Names of internal values

    // General constants
    double Grav;     ///< Gravity. Ex.: 9.81 m/s2
    double GamW;     ///< Unit weight of water. Ex.: 9.81 kN/m3
    double Rho;      ///< Density of the material. Ex.: 2.0 Mg/m3
    double RhoS;     ///< Density of solids. Ex.: 2.0 Mg/m3
    double Por;      ///< Porosity of material. Ex.: 0.2
    double GamNat;   ///< Natural unit weight of geo-material: Ex.: 20.0 kN/m3
    double GamSat;   ///< Saturated unit weight of geo-material: Ex.: 21.0 kN/m3

#define STRESSUPDATE_DECLARE
    #include <mechsys/models/stressupdate.h>
    mutable StressUpdate SUp;
#undef STRESSUPDATE_DECLARE
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Model::Model (int TheNDim, SDPair const & ThePrms, size_t NIv, char const * TheName)
    : NDim(TheNDim), Prms(ThePrms), GTy(SDPairToGType(ThePrms,(TheNDim==3?"d3d":"d2d"))), 
      Name(TheName), NCps(2*NDim),  NIvs(NIv)
{
    // constants
    Grav   = (Prms.HasKey("grav")   ? Prms("grav")   : -1);
    GamW   = (Prms.HasKey("gamW")   ? Prms("gamW")   : -1);
    Rho    = (Prms.HasKey("rho")    ? Prms("rho")    : -1);
    RhoS   = (Prms.HasKey("rhoS")   ? Prms("rhoS")   : -1);
    Por    = (Prms.HasKey("por")    ? Prms("por")    : -1);
    GamNat = (Prms.HasKey("gamNat") ? Prms("gamNat") : -1);
    GamSat = (Prms.HasKey("gamSat") ? Prms("gamSat") : -1);
    if (Prms.HasKey("gamNat") && !Prms.HasKey("gamSat")) GamSat = GamNat;
    if (Prms.HasKey("gamSat") && !Prms.HasKey("gamNat")) GamNat = GamSat;
    if (Prms.HasKey("rhoS"))
    {
        double rho_w = Prms("gamW") / Prms("grav");
        double n     = Prms("por");
        Rho = n*rho_w + (1.0-n)*Prms("rhoS");
        if (!Prms.HasKey("gamNat")) GamNat = Rho * Grav;
        if (!Prms.HasKey("gamSat")) GamSat = GamNat;
    }

    // stress update
    SUp.SetModel   (this);
    IvNames.Resize (NIv);
}

std::ostream & operator<< (std::ostream & os, Model const & D)
{
    os << D.Name << " " << D.NDim << "D " << GTypeToStr(D.GTy) << std::endl;
    os << "  Prms = {" << D.Prms << "}\n";
    os << "  Gravity                               : grav   = " << D.Grav   << std::endl;
    os << "  Unit weight of water                  : gamW   = " << D.GamW   << std::endl;
    os << "  Density of the material               : rho    = " << D.Rho    << std::endl;
    os << "  Density of solids                     : rhoS   = " << D.RhoS   << std::endl;
    os << "  Porosity of material                  : por    = " << D.Por    << std::endl;
    os << "  Natural unit weight of geo-material   : gamNat = " << D.GamNat << std::endl;
    os << "  Saturated unit weight of geo-material : gamSat = " << D.GamSat << std::endl;
    return os;
}


#define STRESSUPDATE_IMPLEMENT
  #include <mechsys/models/stressupdate.h>
#undef STRESSUPDATE_IMPLEMENT


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


SDPair MODEL; ///< Model name

typedef std::map<String,Array<String> > Str2ArrayStr_t; ///< Map string to array of strings

Str2ArrayStr_t MODEL_PRM_NAMES; ///< Model name => Parameters names
Str2ArrayStr_t MODEL_IVS_NAMES; ///< Model name => Initial values names

typedef Model * (*ModelMakerPtr)(int NDim, SDPair const & Prms);

typedef std::map<String, ModelMakerPtr> ModelFactory_t;

ModelFactory_t ModelFactory;

Model * AllocModel(String const & Name, int NDim, SDPair const & Prms)
{
    ModelFactory_t::iterator it = ModelFactory.find(Name);
    if (it==ModelFactory.end()) throw new Fatal("AllocModel: '%s' is not available", Name.CStr());

    Model * ptr = (*it->second)(NDim,Prms);

    return ptr;
}

Model * AllocModel(double IDinMODEL, int NDim, SDPair const & Prms)
{
    String model_name;
    MODEL.Val2Key (IDinMODEL, model_name);
    return AllocModel (model_name, NDim, Prms);
}

Array<String> MODEL_CTE_NAMES; ///< Constants names

int ModelRegister()
{
    MODEL_CTE_NAMES.Resize(7);
    MODEL_CTE_NAMES = "grav", "gamW", "rho", "rhoS", "por", "gamNat", "gamSat";
    return 0;
}

int __Model_dummy_int = ModelRegister();


//////////////////////////////////////////////////////////////////////////////////////////////// Functions /////


#ifdef USE_BOOST_PYTHON
double PyMODEL (BPy::str const & Key) { return MODEL(BPy::extract<char const *>(Key)()); }
#endif

#endif // MECHSYS_MODEL_H
