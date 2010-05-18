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

// MechSys
#include <mechsys/geomtype.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>

class State
{
public:
    State (int NDim) {}
    virtual void Init    (SDPair const & Ini, size_t NIvs=0) =0;
    virtual void Backup  () =0;
    virtual void Restore () =0;
};

class Model
{
public:
    // Constructor
    Model (int NDim, SDPair const & Prms, char const * Name="__unnamed_model__"); ///< NDim:space dimension, Prms:parameters

    // Methods
    virtual void   InitIvs      (SDPair const & Ini, State * Sta)                              const =0;
    virtual void   Stiffness    (State const * Sta, Mat_t & D, Vec_t * h=NULL, Vec_t * d=NULL) const {}
    virtual bool   LoadCond     (State const * Sta, Vec_t const & DEps, double & alpInt)       const { alpInt=-1.0; return false; }
    virtual void   CorrectDrift (State * Sta)                                                  const {}
    virtual double CalcDEz      (State const * Sta, Vec_t const & DSig)                        const { throw new Fatal("Model::CalcDEz: This method is not available yet"); return 0; }

    // Data
    int      NDim; ///< Space dimension: 2 or 3
    SDPair   Prms; ///< Parameters
    GeomType GTy;  ///< Geometry type
    String   Name; ///< Model name. Ex: LinElastic
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Model::Model (int TheNDim, SDPair const & ThePrms, char const * TheName)
    : NDim(TheNDim), Prms(ThePrms), GTy(SDPairToGType(ThePrms,(TheNDim==3?"d3d":"d2d"))), Name(TheName)
{
}

std::ostream & operator<< (std::ostream & os, Model const & D)
{
    os << D.Name << " " << D.NDim << "D " << GTypeToStr(D.GTy) << " " << D.Prms;
    return os;
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


SDPair MODEL;

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

#ifdef USE_BOOST_PYTHON
double PyMODEL (BPy::str const & Key) { return MODEL(BPy::extract<char const *>(Key)()); }
#endif

#endif // MECHSYS_MODEL_H
