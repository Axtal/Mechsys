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
    Model (int NDim, SDPair const & Prms, char const * Name="__unnamed_model__"); ///< NDim:space dimension, Prms:parameters
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

#define STRESSUPDATE_DECLARE
    #include <mechsys/models/stressupdate.h>
    mutable StressUpdate SUp;
#undef STRESSUPDATE_DECLARE
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Model::Model (int TheNDim, SDPair const & ThePrms, char const * TheName)
    : NDim(TheNDim), Prms(ThePrms), GTy(SDPairToGType(ThePrms,(TheNDim==3?"d3d":"d2d"))), 
      Name(TheName), NCps(2*NDim),  NIvs(0)
{
}

std::ostream & operator<< (std::ostream & os, Model const & D)
{
    os << D.Name << " " << D.NDim << "D " << GTypeToStr(D.GTy) << " " << D.Prms;
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


//////////////////////////////////////////////////////////////////////////////////////////////// Functions /////


#ifdef USE_BOOST_PYTHON
double PyMODEL (BPy::str const & Key) { return MODEL(BPy::extract<char const *>(Key)()); }
#endif

inline void ReadMaterial (int Tag, int MatID, const char * FileName, String & ModelName, Dict & Prms, Dict & Inis, bool WithMODEL=true)
{
    // parse materials file
    ModelName = "__empty__";
    size_t idxprm   = 0;
    size_t idxini   = 0;
    size_t nprms    = 0;
    size_t ninis    = 0;
    size_t line_num = 1;
    std::fstream mat_file(FileName, std::ios::in);
    if (!mat_file.is_open()) throw new Fatal("ReadMaterial: Could not open file <%s>",FileName);
    bool reading_model = false;
    bool reading_inis  = false;
    bool model_read    = false;
    while (!mat_file.eof() && !model_read)
    {
        String line,key,equal,strval;
        std::getline (mat_file,line);
        std::istringstream iss(line);
        if (iss >> key >> equal >> strval)
        {
            if (key[0]=='#') { line_num++; continue; }
            if (reading_model)
            {
                if (ModelName=="__empty__")
                {
                    if (key=="name") ModelName = strval;
                    else throw new Fatal("ReadMaterial: Error in file <%s> @ line # %d: 'name' must follow 'ID'. Key==%s is invalid",FileName,line_num,key.CStr());
                    if (WithMODEL) Prms.Set (Tag, "name", MODEL(ModelName));
                }
                else if (nprms==0)
                {
                    if (key=="nprms") nprms = atoi(strval.CStr());
                    else throw new Fatal("ReadMaterial: Error in file <%s> @ line # %d: 'nprms' must follow 'name'. Key==%s is invalid",FileName,line_num,key.CStr());
                }
                else if (key=="ninis") { reading_model=false; ninis=atoi(strval.CStr()); reading_inis=(ninis==0?false:true); model_read=(ninis==0?true:false); }
                else if (idxprm<nprms)
                {
                    Prms.Set (Tag, key.CStr(), atof(strval.CStr()));
                    idxprm++;
                }
                else throw new Fatal("ReadMaterial: Error in file <%s> @ line # %d: there are more parameters than specified by nprms==%d. The reading of parameters finishes when 'ninis' is found. Key==%s is invalid",FileName,line_num,nprms,key.CStr());
            }
            else if (reading_inis)
            {
                if (key=="ID" || key=="name" || key=="nprms" || key=="ninis") throw new Fatal("ReadMaterial: Error in file <%s> @ line # %d: Key==%s found when reading ninis==%d.",FileName,line_num,key.CStr(),ninis);
                if (idxini<ninis)
                {
                    Inis.Set (Tag, key.CStr(), atof(strval.CStr()));
                    idxini++;
                    if (idxini==ninis) { reading_inis=false; model_read=true; }
                }
            }
            else if (key=="ID") { if (atoi(strval.CStr())==MatID) reading_model = true; }
        }
        line_num++;
    }
}

#endif // MECHSYS_MODEL_H
