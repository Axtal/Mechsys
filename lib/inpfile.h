/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

#ifndef MECHSYS_INPFILE_H
#define MECHSYS_INPFILE_H

// Std Lib
#include <fstream>

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

struct PathIncs
{
    double dsx, dsy, dsz, dsxy, dsyz, dszx; // stress increments
    double dex, dey, dez, dexy, deyz, dezx; // strain increments (percentage)
    double lode, dp;                        // path given Lode angle (deg), dpoct and dez (percentage)
    bool   zPath;                           // use lode, dp and dez ?
    int    ninc;                            // number of increments for this path. -1 => use general
    double k;                               // path given Lode, k=dqoct/dpoct, and dez
    bool   kPath;                           // with k=Dq/Dp
    PathIncs () : dsx(0.),dsy(0.),dsz(0.),dsxy(0.),dsyz(0.),dszx(0.), 
                  dex(0.),dey(0.),dez(0.),dexy(0.),deyz(0.),dezx(0.),
                  lode(0.),dp(0.),zPath(false),ninc(-1),k(0.),kPath(false) {}
};

class InpFile
{
public:
    // Constructor
    InpFile ();

    // Methods
    void Read (char const * FileName);

    // Data
    int    MatID;  // material ID
    double pCam0;  // pCam
    size_t NInc;   // general number of increments (for all load-unload paths)
    bool   CDrift; // correct YS drift
    double STOL;   // local error tolerance
    bool   FEM;    // use one Hex8 FEM element instead of point integration
    bool   SSOut;  // output substeps ?

    // path increments
    Array<PathIncs> Path;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline InpFile::InpFile ()
    : MatID  (0),
      pCam0  (100.0),
      NInc   (10),
      CDrift (true),
      STOL   (1.0e-5),
      FEM    (false),
      SSOut  (true)
{
}

inline void InpFile::Read (char const * FileName)
{
    // parse input file
    std::fstream inp_file(FileName, std::ios::in);
    if (!inp_file.is_open()) throw new Fatal("InpFile::Read: Could not open file <%s>",FileName);
    bool   reading_path = false;
    int    ndat         = -1;
    size_t line_num     = 1;
    int    idxdat       = 0;
    size_t idxpath      = 0;
    while (!inp_file.eof())
    {
        String line,key,equal;
        double val;
        std::getline (inp_file,line);
        std::istringstream iss(line);
        if (iss >> key >> equal >> val)
        {
            if (key[0]=='#') { line_num++; continue; }
            if (reading_path)
            {
                if      (key=="ndat") ndat = (int)val;
                else if (ndat<0) throw new Fatal("InpFile::Read: Error in file <%s> @ line # %d: key 'ndat' must come before data. Key==%s is in the wrong position",FileName,line_num,key.CStr());
                else if (key=="kcam")  { Path[idxpath].k    = Util::SQ2*val/3.; Path[idxpath].kPath=true; Path[idxpath].zPath=false; idxdat++; }
                else if (key=="dpcam") { Path[idxpath].dp   = val*Util::SQ3;    Path[idxpath].zPath=true; Path[idxpath].kPath=false; idxdat++; }
                else if (key=="lode")  { Path[idxpath].lode = val;       idxdat++; if (val<30. || val>90.) throw new Fatal("Lode angle alpha must be inside [30,90]. Alpha==%g is invalid",val); }
                else if (key=="dex")   { Path[idxpath].dex  = val/100.;  idxdat++; }
                else if (key=="dey")   { Path[idxpath].dey  = val/100.;  idxdat++; }
                else if (key=="dez")   { Path[idxpath].dez  = val/100.;  idxdat++; }
                else if (key=="dexy")  { Path[idxpath].dexy = val/100.;  idxdat++; }
                else if (key=="deyz")  { Path[idxpath].deyz = val/100.;  idxdat++; }
                else if (key=="dezx")  { Path[idxpath].dezx = val/100.;  idxdat++; }
                else if (key=="dsx")   { Path[idxpath].dsx  = val;       idxdat++; }
                else if (key=="dsy")   { Path[idxpath].dsy  = val;       idxdat++; }
                else if (key=="dsz")   { Path[idxpath].dsz  = val;       idxdat++; }
                else if (key=="dsxy")  { Path[idxpath].dsxy = val;       idxdat++; }
                else if (key=="dsyz")  { Path[idxpath].dsyz = val;       idxdat++; }
                else if (key=="dszx")  { Path[idxpath].dszx = val;       idxdat++; }
                else if (key=="ninc")  { Path[idxpath].ninc = (int)val;  idxdat++; }
                else throw new Fatal("InpFile::Read: Error in file<%s> @ line # %d: reading data of Path # %d. Key==%s is invalid",FileName,line_num,idxpath,key.CStr());
                if (idxdat==ndat)
                {
                    ndat   = -1;
                    idxdat = 0;
                    idxpath++;
                    if (idxpath==Path.Size()) break;
                }
            }
            else
            {
                if      (key=="matid")  MatID  = (int)val;
                else if (key=="pcam0")  pCam0  = val;
                else if (key=="ninc")   NInc   = (int)val;
                else if (key=="cdrift") CDrift = (bool)val;
                else if (key=="stol")   STOL   = val;
                else if (key=="fem")    FEM    = val;
                else if (key=="ssout")  SSOut  = (bool)val;
                else if (key=="npath")
                {
                    Path.Resize ((size_t)val);
                    reading_path = true;
                }
                else throw new Fatal("InpFile::Read: Error in file <%s> @ line # %d: Key==%s in invalid",FileName,line_num,key.CStr());
            }
        }
        line_num++;
    }
    if (idxpath!=Path.Size()) throw new Fatal("InpFile::Read: Error in file <%s>: Not all Path data could be read for npath==%d",FileName,Path.Size());
}

std::ostream & operator<< (std::ostream & os, PathIncs const & P)
{
    os << "ninc="<<P.ninc << " ";
    if      (P.kPath) { os << "lode="<<P.lode << " kcam="<<P.k*3./Util::SQ2 << " dez="<<P.dez << std::endl; }
    else if (P.zPath) { os << "lode="<<P.lode << " dpcam="<<P.dp/Util::SQ3  << " dez="<<P.dez << std::endl; }
    else
    {
        os << "ds=["<<P.dsx << " "<<P.dsy << " "<<P.dsz << " "<<P.dsxy << " "<<P.dsyz << " "<<P.dszx << "] ";
        os << "de=["<<P.dex << " "<<P.dey << " "<<P.dez << " "<<P.dexy << " "<<P.deyz << " "<<P.dezx << "]\n";
    }
    return os;
}

std::ostream & operator<< (std::ostream & os, Array<PathIncs> const & A)
{
    for (size_t i=0; i<A.Size(); ++i) os << "  " << A[i];
    return os;
}

std::ostream & operator<< (std::ostream & os, InpFile const & IF)
{
    os << "Input data:\n";
    os << "  matid  = " <<  IF.MatID   << "\n";
    os << "  pcam0  = " <<  IF.pCam0   << "\n";
    os << "  ninc   = " <<  IF.NInc    << "\n";
    os << "  cdrift = " <<  IF.CDrift  << "\n";
    os << "  stol   = " <<  IF.STOL    << "\n";
    os << "  fem    = " <<  IF.FEM     << "\n";
    os << "  ssout  = " <<  IF.SSOut   << "\n";
    return os;
}

#endif // MECHSYS_INPFILE_H
