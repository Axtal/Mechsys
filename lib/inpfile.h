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
    double dpw, dSw;                        // increment of pore-water pressure and water saturation
    bool   HasDpw, HasDSw;                  // has dpw or dSw ?
    PathIncs () : dsx(0.),dsy(0.),dsz(0.),dsxy(0.),dsyz(0.),dszx(0.), 
                  dex(0.),dey(0.),dez(0.),dexy(0.),deyz(0.),dezx(0.),
                  lode(0.),dp(0.),zPath(false),ninc(-1),k(0.),kPath(false),
                  dpw(0.),dSw(0.),HasDpw(false),HasDSw(false) {}
};

class InpFile
{
public:
    // Constructor
    InpFile ();

    // Methods
    void Read (char const * FileName);

    // Data
    int    matid;       ///<  1 material ID
    int    flwid;       ///<  2 flow material ID
    int    ninc;        ///<  3 general number of increments (for all load-unload paths)
    bool   cdrift;      ///<  4 correct YS drift
    double stol;        ///<  5 local error tolerance
    bool   ssout;       ///<  6 output substeps ?
    bool   ctetg;       ///<  7 Constant stiffness (linear) ?
    bool   fem;         ///<  8 use one Hex8 FEM element instead of point integration
    bool   dyn;         ///<  9 dynamic analysis ?
    bool   hm;          ///< 10 HydroMech ?
    double tf;          ///< 11 final time
    double dt;          ///< 12 time step
    double dtout;       ///< 13 output time step
    double tsw;         ///< 14 switch time (for dynamic simulation)
    int    ndiv;        ///< 15 mesh number of divisions
    int    nip;         ///< 16 Number of integration points in element
    bool   o2;          ///< 17 quadratic elements ?
    bool   ray;         ///< 18 Rayleigh damping ?
    double am;          ///< 19 Damping Am
    double ak;          ///< 20 Damping Ak
    bool   rk;          ///< 21 Runge-Kutta instead of GN22 ?
    String rkscheme;    ///< 22 Runge-Kutta scheme 
    double rkstol;      ///< 23 Runge-Kutta tolerance
    String refdat;      ///< 24 reference data file
    String refsim;      ///< 25 reference simulation file
    String refana;      ///< 26 reference analytical solution file
    int    idxvert1;    ///< 27 index of vertex # 1 for output
    int    idxvert2;    ///< 28 index of vertex # 1 for output
    int    idxvert3;    ///< 29 index of vertex # 1 for output
    double optdbl1;     ///< 30 optional double 1
    double optdbl2;     ///< 31 optional double 2
    double optdbl3;     ///< 32 optional double 3
    bool   hasoptdbl1, hasoptdbl2, hasoptdbl3;
    int    nldt_nsml;   ///< 33 nonlinear timesteps Nsml
    int    nldt_nn;     ///< 34 nonlinear timesteps N
    int    nldt_n;      ///< 35 nonlinear timesteps n
    double nldt_ll;     ///< 36 nonlinear timesteps denominator
    int    nldt_sch;    ///< 37 nonlinear timesteps scheme (0 or 1)
    double nldt_m;      ///< 38 nonlinear timesteps multiplier for larger timesteps in sch==1
    int    maxit;       ///< 39 max num of iterations
    double tolr;        ///< 40 tolerance for residual

    // path increments
    Array<PathIncs> Path;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline InpFile::InpFile ()
    : matid     (-1),     //  1
      flwid     (-1),     //  2
      ninc      (-1),     //  3
      cdrift    (false),  //  4
      stol      (-1),     //  5
      ssout     (false),  //  6
      ctetg     (false),  //  7
      fem       (false),  //  8
      dyn       (false),  //  9
      hm        (false),  // 10
      tf        (-1),     // 11
      dt        (-1),     // 12
      dtout     (-1),     // 13
      tsw       (-1),     // 14
      ndiv      (-1),     // 15
      nip       (-1),     // 16
      o2        (false),  // 17
      ray       (false),  // 18
      am        (-1),     // 19
      ak        (-1),     // 20
      rk        (false),  // 21
      rkscheme  (""),     // 22
      rkstol    (-1),     // 23
      refdat    (""),     // 24
      refsim    (""),     // 25
      refana    (""),     // 26
      idxvert1  (-1),     // 27
      idxvert2  (-1),     // 28
      idxvert3  (-1),     // 29
      optdbl1   (0),      // 30
      optdbl2   (0),      // 31
      optdbl3   (0),      // 32
      hasoptdbl1(false), hasoptdbl2(false), hasoptdbl3(false),
      nldt_nsml (-1),     // 33
      nldt_nn   (-1),     // 34
      nldt_n    (-1),     // 35
      nldt_ll   (-1),     // 36
      nldt_sch  (-1),     // 37
      nldt_m    (-1),     // 38
      maxit     (-1),     // 39
      tolr      (-1)      // 40
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
        String line,key,equal,str_val;
        double val;
        std::getline (inp_file,line);
        std::istringstream iss(line);
        if (iss >> key >> equal >> str_val)
        {
            val = atof(str_val.CStr());
            if (key[0]=='#') { line_num++; continue; }
            if (reading_path)
            {
                if      (key=="ndat") ndat = atoi(str_val.CStr());
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
                else if (key=="dpw")   { Path[idxpath].dpw  = val;  Path[idxpath].HasDpw=true;  idxdat++; }
                else if (key=="dSw")   { Path[idxpath].dSw  = val;  Path[idxpath].HasDSw=true;  idxdat++; }
                else if (key=="ninc")  { Path[idxpath].ninc = atoi(str_val.CStr());  idxdat++; }
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
                if      (key=="matid")      matid     = atoi(str_val.CStr());                    //  1
                else if (key=="flwid")      flwid     = atoi(str_val.CStr());                    //  2
                else if (key=="ninc")       ninc      = atoi(str_val.CStr());                    //  3
                else if (key=="cdrift")     cdrift    = static_cast<bool>(atoi(str_val.CStr())); //  4
                else if (key=="stol")       stol      = val;                                     //  5
                else if (key=="ssout")      ssout     = static_cast<bool>(atoi(str_val.CStr())); //  6
                else if (key=="ctetg")      ctetg     = static_cast<bool>(atoi(str_val.CStr())); //  7
                else if (key=="fem")        fem       = val;                                     //  8
                else if (key=="dyn")        dyn       = static_cast<bool>(atoi(str_val.CStr())); //  9
                else if (key=="hm")         hm        = static_cast<bool>(atoi(str_val.CStr())); // 10
                else if (key=="tf")         tf        = val;                                     // 11
                else if (key=="dt")         dt        = val;                                     // 12
                else if (key=="dtout")      dtout     = val;                                     // 13
                else if (key=="tsw")        tsw       = val;                                     // 14
                else if (key=="ndiv")       ndiv      = atoi(str_val.CStr());                    // 15
                else if (key=="nip")        nip       = atoi(str_val.CStr());                    // 16
                else if (key=="o2")         o2        = static_cast<bool>(atoi(str_val.CStr())); // 17
                else if (key=="ray")        ray       = static_cast<bool>(atoi(str_val.CStr())); // 18
                else if (key=="am")         am        = val;                                     // 19
                else if (key=="ak")         ak        = val;                                     // 20
                else if (key=="rk")         rk        = static_cast<bool>(atoi(str_val.CStr())); // 21
                else if (key=="rkscheme")   rkscheme  = str_val;                                 // 22
                else if (key=="rkstol")     rkstol    = val;                                     // 23
                else if (key=="refdat")     refdat    = str_val;                                 // 24
                else if (key=="refsim")     refsim    = str_val;                                 // 25
                else if (key=="refana")     refana    = str_val;                                 // 26
                else if (key=="idxvert1")   idxvert1  = atoi(str_val.CStr());                    // 27
                else if (key=="idxvert2")   idxvert2  = atoi(str_val.CStr());                    // 28
                else if (key=="idxvert3")   idxvert3  = atoi(str_val.CStr());                    // 29
                else if (key=="optdbl1")  { optdbl1   = val;   hasoptdbl1=true; }                // 30
                else if (key=="optdbl2")  { optdbl2   = val;   hasoptdbl2=true; }                // 31
                else if (key=="optdbl3")  { optdbl3   = val;   hasoptdbl3=true; }                // 32
                else if (key=="nldt_nsml")  nldt_nsml = atoi(str_val.CStr());                    // 33
                else if (key=="nldt_nn")    nldt_nn   = atoi(str_val.CStr());                    // 34
                else if (key=="nldt_n")     nldt_n    = atoi(str_val.CStr());                    // 35
                else if (key=="nldt_ll")    nldt_ll   = val;                                     // 36
                else if (key=="nldt_sch")   nldt_sch  = atoi(str_val.CStr());                    // 37
                else if (key=="nldt_m")     nldt_m    = val;                                     // 38
                else if (key=="maxit")      maxit     = atoi(str_val.CStr());                    // 39
                else if (key=="tolr")       tolr      = val;                                     // 40
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
        os << "de=["<<P.dex << " "<<P.dey << " "<<P.dez << " "<<P.dexy << " "<<P.deyz << " "<<P.dezx << "]";
    }
    if (P.HasDpw) os << " dpw=" << P.dpw << " ";
    if (P.HasDSw) os << " dSw=" << P.dSw << " ";
    os << "\n";
    return os;
}

std::ostream & operator<< (std::ostream & os, Array<PathIncs> const & A)
{
    for (size_t i=0; i<A.Size(); ++i) os << "  " << A[i];
    return os;
}

std::ostream & operator<< (std::ostream & os, InpFile const & IF)
{
    if (IF.matid    >=0) os << "matid     = " << IF.matid     << "\n"; //   1
    if (IF.flwid    >=0) os << "flwid     = " << IF.flwid     << "\n"; //   2
    if (IF.ninc     >=0) os << "ninc      = " << IF.ninc      << "\n"; //   3
    if (IF.cdrift      ) os << "cdrift    = " << IF.cdrift    << "\n"; //   4
    if (IF.stol     >=0) os << "stol      = " << IF.stol      << "\n"; //   5
    if (IF.ssout       ) os << "ssout     = " << IF.ssout     << "\n"; //   6
    if (IF.ctetg       ) os << "ctetg     = " << IF.ctetg     << "\n"; //   7
    if (IF.fem         ) os << "fem       = " << IF.fem       << "\n"; //   8
    if (IF.dyn         ) os << "dyn       = " << IF.dyn       << "\n"; //   9
    if (IF.hm          ) os << "hm        = " << IF.hm        << "\n"; //  10
    if (IF.tf       >=0) os << "tf        = " << IF.tf        << "\n"; //  11
    if (IF.dt       >=0) os << "dt        = " << IF.dt        << "\n"; //  12
    if (IF.dtout    >=0) os << "dtout     = " << IF.dtout     << "\n"; //  13
    if (IF.tsw      >=0) os << "tsw       = " << IF.tsw       << "\n"; //  14
    if (IF.ndiv     >=0) os << "ndiv      = " << IF.ndiv      << "\n"; //  15
    if (IF.nip      >=0) os << "nip       = " << IF.nip       << "\n"; //  16
    if (IF.o2          ) os << "o2        = " << IF.o2        << "\n"; //  17
    if (IF.ray         ) os << "ray       = " << IF.ray       << "\n"; //  18
    if (IF.am       >=0) os << "am        = " << IF.am        << "\n"; //  19
    if (IF.ak       >=0) os << "ak        = " << IF.ak        << "\n"; //  20
    if (IF.rk          ) os << "rk        = " << IF.rk        << "\n"; //  21
    if (IF.rkscheme!="") os << "rkscheme  = " << IF.rkscheme  << "\n"; //  22
    if (IF.rkstol   >=0) os << "rkstol    = " << IF.rkstol    << "\n"; //  23
    if (IF.refdat  !="") os << "refdat    = " << IF.refdat    << "\n"; //  24
    if (IF.refsim  !="") os << "refsim    = " << IF.refsim    << "\n"; //  25
    if (IF.refana  !="") os << "refana    = " << IF.refana    << "\n"; //  26
    if (IF.idxvert1 >=0) os << "idxvert1  = " << IF.idxvert1  << "\n"; //  27
    if (IF.idxvert2 >=0) os << "idxvert2  = " << IF.idxvert2  << "\n"; //  28
    if (IF.idxvert3 >=0) os << "idxvert3  = " << IF.idxvert3  << "\n"; //  29
    if (IF.hasoptdbl1  ) os << "optdbl1   = " << IF.optdbl1   << "\n"; //  30
    if (IF.hasoptdbl2  ) os << "optdbl2   = " << IF.optdbl2   << "\n"; //  31
    if (IF.hasoptdbl3  ) os << "optdbl3   = " << IF.optdbl3   << "\n"; //  32
    if (IF.nldt_nsml >0) os << "nldt_nsml = " << IF.nldt_nsml << "\n"; //  33
    if (IF.nldt_nn   >0) os << "nldt_nn   = " << IF.nldt_nn   << "\n"; //  34
    if (IF.nldt_n    >0) os << "nldt_n    = " << IF.nldt_n    << "\n"; //  35
    if (IF.nldt_ll   >0) os << "nldt_ll   = " << IF.nldt_ll   << "\n"; //  36
    if (IF.nldt_sch  >0) os << "nldt_sch  = " << IF.nldt_sch  << "\n"; //  37
    if (IF.nldt_m    >0) os << "nldt_m    = " << IF.nldt_m    << "\n"; //  38
    if (IF.maxit     >0) os << "maxit     = " << IF.maxit     << "\n"; //  39
    if (IF.tolr     >=0) os << "tolr      = " << IF.tolr      << "\n"; //  40
    return os;
}

#endif // MECHSYS_INPFILE_H
