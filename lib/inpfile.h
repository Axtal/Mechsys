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

// wxWidgets
#ifdef USE_WXWIDGETS
  #include <mechsys/gui/common.h>
  #include <mechsys/gui/wxdict.h>
  #include <mechsys/gui/wxsipair.h>
  #include <mechsys/gui/wxarrayint.h>
#endif

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/fem/element.h>   // for PROB
#include <mechsys/fem/geomelem.h>  // for GEOM
#include <mechsys/fem/solver.h>
#include <mechsys/matfile.h>
#include <mechsys/models/model.h>

/*
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
*/

#ifdef USE_WXWIDGETS
class InpFile : public wxWindow
#else
class InpFile
#endif
{
public:
    // Constructor & Destructor
#ifdef USE_WXWIDGETS
     InpFile (wxFrame * Parent);
#else
     InpFile ();
#endif
    ~InpFile ();

    // Methods
    void Defaults    ();
    void Read        (char const * FileName);
    void SetPrmsInis (MatFile const & Mat, bool ForceGTY=false);
    void GetIncs     (int PathKey, double Div, Vec_t & dsig, Vec_t & deps, Array<bool> & PrescDeps, double & dpw, double & dSw, bool & PrescDpw, bool & PrescDSw) const;
    void SetSolver   (FEM::Solver & Sol) const;
    void SetSUp      (Model const * Mdl, Model::StressUpdate::pDbgFun pFun=NULL, void * UserData=NULL) const;

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
    String fnkey;       ///< 41 filename key
    double pcam0;       ///< 42 pcam0
    bool   haspcam0;    ///< has pcam0 ?
    String scheme;      ///< 43 solver scheme
    bool   vtufile;     ///< 44 write vtu file ?
    String suscheme;    ///< 45 stress-update scheme
    double sustol;      ///< 46 stress-update STOL
    String surkscheme;  ///< 47 stress-update RK scheme
    size_t dcmaxit;     ///< 48 drift correction max iterations
    double dcftol;      ///< 49 drift correction f tolerance

    // Additional data
    Dict * Prms; ///< parameters (set by SetMat)
    Dict * Inis; ///< initial values (set by SetMat)

#ifdef USE_WXWIDGETS
    // Methods
    void Sync (bool Dat2Ctrl=false) { if (Dat2Ctrl) TransferDataToWindow(); else TransferDataFromWindow(); } ///< Synchronise (validate/transfer) data in controls

    // Data
    wxAuiManager       Aui;    ///< Aui manager
    wxString           LstDir; ///< Last accessed directory
    String             FName;  ///< Input file (.inp) filename
    GUI::WxDictTable * Path;   ///< Path increments
    GUI::WxDict      * GPath;  ///< Grid for editing the path

    // Additional data
    GUI::WxDictTable     * Prps;       ///< elements properties
    GUI::WxArrayIntTable * OutNods;    ///< output nodes
    GUI::WxDict          * GPrps;      ///< grid: elements properties
    GUI::WxArrayInt      * GOutNods;   ///< grid: output nodes
    GUI::WxSIPairTable   * MatId2Tag;  ///< material ID to element tag
    GUI::WxSIPairTable   * XMaId2Tag;  ///< material ID to element tag
    GUI::WxSIPair        * GMatId2Tag; ///< material ID to element tag

    // Boundary conditions
    Array<GUI::WxDictTable*> Stages;  ///< boundary conditions
    //Array<GUI::WxDict*>      GStages; ///< grid: boundary conditions

    // Events
    void OnLoad (wxCommandEvent & Event);
    void OnSave (wxCommandEvent & Event);
    DECLARE_EVENT_TABLE();

#else
    // Data
    Dict * Path; ///< Path increments

    // Additional data
    Dict       * Prps;      ///< elements properties
    Array<int> * OutNods;   ///< output nodes
    SIPair     * MatId2Tag; ///< material ID to tag
    SIPair     * XMaId2Tag; ///< X material ID to tag

    // Boundary conditions
    Array<Dict*> Stages; ///< boundary conditions
#endif

#ifdef USE_BOOST_PYTHON
    BPy::tuple PyReadPrmIni (char const * MatFN, char const * InpFN, int Tag)
    {
        BPy::dict prm, ini;
        MatFile   mat;
        mat.Read    (MatFN);
        Read        (InpFN);
        SetPrmsInis (mat);
        Prms->PyGet (Tag, prm);
        Inis->PyGet (Tag, ini);
        return BPy::make_tuple (prm, ini);
    }
#endif
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


#ifdef USE_WXWIDGETS

inline InpFile::~InpFile ()
{
    Aui.UnInit ();
    delete Prms;
    delete Inis;
}

#else

inline InpFile::InpFile ()
{
    Defaults ();
    Path      = new Dict;
    Prps      = new Dict;
    OutNods   = new Array<int>;
    MatId2Tag = new SIPair;
    XMaId2Tag = new SIPair;
    Prms      = new Dict;
    Inis      = new Dict;
}

inline InpFile::~InpFile ()
{
    delete Path;
    delete Prps;
    delete OutNods;
    delete MatId2Tag;
    delete XMaId2Tag;
    delete Prms;
    delete Inis;
    for (size_t i=0; i<Stages.Size(); ++i) delete Stages[i];
}

#endif

inline void InpFile::Defaults ()
{
    matid      = -1;     //  1
    flwid      = -1;     //  2
    ninc       = -1;     //  3
    cdrift     = false;  //  4
    stol       = -1;     //  5
    ssout      = false;  //  6
    ctetg      = false;  //  7
    fem        = false;  //  8
    dyn        = false;  //  9
    hm         = false;  // 10
    tf         = -1;     // 11
    dt         = -1;     // 12
    dtout      = -1;     // 13
    tsw        = -1;     // 14
    ndiv       = -1;     // 15
    nip        = -1;     // 16
    o2         = false;  // 17
    ray        = false;  // 18
    am         = -1;     // 19
    ak         = -1;     // 20
    rk         = false;  // 21
    rkscheme   = "";     // 22
    rkstol     = -1;     // 23
    refdat     = "";     // 24
    refsim     = "";     // 25
    refana     = "";     // 26
    idxvert1   = -1;     // 27
    idxvert2   = -1;     // 28
    idxvert3   = -1;     // 29
    optdbl1    = 0;      // 30
    optdbl2    = 0;      // 31
    optdbl3    = 0;      // 32
    hasoptdbl1 = false; hasoptdbl2=false; hasoptdbl3=false;
    nldt_nsml  = -1;     // 33
    nldt_nn    = -1;     // 34
    nldt_n     = -1;     // 35
    nldt_ll    = -1;     // 36
    nldt_sch   = -1;     // 37
    nldt_m     = -1;     // 38
    maxit      = -1;     // 39
    tolr       = -1;     // 40
    fnkey      = "";     // 41
    pcam0      = 0;      // 42
    haspcam0   = false;
    scheme     = "";     // 43
    vtufile    = false;  // 44
    suscheme   = "";     // 45
    sustol     = -1;     // 46
    surkscheme = "";     // 47
    dcmaxit    = 0;      // 48
    dcftol     = -1;     // 49
}

inline void InpFile::Read (char const * FileName)
{
    // parse input file
    std::fstream inp_file(FileName, std::ios::in);
    if (!inp_file.is_open()) throw new Fatal("InpFile::Read: Could not open file <%s>",FileName);
    bool   reading_path   = false;
    bool   reading_eprps  = false;
    bool   reading_stages = false;
    bool   reading_bcs    = false;
    int    npath          = 0;
    int    nelemprps      = 0;
    int    nstages        = 0;
    int    nbcs           = -1;
    int    ndat           = -1;
    size_t line_num       = 1;
    int    idxdat         = 0;
    int    idxbcs         = 0;
    int    idxpath        = 0;
    int    idxeprps       = 0;
    int    idxstage       = 0;
    int    elemtag        = -1;
    int    bcstag         = 0;
    SDPair elemprps;
    SDPair bcs;
    Path      -> clear();
    Prps      -> clear();
    OutNods   -> Resize(0);
    MatId2Tag -> clear();
    XMaId2Tag -> clear();
    Prms      -> clear();
    Inis      -> clear();
    for (size_t i=0; i<Stages.Size(); ++i) delete Stages[i];
    Stages.Resize(0);
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
                else if (ndat<0) throw new Fatal("InpFile::Read: Reading path. Error in file <%s> at line # %d: key 'ndat' must come before data. '%s' is in the wrong place",FileName,line_num,key.CStr());
                else if (key=="kcam")  { Path->Set(idxpath, "kcam" , val     ); idxdat++; }
                else if (key=="dpcam") { Path->Set(idxpath, "dpcam", val     ); idxdat++; }
                else if (key=="lode")  { Path->Set(idxpath, "lode" , val     ); idxdat++; if (val<30. || val>90.) throw new Fatal("InpFile::Read: Error in file <%s> at line # %d: Lode angle alpha must be inside [30,90]. Alpha==%g is invalid",FileName,line_num,val); }
                else if (key=="dex")   { Path->Set(idxpath, "dex"  , val/100.); idxdat++; }
                else if (key=="dey")   { Path->Set(idxpath, "dey"  , val/100.); idxdat++; }
                else if (key=="dez")   { Path->Set(idxpath, "dez"  , val/100.); idxdat++; }
                else if (key=="dexy")  { Path->Set(idxpath, "dexy" , val/100.); idxdat++; }
                else if (key=="deyz")  { Path->Set(idxpath, "deyz" , val/100.); idxdat++; }
                else if (key=="dezx")  { Path->Set(idxpath, "dezx" , val/100.); idxdat++; }
                else if (key=="dsx")   { Path->Set(idxpath, "dsx"  , val     ); idxdat++; }
                else if (key=="dsy")   { Path->Set(idxpath, "dsy"  , val     ); idxdat++; }
                else if (key=="dsz")   { Path->Set(idxpath, "dsz"  , val     ); idxdat++; }
                else if (key=="dsxy")  { Path->Set(idxpath, "dsxy" , val     ); idxdat++; }
                else if (key=="dsyz")  { Path->Set(idxpath, "dsyz" , val     ); idxdat++; }
                else if (key=="dszx")  { Path->Set(idxpath, "dszx" , val     ); idxdat++; }
                else if (key=="dpw")   { Path->Set(idxpath, "dpw"  , val     ); idxdat++; }
                else if (key=="dSw")   { Path->Set(idxpath, "dSw"  , val     ); idxdat++; }
                else if (key=="ninc")  { Path->Set(idxpath, "ninc" , val     ); idxdat++; }
                else throw new Fatal("InpFile::Read: Reading path. Error in file <%s> at line # %d when reading data of Path # %d. Key==%s is invalid or in the wrong place",FileName,line_num,idxpath,key.CStr());
                if (idxdat==ndat)
                {
                    ndat   = -1;
                    idxdat = 0;
                    idxpath++;
                    if (idxpath==npath) reading_path = false;
                }
            }
            else if (reading_eprps)
            {
                if      (key=="ndat") ndat = atoi(str_val.CStr());
                else if (ndat<0) throw new Fatal("InpFile::Read: Reading elements properties. Error in file <%s> at line # %d: key 'ndat' must come before data. '%s' is in the wrong place",FileName,line_num,key.CStr());
                else if (key=="elemtag") { elemtag = atoi(str_val.CStr());                       idxdat++; }
                else if (key=="prob")    { elemprps.Set (key.CStr(), FEM::PROB(str_val.CStr())); idxdat++; }
                else if (key=="geom")    { elemprps.Set (key.CStr(), FEM::GEOM(str_val.CStr())); idxdat++; }
                else if (key=="psa")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="pse")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="fra")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="d2d")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="d3d")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="rho")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="geosta")  { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="pospw")   { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="K0")      { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="surf")    { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="water")   { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else throw new Fatal("InpFile::Read: Reading elements properties. Error in file <%s> at line # %d when reading data of Properties # %d. Key==%s is invalid or in the wrong place",FileName,line_num,idxeprps,key.CStr());
                if (idxdat==ndat)
                {
                    Prps->Set (elemtag, elemprps);
                    ndat    = -1;
                    idxdat  = 0;
                    elemtag = -1;
                    elemprps.clear();
                    idxeprps++;
                    if (idxeprps==nelemprps) reading_eprps = false;
                }
            }
            else if (reading_stages && !reading_bcs)
            {
#ifdef USE_WXWIDGETS
                if (key=="nbcs") { nbcs = atoi(str_val.CStr());  Stages.Push(new GUI::WxDictTable);  reading_bcs=true; }
#else
                if (key=="nbcs") { nbcs = atoi(str_val.CStr());  Stages.Push(new Dict);  reading_bcs=true; }
#endif
                else throw new Fatal("InpFile::Read: Reading boundary conditions (stages). Error in file <%s> at line # %d: key '%s' is in the wrong place",FileName,line_num,key.CStr());
            }
            else if (reading_bcs)
            {
                if      (key=="ndat") ndat = atoi(str_val.CStr());
                else if (ndat<0) throw new Fatal("InpFile::Read: Reading boundary conditions (stages). Error in file <%s> at line # %d: key 'ndat' must come after 'nbcs' and before data. '%s' is in the wrong place",FileName,line_num,key.CStr());
                else if (key=="tag") { bcstag = atoi(str_val.CStr());              idxdat++; }
                else if (key=="ux")  { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else if (key=="uy")  { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else if (key=="uz")  { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else if (key=="fx")  { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else if (key=="fy")  { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else if (key=="fz")  { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else if (key=="qn")  { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else throw new Fatal("InpFile::Read: Reading boundary conditions (stages). Error in file <%s> at line # %d when reading data of Stage # %d. Key==%s is invalid or in the wrong place",FileName,line_num,idxstage,key.CStr());
                if (idxdat==ndat)
                {
                    Stages[idxstage]->Set (bcstag, bcs);
                    ndat   = -1;
                    idxdat = 0;
                    bcstag = 0;
                    bcs.clear ();
                    idxbcs++;
                    if (idxbcs==nbcs)
                    {
                        reading_bcs = false;
                        nbcs        = -1;
                        idxbcs      = 0;
                        idxstage++;
                        if (idxstage==nstages) reading_stages = false;
                    }
                }
            }
            else
            {
                if      (key=="matid")      matid     = atoi(str_val.CStr());       //  1
                else if (key=="flwid")      flwid     = atoi(str_val.CStr());       //  2
                else if (key=="ninc")       ninc      = atoi(str_val.CStr());       //  3
                else if (key=="cdrift")     cdrift    = (bool)atoi(str_val.CStr()); //  4
                else if (key=="stol")       stol      = val;                        //  5
                else if (key=="ssout")      ssout     = (bool)atoi(str_val.CStr()); //  6
                else if (key=="ctetg")      ctetg     = (bool)atoi(str_val.CStr()); //  7
                else if (key=="fem")        fem       = val;                        //  8
                else if (key=="dyn")        dyn       = (bool)atoi(str_val.CStr()); //  9
                else if (key=="hm")         hm        = (bool)atoi(str_val.CStr()); // 10
                else if (key=="tf")         tf        = val;                        // 11
                else if (key=="dt")         dt        = val;                        // 12
                else if (key=="dtout")      dtout     = val;                        // 13
                else if (key=="tsw")        tsw       = val;                        // 14
                else if (key=="ndiv")       ndiv      = atoi(str_val.CStr());       // 15
                else if (key=="nip")        nip       = atoi(str_val.CStr());       // 16
                else if (key=="o2")         o2        = (bool)atoi(str_val.CStr()); // 17
                else if (key=="ray")        ray       = (bool)atoi(str_val.CStr()); // 18
                else if (key=="am")         am        = val;                        // 19
                else if (key=="ak")         ak        = val;                        // 20
                else if (key=="rk")         rk        = (bool)atoi(str_val.CStr()); // 21
                else if (key=="rkscheme")   rkscheme  = str_val;                    // 22
                else if (key=="rkstol")     rkstol    = val;                        // 23
                else if (key=="refdat")     refdat    = str_val;                    // 24
                else if (key=="refsim")     refsim    = str_val;                    // 25
                else if (key=="refana")     refana    = str_val;                    // 26
                else if (key=="idxvert1")   idxvert1  = atoi(str_val.CStr());       // 27
                else if (key=="idxvert2")   idxvert2  = atoi(str_val.CStr());       // 28
                else if (key=="idxvert3")   idxvert3  = atoi(str_val.CStr());       // 29
                else if (key=="optdbl1")  { optdbl1   = val;   hasoptdbl1=true; }   // 30
                else if (key=="optdbl2")  { optdbl2   = val;   hasoptdbl2=true; }   // 31
                else if (key=="optdbl3")  { optdbl3   = val;   hasoptdbl3=true; }   // 32
                else if (key=="nldt_nsml")  nldt_nsml = atoi(str_val.CStr());       // 33
                else if (key=="nldt_nn")    nldt_nn   = atoi(str_val.CStr());       // 34
                else if (key=="nldt_n")     nldt_n    = atoi(str_val.CStr());       // 35
                else if (key=="nldt_ll")    nldt_ll   = val;                        // 36
                else if (key=="nldt_sch")   nldt_sch  = atoi(str_val.CStr());       // 37
                else if (key=="nldt_m")     nldt_m    = val;                        // 38
                else if (key=="maxit")      maxit     = atoi(str_val.CStr());       // 39
                else if (key=="tolr")       tolr      = val;                        // 40
                else if (key=="fnkey")      fnkey     = str_val;                    // 41
                else if (key=="pcam0")    { pcam0     = val;     haspcam0 = true; } // 42
                else if (key=="scheme")     scheme    = str_val;                    // 43
                else if (key=="vtufile")    vtufile   = (int)val;                   // 44
                else if (key=="suscheme")   suscheme  = str_val;                    // 45
                else if (key=="sustol")     sustol    = val;                        // 46
                else if (key=="surkscheme") surkscheme= str_val;                    // 47
                else if (key=="dcmaxit")    dcmaxit   = (int)val;                   // 48
                else if (key=="dcftol")     dcftol    = val;                        // 49
                else if (key=="npath")    { npath     = (int)val;  reading_path   = true; }
                else if (key=="nelemprps"){ nelemprps = (int)val;  reading_eprps  = true; }
                else if (key=="nstages")  { nstages   = (int)val;  reading_stages = true; }
                else if (key=="matids" || key=="xmatids")
                {
                    String left, right, str_id, str_tag;
                    line.Split (left, right, "=");
                    std::istringstream subiss(right);
                    while (subiss >> str_id >> str_tag)
                    {
                        int id  = atoi(str_id.CStr());
                        int tag = atoi(str_tag.CStr());
                        if (id<0 || tag>=0) throw new Fatal("InpFile::Read: Error in file <%s> @ line # %d with Key==%s. Material ids must be zero or positive and element tags must be negative. Ex.: matids = 0 -1  1 -2  2 -3",FileName,line_num,key.CStr());
                        if (key=="xmatids") XMaId2Tag->Set (str_id.CStr(), tag);
                        else                MatId2Tag->Set (str_id.CStr(), tag);
                    }
                }
                else if (key=="outnods")
                {
                    String left, right, nod;
                    line.Split (left, right, "=");
                    std::istringstream subiss(right);
                    while (subiss >> nod) OutNods->Push (atoi(nod.CStr()));
                }
                else throw new Fatal("InpFile::Read: Error in file <%s> @ line # %d: Key==%s in invalid (as general data)",FileName,line_num,key.CStr());
            }
        }
        line_num++;
    }
    if ((size_t)idxpath!=Path->Keys.Size()) throw new Fatal("InpFile::Read: Error in file <%s>: Not all Path data could be read for npath==%zd",FileName,npath);

    // filename key
    if (fnkey=="")
    {
        String buf(FileName);
        buf.GetFNKey (fnkey);
    }

#ifdef USE_WXWIDGETS
    Sync (/*Dat2Ctrl*/true);
    GPath      -> ReBuild ();
    GPrps      -> ReBuild ();
    GOutNods   -> ReBuild ();
    GMatId2Tag -> ReBuild ();
#endif
}

inline void InpFile::SetPrmsInis (MatFile const & Mat, bool ForceGTY)
{
    for (size_t i=0; i<MatId2Tag->Keys.Size(); ++i)
    {
        int            id   = atoi(MatId2Tag->Keys[i].CStr());
        int            tag  = (*MatId2Tag)(MatId2Tag->Keys[i]);
        SDPair const & prms = (*Mat.ID2Prms)(id);
        SDPair const & inis = (*Mat.ID2Inis)(id);
        Prms->Set (tag, prms);
        Inis->Set (tag, inis);
        if (ForceGTY)
        {
            if (Prps->Keys.Size()==0) throw new Fatal("InpFile::SetPrmsInis: Dictionary of properties (Prps) is empty => ForceGTY (geometry type) cannot be applied");
            SDPair const & prps = (*Prps)(tag);
            Array<String> gtypes("d3d", "d2d", "psa", "pse", "fra");
            for (size_t i=0; i<gtypes.Size(); ++i)
            {
                if (prps.HasKey(gtypes[i])) (*Prms)(tag).Set (gtypes[i].CStr(), 1.0);
            }
        }
        if (haspcam0) (*Inis)(tag).Set ("sx sy sz", -pcam0, -pcam0, -pcam0);
    }
    for (size_t i=0; i<XMaId2Tag->Keys.Size(); ++i)
    {
        int            id   = atoi(XMaId2Tag->Keys[i].CStr());
        int            tag  = (*XMaId2Tag)(XMaId2Tag->Keys[i]);
        SDPair         prms = (*Mat.ID2Prms)(id);
        SDPair const & inis = (*Mat.ID2Inis)(id);
        double nam = prms("name");
        prms.Del("name");
        prms.Set("xname", nam);
        Prms->Set (tag, prms);
        Inis->Set (tag, inis);
    }
}

inline void InpFile::GetIncs (int PathKey, double Div, Vec_t & dsig, Vec_t & deps, Array<bool> & PrescDeps, double & dpw, double & dSw, bool & PrescDpw, bool & PrescDSw) const
{
    // number of increments
    SDPair const & path = (*Path)(Path->Keys[PathKey]);

    // set prescribed increments
    set_to_zero (dsig);
    set_to_zero (deps);
    PrescDeps = false,false,false, false,false,false;
    dpw       = 0.0;
    dSw       = 0.0;
    PrescDpw  = false;
    PrescDSw  = false;
    if (path.HasKey("dsx")) { dsig(0) = path("dsx")/Div;                       }
    if (path.HasKey("dsy")) { dsig(1) = path("dsy")/Div;                       }
    if (path.HasKey("dsz")) { dsig(2) = path("dsz")/Div;                       }
    if (path.HasKey("dex")) { deps(0) = path("dex")/Div;                       }
    if (path.HasKey("dey")) { deps(1) = path("dey")/Div;  PrescDeps[1] = true; }
    if (path.HasKey("dez")) { deps(2) = path("dez")/Div;  PrescDeps[2] = true; }
    if (path.HasKey("dpw")) { dpw     = path("dpw")/Div;  PrescDpw     = true; }
    if (path.HasKey("dSw")) { dSw     = path("dSw")/Div;  PrescDSw     = true; }
}

inline void InpFile::SetSolver (FEM::Solver & Sol) const
{
    if (ssout    ) Sol.SSOut  = ssout;
    if (ctetg    ) Sol.CteTg  = ctetg;
    if (hm       ) Sol.DampTy = FEM::Solver::HMCoup_t;
    if (ray      ) Sol.DampTy = FEM::Solver::Rayleigh_t;
    if (am    >=0) Sol.DampAm = am;
    if (ak    >=0) Sol.DampAk = ak;
    if (rk       )
    {
        Sol.DScheme = FEM::Solver::RK_t;
        if (rkscheme!="") Sol.RKScheme = rkscheme;
        if (rkstol   >=0) Sol.RKSTOL   = rkstol;
    }
    if (maxit  >0 ) Sol.MaxIt = maxit;
    if (tolr   >=0) Sol.TolR  = tolr;
    if (scheme!="") Sol.SetScheme (scheme.CStr());
}

inline void InpFile::SetSUp (Model const * Mdl, Model::StressUpdate::pDbgFun pFun, void * UserData) const
{
    Mdl->SUp.CDrift = cdrift;
    Mdl->SUp.DbgFun = pFun;
    Mdl->SUp.DbgDat = UserData;
    if (suscheme  !="") Mdl->SUp.SetScheme (suscheme);
    if (surkscheme!="") Mdl->SUp.RKScheme = surkscheme;
    if (sustol      >0) Mdl->SUp.STOL     = sustol;
}

std::ostream & operator<< (std::ostream & os, InpFile const & IF)
{
    if (IF.matid       >=0) os << "matid     = " << IF.matid     << "\n"; //   1
    if (IF.flwid       >=0) os << "flwid     = " << IF.flwid     << "\n"; //   2
    if (IF.ninc        >=0) os << "ninc      = " << IF.ninc      << "\n"; //   3
    if (IF.cdrift         ) os << "cdrift    = " << IF.cdrift    << "\n"; //   4
    if (IF.stol        >=0) os << "stol      = " << IF.stol      << "\n"; //   5
    if (IF.ssout          ) os << "ssout     = " << IF.ssout     << "\n"; //   6
    if (IF.ctetg          ) os << "ctetg     = " << IF.ctetg     << "\n"; //   7
    if (IF.fem            ) os << "fem       = " << IF.fem       << "\n"; //   8
    if (IF.dyn            ) os << "dyn       = " << IF.dyn       << "\n"; //   9
    if (IF.hm             ) os << "hm        = " << IF.hm        << "\n"; //  10
    if (IF.tf          >=0) os << "tf        = " << IF.tf        << "\n"; //  11
    if (IF.dt          >=0) os << "dt        = " << IF.dt        << "\n"; //  12
    if (IF.dtout       >=0) os << "dtout     = " << IF.dtout     << "\n"; //  13
    if (IF.tsw         >=0) os << "tsw       = " << IF.tsw       << "\n"; //  14
    if (IF.ndiv        >=0) os << "ndiv      = " << IF.ndiv      << "\n"; //  15
    if (IF.nip         >=0) os << "nip       = " << IF.nip       << "\n"; //  16
    if (IF.o2             ) os << "o2        = " << IF.o2        << "\n"; //  17
    if (IF.ray            ) os << "ray       = " << IF.ray       << "\n"; //  18
    if (IF.am          >=0) os << "am        = " << IF.am        << "\n"; //  19
    if (IF.ak          >=0) os << "ak        = " << IF.ak        << "\n"; //  20
    if (IF.rk             ) os << "rk        = " << IF.rk        << "\n"; //  21
    if (IF.rkscheme   !="") os << "rkscheme  = " << IF.rkscheme  << "\n"; //  22
    if (IF.rkstol      >=0) os << "rkstol    = " << IF.rkstol    << "\n"; //  23
    if (IF.refdat     !="") os << "refdat    = " << IF.refdat    << "\n"; //  24
    if (IF.refsim     !="") os << "refsim    = " << IF.refsim    << "\n"; //  25
    if (IF.refana     !="") os << "refana    = " << IF.refana    << "\n"; //  26
    if (IF.idxvert1    >=0) os << "idxvert1  = " << IF.idxvert1  << "\n"; //  27
    if (IF.idxvert2    >=0) os << "idxvert2  = " << IF.idxvert2  << "\n"; //  28
    if (IF.idxvert3    >=0) os << "idxvert3  = " << IF.idxvert3  << "\n"; //  29
    if (IF.hasoptdbl1     ) os << "optdbl1   = " << IF.optdbl1   << "\n"; //  30
    if (IF.hasoptdbl2     ) os << "optdbl2   = " << IF.optdbl2   << "\n"; //  31
    if (IF.hasoptdbl3     ) os << "optdbl3   = " << IF.optdbl3   << "\n"; //  32
    if (IF.nldt_nsml   >0 ) os << "nldt_nsml = " << IF.nldt_nsml << "\n"; //  33
    if (IF.nldt_nn     >0 ) os << "nldt_nn   = " << IF.nldt_nn   << "\n"; //  34
    if (IF.nldt_n      >0 ) os << "nldt_n    = " << IF.nldt_n    << "\n"; //  35
    if (IF.nldt_ll     >0 ) os << "nldt_ll   = " << IF.nldt_ll   << "\n"; //  36
    if (IF.nldt_sch    >0 ) os << "nldt_sch  = " << IF.nldt_sch  << "\n"; //  37
    if (IF.nldt_m      >0 ) os << "nldt_m    = " << IF.nldt_m    << "\n"; //  38
    if (IF.maxit       >0 ) os << "maxit     = " << IF.maxit     << "\n"; //  39
    if (IF.tolr        >=0) os << "tolr      = " << IF.tolr      << "\n"; //  40
    if (IF.fnkey.size()>0 ) os << "fnkey     = " << IF.fnkey     << "\n"; //  41
    if (IF.haspcam0       ) os << "pcam0     = " << IF.pcam0     << "\n"; //  42
    if (IF.scheme     !="") os << "scheme    = " << IF.scheme    << "\n"; //  43
    if (IF.vtufile        ) os << "vtufile   = " << IF.vtufile   << "\n"; //  44
    if (IF.suscheme   !="") os << "suscheme  = " << IF.suscheme  << "\n"; //  45
    if (IF.sustol      >=0) os << "sustol    = " << IF.sustol    << "\n"; //  46
    if (IF.surkscheme !="") os << "surkscheme= " << IF.surkscheme<< "\n"; //  47
    if (IF.dcmaxit      >0) os << "dcmaxit   = " << IF.dcmaxit   << "\n"; //  48
    if (IF.dcftol       >0) os << "dcftol    = " << IF.dcftol    << "\n"; //  49
    os << "\nPath:\n"                        << (*IF.Path)      << "\n";
    os << "\nElements properties:\n"         << (*IF.Prps)      << "\n";
    os << "\nOutput nodes:\n"                << (*IF.OutNods)   << "\n";
    os << "\nMaterial IDs => Tags:\n"        << (*IF.MatId2Tag) << "\n";
    os << "\nE(x)tra Material IDs => Tags:\n"<< (*IF.XMaId2Tag) << "\n";
    os << "\nParameters:\n"                  << (*IF.Prms)      << "\n";
    os << "\nInitial values:\n"              << (*IF.Inis)      << "\n";
    os << "\nBoundary conditions (stages):\n";
    for (size_t i=0; i<IF.Stages.Size(); ++i)
    {
        os << "Stage # " << i << ":\n";
        os << (*IF.Stages[i]) << "\n";
    }
    return os;
}

#ifdef USE_WXWIDGETS

enum
{
    ID_INPFILE_LOAD = wxID_HIGHEST+1000,
    ID_INPFILE_SAVE,
};

BEGIN_EVENT_TABLE(InpFile, wxWindow)
    EVT_BUTTON (ID_INPFILE_LOAD, InpFile::OnLoad)
    EVT_BUTTON (ID_INPFILE_SAVE, InpFile::OnSave)
END_EVENT_TABLE()

inline InpFile::InpFile (wxFrame * Parent)
    : wxWindow (Parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE)
{
    // default values
    Defaults();

    // force validation of child controls
    SetExtraStyle (wxWS_EX_VALIDATE_RECURSIVELY);

    // tell wxAuiManager to manage this window
    Aui.SetManagedWindow (this);

    // control panel
    ADD_WXPANEL  (pnl, szt, szr, 1, 2);
    ADD_WXBUTTON (pnl, szr, ID_INPFILE_LOAD, c0, "Load");
    ADD_WXBUTTON (pnl, szr, ID_INPFILE_SAVE, c1, "Save");

    // main
    ADD_WXPANEL     (p_mai, sz_mai_t, sz_mai, 7, 2);
    ADD_WXNUMINPUT2 (p_mai, sz_mai, wxID_ANY, c_matid     , "matid    ", matid    ); //   1
    ADD_WXNUMINPUT2 (p_mai, sz_mai, wxID_ANY, c_flwid     , "flwid    ", flwid    ); //   2
    ADD_WXNUMINPUT2 (p_mai, sz_mai, wxID_ANY, c_ninc      , "ninc     ", ninc     ); //   3
    ADD_WXCHECKBOX2 (p_mai, sz_mai, wxID_ANY, c_fem       , "fem      ", fem      ); //   8
    ADD_WXCHECKBOX2 (p_mai, sz_mai, wxID_ANY, c_dyn       , "dyn      ", dyn      ); //   9
    ADD_WXCHECKBOX2 (p_mai, sz_mai, wxID_ANY, c_hm        , "hm       ", hm       ); //  10
    ADD_WXNUMINPUT2 (p_mai, sz_mai, wxID_ANY, c_pcam0     , "pcam0    ", pcam0    ); //  42
                                                                                            
    // local integration                                                                    
    ADD_WXPANEL     (p_loc, sz_loc_t, sz_loc, 2, 2);                                        
    ADD_WXCHECKBOX2 (p_loc, sz_loc, wxID_ANY, c_cdrift    , "cdrift   ", cdrift   ); //   4
    ADD_WXNUMINPUT2 (p_loc, sz_loc, wxID_ANY, c_stol      , "stol     ", stol     ); //   5
                                                                                            
    // fem solution                                                                         
    ADD_WXPANEL     (p_fem, sz_fem_t, sz_fem, 18, 2);                             
    ADD_WXCHECKBOX2 (p_fem, sz_fem, wxID_ANY, c_ssout     , "ssout    ", ssout    ); //   6
    ADD_WXCHECKBOX2 (p_fem, sz_fem, wxID_ANY, c_ctetg     , "ctetg    ", ctetg    ); //   7
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_tf        , "tf       ", tf       ); //  11
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_dt        , "dt       ", dt       ); //  12
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_dtout     , "dtout    ", dtout    ); //  13
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_tsw       , "tsw      ", tsw      ); //  14
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_ndiv      , "ndiv     ", ndiv     ); //  15
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_nip       , "nip      ", nip      ); //  16
    ADD_WXCHECKBOX2 (p_fem, sz_fem, wxID_ANY, c_o2        , "o2       ", o2       ); //  17
    ADD_WXCHECKBOX2 (p_fem, sz_fem, wxID_ANY, c_ray       , "ray      ", ray      ); //  18
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_am        , "am       ", am       ); //  19
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_ak        , "ak       ", ak       ); //  20
    ADD_WXCHECKBOX2 (p_fem, sz_fem, wxID_ANY, c_rk        , "rk       ", rk       ); //  21
    ADD_WXTEXTCTRL2 (p_fem, sz_fem, wxID_ANY, c_rkscheme  , "rkscheme ", rkscheme ); //  22
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_rkstol    , "rkstol   ", rkstol   ); //  23
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_maxit     , "maxit    ", maxit    ); //  39
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_tolr      , "tolr     ", tolr     ); //  40
    ADD_WXTEXTCTRL2 (p_fem, sz_fem, wxID_ANY, c_scheme    , "scheme   ", scheme   ); //  43
    ADD_WXCHECKBOX2 (p_fem, sz_fem, wxID_ANY, c_vtufile   , "vtufile  ", vtufile  ); //  44
                                                                                            
    // reference files                                                                      
    ADD_WXPANEL     (p_rfi, sz_rfi_t, sz_rfi, 3, 2);                                        
    ADD_WXTEXTCTRL2 (p_rfi, sz_rfi, wxID_ANY, c_refdat    , "refdat   ", refdat   ); //  24
    ADD_WXTEXTCTRL2 (p_rfi, sz_rfi, wxID_ANY, c_refsim    , "refsim   ", refsim   ); //  25
    ADD_WXTEXTCTRL2 (p_rfi, sz_rfi, wxID_ANY, c_refana    , "refana   ", refana   ); //  26
    c_refdat->SetMinSize (wxSize(200,20));
    c_refsim->SetMinSize (wxSize(200,20));
    c_refana->SetMinSize (wxSize(200,20));
                                                                                            
    // nonlinear steps                                                                      
    ADD_WXPANEL     (p_nls, sz_nls_t, sz_nls, 6, 2);                                        
    ADD_WXNUMINPUT2 (p_nls, sz_nls, wxID_ANY, c_nldt_nsml , "nldt_nsml", nldt_nsml); //  33
    ADD_WXNUMINPUT2 (p_nls, sz_nls, wxID_ANY, c_nldt_nn   , "nldt_nn  ", nldt_nn  ); //  34
    ADD_WXNUMINPUT2 (p_nls, sz_nls, wxID_ANY, c_nldt_n    , "nldt_n   ", nldt_n   ); //  35
    ADD_WXNUMINPUT2 (p_nls, sz_nls, wxID_ANY, c_nldt_ll   , "nldt_ll  ", nldt_ll  ); //  36
    ADD_WXNUMINPUT2 (p_nls, sz_nls, wxID_ANY, c_nldt_sch  , "nldt_sch ", nldt_sch ); //  37
    ADD_WXNUMINPUT2 (p_nls, sz_nls, wxID_ANY, c_nldt_m    , "nldt_m   ", nldt_m   ); //  38
                                                                                            
    // others                                                                               
    ADD_WXPANEL     (p_oth, sz_oth_t, sz_oth, 6, 2);                                        
    ADD_WXNUMINPUT2 (p_oth, sz_oth, wxID_ANY, c_idxvert1  , "idxvert1 ", idxvert1 ); //  27
    ADD_WXNUMINPUT2 (p_oth, sz_oth, wxID_ANY, c_idxvert2  , "idxvert2 ", idxvert2 ); //  28
    ADD_WXNUMINPUT2 (p_oth, sz_oth, wxID_ANY, c_idxvert3  , "idxvert3 ", idxvert3 ); //  29
    ADD_WXNUMINPUT2 (p_oth, sz_oth, wxID_ANY, c_optdbl1   , "optdbl1  ", optdbl1  ); //  30
    ADD_WXNUMINPUT2 (p_oth, sz_oth, wxID_ANY, c_optdbl2   , "optdbl2  ", optdbl2  ); //  31
    ADD_WXNUMINPUT2 (p_oth, sz_oth, wxID_ANY, c_optdbl3   , "optdbl3  ", optdbl3  ); //  32

    // path
    Path  = new GUI::WxDictTable;         Path->Transposed = false;
    GPath = new GUI::WxDict (this, Path); GPath->FitCol    = true;

    // additional data
    Prps       = new GUI::WxDictTable;
    OutNods    = new GUI::WxArrayIntTable;
    MatId2Tag  = new GUI::WxSIPairTable;
    //XMaId2Tag  = new GUI::WxSIPairTable;
    GPrps      = new GUI::WxDict     (this, Prps);
    GOutNods   = new GUI::WxArrayInt (this, OutNods);
    GMatId2Tag = new GUI::WxSIPair   (this, MatId2Tag);
    Prms       = new Dict;
    Inis       = new Dict;

    // notebook
    ADD_WXNOTEBOOK (this, nbk0);
    ADD_WXNOTEBOOK (this, nbk1);
    ADD_WXNOTEBOOK (this, nbk2);
    nbk0->AddPage  (p_mai,      "Main",                 false);
    nbk0->AddPage  (p_loc,      "Local Integration",    false);
    nbk0->AddPage  (p_oth,      "Others",               false);
    nbk2->AddPage  (p_fem,      "FEM Solution",         false);
    nbk2->AddPage  (p_nls,      "Nonlinear Steps",      false);
    nbk0->AddPage  (p_rfi,      "Reference Files",      false);
    nbk1->AddPage  (GPath,      "Path",                 false);
    nbk2->AddPage  (GPrps,      "Elements Properties",  false);
    nbk2->AddPage  (GOutNods,   "Output Nodes",         false);
    nbk2->AddPage  (GMatId2Tag, "Material IDs => Tags", false);

    // commit all changes to wxAuiManager
    Aui.AddPane (pnl,  wxAuiPaneInfo().Name("cpnl").Caption("cpnl").Top().MinSize(wxSize(100,50)).DestroyOnClose(false).CaptionVisible(false) .CloseButton(false));
    Aui.AddPane (nbk0, wxAuiPaneInfo().Name("nbk0").Caption("General Input Data").Centre().Position(0).DestroyOnClose(false).CaptionVisible(true).CloseButton(false));
    Aui.AddPane (nbk1, wxAuiPaneInfo().Name("nbk1").Caption("Stress/Strain Path").Centre().Position(1).DestroyOnClose(false).CaptionVisible(true).CloseButton(false));
    Aui.AddPane (nbk2, wxAuiPaneInfo().Name("nbk2").Caption("FEM Input Data")    .Centre().Position(2).DestroyOnClose(false).CaptionVisible(true).CloseButton(false));
    Aui.Update  ();
}

inline void InpFile::OnLoad (wxCommandEvent & Event)
{
    wxFileDialog fd(this, "Load input (.inp) file", LstDir, "", "*.inp");
    if (fd.ShowModal()==wxID_OK)
    {
        LstDir = fd.GetDirectory ();
        try { Read (fd.GetPath().ToStdString().c_str()); }
        catch (Fatal * e) { WxError(e->Msg().CStr()); }
    }
}

inline void InpFile::OnSave (wxCommandEvent & Event)
{
    Sync ();
    wxFileDialog fd(this, "Save input (.inp) file", LstDir, "", "*.inp", wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
    if (fd.ShowModal()==wxID_OK)
    {
        std::cout << "fem = " << fem << std::endl;
        std::fstream of(fd.GetPath().ToStdString().c_str(), std::ios::out);
        of << "matid     = " << matid     << std::endl; //  1
        of << "flwid     = " << flwid     << std::endl; //  2
        of << "ninc      = " << ninc      << std::endl; //  3
        of << "cdrift    = " << cdrift    << std::endl; //  4
        of << "stol      = " << stol      << std::endl; //  5
        of << "ssout     = " << ssout     << std::endl; //  6
        of << "ctetg     = " << ctetg     << std::endl; //  7
        of << "fem       = " << fem       << std::endl; //  8
        of << "dyn       = " << dyn       << std::endl; //  9
        of << "hm        = " << hm        << std::endl; // 10
        of << "tf        = " << tf        << std::endl; // 11
        of << "dt        = " << dt        << std::endl; // 12
        of << "dtout     = " << dtout     << std::endl; // 13
        of << "tsw       = " << tsw       << std::endl; // 14
        of << "ndiv      = " << ndiv      << std::endl; // 15
        of << "nip       = " << nip       << std::endl; // 16
        of << "o2        = " << o2        << std::endl; // 17
        of << "ray       = " << ray       << std::endl; // 18
        of << "am        = " << am        << std::endl; // 19
        of << "ak        = " << ak        << std::endl; // 20
        of << "rk        = " << rk        << std::endl; // 21
        of << "rkscheme  = " << rkscheme  << std::endl; // 22
        of << "rkstol    = " << rkstol    << std::endl; // 23
        of << "refdat    = " << refdat    << std::endl; // 24
        of << "refsim    = " << refsim    << std::endl; // 25
        of << "refana    = " << refana    << std::endl; // 26
        of << "idxvert1  = " << idxvert1  << std::endl; // 27
        of << "idxvert2  = " << idxvert2  << std::endl; // 28
        of << "idxvert3  = " << idxvert3  << std::endl; // 29
        of << "optdbl1   = " << optdbl1   << std::endl; // 30
        of << "optdbl2   = " << optdbl2   << std::endl; // 31
        of << "optdbl3   = " << optdbl3   << std::endl; // 32
        of << "nldt_nsml = " << nldt_nsml << std::endl; // 33
        of << "nldt_nn   = " << nldt_nn   << std::endl; // 34
        of << "nldt_n    = " << nldt_n    << std::endl; // 35
        of << "nldt_ll   = " << nldt_ll   << std::endl; // 36
        of << "nldt_sch  = " << nldt_sch  << std::endl; // 37
        of << "nldt_m    = " << nldt_m    << std::endl; // 38
        of << "maxit     = " << maxit     << std::endl; // 39
        of << "tolr      = " << tolr      << std::endl; // 40
        of << "npath     = " << Path->Keys.Size() << std::endl;

        String buf;
        for (size_t i=0; i<Path->Keys.Size(); ++i)
        {
            SDPair const & path = (*Path)(Path->Keys[i]);
            of << std::endl;
            of << "ndat = " << path.Keys.Size() << std::endl;
            for (size_t j=0; j<path.Keys.Size(); ++j)
            {
                buf.Printf("  %-5s = %g\n", path.Keys[j].CStr(), path(path.Keys[j]));
                of << buf;
            }
        }
        of.close();
    }
}

#endif


/////////////////////////////////////////////////////////////////////////////////////////// Macros ///////


#define INIT_MAT_INP(argc, argv, inpfn, matfn, verbose, forcegty,   MAT, INP) \
    if (argc>1) inpfn   =      argv[1];                                       \
    if (argc>2) matfn   =      argv[2];                                       \
    if (argc>3) verbose = atoi(argv[3]);                                      \
    MatFile MAT;                                                              \
    InpFile INP;                                                              \
    MAT.Read        (matfn.CStr());                                           \
    INP.Read        (inpfn.CStr());                                           \
    INP.SetPrmsInis (MAT, forcegty);                                          \
    if (verbose)                                                              \
    {                                                                         \
        printf("\n%s--- Materials data <%s> ----------------------------------------------------%s\n",TERM_CLR1,matfn.CStr(),TERM_RST); \
        std::cout << MAT << std::endl;                                        \
        printf("\n%s--- Input data <%s> --------------------------------------------------------%s\n",TERM_CLR1,inpfn.CStr(),TERM_RST); \
        std::cout << INP << std::endl;                                        \
    }

#define INIT_MAT_INP_(argc, argv,   MAT, INP) \
    String inpfn("input.inp");                \
    String matfn("materials.mat");            \
    bool   verbose  = true;                   \
    bool   forcegty = false;                  \
    INIT_MAT_INP(argc, argv, inpfn, matfn, verbose, forcegty,   MAT, INP);


#endif // MECHSYS_INPFILE_H
