/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso                                *
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

/* MPoints2D - Copyright (C) 2007 Dorival de Moraes Pedroso */

#ifndef MPM_MPOINTS2D_H
#define MPM_MPOINTS2D_H

// STL
#include <cmath>  // for sqrt, pow, etc.
#include <cfloat> // for DBL_EPSILON

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/models/model.h>
#include <mechsys/models/equilibstate.h>

// Local
#include <mechsys/mpm/defs.h>
#include <mechsys/mpm/grid2d.h>
#include <mechsys/mpm/tensors.h>

/** MPM: Material Point Method. */
namespace MPM {

class MPoints2D
{
public:
    static const int NDIM = 2;
    typedef std::map<int,Model*> Models_t; ///< Map tag to model pointer

    /** Alternative constructor. nPCell=number of points per cell */
    MPoints2D (size_t nPcell, Grid2D * G,
               ptIsPointInGeom pIsPointInGeom,
               ptDensity       pDensity,
               ptModelData     pModelData,
               ptIniVelocity   pIniVelocity,
               ptHasTraction   pHasTraction,
               ptHasAppDisp    pHasAppDisp,
               ptAppDisp       pAppDisp=NULL);

    // Destructor
    ~MPoints2D ();

    // Methods
    size_t nPoints      () const { return _p.Size();     } ///< Return the number of material points
    size_t nPointsOnBry () const { return _ipbry.Size(); } ///< Return the number of material points on boundary
    void   SetGrid      (Grid2D * G) { _g = G; }           ///< Set access (read/write) to the grid
    void   ReInit       () { _initialize_points(); }       ///< Reinitialize points (during initial refinement)

    // Update method
    bool TimeUpdate (bool FwEuler, double Dt, ptB pB, ptLdM pLdM, double & t, double & DsE); ///< Advance one step in time, considering new b-body forces. Returns DsE, the increment of strain Energy.

    // Set Methods
    void SetCPDI (bool UseCPDI=true) { _cpdi = UseCPDI; } ///< Use CPDI ?
    void SetMPM  (bool UseMPM=true)  { _mpm  = UseMPM;  } ///< Use (original) MPM instead of the generalized MPM (GMPM) ?
    void SetUSF  (bool DoUSF =true)  { _usf  = DoUSF;   } ///< Do update stress first (USF) instead of updating stress last (USL) ?

    // Access methods
    Array<ATensor2> const & Fs   ()         const { return _F;             } ///< Returns access to all deformation gradients
    Array<Vector3D> const & Ps   ()         const { return _p;             } ///< Returns access to all mat points positions
    Array<Vector3D> const & v    ()         const { return _v;             } ///< Returns access to all mat points velocities
    Array<Vector3D> const & u    ()         const { return _u;             } ///< Returns access to all mat points displacements
    Vector3D        const & P    (size_t i) const { return _p[i];          } ///< Returns access to a material point (read)
    Vector3D        const & B    (size_t i) const { return _p[_ipbry[i]];  } ///< Returns access to a material point which is on boundary (read)
    double          const & m    (size_t i) const { return _m[i];          } ///< Mass of the point/particle
    Vector3D        const & f    (size_t i) const { return _f[i];          } ///< Force on the point/particle
    Vector3D        const & t    (size_t i) const { return _t[i];          } ///< Tractions on point/particle
    Vector3D        const & fu   (size_t i) const { return _fu[i];         } ///< Applied displacements on point/particle
    Vector3D        const & u    (size_t i) const { return _u[i];          } ///< Displacement of the point/particle
    Vector3D        const & v    (size_t i) const { return _v[i];          } ///< Velocity of the point/particle
    Vec_t           const & s    (size_t i) const { return _sta[i]->Sig;   } ///< Stress on the point/particle
    Vec_t           const & e    (size_t i) const { return _sta[i]->Eps;   } ///< Strain on the point/particle
    bool                    HasT (size_t i) const { return _has_t[i];      } ///< Is mat point on boundary with applied tractions ?
    bool                    HasU (size_t i) const { return _has_fu[i];     } ///< Is mat point on boundary with applied displacements ?
    int                     Clr  (size_t i) const { return _clr[i];        } ///< Color index
    Array<EquilibState*> const & Sta () const { return _sta; }

    // Access corners
    mutable Array<Vector3D> C; ///< 4 corners
    void CalcC (size_t p) const ///< calculate corners (must be called before accessing C)
    {
        if (NRVEC==2)
        {
            C[0] = _p[p] - _sg[p].R[0] - _sg[p].R[1];
            C[1] = _p[p] + _sg[p].R[0] - _sg[p].R[1];
            C[2] = _p[p] + _sg[p].R[0] + _sg[p].R[1];
            C[3] = _p[p] - _sg[p].R[0] + _sg[p].R[1];
        }
        else for (size_t i=0; i<4; ++i) C[i] = _p[p] + _sg[p].R[i];
    }

    // Limit methods
    void MinMaxS    (int iComp, double & MinS, double & MaxS) const; ///< limits: stresses
    void MinMaxE    (int iComp, double & MinE, double & MaxE) const; ///< limits: strains
    void MinMaxCamp (double & Minp, double & Maxp)            const; ///< limits: mean stress
    void MinMaxCamq (double & Minq, double & Maxq)            const; ///< limits: deviatoric stress
    void MinMaxVelN (double & MinNv, double & MaxNv)          const; ///< limits: norm of velocity

    // Constants
    double MINMASS; ///< Minimum mass of the point/particle

    // Data
    Models_t Mdls; ///< All models

private:
    // Data
    bool              _cpdi;             ///< Use CPDI
    bool              _mpm;              ///< Use original MPM shape functions
    bool              _usf;              ///< Update stress first ?
    size_t            _npcell;           ///< Number of points per cell (initial only)
    Grid2D          * _g;                ///< Access to the grid
    ptIsPointInGeom   _is_point_in_geom; ///< Is point in geometry ?
    ptDensity         _density;          ///< Density function
    ptModelData       _model_data;       ///< Get model data
    ptIniVelocity     _ini_velocity;     ///< Initial velocity
    ptHasTraction     _has_traction;     ///< Has traction BC ?
    ptHasAppDisp      _has_app_disp;     ///< Has applied displacements BC ?
    ptAppDisp         _app_disp;         ///< Applied displacements

    // Geometry
    Array<Vector3D> _p;     ///< Mat points current positions (x,y coords)
    Array<int>      _ipbry; ///< Index of the mat points on boundary

    // Data (size == _np)
    Array<Vector3D>        _p0;  ///< Mat points initial positions (x,y coords)
    Array<double>          _V0;  ///< Initial volume
    Array<double>          _V;   ///< Volume
    Array<double>          _m;   ///< Masses
    Array<Vector3D>        _f;   ///< Force (x,y comps)
    Array<Vector3D>        _u;   ///< Displacements (x,y comps)
    Array<Vector3D>        _v;   ///< velocities (x,y comps)
    Array<ATensor2>        _F;   ///< Def grad
    Array<EquilibState*>   _sta; ///< Constitutive models for each material point
    Array<ShapeAndGrads2D> _sg;  ///< Shape and gradient functions
    Array<Model*>          _mdl; ///< Pointer to models
    Array<int>             _clr; ///< Color index

    // Boundary conditions (size == np)
    Array<bool>     _has_t;  ///< Is mat point on boundary with applied traction ?
    Array<bool>     _has_fu; ///< Is mat point on boundary with applied displacements ?
    Array<Vector3D> _t;      ///< Mat points tractions
    Array<Vector3D> _fu;     ///< Mat points force due to applied displacements

    // Scratch pad
    mutable STensor2 _sp0; ///< initial stress (before update). to calculate strain energy
    mutable Vector3D _vn;  ///< node velocity. to calculate _Lp
    mutable ATensor2 _Lp;  ///< velocity gradient at point p
    mutable ATensor2 _DF;  ///< increment of F

    // Update functions
    void _initialize_points ();                                                                                ///< Initialize/generate mat points
    void _calc_veloc_grad   (size_t p);                                                                        ///< Calculate _Lp
    void _stress_update     (double Dt, double & DsE);                                                         ///< Update stresses; returns DsE (increment of strain Energy)
    bool _explicit_update   (double Dt, Vector3D const & B, double LdM, Array<Vector3D> * DqDt, double & DsE); ///< Explicit/Forward-Euler update

    // Interpolation functions
    double _mpm_snp  (int n, int p, int icomp)    const; ///< 1D Shape functions (MPM)
    double _mpm_gnp  (int n, int p, int icomp)    const; ///< 1D Derivative of shape functions (MPM)
    double _gmpm_snp (int n, int p, int icomp)    const; ///< 1D Shape functions (GMPM)
    double _gmpm_gnp (int n, int p, int icomp)    const; ///< 1D Derivative of shape functions (GMPM)
    void   _Nnp      (int n, int p, double & N)   const; ///< 2D shape functions (C0)
    double _S        (int n, Vector3D const & X)  const; ///< 2D shape (hat-function)
    void   _Snp      (int n, int p, double & S)   const; ///< 2D shape functions (C1)
    void   _Gnp      (int n, int p, Vector3D & G) const; ///< 2D Derivative of shape functions
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline MPoints2D::MPoints2D(size_t nPcell, Grid2D * G, ptIsPointInGeom pIsPointInGeom, ptDensity pDensity, ptModelData pModelData, ptIniVelocity pIniVelocity, ptHasTraction pHasTraction, ptHasAppDisp pHasAppDisp, ptAppDisp pAppDisp)
    : MINMASS           (DBL_EPSILON*10.0),
      _cpdi             (true),
      _mpm              (false),
      _usf              (true),
      _npcell           (nPcell),
      _g                (G),
      _is_point_in_geom (pIsPointInGeom),
      _density          (pDensity),
      _model_data       (pModelData),
      _ini_velocity     (pIniVelocity),
      _has_traction     (pHasTraction),
      _has_app_disp     (pHasAppDisp),
      _app_disp         (pAppDisp)
{
    C.Resize(4);
    _initialize_points();
}

inline MPoints2D::~MPoints2D ()
{
    for (size_t i=0; i<_sta.Size(); ++i) if (_sta[i]!=NULL) delete _sta[i];
    for (Models_t::iterator it=Mdls.begin(); it!=Mdls.end(); ++it) delete it->second;
}

inline bool MPoints2D::TimeUpdate (bool FwEuler, double Dt, ptB pB, ptLdM pM, double & t, double & DsE)
{
    // Check if access to grid was set
    if (_g==NULL) return false; // failed

    // Variables
    double   dsek;
    Vector3D bk;
    double   Mk;

    // Update
    if (FwEuler)
    {
        double nincs = 1;        // Number of sub-increments
        double dtk   = Dt/nincs; // Sub-divide time-step
        for (int k=0; k<nincs; ++k)
        {
            // Increments
            (*pB) (t+dtk, bk); // Body forces
            (*pM) (t+dtk, Mk); // Multiplier for external forces

            // Update points
            bool ok = _explicit_update (dtk, bk, Mk, &_g->dqdt(), dsek);

            // Update time and strain energy
            t   += dtk;
            DsE += dsek;

            // Check
            if (!ok) return false;
        }
        return true;
    }
    else return false; // not implemented yet
}

inline void MPoints2D::MinMaxS (int iComp, double & MinS, double & MaxS) const
{
    if (_sta.Size()<1) { MinS=0; MaxS=0; return; }
    MinS = _sta[0]->Sig(iComp);
    MaxS = _sta[0]->Sig(iComp);
    for (size_t i=1; i<_p.Size(); ++i)
    {
        if (_sta[i]->Sig(iComp)<MinS) MinS = _sta[i]->Sig(iComp);
        if (_sta[i]->Sig(iComp)>MaxS) MaxS = _sta[i]->Sig(iComp);
    }
}

inline void MPoints2D::MinMaxE (int iComp, double & MinE, double & MaxE) const
{
    if (_sta.Size()<1) { MinE=0; MaxE=0; return; }
    MinE = _sta[0]->Eps(iComp);
    MaxE = _sta[0]->Eps(iComp);
    for (size_t i=1; i<_p.Size(); ++i)
    {
        if (_sta[i]->Eps(iComp)<MinE) MinE = _sta[i]->Eps(iComp);
        if (_sta[i]->Eps(iComp)>MaxE) MaxE = _sta[i]->Eps(iComp);
    }
}

inline void MPoints2D::MinMaxCamp (double & Minp, double & Maxp) const
{
    if (_sta.Size()<1) { Minp=0; Maxp=0; return; }
    Minp = Calc_pcam(_sta[0]->Sig);
    Maxp = Minp;
    for (size_t i=1; i<_p.Size(); ++i)
    {
        double p = Calc_pcam(_sta[i]->Sig);
        if (p<Minp) Minp = p;
        if (p>Maxp) Maxp = p;
    }
}

inline void MPoints2D::MinMaxCamq (double & Minq, double & Maxq) const
{
    if (_sta.Size()<1) { Minq=0; Maxq=0; return; }
    Minq = Calc_qcam(_sta[0]->Sig);
    Maxq = Minq;
    for (size_t i=1; i<_p.Size(); ++i)
    {
        double q = Calc_qcam(_sta[i]->Sig);
        if (q<Minq) Minq = q;
        if (q>Maxq) Maxq = q;
    }
}

inline void MPoints2D::MinMaxVelN (double & MinNv, double & MaxNv) const
{
    if (_v.Size()<1) { MinNv=0; MaxNv=0; return; }
    MinNv = MPM::NormV(_v[0]);
    MaxNv = MinNv;
    for (size_t i=1; i<_v.Size(); ++i)
    {
        double nv = MPM::NormV(_v[i]);
        if (nv<MinNv) MinNv = nv;
        if (nv>MaxNv) MaxNv = nv;
    }
}


/* update functions */

inline void MPoints2D::_initialize_points ()
{
    // Check input
    if (!(_npcell==1 || _npcell==4 || _npcell==9 || _npcell==16)) throw new Fatal("MPoints2D::_initialize_points: Number of points per cell (NPcell=%d) must be 1, 4, 9 or 16.",_npcell);

    /*   x  sqrt(x)  2*sqrt(x)  log2(2*sqrt(x))  1/sqrt(x)  1-1/sqrt(x)
     *   1        1          2                1       1            0.0
     *   4        2          4                2       0.5          0.5
     *  16        4          8                3       0.25         0.75
     *
     *        1                    4                    9                    16
     *  1_ _______           1_ _______          1_ ________            1_ _______ 
     *    |       |        0.5_|   | _ |_3L/4      |  |  |  |        0.75>|_|_|_|_|<7L/8
     *  0_|   _   |_L/2      0_|___|___|         0_|--|--|--|_L/2       0_|_|_|_|_|
     *    |       |            |   | _ |_L/4       |--|--|--|             |_|_|_|_|
     * -1_|_______|         -1_|___|___|         1_|__|__|__|          -1_|_|_|_|_|<L/8    */

    // Gauss points
    double sqnpc   = sqrt(static_cast<double>(_npcell));
    double gv      = 1.0-1.0/sqnpc; // Gauss point value
    double ksi[ 4] = {-1,1,1,-1};   // Natural coordinates
    double eta[ 4] = {-1,-1,1,1};   // Natural coordinates
    double kgp[16] = {-gv, gv,gv,-gv,  -gv,-gv,       -gv/3.,-gv/3.,-gv/3.,-gv/3.,   gv/3.,gv/3.,gv/3.,gv/3.,   gv,gv       }; // ksi Gauss point
    double egp[16] = {-gv,-gv,gv, gv,  -gv/3.,gv/3.,  -gv,-gv/3.,gv/3.,gv,          -gv,-gv/3.,gv/3.,gv,       -gv/3.,gv/3. }; // eta Gauss point
    if (_npcell==9)
    {
        gv = (1.0-(-1.0))/3.0; // one third of cell length in nat cooords
        kgp[0]=-gv;   kgp[1]=0.0;   kgp[2]=gv;   kgp[3]=-gv;   kgp[4]=0.0;   kgp[5]=gv;   kgp[6]=-gv;   kgp[7]=0.0;   kgp[8]=gv;
        egp[0]=-gv;   egp[1]=-gv;   egp[2]=-gv;  egp[3]=0.0;   egp[4]=0.0;   egp[5]=0.0;  egp[6]=gv;    egp[7]=gv;    egp[8]=gv;
    }

    // Loop over cells
    _p     .Resize(0);
    _clr   .Resize(0);
    _t     .Resize(0);
    _has_t .Resize(0);
    _fu    .Resize(0);
    _has_fu.Resize(0);
    for (size_t c=0; c<_g->nCells(); ++c)
    {
        // For each cell c, generate _npcell mat points
        for (size_t j=0; j<_npcell; ++j)
        {
            // For each node j, sum the contribution of node c
            Vector3D p; p = 0.0;
            Array<Vector3D> nodes(4);
            for (size_t l=0; l<4; ++l)
            {
                nodes[l] = _g->X(_g->C(c)(l));
                double shape = (1.0+ksi[l]*kgp[j])*(1.0+eta[l]*egp[j])/4.0;
                p(0) += shape*nodes[l](0);
                p(1) += shape*nodes[l](1);
            }

            // Check if the mat point is inside the geometry
            int tag = -1;
            if ((*_is_point_in_geom)(p, tag))
            {
                // Coords
                _p  .Push (p);   // coordinates
                _clr.Push (tag); // clr index

                // Boundary conds - tractions
                _t    .Push(0.0);
                _has_t.Push((*_has_traction) (_p[_p.Size()-1], nodes, _g->L(), _npcell, _t[_p.Size()-1]));

                // Boundary conds - displacements
                _fu    .Push(0.0);
                _has_fu.Push((*_has_app_disp) (_p[_p.Size()-1], nodes, _g->L(), _npcell, _fu[_p.Size()-1]));
            }
        }
    }

    // Delete previous EquilibState
    for (size_t i=0; i<_sta.Size(); ++i) if (_sta[i]!=NULL) delete _sta[i];

    // Delete previous Models
    for (Models_t::iterator it=Mdls.begin(); it!=Mdls.end(); ++it) delete it->second;
    Mdls.clear();

    // Data for models
    String name, scheme;
    SDPair prms, inis;

    // Initialize arrays
    _p0 .Resize(_p.Size());
    _V0 .Resize(_p.Size());
    _V  .Resize(_p.Size());
    _m  .Resize(_p.Size());
    _f  .Resize(_p.Size());
    _u  .Resize(_p.Size());
    _v  .Resize(_p.Size());
    _F  .Resize(_p.Size());
    _sta.Resize(_p.Size());
    _sg .Resize(_p.Size());
    _mdl.Resize(_p.Size());
    double mt = 0.0; // total mass
    for (size_t p=0; p<_p.Size(); ++p)
    {
        double lx = 0.5*_g->L(0)/sqnpc; // half-length
        double ly = 0.5*_g->L(1)/sqnpc; // half-length
        _p0  [p] = _p[p];
        _V0  [p] = 4.0*lx*ly;
        _V   [p] = _V0[p];
        _m   [p] = _V[p] * (*_density) (_p[p]);
        _f   [p] = 0.0, 0.0, 0.0;
        _u   [p] = 0.0, 0.0, 0.0;
        _v   [p] = 0.0, 0.0, 0.0;
        (*_ini_velocity) (_p[p], _v[p]);
        _F   [p] = 1.0,1.0,1.0, 0.0,0.0,0.0, 0.0,0.0,0.0;
        if (NRVEC==2)
        {
            _sg[p].R0[0] = lx, 0., 0.;
            _sg[p].R0[1] = 0., ly, 0.;
        }
        else
        {
            _sg[p].R0[0] = -lx, -ly, 0.0;
            _sg[p].R0[1] =  lx, -ly, 0.0;
            _sg[p].R0[2] =  lx,  ly, 0.0;
            _sg[p].R0[3] = -lx,  ly, 0.0;
        }
        for (size_t k=0; k<NRVEC; ++k) _sg[p].R[k] = _sg[p].R0[k];
        mt += _m[p];

        // allocate models and set initial state
        int tag = -1;
        (*_model_data) (_p[p], name, tag, NDIM, prms, inis, scheme);
        Models_t::iterator it = Mdls.find(tag);
        if (it==Mdls.end())
        {
            Mdls[tag] = AllocModel (name, NDIM, prms, NULL);
            it = Mdls.find(tag);
            if (scheme!="") it->second->SUp.SetScheme (scheme);
        }
        _mdl[p] = it->second;
        _sta[p] = new EquilibState (NDIM);
        _mdl[p]->InitIvs (inis, _sta[p]);
    }
    std::cout << "MPoints2D::_initialize_points: total mass = " << mt << std::endl;
}

inline void MPoints2D::_calc_veloc_grad (size_t p)
{
    _Lp = 0.,0.,0., 0.,0.,0., 0.,0.,0.;
    for (n2i_it it=_sg[p].n2i.begin(); it!=_sg[p].n2i.end(); ++it)
    {
        int n = it->first;                             // grid node number
        int i = it->second;                            // index in N,S,G
        _vn = 0., 0., 0.;                              // node velocity
        if (_g->m(n)>MINMASS) _vn = _g->q(n)/_g->m(n); // node velocity
        if (_g->IsFixed(n,0)) _vn(0) = 0.0;            // Fix nodes: x-direction
        if (_g->IsFixed(n,1)) _vn(1) = 0.0;            // Fix nodes: y-direction
        DyadUp (_vn, _sg[p].G[i],  _Lp);               // Lp += vp dyad G
    }
}

inline void MPoints2D::_stress_update (double Dt, double & DsE)
{
    DsE = 0.0;
    for (size_t p=0; p<_p.Size(); ++p)
    {
        // Compute velocity gradient: _Lp
        _calc_veloc_grad (p);

        // Update F and stress state
        DsE += _mdl[p]->SUp.Update (Dt, _Lp, _F[p], _sta[p]) * _V[p];

        // Update volume
        _V[p] = Det(_F[p]) * _V0[p];

        // Update centre-to-corner vectors
        for (size_t i=0; i<NRVEC; ++i) Dot (_F[p], _sg[p].R0[i], _sg[p].R[i]);
    }
}

inline bool MPoints2D::_explicit_update (double Dt, Vector3D const & B, double LdM, Array<Vector3D> * DqDt, double & DsE)
{
    // 1) Discard previous grid
    _g->ClearState();

    // 2) Compute interpolation values
    if (_cpdi)
    {
        for (size_t p=0; p<_p.Size(); ++p)
        {
            CalcC(p); // calculate corners
            _sg[p].n2i.clear();
            int k = 0.0;
            for (size_t i=0; i<4; ++i) // for each corner i
            {
                int l   = static_cast<int>((C[i](0)-_g->xMin())/_g->L(0)); // left-grid-index
                int b   = static_cast<int>((C[i](1)-_g->yMin())/_g->L(1)); // bottom-grid-index
                int r   = l+1;                                             // right-grid-index
                int t   = b+1;                                             // top-grid-index
                int nlb = l+b*_g->nCol();
                int nrb = r+b*_g->nCol();
                int nlt = l+t*_g->nCol();
                int nrt = r+t*_g->nCol();
                _sg[p].n2i[nlb]=k;  _Nnp(nlb,p,_sg[p].N[k]);  _Snp(nlb,p,_sg[p].S[k]);  _Gnp(nlb,p,_sg[p].G[k]);  k++;
                _sg[p].n2i[nrb]=k;  _Nnp(nrb,p,_sg[p].N[k]);  _Snp(nrb,p,_sg[p].S[k]);  _Gnp(nrb,p,_sg[p].G[k]);  k++;
                _sg[p].n2i[nlt]=k;  _Nnp(nlt,p,_sg[p].N[k]);  _Snp(nlt,p,_sg[p].S[k]);  _Gnp(nlt,p,_sg[p].G[k]);  k++;
                _sg[p].n2i[nrt]=k;  _Nnp(nrt,p,_sg[p].N[k]);  _Snp(nrt,p,_sg[p].S[k]);  _Gnp(nrt,p,_sg[p].G[k]);  k++;
            }
        }
    }
    else
    {
        int nctb = (_mpm ? 2 : 4); // number of contributions
        for (size_t p=0; p<_p.Size(); ++p)
        {
            _sg[p].n2i.clear();
            int k   = 0;
            int lbn = _g->LBN(_p[p], _mpm); // Compute reference node
            for (int i=0; i<nctb; ++i)
            for (int j=0; j<nctb; ++j)
            {
                int n = lbn+i+j*_g->nCol(); // grid node number
                _sg[p].n2i[n] = k;          // set map: node => index in S,N,G
                _Nnp (n,p, _sg[p].N[k]);    // calc interpolation function (C0)
                _Snp (n,p, _sg[p].S[k]);    // calc interpolation function (C1)
                _Gnp (n,p, _sg[p].G[k]);    // calc deriv. interp. function
                k++;
            }
        }
    }

    // 3) Initialize grid state (mass and momentum)
    for (size_t p=0; p<_p.Size(); ++p)
    {
        for (n2i_it it=_sg[p].n2i.begin(); it!=_sg[p].n2i.end(); ++it)
        {
            int n = it->first;                   // grid node number
            int i = it->second;                  // index in N,S,G
            _g->m(n) += _sg[p].S[i]*_m[p];       // node mass
            _g->q(n) += _sg[p].S[i]*_m[p]*_v[p]; // node momentum
        }
    }

    // 4) (USF) Update strain and stress
    if (_usf) _stress_update (Dt, DsE);

    // 5) Compute internal and external forces
    for (size_t p=0; p<_p.Size(); ++p)
    {
        for (n2i_it it=_sg[p].n2i.begin(); it!=_sg[p].n2i.end(); ++it)
        {
            int n = it->first;                                   // grid node number
            int i = it->second;                                  // index in N,S,G
            ScDotUp (_V[p],_sta[p]->Sig,_sg[p].G[i], _g->fi(n)); // fi += V_p*Sig.G
            _g->fe(n) += (_sg[p].S[i]*_m[p])*B;                  // fe += mp*b*S
            if (_has_t[p]) _g->fe(n) += LdM*_sg[p].N[i]*_t[p];   // fe += int(S*t)dA
        }
    }

    // 6) Compute rate of momentum and momentum
    for (size_t n=0; n<_g->nNodes(); ++n)
    {
        (*DqDt)[n] = _g->fe(n) - _g->fi(n); // Nodes rate of momentum
        _g->q(n)  += (*DqDt)[n]*Dt;         // Update nodes momentum
        if (_g->IsFixed(n,0)) { (*DqDt)[n](0)=0.0;  _g->q(n)(0)=0.0; } // Fix nodes: x-direction
        if (_g->IsFixed(n,1)) { (*DqDt)[n](1)=0.0;  _g->q(n)(1)=0.0; } // Fix nodes: y-direction
    }

    // 7) Update material points
    for (size_t p=0; p<_p.Size(); ++p)
    {
        // Update position and velocity
        for (n2i_it it=_sg[p].n2i.begin(); it!=_sg[p].n2i.end(); ++it)
        {
            int n = it->first;  // grid node number
            int i = it->second; // index in N,S,G
            if (_g->m(n)>MINMASS)
            {
                _p[p] += (Dt*_sg[p].S[i]/_g->m(n))*_g->q(n);   // point position
                _v[p] += (Dt*_sg[p].S[i]/_g->m(n))*(*DqDt)[n]; // point velocity
                _f[p] +=     _sg[p].S[i]*_g->fe(n);            // point force
            }
        }

        // Update displacements
        _u[p] = _p[p] - _p0[p];

        // Check
        if (!_g->IsInside(_p[p], /*WithBorder*/true))
        {
            std::cout << "  [1;31mError:[0m Material point moved outside the valid region inside the grid (with borders)" << std::endl;
            return false;
        }
    }

    // 8) (USL) Update strain and stress
    if (!_usf) _stress_update (Dt, DsE);

    // end
    return true;
}


/* interpolation functions */

inline double MPoints2D::_mpm_snp (int n, int p, int icomp) const
{
    double d = _p[p](icomp)-_g->x(n)(icomp);
    double L = _g->L(icomp);
    return (d<=-L  ? 0.0     :
           (d>  L  ? 0.0     :
           (d<=0.0 ? 1.0+d/L :
                     1.0-d/L )));
}

inline double MPoints2D::_mpm_gnp (int n, int p, int icomp) const
{
    double d = _p[p](icomp)-_g->x(n)(icomp);
    double L = _g->L(icomp);
    return (d<=-L   ?  0.0   :
           (d>  L   ?  0.0   :
           (d<= 0.0 ?  1.0/L :
                      -1.0/L )));
}

inline double MPoints2D::_gmpm_snp (int n, int p, int icomp) const
{
    double d  = _p[p](icomp)-_g->x(n)(icomp);
    double L  = _g->L(icomp);
    double lp = fabs(_sg[p].R0[0](icomp));
    return (d<=-L-lp ? 0.0                        :
           (d>  L+lp ? 0.0                        :
           (d<=-L+lp ? pow(L+lp+d,2.0)/(4.0*L*lp) :
           (d<=  -lp ? 1.0+d/L                    :
           (d<=   lp ? 1.0-(d*d+lp*lp)/(2.0*L*lp) :
           (d<= L-lp ? 1.0-d/L                    :
                       pow(L+lp-d,2.0)/(4.0*L*lp) ))))));
}

inline double MPoints2D::_gmpm_gnp (int n, int p, int icomp) const
{
    double d  = _p[p](icomp)-_g->x(n)(icomp);
    double L  = _g->L(icomp);
    double lp = fabs(_sg[p].R0[0](icomp));
    return (d<=-L-lp ? 0.0                  :
           (d>  L+lp ? 0.0                  :
           (d<=-L+lp ? (L+d+lp)/(2.0*L*lp)  :
           (d<=  -lp ? 1.0/L                :
           (d<=   lp ? -d/(L*lp)            :
           (d<= L-lp ? -1.0/L               :
                       (-L+d-lp)/(2.0*L*lp) ))))));
}

inline void MPoints2D::_Nnp (int n, int p, double & N) const
{
    N = _mpm_snp(n,p,0)*_mpm_snp(n,p,1);
}

inline double MPoints2D::_S (int n, Vector3D const & X) const
{
    double dx = X(0)-_g->x(n)(0);
    double dy = X(1)-_g->x(n)(1);
    double sx = ( dx<-_g->L(0) || dx>_g->L(0) ? 0.0 : (dx<0.0 ? 1.0+dx/_g->L(0) : 1.0-dx/_g->L(0)) );
    double sy = ( dy<-_g->L(1) || dy>_g->L(1) ? 0.0 : (dy<0.0 ? 1.0+dy/_g->L(1) : 1.0-dy/_g->L(1)) );
    return sx*sy;
}

inline void MPoints2D::_Snp (int n, int p, double & S) const
{
    if (_cpdi)
    {
        // remember to calculate corners first !!!
        S = (_S(n,C[0])+_S(n,C[1])+_S(n,C[2])+_S(n,C[3]))/4.0;
    }
    else
    {
        if (_mpm) S = _mpm_snp (n,p,0)*_mpm_snp (n,p,1);
        else      S = _gmpm_snp(n,p,0)*_gmpm_snp(n,p,1);
    }
}

inline void MPoints2D::_Gnp (int n, int p, Vector3D & G) const 
{
    if (_cpdi)
    {
        // remember to calculate corners first !!!
        double a = _S(n,C[0]) - _S(n,C[2]);
        double b = _S(n,C[1]) - _S(n,C[3]);
        G = (a*(_sg[p].R[0](1)-_sg[p].R[1](1))+b*(_sg[p].R[0](1)+_sg[p].R[1](1)))/_V[p],
            (a*(_sg[p].R[1](0)-_sg[p].R[0](0))-b*(_sg[p].R[1](0)+_sg[p].R[0](0)))/_V[p],  0.0;
    }
    else
    {
        if (_mpm) G = _mpm_snp (n,p,1)*_mpm_gnp (n,p,0),  _mpm_snp (n,p,0)*_mpm_gnp (n,p,1),  0.0;
        else      G = _gmpm_snp(n,p,1)*_gmpm_gnp(n,p,0),  _gmpm_snp(n,p,0)*_gmpm_gnp(n,p,1),  0.0;
    }
}

}; // namespace MPM

#endif // MPM_MPOINTS2D_H
