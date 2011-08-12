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

#ifndef MECHSYS_FEM_NODE_H
#define MECHSYS_FEM_NODE_H

// Std Lib
#include <sstream> // for istringstream, ostringstream
#include <map>

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/linalg/sparse_triplet.h>

namespace FEM
{

/** Boundary conditions functions. */
class BCFuncs
{
public:
    virtual ~BCFuncs () {}
    virtual double u  (double Time) { return 0; } ///< Value of U component at Time
    virtual double v  (double Time) { return 0; } ///< Value of V component at Time
    virtual double a  (double Time) { return 0; } ///< Value of A component at Time
    virtual double fm (double Time) { return 0; } ///< Value of f component multiplier at Time
};

class Node
{
public:
    // Constructor
    Node (Mesh::Vertex const & TheVert) : Vert(TheVert), NShares(0), _has_incsup(false) {}

    // Methods
    void           AddDOF  (char const * StrU, char const * StrF);                      ///< Add DOF if it doesn't exist. Ex.: StrU="ux uy", StrF="fx fy"
    size_t         NDOF    ()                      const { return _eq.Size();         } ///< Get the number of DOFs
    double         U       (char   const * UKey  ) const { return _U(UKey);           } ///< Get U value
    double         F       (char   const * FKey  ) const { return _F(FKey);           } ///< Get F value
    double       & U       (String const & UKey  )       { return _U(UKey);           } ///< Get/Set U value
    double       & F       (String const & FKey  )       { return _F(FKey);           } ///< Get/Set F value
    double       & U       (char   const * UKey  )       { return _U(UKey);           } ///< Get/Set U value
    double       & F       (char   const * FKey  )       { return _F(FKey);           } ///< Get/Set F value
    double       & V       (size_t         IdxDOF)       { return _V(UKey(IdxDOF));   } ///< Get/Set V value given IdxDOF
    double       & U       (size_t         IdxDOF)       { return _U(UKey(IdxDOF));   } ///< Get/Set U value given IdxDOF
    double       & F       (size_t         IdxDOF)       { return _F(FKey(IdxDOF));   } ///< Get/Set F value given IdxDOF
    double         UOrZero (String const & UKey  ) const { return _U.ValOrZero(UKey); } ///< U value or zero if U key does not exist
    double         FOrZero (String const & FKey  ) const { return _F.ValOrZero(FKey); } ///< F value or zero if F key does not exist
    void           TrySetU (String const & UKey, double Val) { if (_U.HasKey(UKey)) _U(UKey)=Val; } ///< Try to set U if key exists
    void           TrySetF (String const & FKey, double Val) { if (_F.HasKey(FKey)) _F(FKey)=Val; } ///< Try to set F if key exists
    int            Eq      (char   const * UKey  ) const { return _eq[_U2IDOF(UKey)]; } ///< Get Equation number
    int            Eq      (String const & UKey  ) const { return _eq[_U2IDOF(UKey)]; } ///< Get Equation number
    int            Eq      (size_t         IdxDOF) const { return _eq[IdxDOF];        } ///< Get Equation number given index of DOF
    int          & Eq      (size_t         IdxDOF)       { return _eq[IdxDOF];        } ///< Get/Set Equation number given index of DOF
    String const & UKey    (size_t         IdxDOF) const { return _U.Keys[IdxDOF];    } ///< Get U key
    String const & FKey    (size_t         IdxDOF) const { return _F.Keys[IdxDOF];    } ///< Get F key
    bool           pU      (size_t         IdxDOF) const { return _pU[IdxDOF];        } ///< Prescribed U ?
    void           SetUF   (Vec_t const & UVec, Vec_t const & FVec);                    ///< Set internal U and F values given global vectors
    void           GetUF   (Vec_t       & UVec, Vec_t       & FVec);                    ///< Get U and F values and store into UVec and FVec vectors
    void           Clear   ();                                                          ///< Clear all values but keep structure of DOFs
    void           ClearU  () { _U.SetValues (0.0); }                                   ///< Clear U values only

    // Methods to set/access prescribed U values
    void           SetPU     (String const & UKey, double Val, BCFuncs * BCF=NULL);                      ///< Set prescribed U
    void           SetMPU    (size_t IdxPU, BCFuncs * BCF=NULL){ _MPU[IdxPU] = BCF; }                    ///< Set multiplier to PU
    size_t         NPU       ()                          const { return _PU.Keys.Size(); }               ///< Get number of prescribed U
    int            EqPU      (size_t IdxPU)              const { return _eq[_U2IDOF(_PU.Keys[IdxPU])]; } ///< Get Eq number corresponding to prescribed U
    void           DelPUs    ()                                { _PU.clear(); _MPU.Resize(0); }          ///< Delete PU structure
    String const & PUKey     (size_t IdxPU)              const { return _PU.Keys[IdxPU]; }               ///< Get PU key
    double         PU        (size_t IdxPU, double Time) const;                                          ///< Get prescribed U value given index to PU
    double         PV        (size_t IdxPU, double Time) const;                                          ///< Get prescribed V=dUdt value given index to PU
    double         PA        (size_t IdxPU, double Time) const;                                          ///< Get prescribed A=d2Udt2 value given index to PU
    double         PUIdxDOF  (size_t IdxDOF, double Time) const { return PU(_PU.IdxKey(UKey(IdxDOF)),Time); } ///< Get prescribed U given index of DOF
    double         PVIdxDOF  (size_t IdxDOF, double Time) const { return PV(_PU.IdxKey(UKey(IdxDOF)),Time); } ///< Get prescribed V given index of DOF
    double         PAIdxDOF  (size_t IdxDOF, double Time) const { return PA(_PU.IdxKey(UKey(IdxDOF)),Time); } ///< Get prescribed A given index of DOF

    // Methods to set/access prescribed F values
    void           AddToPF   (String const & FKey, double Val, BCFuncs * BCF=NULL);                      ///< Add value to prescribed F
    void           SetMPF    (size_t IdxPF, BCFuncs * BCF=NULL){ _MPF[IdxPF] = BCF; }                    ///< Set multiplier to PF
    size_t         NPF       ()                          const { return _PF.Keys.Size(); }               ///< Get number of prescribed F
    int            EqPF      (size_t IdxPF)              const { return _eq[_F2IDOF(_PF.Keys[IdxPF])]; } ///< Get Eq number corresponding to prescribed F
    void           DelPFs    ()                                { _PF.clear(); _MPF.Resize(0); }          ///< Delete PF structure
    String const & PFKey     (size_t IdxPF)              const { return _PF.Keys[IdxPF]; }               ///< Get PF key
    double         PF        (size_t IdxPF, double Time) const;                                          ///< Get prescribed F value given index to PF
    double         PFOrZero  (size_t IdxDOF, double Time) const { int idxpf=_PF.IdxKey(FKey(IdxDOF)); return (idxpf<0 ? 0.0 : PF(idxpf,Time)); } ///< Get prescribed F given index of DOF
    void           AccumPF   ();                                                                         ///< Accumulate prescribed F (to calculate Reactions later). Copy from _PF into _aPF.
    void           Reactions (std::map<String,double> & R) const;                                        ///< Calculate reactions

    // Pins
    void SetPin    (Array<Node*> const & ConnectedNodes) { _pin_nodes=ConnectedNodes; } ///< Set Pin
    void SetLagPin (int & EqLag, Sparse::Triplet<double,int> & A) const;                ///< Set A matrix with the equations for the Lagrangian multipliers corresponding to pins. Will increment EqLag
    void ClrRPin   (int & EqLag, Vec_t & R)                       const;                ///< Clear R components corresponding pins. Will inclined EqLag

    // Inclined supports: 2D
    void SetIncSup    (double Alpha);                                            ///< Set inclined support
    void DelIncSup    ()             { _has_incsup=false; }                      ///< Delete inclined support
    bool HasIncSup    () const       { return _has_incsup; }                     ///< Has inclined support ?
    void SetLagIncSup (int & EqLag, Sparse::Triplet<double,int> & A) const;      ///< Set A matrix with the equations for the Lagrangian multipliers corresponding to inclined supports. Will increment EqLag
    void ClrRIncSup   (int & EqLag, Vec_t & R)                       const;      ///< Clear R components corresponding to inclined supports. Will inclined EqLag

    // Data
    Mesh::Vertex const & Vert;    ///< Geometric information: ID, Tag, coordinates
    long                 NShares; ///< Number of active elements sharing this node

private:
    // Data at every node
    SDPair     _V;      ///< V = dUdt
    SDPair     _U;      ///< U values mapped by keys: "ux", "uy", etc.
    SDPair     _F;      ///< F values mapped by keys: "fx", "fy", etc.
    SIPair     _U2IDOF; ///< Maps: UKey to index of DOF
    SIPair     _F2IDOF; ///< Maps: FKey to index of DOF
    Array<int> _eq;     ///< Equation numbers (of each DOF)
    Array<bool> _pU;    ///< Prescribed U ?

    // Data at nodes with prescribed values
    SDPair          _PU;  ///< Prescribed U
    SDPair          _PF;  ///< Prescribed F
    SDPair          _aPF; ///< Accumulated prescribed F (not erased in DelPFs)
    Array<BCFuncs*> _MPU; ///< Multipliers of PU
    Array<BCFuncs*> _MPF; ///< Multipliers of PF

    // Data for pins
    Array<Node*> _pin_nodes; ///< Other nodes connected to this pin

    // Data at nodes with inclined supports
    double _incsup_alpha; ///< Inclined support alpha's
    bool   _has_incsup;   ///< Has inclined support ?

friend std::ostream & operator<< (std::ostream & os, Node const & N);
};


///////////////////////////////////////////////////////////////////////////////////// Node: Implementation /////


inline void Node::AddDOF (char const * StrU, char const * StrF)
{
    std::istringstream u_iss(StrU);
    std::istringstream f_iss(StrF);
    String u_key, f_key;
    while ((u_iss>>u_key) && (f_iss>>f_key))
    {
        if (!_U.HasKey(u_key)) // add only if not found
        {
            _V.Set      (u_key.CStr(), 0.0);
            _U.Set      (u_key.CStr(), 0.0);
            _F.Set      (f_key.CStr(), 0.0);
            _U2IDOF.Set (u_key.CStr(), _eq.Size());
            _F2IDOF.Set (f_key.CStr(), _eq.Size());
            _eq.Push    (-1);
            _pU.Push    (false);
        }
        else
        {
            if (!_F.HasKey(f_key)) throw new Fatal("FEM::Node: __internal_error__ Node has UKey=%s but not FKey=%s",u_key.CStr(),f_key.CStr());
        }
    }
}

inline void Node::Clear ()
{
    _U .SetValues (0.0);
    _F .SetValues (0.0);
    _eq.SetValues (-1);
    _pU.SetValues (false);
    DelPUs ();
    DelPFs ();
    _has_incsup = false;
}

inline void Node::SetUF (Vec_t const & UVec, Vec_t const & FVec)
{
    for (size_t idof=0; idof<NDOF(); ++idof)
    {
        int eq = Eq(idof);
        _U(UKey(idof)) = UVec(eq);
        _F(FKey(idof)) = FVec(eq);
    }
}

inline void Node::GetUF (Vec_t & UVec, Vec_t & FVec)
{
    for (size_t idof=0; idof<NDOF(); ++idof)
    {
        int eq = Eq(idof);
        UVec(eq) = _U(UKey(idof));
        FVec(eq) = _F(FKey(idof));
    }
}

inline void Node::SetPU (String const & UKey, double Val, BCFuncs * BCF)
{
    size_t idx_pu = _PU.ReSet (UKey.CStr(), Val);
    if (idx_pu==_MPU.Size()) _MPU.Push (BCF);
    else                     _MPU[idx_pu] = BCF;
    _has_incsup = false;
    _pU[_U2IDOF(UKey)] = true;
}

inline void Node::AddToPF (String const & FKey, double Val, BCFuncs * BCF)
{
    size_t idx_pf = _PF.AddVal (FKey.CStr(), Val);
    if (idx_pf==_MPF.Size()) _MPF.Push (BCF);
    else                     _MPF[idx_pf] = BCF;
}

inline double Node::PU (size_t IdxPU, double Time) const
{
    if (_MPU[IdxPU]==NULL) return _PU(_PU.Keys[IdxPU]);
    else                   return _MPU[IdxPU]->u(Time);
}

inline double Node::PV (size_t IdxPU, double Time) const
{
    if (_MPU[IdxPU]==NULL) return 0.0;
    else                   return _MPU[IdxPU]->v(Time);
}

inline double Node::PA (size_t IdxPU, double Time) const
{
    if (_MPU[IdxPU]==NULL) return 0.0;
    else                   return _MPU[IdxPU]->a(Time);
}

inline double Node::PF (size_t IdxPF, double Time) const
{
    if (_MPF[IdxPF]==NULL) return _PF(_PF.Keys[IdxPF]);
    else                   return _PF(_PF.Keys[IdxPF]) * _MPF[IdxPF]->fm(Time);
}

inline void Node::AccumPF ()
{
    for (StrDbl_t::iterator it=_PF.begin(); it!=_PF.end(); ++it)
        _aPF.AddVal (it->first.CStr(), it->second);
}

inline void Node::Reactions (std::map<String,double> & R) const
{
    for (size_t i=0; i<_PU.Keys.Size(); ++i)
    {
        String const & ukey = _PU.Keys[i];
        int            idof = _U2IDOF(ukey);
        String const & fkey = _F.Keys[idof];
        double         reac = _F(fkey);
        if (_aPF.HasKey(fkey)) reac -= _aPF(fkey);
        R[ukey] = reac;
    }
}

inline void Node::SetLagPin (int & EqLag, Sparse::Triplet<double,int> & A) const
{
    for (size_t i=0; i<_pin_nodes.Size(); ++i)
    {
        for (size_t j=0; j<_U.Keys.Size(); ++j)
        {
            if (_U.Keys[j]=="ux" || _U.Keys[j]=="uy" || _U.Keys[j]=="uz")
            {
                int eq0 =                Eq(_U.Keys[j]);
                int eq1 = _pin_nodes[i]->Eq(_U.Keys[j]);
                A.PushEntry (eq0,EqLag,1.0);   A.PushEntry (eq1,EqLag,-1.0);
                A.PushEntry (EqLag,eq0,1.0);   A.PushEntry (EqLag,eq1,-1.0);
                EqLag++; // each pin adds one equation per DOF
            }
        }
    }
}

inline void Node::ClrRPin (int & EqLag, Vec_t & R) const
{
    for (size_t i=0; i<_pin_nodes.Size(); ++i)
    {
        for (size_t j=0; j<_U.Keys.Size(); ++j)
        {
            if (_U.Keys[j]=="ux" || _U.Keys[j]=="uy" || _U.Keys[j]=="uz")
            {
                int eq0 =                Eq(_U.Keys[j]);
                int eq1 = _pin_nodes[i]->Eq(_U.Keys[j]);
                R(eq0) = 0.0;
                R(eq1) = 0.0;
                EqLag++; // each pin adds one equation per DOF
            }
        }
    }
}

inline void Node::SetIncSup (double Alpha)
{
    if (NPU()>0) return; // skip if this node has prescribed U already
    _incsup_alpha = Alpha*Util::PI/180.0;
    _has_incsup   = true; 
}

inline void Node::SetLagIncSup (int & EqLag, Sparse::Triplet<double,int> & A) const
{
    double s   = sin(_incsup_alpha);
    double c   = cos(_incsup_alpha);
    int    eq0 = Eq("ux");
    int    eq1 = Eq("uy");
    A.PushEntry (EqLag, eq0,  s);
    A.PushEntry (EqLag, eq1, -c);
    A.PushEntry (eq0, EqLag,  s);
    A.PushEntry (eq1, EqLag, -c);
    EqLag++; // each inclined support adds one equation
}

inline void Node::ClrRIncSup (int & EqLag, Vec_t & R) const
{
    int eq0 = Eq("ux");
    int eq1 = Eq("uy");
    R(eq0) = 0.0;
    R(eq1) = 0.0;
    EqLag++; // each inclined support adds one equation
}

std::ostream & operator<< (std::ostream & os, Node const & N)
{
    os << Util::_4 << N.Vert.ID << " " << Util::_4 << N.Vert.Tag << " ";
    for (size_t i=0; i<N.NDOF(); ++i)
    {
        os << N.UKey(i) << " ";
        os << N.FKey(i) << " ";
        if (i!=N.NDOF()-1) os << " ";
    }
    os << " NShares=" << (N.NShares>0?TERM_GREEN:TERM_RED) << N.NShares << TERM_RST;
    os << " EQ=[";
    for (size_t i=0; i<N.NDOF(); ++i)
    {
        os << N.Eq(i);
        if (i!=N.NDOF()-1) os << ",";
        else               os << "]";
    }
    os << " (";
    for (int j=0; j<3; ++j)
    {
        os << Util::_6_3 << N.Vert.C[j];
        if (j==2) os << ")";
        else      os << ", ";
    }
    if (N._pin_nodes.Size()>0)
    {
        os << TERM_GREEN <<" PIN:[";
        for (size_t i=0; i<N._pin_nodes.Size(); ++i) os << N._pin_nodes[i]->Vert.ID << (i==N._pin_nodes.Size()-1?"]":",");
        os << TERM_RST;
    }
    if (N._has_incsup)         os << " INCSUP";
    os << TERM_CLR4 << Util::_reset << " PU=[";
    for (size_t i=0; i<N.NPU(); ++i)
    {
        os << N._PU.Keys[i] << "=" << N._PU(N._PU.Keys[i]);
        if (N._MPU[i]!=NULL) os << "*";
        if (i!=N.NPU()-1) os << " ";
    }
    os << "]" << TERM_CLR5 << " PF=[";
    for (size_t i=0; i<N.NPF(); ++i)
    {
        os << N._PF.Keys[i] << "=" << N._PF(N._PF.Keys[i]);
        if (N._MPF[i]!=NULL) os << "*";
        if (i!=N.NPF()-1) os << " ";
    }
    os << "]" << TERM_RST;
    return os;
}

}; //namespace FEM

#endif // MECHSYS_FEM_NODE_H
