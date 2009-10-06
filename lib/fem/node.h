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
#include "mesh/mesh.h"
#include "util/maps.h"
#include "util/fatal.h"
#include "util/numstreams.h"
#include "linalg/matvec.h"

namespace FEM
{

class Node
{
public:
    // Constructor
    Node (Mesh::Vertex const & TheVert) : Vert(TheVert) {}

    // Methods
    void   AddDOF (char const * StrU, char const * StrF);
    void   SetBCs (SDPair const & BCs);
    void   ClrBCs ();
    size_t nDOF   () const { return UMap.Keys.Size(); }
    double Fa     (size_t iDOF, double Time) const { return Time*DF[iDOF]; } ///< F(applied)

    // Data
    Mesh::Vertex const & Vert;  ///< Geometric information: ID, Tag, coordinates
    SIPair               UMap;  ///< U keys ("ux", "uy", ...) to index (0, 1, ...) map
    SIPair               FMap;  ///< F keys ("fx", "fy", ...) to index (0, 1, ...) map
    Array<double>        DU;    ///< Delta U (if prescribed U) to be applied in one stage
    Array<double>        DF;    ///< Delta F (if NOT pU) to be applied in one stage
    Array<bool>          pU;    ///< Prescribed U (essential) ?
    Array<double>        U;     ///< Current U value. Size = num DOF
    Array<double>        F;     ///< Current F value. Size = num DOF
    Array<long>          EQ;    ///< Equation numbers (of each DOF)
};


///////////////////////////////////////////////////////////////////////////////////// operator << //////////////


std::ostream & operator<< (std::ostream & os, Node const & N)
{
    os << Util::_6 << N.Vert.ID << " ";
    for (size_t i=0; i<N.nDOF(); ++i)
    {
        os << "[";
        os << "d" << N.UMap.Keys[i] << "=" << Util::_6_3 << N.DU[i] << ", ";
        os << "d" << N.FMap.Keys[i] << "=" << Util::_6_3 << N.DF[i] << ", ";
        os << "p" << N.UMap.Keys[i] << "=" << (N.pU[i]?"[1;31mT[0m":"F") << "]";
        if (i!=N.nDOF()-1) os << " ";
    }
    os << " EQ=[";
    for (size_t i=0; i<N.EQ.Size(); ++i)
    {
        os << N.EQ[i];
        if (i!=N.EQ.Size()-1) os << ",";
        else                  os << "]";
    }
    /*
    os << " U=[";
    for (size_t i=0; i<N.U.Size(); ++i)
    {
        os << N.U[i];
        if (i!=N.U.Size()-1) os << ",";
        else                 os << "]";
    }
    os << " F=[";
    for (size_t i=0; i<N.F.Size(); ++i)
    {
        os << N.F[i];
        if (i!=N.F.Size()-1) os << ",";
        else                 os << "]";
    }
    */
    os << " (";
    for (int j=0; j<3; ++j)
    {
        os << Util::_6_3 << N.Vert.C[j];
        if (j==2) os << ")";
        else      os << ", ";
    }
    os << Util::_4 << N.Vert.Tag;
    return os;
}


///////////////////////////////////////////////////////////////////////////////////// Node: Implementation /////


inline void Node::AddDOF (char const * StrU, char const * StrF)
{
    std::istringstream u_iss(StrU);
    std::istringstream f_iss(StrF);
    String u_key, f_key;
    while ((u_iss>>u_key) && (f_iss>>f_key))
    {
        if (!UMap.HasKey(u_key)) // add only if not found
        {
            int idx = static_cast<int>(nDOF());
            UMap.Set   (u_key.CStr(), idx);
            FMap.Set   (f_key.CStr(), idx);
            DU  .Push  (0.0);
            DF  .Push  (0.0);
            pU  .Push  (false);
            U   .Push  (0.0);
            F   .Push  (0.0);
            EQ  .Push  (-1);
        }
    }
}

inline void Node::SetBCs (SDPair const & BCs)
{
    for (size_t i=0; i<BCs.Keys.Size(); ++i)
    {
        String key = BCs.Keys[i];
        if (UMap.HasKey(key)) // essential
        {
            int idx = UMap(key);
            DU[idx] = BCs(key); // set
            //DF[idx] = 0.0;
            pU[idx] = true;
        }
        else if (FMap.HasKey(key)) // natural
        {
            int idx  = FMap(key);
            DU[idx]  = 0.0;
            DF[idx] += BCs(key); // accumulate
            pU[idx]  = false;
        }
        else
        {
            std::ostringstream oss;
            oss << (*this);
            throw new Fatal("Node::SetBCs: Node does not have key=%s.  %s", key.CStr(), oss.str().c_str());
        }
    }
}

inline void Node::ClrBCs ()
{
    for (size_t i=0; i<nDOF(); ++i)
    {
        DU[i] = 0.0;
        DF[i] = 0.0;
        pU[i] = false;
    }
}

}; //namespace FEM

#endif // MECHSYS_FEM_NODE_H
