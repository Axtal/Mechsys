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

#ifndef __FACE3D_H
#define __FACE3D_H

#include <vector>
#include <map>
#include "Edge3d.h"

class CFace3d {
    protected:
        int m_Ns;
        vector<CEdge3d> m_sides;
        //CEdge3d m_sides[4];
    public:
        CFace3d() {};
        CFace3d(Vec3 *,int N);
        CFace3d(int N) {m_Ns=N;};
        void translation(Vec3 & a);
        void rotation(Quaternion &,Vec3 &);
        friend CFace3d rotation_face(const CFace3d &,Quaternion & Q,Vec3 & v);
        friend CEdge3d distance_point_Face3d(Vec3 &,CFace3d &);
        friend class CGraph;
        friend class CParticle;
        friend class CEdge3d;
};


#endif
