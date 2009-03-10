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

#ifndef MECHSYS_LBM_LATTICE_H
#define MECHSYS_LBM_LATTICE_H

// MechSys
#include "fem/node.h"

typedef blitz::TinyVec<double,3> Vec3_t;

namespace LBM
{

class Lattice
{
public:

protected:

}; // class Lattice

const Vec3_t LOCAL_VELOC2D[ 9] = {};
const Vec3_t LOCAL_VELOC3D[27] = {};

class Lattice2D : public Lattice
{
public:
	Lattice2D() { _w = WEIGHTS2D; }
private:
}; // class Lattice2D


class Lattice3D : public Lattice
{
public:
	Lattice3D() { _w = WEIGHTS3D; }
private:
}; // class Lattice3D


}; // namespace LBM

#endif // MECHSYS_LBM_LATTICE_H
