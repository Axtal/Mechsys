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
 * You should have received A copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef DEM3D_QUATERNION_H
#define DEM3D_QUATERNION_H

#include <math.h>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include <linalg/matvec.h>

typedef blitz::TinyVector<double,4> Quaternion_t;

void NormalizeRotation(double Theta, Vec3_t const & Axis, Quaternion_t & C)
{
    Vec3_t A = Axis/norm(Axis);
    C(0)     = cos(Theta/2.0);
    C(1)     = A(0)*sin(Theta/2.0);
    C(2)     = A(1)*sin(Theta/2.0);
    C(3)     = A(2)*sin(Theta/2.0);
}

void Conjugate(Quaternion_t const & A, Quaternion_t & C)
{
    C(0) =  A(0);
    C(1) = -A(1);
    C(2) = -A(2);
    C(3) = -A(3);
}

void GetVector(Quaternion_t const & A, Vec3_t & C)
{
    C(0) = A(1);
    C(1) = A(2);
    C(2) = A(3);
}

void SetQuaternion(double Scalar, Vec3_t const & A, Quaternion_t & C)
{
    C(0) = Scalar;
    C(1) = A(0);
    C(2) = A(1);
    C(3) = A(2);
}

void QuaternionProduct(Quaternion_t const & A, Quaternion_t const & B, Quaternion_t & C)
{
    Vec3_t t1,t2;
    GetVector(A,t1);
    GetVector(B,t2);
    double scalar = A(0)*B(0)-dot(t1,t2);
    Vec3_t vect = A(0)*t2 + B(0)*t1 + cross(t1,t2);
    SetQuaternion(scalar,vect,C);
}

void Rotation(Vec3_t const & A, Quaternion_t const & B, Vec3_t & C)
{
    Quaternion_t t1,t2,t3;
    SetQuaternion     (0.0,A,t1);
    QuaternionProduct (B,t1,t2);
    Conjugate         (B,t3);
    QuaternionProduct (t2,t3,t1);
    GetVector         (t1,C);
}

#endif // DEM3D_QUATERNION_H
