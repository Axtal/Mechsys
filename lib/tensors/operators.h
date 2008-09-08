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

#ifndef MECHSYS_TENSORS_OPERATORS_H
#define MECHSYS_TENSORS_OPERATORS_H

// MechSys
#include "tensors/tensors.h"
#include "util/util.h"

using Util::SQ2;

/* Operators:
 * sc: scalar
 * T2: 2nd order tensors;  6 components (Mandel representation)
 * T4: 4th order tensors; 36 components (Mandel representation)
 *
 * + add:
 *    T2 = T2 + T2
 *    T4 = T4 + T4
 * - sub:
 *    T2 = T2 - T2
 *    T4 = T4 - T4
 * * dot:
 *    sc = T2 * T2;     T2 = T2 * T4;
 *    T2 = T4 * T2;     T4 = T4 * T4;
 * & dyadic:
 *    T4 = (T2 & T2);
 * 
 */

namespace Tensors
{

// 0) v = A * u
/// \f$ \{v\}\gets[A]\bullet\{u\} \quad\equiv\quad \TeSe{v}\gets\TeFo{A}:\TeSe{u} \f$
void Dot (Tensor2 const & A,
          Tensor1 const & u,
          Tensor1       & v);

// 1) sc = x * y
/// \f$ sc\gets\{x\}\bullet\{y\} \quad\equiv\quad sc\gets\TeSe{x}:\TeSe{y} \f$
double Dot (Tensor2 const & x,
            Tensor2 const & y);

// 2) y = x * A
/// \f$ \{y\}\gets\{x\}\bullet[A] \quad\equiv\quad \TeSe{y}\gets\TeSe{x}:\TeFo{A} \f$
void Dot (Tensor2 const & x,
          Tensor4 const & A,
          Tensor2       & y);

// 3) y = A * x
/// \f$ \{y\}\gets[A]\bullet\{x\} \quad\equiv\quad \TeSe{y}\gets\TeFo{A}:\TeSe{x} \f$
void Dot (Tensor4 const & A,
          Tensor2 const & x,
          Tensor2       & y);

// 4) C = A * B
/// \f$ [C]\gets[A]\bullet[B] \quad\equiv\quad \TeFo{C}\gets\TeFo{A}:\TeFo{B} \f$
void Dot (Tensor4 const & A,
          Tensor4 const & B,
          Tensor4       & C);

// 5) A = x dyadic y
/// \f$ [A]=\{x\}\otimes\{y\} \quad\equiv\quad \TeFo{A}\gets\TeSe{x}\otimes\TeSe{y} \f$
void Dyad (Tensor2 const & x,
           Tensor2 const & y,
           Tensor4       & A);

// 6) Z = a*X + b*Y
/// Add scaled tensors: \f$ [Z] \gets a[X]+b[Y] \quad\equiv\quad \TeFo{Z}\gets a\TeFo{X}+b\TeFo{Y} \f$
void AddScaled (double  const & a,
                Tensor4 const & X,
                double  const & b,
                Tensor4 const & Y,
                Tensor4       & Z);

// 7) A = a * (x dyadic y) + A
/// \f$ [A]=\alpha\{x\}\otimes\{y\}+[A] \quad\equiv\quad \TeFo{A}\gets\alpha\TeSe{x}\otimes\TeSe{y}+\TeFo{A} \f$
void Ger (double  const & a,
          Tensor2 const & x,
          Tensor2 const & y,
          Tensor4       & A);

// 8) D = a * (A*x) dyadic (y*B) + C
/// \f$ [D]=\alpha([A]\bullet\{x\})\otimes(\{y\}\bullet[B])+[C] \quad\equiv\quad \TeFo{D}\gets\alpha(\TeFo{A}:\TeSe{x})\otimes(\TeSe{y}:\TeFo{B})+\TeFo{C} \f$
void GerX (double  const & a,
           Tensor4 const & A,
           Tensor2 const & x,
           Tensor2 const & y,
           Tensor4 const & B,
           Tensor4 const & C, 
           Tensor4       & D);

// 9) B = a * B
/// \f$ [B]=\alpha([B]) \quad\equiv\quad \TeFo{B}\gets\alpha\TeFo{B} \f$
void Scale (double  const & a,
            Tensor4       & B);

// 10) B = a * A
/// \f$ [B]=\alpha([A]) \quad\equiv\quad \TeFo{B}\gets\alpha\TeFo{A} \f$
void CopyScale (double  const & a,
                Tensor4 const & A,
                Tensor4       & B);

// 11) s = xt * A * y
/// \f$ s=\{x\}^T[A]\{y\} \quad\equiv\quad s=\Dc{\Dc{\TeSe{x}}{ \TeFo{A}} }{\TeSe{y}} \f$
double Reduce(Tensor2 const & x,
              Tensor4 const & A,
              Tensor2 const & y);


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

// 0) v = A * u
inline void Dot(Tensor2 const & A, Tensor1 const & u, Tensor1 & v)
{
	v(0) = A(0)*u(0) + A(3)*u(1) + A(5)*u(2);
	v(1) = A(3)*u(0) + A(1)*u(1) + A(4)*u(2);
	v(2) = A(5)*u(0) + A(4)*u(1) + A(2)*u(2);
}

// 1) sc = x * y
inline double Dot (Tensor2 const & x, Tensor2 const & y)
{
	return blitz::dot(x,y);
}

// 2) y = x * A
inline void Dot (Tensor2 const & x, Tensor4 const & A, Tensor2  & y)
{
	y(0) = x(0)*A(0,0) + x(1)*A(1,0) + x(2)*A(2,0) + x(3)*A(3,0) + x(4)*A(4,0) + x(5)*A(5,0);
	y(1) = x(0)*A(0,1) + x(1)*A(1,1) + x(2)*A(2,1) + x(3)*A(3,1) + x(4)*A(4,1) + x(5)*A(5,1);
	y(2) = x(0)*A(0,2) + x(1)*A(1,2) + x(2)*A(2,2) + x(3)*A(3,2) + x(4)*A(4,2) + x(5)*A(5,2);
	y(3) = x(0)*A(0,3) + x(1)*A(1,3) + x(2)*A(2,3) + x(3)*A(3,3) + x(4)*A(4,3) + x(5)*A(5,3);
	y(4) = x(0)*A(0,4) + x(1)*A(1,4) + x(2)*A(2,4) + x(3)*A(3,4) + x(4)*A(4,4) + x(5)*A(5,4);
	y(5) = x(0)*A(0,5) + x(1)*A(1,5) + x(2)*A(2,5) + x(3)*A(3,5) + x(4)*A(4,5) + x(5)*A(5,5);
}

// 2a) y = s * x * A
inline void DotScaled (double s, Tensor2 const & x, Tensor4 const & A, Tensor2  & y)
{
	y(0) = s*x(0)*A(0,0) + s*x(1)*A(1,0) + s*x(2)*A(2,0) + s*x(3)*A(3,0) + s*x(4)*A(4,0) + s*x(5)*A(5,0);
	y(1) = s*x(0)*A(0,1) + s*x(1)*A(1,1) + s*x(2)*A(2,1) + s*x(3)*A(3,1) + s*x(4)*A(4,1) + s*x(5)*A(5,1);
	y(2) = s*x(0)*A(0,2) + s*x(1)*A(1,2) + s*x(2)*A(2,2) + s*x(3)*A(3,2) + s*x(4)*A(4,2) + s*x(5)*A(5,2);
	y(3) = s*x(0)*A(0,3) + s*x(1)*A(1,3) + s*x(2)*A(2,3) + s*x(3)*A(3,3) + s*x(4)*A(4,3) + s*x(5)*A(5,3);
	y(4) = s*x(0)*A(0,4) + s*x(1)*A(1,4) + s*x(2)*A(2,4) + s*x(3)*A(3,4) + s*x(4)*A(4,4) + s*x(5)*A(5,4);
	y(5) = s*x(0)*A(0,5) + s*x(1)*A(1,5) + s*x(2)*A(2,5) + s*x(3)*A(3,5) + s*x(4)*A(4,5) + s*x(5)*A(5,5);
}

// 3) y = A * x
inline void Dot (Tensor4 const & A, Tensor2 const & x, Tensor2 & y)
{
	y(0) = A(0,0)*x(0) + A(0,1)*x(1) + A(0,2)*x(2) + A(0,3)*x(3) + A(0,4)*x(4) + A(0,5)*x(5);
	y(1) = A(1,0)*x(0) + A(1,1)*x(1) + A(1,2)*x(2) + A(1,3)*x(3) + A(1,4)*x(4) + A(1,5)*x(5);
	y(2) = A(2,0)*x(0) + A(2,1)*x(1) + A(2,2)*x(2) + A(2,3)*x(3) + A(2,4)*x(4) + A(2,5)*x(5);
	y(3) = A(3,0)*x(0) + A(3,1)*x(1) + A(3,2)*x(2) + A(3,3)*x(3) + A(3,4)*x(4) + A(3,5)*x(5);
	y(4) = A(4,0)*x(0) + A(4,1)*x(1) + A(4,2)*x(2) + A(4,3)*x(3) + A(4,4)*x(4) + A(4,5)*x(5);
	y(5) = A(5,0)*x(0) + A(5,1)*x(1) + A(5,2)*x(2) + A(5,3)*x(3) + A(5,4)*x(4) + A(5,5)*x(5);
}

// 4) C = A * B
inline void Dot (Tensor4 const & A, Tensor4 const & B, Tensor4 & C)
{
	C = blitz::product(A,B);
}

// 5) A = x dyadic y
inline void Dyad (Tensor2 const & x, Tensor2 const & y, Tensor4  & A)
{
	A(0,0)=x(0)*y(0); A(0,1)=x(0)*y(1); A(0,2)=x(0)*y(2); A(0,3)=x(0)*y(3); A(0,4)=x(0)*y(4); A(0,5)=x(0)*y(5);
	A(1,0)=x(1)*y(0); A(1,1)=x(1)*y(1); A(1,2)=x(1)*y(2); A(1,3)=x(1)*y(3); A(1,4)=x(1)*y(4); A(1,5)=x(1)*y(5);
	A(2,0)=x(2)*y(0); A(2,1)=x(2)*y(1); A(2,2)=x(2)*y(2); A(2,3)=x(2)*y(3); A(2,4)=x(2)*y(4); A(2,5)=x(2)*y(5);
	A(3,0)=x(3)*y(0); A(3,1)=x(3)*y(1); A(3,2)=x(3)*y(2); A(3,3)=x(3)*y(3); A(3,4)=x(3)*y(4); A(3,5)=x(3)*y(5);
	A(4,0)=x(4)*y(0); A(4,1)=x(4)*y(1); A(4,2)=x(4)*y(2); A(4,3)=x(4)*y(3); A(4,4)=x(4)*y(4); A(4,5)=x(4)*y(5);
	A(5,0)=x(5)*y(0); A(5,1)=x(5)*y(1); A(5,2)=x(5)*y(2); A(5,3)=x(5)*y(3); A(5,4)=x(5)*y(4); A(5,5)=x(5)*y(5);
}

// 6aa) Z += a*X
inline void AddScaledUp (double const & a, Tensor4 const & X, Tensor4 & Z)
{
	Z(0,0)+=a*X(0,0); Z(0,1)+=a*X(0,1); Z(0,2)+=a*X(0,2); Z(0,3)+=a*X(0,3); Z(0,4)+=a*X(0,4); Z(0,5)+=a*X(0,5);
	Z(1,0)+=a*X(1,0); Z(1,1)+=a*X(1,1); Z(1,2)+=a*X(1,2); Z(1,3)+=a*X(1,3); Z(1,4)+=a*X(1,4); Z(1,5)+=a*X(1,5);
	Z(2,0)+=a*X(2,0); Z(2,1)+=a*X(2,1); Z(2,2)+=a*X(2,2); Z(2,3)+=a*X(2,3); Z(2,4)+=a*X(2,4); Z(2,5)+=a*X(2,5);
	Z(3,0)+=a*X(3,0); Z(3,1)+=a*X(3,1); Z(3,2)+=a*X(3,2); Z(3,3)+=a*X(3,3); Z(3,4)+=a*X(3,4); Z(3,5)+=a*X(3,5);
	Z(4,0)+=a*X(4,0); Z(4,1)+=a*X(4,1); Z(4,2)+=a*X(4,2); Z(4,3)+=a*X(4,3); Z(4,4)+=a*X(4,4); Z(4,5)+=a*X(4,5);
	Z(5,0)+=a*X(5,0); Z(5,1)+=a*X(5,1); Z(5,2)+=a*X(5,2); Z(5,3)+=a*X(5,3); Z(5,4)+=a*X(5,4); Z(5,5)+=a*X(5,5);
}

// 6a) Z = a*X + b*Y 
inline void AddScaled (double const & a, Tensor4 const & X, double const & b, Tensor4 const & Y, Tensor4 & Z)
{
	Z(0,0)=a*X(0,0)+b*Y(0,0); Z(0,1)=a*X(0,1)+b*Y(0,1); Z(0,2)=a*X(0,2)+b*Y(0,2); Z(0,3)=a*X(0,3)+b*Y(0,3); Z(0,4)=a*X(0,4)+b*Y(0,4); Z(0,5)=a*X(0,5)+b*Y(0,5);
	Z(1,0)=a*X(1,0)+b*Y(1,0); Z(1,1)=a*X(1,1)+b*Y(1,1); Z(1,2)=a*X(1,2)+b*Y(1,2); Z(1,3)=a*X(1,3)+b*Y(1,3); Z(1,4)=a*X(1,4)+b*Y(1,4); Z(1,5)=a*X(1,5)+b*Y(1,5);
	Z(2,0)=a*X(2,0)+b*Y(2,0); Z(2,1)=a*X(2,1)+b*Y(2,1); Z(2,2)=a*X(2,2)+b*Y(2,2); Z(2,3)=a*X(2,3)+b*Y(2,3); Z(2,4)=a*X(2,4)+b*Y(2,4); Z(2,5)=a*X(2,5)+b*Y(2,5);
	Z(3,0)=a*X(3,0)+b*Y(3,0); Z(3,1)=a*X(3,1)+b*Y(3,1); Z(3,2)=a*X(3,2)+b*Y(3,2); Z(3,3)=a*X(3,3)+b*Y(3,3); Z(3,4)=a*X(3,4)+b*Y(3,4); Z(3,5)=a*X(3,5)+b*Y(3,5);
	Z(4,0)=a*X(4,0)+b*Y(4,0); Z(4,1)=a*X(4,1)+b*Y(4,1); Z(4,2)=a*X(4,2)+b*Y(4,2); Z(4,3)=a*X(4,3)+b*Y(4,3); Z(4,4)=a*X(4,4)+b*Y(4,4); Z(4,5)=a*X(4,5)+b*Y(4,5);
	Z(5,0)=a*X(5,0)+b*Y(5,0); Z(5,1)=a*X(5,1)+b*Y(5,1); Z(5,2)=a*X(5,2)+b*Y(5,2); Z(5,3)=a*X(5,3)+b*Y(5,3); Z(5,4)=a*X(5,4)+b*Y(5,4); Z(5,5)=a*X(5,5)+b*Y(5,5);
}

// 6b) Z += a*X + b*Y 
inline void AddScaledUp (double const & a, Tensor4 const & X, double const & b, Tensor4 const & Y, Tensor4 & Z)
{
	Z(0,0)+=a*X(0,0)+b*Y(0,0); Z(0,1)+=a*X(0,1)+b*Y(0,1); Z(0,2)+=a*X(0,2)+b*Y(0,2); Z(0,3)+=a*X(0,3)+b*Y(0,3); Z(0,4)+=a*X(0,4)+b*Y(0,4); Z(0,5)+=a*X(0,5)+b*Y(0,5);
	Z(1,0)+=a*X(1,0)+b*Y(1,0); Z(1,1)+=a*X(1,1)+b*Y(1,1); Z(1,2)+=a*X(1,2)+b*Y(1,2); Z(1,3)+=a*X(1,3)+b*Y(1,3); Z(1,4)+=a*X(1,4)+b*Y(1,4); Z(1,5)+=a*X(1,5)+b*Y(1,5);
	Z(2,0)+=a*X(2,0)+b*Y(2,0); Z(2,1)+=a*X(2,1)+b*Y(2,1); Z(2,2)+=a*X(2,2)+b*Y(2,2); Z(2,3)+=a*X(2,3)+b*Y(2,3); Z(2,4)+=a*X(2,4)+b*Y(2,4); Z(2,5)+=a*X(2,5)+b*Y(2,5);
	Z(3,0)+=a*X(3,0)+b*Y(3,0); Z(3,1)+=a*X(3,1)+b*Y(3,1); Z(3,2)+=a*X(3,2)+b*Y(3,2); Z(3,3)+=a*X(3,3)+b*Y(3,3); Z(3,4)+=a*X(3,4)+b*Y(3,4); Z(3,5)+=a*X(3,5)+b*Y(3,5);
	Z(4,0)+=a*X(4,0)+b*Y(4,0); Z(4,1)+=a*X(4,1)+b*Y(4,1); Z(4,2)+=a*X(4,2)+b*Y(4,2); Z(4,3)+=a*X(4,3)+b*Y(4,3); Z(4,4)+=a*X(4,4)+b*Y(4,4); Z(4,5)+=a*X(4,5)+b*Y(4,5);
	Z(5,0)+=a*X(5,0)+b*Y(5,0); Z(5,1)+=a*X(5,1)+b*Y(5,1); Z(5,2)+=a*X(5,2)+b*Y(5,2); Z(5,3)+=a*X(5,3)+b*Y(5,3); Z(5,4)+=a*X(5,4)+b*Y(5,4); Z(5,5)+=a*X(5,5)+b*Y(5,5);
}

// 6c) W += a*X + b*Y +c*Z
inline void AddScaledUp (double const & a, Tensor4 const & X, double const & b, Tensor4 const & Y, double const & c, Tensor4 const & Z, Tensor4 & W)
{
	W(0,0)+=a*X(0,0)+b*Y(0,0)+c*Z(0,0); W(0,1)+=a*X(0,1)+b*Y(0,1)+c*Z(0,1); W(0,2)+=a*X(0,2)+b*Y(0,2)+c*Z(0,2); W(0,3)+=a*X(0,3)+b*Y(0,3)+c*Z(0,3); W(0,4)+=a*X(0,4)+b*Y(0,4)+c*Z(0,4); W(0,5)+=a*X(0,5)+b*Y(0,5)+c*Z(0,5);
	W(1,0)+=a*X(1,0)+b*Y(1,0)+c*Z(1,0); W(1,1)+=a*X(1,1)+b*Y(1,1)+c*Z(1,1); W(1,2)+=a*X(1,2)+b*Y(1,2)+c*Z(1,2); W(1,3)+=a*X(1,3)+b*Y(1,3)+c*Z(1,3); W(1,4)+=a*X(1,4)+b*Y(1,4)+c*Z(1,4); W(1,5)+=a*X(1,5)+b*Y(1,5)+c*Z(1,5);
	W(2,0)+=a*X(2,0)+b*Y(2,0)+c*Z(2,0); W(2,1)+=a*X(2,1)+b*Y(2,1)+c*Z(2,1); W(2,2)+=a*X(2,2)+b*Y(2,2)+c*Z(2,2); W(2,3)+=a*X(2,3)+b*Y(2,3)+c*Z(2,3); W(2,4)+=a*X(2,4)+b*Y(2,4)+c*Z(2,4); W(2,5)+=a*X(2,5)+b*Y(2,5)+c*Z(2,5);
	W(3,0)+=a*X(3,0)+b*Y(3,0)+c*Z(3,0); W(3,1)+=a*X(3,1)+b*Y(3,1)+c*Z(3,1); W(3,2)+=a*X(3,2)+b*Y(3,2)+c*Z(3,2); W(3,3)+=a*X(3,3)+b*Y(3,3)+c*Z(3,3); W(3,4)+=a*X(3,4)+b*Y(3,4)+c*Z(3,4); W(3,5)+=a*X(3,5)+b*Y(3,5)+c*Z(3,5);
	W(4,0)+=a*X(4,0)+b*Y(4,0)+c*Z(4,0); W(4,1)+=a*X(4,1)+b*Y(4,1)+c*Z(4,1); W(4,2)+=a*X(4,2)+b*Y(4,2)+c*Z(4,2); W(4,3)+=a*X(4,3)+b*Y(4,3)+c*Z(4,3); W(4,4)+=a*X(4,4)+b*Y(4,4)+c*Z(4,4); W(4,5)+=a*X(4,5)+b*Y(4,5)+c*Z(4,5);
	W(5,0)+=a*X(5,0)+b*Y(5,0)+c*Z(5,0); W(5,1)+=a*X(5,1)+b*Y(5,1)+c*Z(5,1); W(5,2)+=a*X(5,2)+b*Y(5,2)+c*Z(5,2); W(5,3)+=a*X(5,3)+b*Y(5,3)+c*Z(5,3); W(5,4)+=a*X(5,4)+b*Y(5,4)+c*Z(5,4); W(5,5)+=a*X(5,5)+b*Y(5,5)+c*Z(5,5);
}

// 6d) W = a*X + b*Y +c*Z + d*D
inline void AddScaled (double const & a, Tensor4 const & X, double const & b, Tensor4 const & Y, double const & c, Tensor4 const & Z, double const & d, Tensor4 const & D,  Tensor4 & W)
{
	W(0,0)=a*X(0,0)+b*Y(0,0)+c*Z(0,0)+d*D(0,0); W(0,1)=a*X(0,1)+b*Y(0,1)+c*Z(0,1)+d*D(0,1); W(0,2)=a*X(0,2)+b*Y(0,2)+c*Z(0,2)+d*D(0,2); W(0,3)=a*X(0,3)+b*Y(0,3)+c*Z(0,3)+d*D(0,3); W(0,4)=a*X(0,4)+b*Y(0,4)+c*Z(0,4)+d*D(0,4); W(0,5)=a*X(0,5)+b*Y(0,5)+c*Z(0,5)+d*D(0,5);
	W(1,0)=a*X(1,0)+b*Y(1,0)+c*Z(1,0)+d*D(1,0); W(1,1)=a*X(1,1)+b*Y(1,1)+c*Z(1,1)+d*D(1,1); W(1,2)=a*X(1,2)+b*Y(1,2)+c*Z(1,2)+d*D(1,2); W(1,3)=a*X(1,3)+b*Y(1,3)+c*Z(1,3)+d*D(1,3); W(1,4)=a*X(1,4)+b*Y(1,4)+c*Z(1,4)+d*D(1,4); W(1,5)=a*X(1,5)+b*Y(1,5)+c*Z(1,5)+d*D(1,5);
	W(2,0)=a*X(2,0)+b*Y(2,0)+c*Z(2,0)+d*D(2,0); W(2,1)=a*X(2,1)+b*Y(2,1)+c*Z(2,1)+d*D(2,1); W(2,2)=a*X(2,2)+b*Y(2,2)+c*Z(2,2)+d*D(2,2); W(2,3)=a*X(2,3)+b*Y(2,3)+c*Z(2,3)+d*D(2,3); W(2,4)=a*X(2,4)+b*Y(2,4)+c*Z(2,4)+d*D(2,4); W(2,5)=a*X(2,5)+b*Y(2,5)+c*Z(2,5)+d*D(2,5);
	W(3,0)=a*X(3,0)+b*Y(3,0)+c*Z(3,0)+d*D(3,0); W(3,1)=a*X(3,1)+b*Y(3,1)+c*Z(3,1)+d*D(3,1); W(3,2)=a*X(3,2)+b*Y(3,2)+c*Z(3,2)+d*D(3,2); W(3,3)=a*X(3,3)+b*Y(3,3)+c*Z(3,3)+d*D(3,3); W(3,4)=a*X(3,4)+b*Y(3,4)+c*Z(3,4)+d*D(3,4); W(3,5)=a*X(3,5)+b*Y(3,5)+c*Z(3,5)+d*D(3,5);
	W(4,0)=a*X(4,0)+b*Y(4,0)+c*Z(4,0)+d*D(4,0); W(4,1)=a*X(4,1)+b*Y(4,1)+c*Z(4,1)+d*D(4,1); W(4,2)=a*X(4,2)+b*Y(4,2)+c*Z(4,2)+d*D(4,2); W(4,3)=a*X(4,3)+b*Y(4,3)+c*Z(4,3)+d*D(4,3); W(4,4)=a*X(4,4)+b*Y(4,4)+c*Z(4,4)+d*D(4,4); W(4,5)=a*X(4,5)+b*Y(4,5)+c*Z(4,5)+d*D(4,5);
	W(5,0)=a*X(5,0)+b*Y(5,0)+c*Z(5,0)+d*D(5,0); W(5,1)=a*X(5,1)+b*Y(5,1)+c*Z(5,1)+d*D(5,1); W(5,2)=a*X(5,2)+b*Y(5,2)+c*Z(5,2)+d*D(5,2); W(5,3)=a*X(5,3)+b*Y(5,3)+c*Z(5,3)+d*D(5,3); W(5,4)=a*X(5,4)+b*Y(5,4)+c*Z(5,4)+d*D(5,4); W(5,5)=a*X(5,5)+b*Y(5,5)+c*Z(5,5)+d*D(5,5);
}

// 7) A = a * (x dyadic y) + A
inline void Ger (double const & a, Tensor2 const & x, Tensor2 const & y, Tensor4  & A)
{
	A(0,0)=a*x(0)*y(0)+A(0,0); A(0,1)=a*x(0)*y(1)+A(0,1); A(0,2)=a*x(0)*y(2)+A(0,2); A(0,3)=a*x(0)*y(3)+A(0,3); A(0,4)=a*x(0)*y(4)+A(0,4); A(0,5)=a*x(0)*y(5)+A(0,5);
	A(1,0)=a*x(1)*y(0)+A(1,0); A(1,1)=a*x(1)*y(1)+A(1,1); A(1,2)=a*x(1)*y(2)+A(1,2); A(1,3)=a*x(1)*y(3)+A(1,3); A(1,4)=a*x(1)*y(4)+A(1,4); A(1,5)=a*x(1)*y(5)+A(1,5);
	A(2,0)=a*x(2)*y(0)+A(2,0); A(2,1)=a*x(2)*y(1)+A(2,1); A(2,2)=a*x(2)*y(2)+A(2,2); A(2,3)=a*x(2)*y(3)+A(2,3); A(2,4)=a*x(2)*y(4)+A(2,4); A(2,5)=a*x(2)*y(5)+A(2,5);
	A(3,0)=a*x(3)*y(0)+A(3,0); A(3,1)=a*x(3)*y(1)+A(3,1); A(3,2)=a*x(3)*y(2)+A(3,2); A(3,3)=a*x(3)*y(3)+A(3,3); A(3,4)=a*x(3)*y(4)+A(3,4); A(3,5)=a*x(3)*y(5)+A(3,5);
	A(4,0)=a*x(4)*y(0)+A(4,0); A(4,1)=a*x(4)*y(1)+A(4,1); A(4,2)=a*x(4)*y(2)+A(4,2); A(4,3)=a*x(4)*y(3)+A(4,3); A(4,4)=a*x(4)*y(4)+A(4,4); A(4,5)=a*x(4)*y(5)+A(4,5);
	A(5,0)=a*x(5)*y(0)+A(5,0); A(5,1)=a*x(5)*y(1)+A(5,1); A(5,2)=a*x(5)*y(2)+A(5,2); A(5,3)=a*x(5)*y(3)+A(5,3); A(5,4)=a*x(5)*y(4)+A(5,4); A(5,5)=a*x(5)*y(5)+A(5,5);
}

// 8) D = a * (A*x) dyadic (y*B) + C
inline void GerX (double const & a, Tensor4 const & A, Tensor2 const & x, Tensor2 const & y, Tensor4 const & B, Tensor4 const & C, Tensor4 & D)
{
	//C(i,j)=a* (A(i,k)*x(k)) * (y(l)*B(l,j)) +D(i,j);

	D(0,0)=a*(A(0,0)*x(0)+A(0,1)*x(1)+A(0,2)*x(2)+A(0,3)*x(3)+A(0,4)*x(4)+A(0,5)*x(5)) * (y(0)*B(0,0)+y(1)*B(1,0)+y(2)*B(2,0)+y(3)*B(3,0)+y(4)*B(4,0)+y(5)*B(5,0)) + C(0,0);
	D(0,1)=a*(A(0,0)*x(0)+A(0,1)*x(1)+A(0,2)*x(2)+A(0,3)*x(3)+A(0,4)*x(4)+A(0,5)*x(5)) * (y(0)*B(0,1)+y(1)*B(1,1)+y(2)*B(2,1)+y(3)*B(3,1)+y(4)*B(4,1)+y(5)*B(5,1)) + C(0,1);
	D(0,2)=a*(A(0,0)*x(0)+A(0,1)*x(1)+A(0,2)*x(2)+A(0,3)*x(3)+A(0,4)*x(4)+A(0,5)*x(5)) * (y(0)*B(0,2)+y(1)*B(1,2)+y(2)*B(2,2)+y(3)*B(3,2)+y(4)*B(4,2)+y(5)*B(5,2)) + C(0,2);
	D(0,3)=a*(A(0,0)*x(0)+A(0,1)*x(1)+A(0,2)*x(2)+A(0,3)*x(3)+A(0,4)*x(4)+A(0,5)*x(5)) * (y(0)*B(0,3)+y(1)*B(1,3)+y(2)*B(2,3)+y(3)*B(3,3)+y(4)*B(4,3)+y(5)*B(5,3)) + C(0,3);
	D(0,4)=a*(A(0,0)*x(0)+A(0,1)*x(1)+A(0,2)*x(2)+A(0,3)*x(3)+A(0,4)*x(4)+A(0,5)*x(5)) * (y(0)*B(0,4)+y(1)*B(1,4)+y(2)*B(2,4)+y(3)*B(3,4)+y(4)*B(4,4)+y(5)*B(5,4)) + C(0,4);
	D(0,5)=a*(A(0,0)*x(0)+A(0,1)*x(1)+A(0,2)*x(2)+A(0,3)*x(3)+A(0,4)*x(4)+A(0,5)*x(5)) * (y(0)*B(0,5)+y(1)*B(1,5)+y(2)*B(2,5)+y(3)*B(3,5)+y(4)*B(4,5)+y(5)*B(5,5)) + C(0,5);
	
	D(1,0)=a*(A(1,0)*x(0)+A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5)) * (y(0)*B(0,0)+y(1)*B(1,0)+y(2)*B(2,0)+y(3)*B(3,0)+y(4)*B(4,0)+y(5)*B(5,0)) + C(1,0);
	D(1,1)=a*(A(1,0)*x(0)+A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5)) * (y(0)*B(0,1)+y(1)*B(1,1)+y(2)*B(2,1)+y(3)*B(3,1)+y(4)*B(4,1)+y(5)*B(5,1)) + C(1,1);
	D(1,2)=a*(A(1,0)*x(0)+A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5)) * (y(0)*B(0,2)+y(1)*B(1,2)+y(2)*B(2,2)+y(3)*B(3,2)+y(4)*B(4,2)+y(5)*B(5,2)) + C(1,2);
	D(1,3)=a*(A(1,0)*x(0)+A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5)) * (y(0)*B(0,3)+y(1)*B(1,3)+y(2)*B(2,3)+y(3)*B(3,3)+y(4)*B(4,3)+y(5)*B(5,3)) + C(1,3);
	D(1,4)=a*(A(1,0)*x(0)+A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5)) * (y(0)*B(0,4)+y(1)*B(1,4)+y(2)*B(2,4)+y(3)*B(3,4)+y(4)*B(4,4)+y(5)*B(5,4)) + C(1,4);
	D(1,5)=a*(A(1,0)*x(0)+A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5)) * (y(0)*B(0,5)+y(1)*B(1,5)+y(2)*B(2,5)+y(3)*B(3,5)+y(4)*B(4,5)+y(5)*B(5,5)) + C(1,5);
	
	D(2,0)=a*(A(2,0)*x(0)+A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5)) * (y(0)*B(0,0)+y(1)*B(1,0)+y(2)*B(2,0)+y(3)*B(3,0)+y(4)*B(4,0)+y(5)*B(5,0)) + C(2,0);
	D(2,1)=a*(A(2,0)*x(0)+A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5)) * (y(0)*B(0,1)+y(1)*B(1,1)+y(2)*B(2,1)+y(3)*B(3,1)+y(4)*B(4,1)+y(5)*B(5,1)) + C(2,1);
	D(2,2)=a*(A(2,0)*x(0)+A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5)) * (y(0)*B(0,2)+y(1)*B(1,2)+y(2)*B(2,2)+y(3)*B(3,2)+y(4)*B(4,2)+y(5)*B(5,2)) + C(2,2);
	D(2,3)=a*(A(2,0)*x(0)+A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5)) * (y(0)*B(0,3)+y(1)*B(1,3)+y(2)*B(2,3)+y(3)*B(3,3)+y(4)*B(4,3)+y(5)*B(5,3)) + C(2,3);
	D(2,4)=a*(A(2,0)*x(0)+A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5)) * (y(0)*B(0,4)+y(1)*B(1,4)+y(2)*B(2,4)+y(3)*B(3,4)+y(4)*B(4,4)+y(5)*B(5,4)) + C(2,4);
	D(2,5)=a*(A(2,0)*x(0)+A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5)) * (y(0)*B(0,5)+y(1)*B(1,5)+y(2)*B(2,5)+y(3)*B(3,5)+y(4)*B(4,5)+y(5)*B(5,5)) + C(2,5);
	
	D(3,0)=a*(A(3,0)*x(0)+A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5)) * (y(0)*B(0,0)+y(1)*B(1,0)+y(2)*B(2,0)+y(3)*B(3,0)+y(4)*B(4,0)+y(5)*B(5,0)) + C(3,0);
	D(3,1)=a*(A(3,0)*x(0)+A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5)) * (y(0)*B(0,1)+y(1)*B(1,1)+y(2)*B(2,1)+y(3)*B(3,1)+y(4)*B(4,1)+y(5)*B(5,1)) + C(3,1);
	D(3,2)=a*(A(3,0)*x(0)+A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5)) * (y(0)*B(0,2)+y(1)*B(1,2)+y(2)*B(2,2)+y(3)*B(3,2)+y(4)*B(4,2)+y(5)*B(5,2)) + C(3,2);
	D(3,3)=a*(A(3,0)*x(0)+A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5)) * (y(0)*B(0,3)+y(1)*B(1,3)+y(2)*B(2,3)+y(3)*B(3,3)+y(4)*B(4,3)+y(5)*B(5,3)) + C(3,3);
	D(3,4)=a*(A(3,0)*x(0)+A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5)) * (y(0)*B(0,4)+y(1)*B(1,4)+y(2)*B(2,4)+y(3)*B(3,4)+y(4)*B(4,4)+y(5)*B(5,4)) + C(3,4);
	D(3,5)=a*(A(3,0)*x(0)+A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5)) * (y(0)*B(0,5)+y(1)*B(1,5)+y(2)*B(2,5)+y(3)*B(3,5)+y(4)*B(4,5)+y(5)*B(5,5)) + C(3,5);
	
	D(4,0)=a*(A(4,0)*x(0)+A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5)) * (y(0)*B(0,0)+y(1)*B(1,0)+y(2)*B(2,0)+y(3)*B(3,0)+y(4)*B(4,0)+y(5)*B(5,0)) + C(4,0);
	D(4,1)=a*(A(4,0)*x(0)+A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5)) * (y(0)*B(0,1)+y(1)*B(1,1)+y(2)*B(2,1)+y(3)*B(3,1)+y(4)*B(4,1)+y(5)*B(5,1)) + C(4,1);
	D(4,2)=a*(A(4,0)*x(0)+A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5)) * (y(0)*B(0,2)+y(1)*B(1,2)+y(2)*B(2,2)+y(3)*B(3,2)+y(4)*B(4,2)+y(5)*B(5,2)) + C(4,2);
	D(4,3)=a*(A(4,0)*x(0)+A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5)) * (y(0)*B(0,3)+y(1)*B(1,3)+y(2)*B(2,3)+y(3)*B(3,3)+y(4)*B(4,3)+y(5)*B(5,3)) + C(4,3);
	D(4,4)=a*(A(4,0)*x(0)+A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5)) * (y(0)*B(0,4)+y(1)*B(1,4)+y(2)*B(2,4)+y(3)*B(3,4)+y(4)*B(4,4)+y(5)*B(5,4)) + C(4,4);
	D(4,5)=a*(A(4,0)*x(0)+A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5)) * (y(0)*B(0,5)+y(1)*B(1,5)+y(2)*B(2,5)+y(3)*B(3,5)+y(4)*B(4,5)+y(5)*B(5,5)) + C(4,5);
	
	D(5,0)=a*(A(5,0)*x(0)+A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5)) * (y(0)*B(0,0)+y(1)*B(1,0)+y(2)*B(2,0)+y(3)*B(3,0)+y(4)*B(4,0)+y(5)*B(5,0)) + C(5,0);
	D(5,1)=a*(A(5,0)*x(0)+A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5)) * (y(0)*B(0,1)+y(1)*B(1,1)+y(2)*B(2,1)+y(3)*B(3,1)+y(4)*B(4,1)+y(5)*B(5,1)) + C(5,1);
	D(5,2)=a*(A(5,0)*x(0)+A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5)) * (y(0)*B(0,2)+y(1)*B(1,2)+y(2)*B(2,2)+y(3)*B(3,2)+y(4)*B(4,2)+y(5)*B(5,2)) + C(5,2);
	D(5,3)=a*(A(5,0)*x(0)+A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5)) * (y(0)*B(0,3)+y(1)*B(1,3)+y(2)*B(2,3)+y(3)*B(3,3)+y(4)*B(4,3)+y(5)*B(5,3)) + C(5,3);
	D(5,4)=a*(A(5,0)*x(0)+A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5)) * (y(0)*B(0,4)+y(1)*B(1,4)+y(2)*B(2,4)+y(3)*B(3,4)+y(4)*B(4,4)+y(5)*B(5,4)) + C(5,4);
	D(5,5)=a*(A(5,0)*x(0)+A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5)) * (y(0)*B(0,5)+y(1)*B(1,5)+y(2)*B(2,5)+y(3)*B(3,5)+y(4)*B(4,5)+y(5)*B(5,5)) + C(5,5);
}

// 9) B = a * B
inline void Scale (double const & a, Tensor4 & B)
{
	B(0,0)=a*B(0,0); B(0,1)=a*B(0,1); B(0,2)=a*B(0,2); B(0,3)=a*B(0,3); B(0,4)=a*B(0,4); B(0,5)=a*B(0,5);
	B(1,0)=a*B(1,0); B(1,1)=a*B(1,1); B(1,2)=a*B(1,2); B(1,3)=a*B(1,3); B(1,4)=a*B(1,4); B(1,5)=a*B(1,5);
	B(2,0)=a*B(2,0); B(2,1)=a*B(2,1); B(2,2)=a*B(2,2); B(2,3)=a*B(2,3); B(2,4)=a*B(2,4); B(2,5)=a*B(2,5);
	B(3,0)=a*B(3,0); B(3,1)=a*B(3,1); B(3,2)=a*B(3,2); B(3,3)=a*B(3,3); B(3,4)=a*B(3,4); B(3,5)=a*B(3,5);
	B(4,0)=a*B(4,0); B(4,1)=a*B(4,1); B(4,2)=a*B(4,2); B(4,3)=a*B(4,3); B(4,4)=a*B(4,4); B(4,5)=a*B(4,5);
	B(5,0)=a*B(5,0); B(5,1)=a*B(5,1); B(5,2)=a*B(5,2); B(5,3)=a*B(5,3); B(5,4)=a*B(5,4); B(5,5)=a*B(5,5);
}

// 10) B = a * A
inline void CopyScale (double const & a, Tensor4 const & A, Tensor4 & B)
{
	B(0,0)=a*A(0,0); B(0,1)=a*A(0,1); B(0,2)=a*A(0,2); B(0,3)=a*A(0,3); B(0,4)=a*A(0,4); B(0,5)=a*A(0,5);
	B(1,0)=a*A(1,0); B(1,1)=a*A(1,1); B(1,2)=a*A(1,2); B(1,3)=a*A(1,3); B(1,4)=a*A(1,4); B(1,5)=a*A(1,5);
	B(2,0)=a*A(2,0); B(2,1)=a*A(2,1); B(2,2)=a*A(2,2); B(2,3)=a*A(2,3); B(2,4)=a*A(2,4); B(2,5)=a*A(2,5);
	B(3,0)=a*A(3,0); B(3,1)=a*A(3,1); B(3,2)=a*A(3,2); B(3,3)=a*A(3,3); B(3,4)=a*A(3,4); B(3,5)=a*A(3,5);
	B(4,0)=a*A(4,0); B(4,1)=a*A(4,1); B(4,2)=a*A(4,2); B(4,3)=a*A(4,3); B(4,4)=a*A(4,4); B(4,5)=a*A(4,5);
	B(5,0)=a*A(5,0); B(5,1)=a*A(5,1); B(5,2)=a*A(5,2); B(5,3)=a*A(5,3); B(5,4)=a*A(5,4); B(5,5)=a*A(5,5);
}

// 11) s = xt * A * y
inline double Reduce(Tensor2 const & x, Tensor4 const & A, Tensor2 const & y)
{
	return x(0)*(A(0,0)*y(0) + A(0,1)*y(1) + A(0,2)*y(2) + A(0,3)*y(3) + A(0,4)*y(4) + A(0,5)*y(5)) +
	       x(1)*(A(1,0)*y(0) + A(1,1)*y(1) + A(1,2)*y(2) + A(1,3)*y(3) + A(1,4)*y(4) + A(1,5)*y(5)) +
	       x(2)*(A(2,0)*y(0) + A(2,1)*y(1) + A(2,2)*y(2) + A(2,3)*y(3) + A(2,4)*y(4) + A(2,5)*y(5)) +
	       x(3)*(A(3,0)*y(0) + A(3,1)*y(1) + A(3,2)*y(2) + A(3,3)*y(3) + A(3,4)*y(4) + A(3,5)*y(5)) +
	       x(4)*(A(4,0)*y(0) + A(4,1)*y(1) + A(4,2)*y(2) + A(4,3)*y(3) + A(4,4)*y(4) + A(4,5)*y(5)) +
	       x(5)*(A(5,0)*y(0) + A(5,1)*y(1) + A(5,2)*y(2) + A(5,3)*y(3) + A(5,4)*y(4) + A(5,5)*y(5));
}

// 12) T = ((A ^ B) + (B ^ A)) % Psd; // symmetric
inline void LeafLeafPsd (Tensor2 const & A, Tensor2 const & B,  Tensor4 & T)
{
	double a3 = A(3)/SQ2;
	double a4 = A(4)/SQ2;
	double a5 = A(5)/SQ2;
	double b3 = B(3)/SQ2;
	double b4 = B(4)/SQ2;
	double b5 = B(5)/SQ2;
	double c = 2.0/3.0;
	double d = 4.0/3.0;
	double e = 1.0/3.0;

	T = d*B(0)*A(0)-c*b3*a3-c*b5*a5                                    , -c*B(0)*A(0)+d*b3*a3-c*b5*a5                                   , -c*B(0)*A(0)-c*b3*a3+d*b5*a5                                   , SQ2*(b3*A(0)+a3*B(0))       , SQ2*(b3*a5+a3*b5)           , SQ2*(B(0)*a5+A(0)*b5)      , 
		d*b3*a3-c*b4*a4-c*B(1)*A(1)                                    , -c*b3*a3-c*b4*a4+d*B(1)*A(1)                                   , -c*b3*a3+d*b4*a4-c*B(1)*A(1)                                   , SQ2*(b3*A(1)+a3*B(1))       , SQ2*(b4*A(1)+a4*B(1))       , SQ2*(b3*a4+a3*b4)          , 
		d*b5*a5-c*b4*a4-c*B(2)*A(2)                                    , -c*b5*a5+d*b4*a4-c*B(2)*A(2)                                   , -c*b5*a5-c*b4*a4+d*B(2)*A(2)                                   , SQ2*(a4*b5+b4*a5)           , SQ2*(A(2)*b4+B(2)*a4)       , SQ2*(B(2)*a5+A(2)*b5)      , 
		SQ2*(-e*a4*b5-e*b3*A(1)+c*b3*A(0)-e*a3*B(1)-e*b4*a5+c*a3*B(0)) , SQ2*(-e*a4*b5+c*b3*A(1)-e*b3*A(0)+c*a3*B(1)-e*b4*a5-e*a3*B(0)) , SQ2*(c*a4*b5-e*b3*A(1)-e*b3*A(0)-e*a3*B(1)+c*b4*a5-e*a3*B(0))  , B(0)*A(1)+A(0)*B(1)+2*b3*a3 , a5*B(1)+b5*A(1)+b3*a4+a3*b4 , b3*a5+A(0)*b4+a3*b5+B(0)*a4, 
		SQ2*(c*b3*a5-e*A(2)*b4-e*B(2)*a4-e*b4*A(1)+c*a3*b5-e*a4*B(1))  , SQ2*(-e*b3*a5-e*A(2)*b4-e*B(2)*a4+c*b4*A(1)-e*a3*b5+c*a4*B(1)) , SQ2*(-e*b3*a5+c*A(2)*b4+c*B(2)*a4-e*b4*A(1)-e*a3*b5-e*a4*B(1)) , a5*B(1)+b5*A(1)+b3*a4+a3*b4 , A(2)*B(1)+2*b4*a4+B(2)*A(1) , a4*b5+a3*B(2)+b3*A(2)+b4*a5, 
		SQ2*(-e*b3*a4-e*B(2)*a5+c*B(0)*a5-e*A(2)*b5-e*a3*b4+c*A(0)*b5) , SQ2*(c*b3*a4-e*B(2)*a5-e*B(0)*a5-e*A(2)*b5+c*a3*b4-e*A(0)*b5)  , SQ2*(-e*b3*a4+c*B(2)*a5-e*B(0)*a5+c*A(2)*b5-e*a3*b4-e*A(0)*b5) , b3*a5+A(0)*b4+a3*b5+B(0)*a4 , a4*b5+a3*B(2)+b3*A(2)+b4*a5 , B(0)*A(2)+2*b5*a5+B(2)*A(0);
}

}; // namespace Tensors

#endif // MECHSYS_TENSORS_OPERATORS_H
