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

/* Jacobi Transformation of a Symmetric Matrix
 *
 * Numerical evaluation of eigenvalues and eigenvectors for 3x3 symmetric
 * matrices using the Jacobi-rotation Method.
 *
 * 2005-09-01 (C): Dorival M. Pedroso (Expanded all loops; Using tensor Tensor2 (Mandel))
 * 2005-05-10 (C): Dorival M. Pedroso (C++)
 * 2004-08-15 (C): Dorival M. Pedroso (Fortran 95)
 * 2002-04-15 (C): Dorival M. Pedroso (Fortran 77)
 *
 * The Jacobi method consists of a sequence of orthogonal similarity transformations.
 * Each transformation (a Jacobi rotation) is just a plane rotation designed to annihilate one of the off-diagonal matrix elements.
 * Successive transformations undo previously set zeros, but the off-diagonal elements nevertheless get smaller and smaller.
 * Accumulating the product of the transformations as you go gives the matrix (Q) of eigenvectors (V), while the elements of the final
 * diagonal matrix (A) are the eigenvalues (L).
 * The Jacobi method is absolutely foolproof for all real symmetric matrices.
 * For matrices of order greater than about 10, say, the algorithm is slower, by a significant constant factor, than the QR method.
 *       _      _         _      _         _      _ 
 *      | Q(0,0) |       | Q(0,1) |       | Q(0,2) |
 *  V0= | Q(1,0) |   V1= | Q(1,1) |   V2= | Q(1,2) |
 *      | Q(2,0) |       | Q(2,1) |       | Q(2,2) |
 *       -      -         -      -         -      - 
 */

#ifndef MECHSYS_TENSORS_JACOBIROT_H
#define MECHSYS_TENSORS_JACOBIROT_H

// STL
#include <cmath> // for sqrt()

// MechSys
#include "tensors/tensors.h"

namespace Tensors
{

/** Jacobi Transformation of a Symmetric Matrix (given as a Tensor2).
 * \f[
 *     V0 = \begin{Bmatrix} Q(0,0) \\ Q(1,0) \\ Q(2,0) \end{Bmatrix}\qquad
 *     V1 = \begin{Bmatrix} Q(0,1) \\ Q(1,1) \\ Q(2,1) \end{Bmatrix}\qquad
 *     V2 = \begin{Bmatrix} Q(0,2) \\ Q(1,2) \\ Q(2,2) \end{Bmatrix}
 * \f]
 * \return the number of iterations
 */
int JacobiRot(Tensors::Tensor2 const & T    , ///< In: Tensor2 (Mandel) corresponding to a the matrix (A) we seek for the eigenvalues (SYMMETRIC and square)
		      double                   V0[3], ///< Out: Eigenvector (array with 3 values)
			  double                   V1[3], ///< Out: Eigenvector (array with 3 values)
			  double                   V2[3], ///< Out: Eigenvector (array with 3 values)
			  double                   L[3]); ///< Out: Eigenvalues (array with 3 values)

/** Jacobi Transformation of a Symmetric Matrix (given as a Tensor2).
 * \return the number of iterations
 */
int JacobiRot(Tensors::Tensor2 const & T , ///< In: Tensor2 (Mandel) corresponding to a the matrix (A) we seek for the eigenvalues (SYMMETRIC and square)
			  Tensors::Tensor1       & L); ///< Out: Eigenvalues (array with 3 values)


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline int JacobiRot(Tensors::Tensor2 const & T, double V0[3], double V1[3], double V2[3], double L[3]) 
{
	const double maxIt  = 30;      // Max number of iterations
	const double errTol = 1.0e-15; // Erro: Tolerance
	const double zero   = 1.0e-10; // Erro: Tolerance
	
	double UT[3];         // Values of the Upper Triangle part of symmetric matrix A
	double th;            // theta = (Aii-Ajj)/2Aij
	double c;             // Cossine
	double s;             // Sine
	double cc;            // Cossine squared
	double ss;            // Sine squared
	double t;             // Tangent
	double Temp;          // Auxiliar variable
	double TM[3];         // Auxiliar array
	double sq2=sqrt(2.0); // Useful for the conversion: Tensor (Mandel) ==> matriz A
	double sumUT;         // Sum of upper triangle of abs(A) that measures the error
	int    it;            // Iteration number
	double h;             // Difference L[i]-L[j]
	
	// Initialize eigenvalues which correnspond to the diagonal part of A
	L[0]=T(0); L[1]=T(1); L[2]=T(2);

	// Initialize Upper Triangle part of A matrix (array[3])
	UT[0]=T(3)/sq2; UT[1]=T(4)/sq2; UT[2]=T(5)/sq2;

	// Initialize eigenvectors
	V0[0]=1.0; V1[0]=0.0; V2[0]=0.0;
	V0[1]=0.0; V1[1]=1.0; V2[1]=0.0;
	V0[2]=0.0; V1[2]=0.0; V2[2]=1.0;

	// Iterations
	for (it=1; it<=maxIt; ++it)
	{
		// Check error
		sumUT = fabs(UT[0])+fabs(UT[1])+fabs(UT[2]);
		if (sumUT<=errTol) return it;
		
		// i=0, j=1 ,r=2 (p=3)
		h=L[0]-L[1];
		if (fabs(h)<zero) t=1.0;
		else
		{
			th=0.5*h/UT[0];
			t=1.0/(fabs(th)+sqrt(th*th+1.0));
			if (th<0.0) t=-t;
		}
		c  = 1.0/sqrt(1.0+t*t);
		s  = c*t;
		cc = c*c;
		ss = s*s;
		// Zeroes term UT[0]
		Temp  = cc*L[0] + 2.0*c*s*UT[0] + ss*L[1];
		L[1]  = ss*L[0] - 2.0*c*s*UT[0] + cc*L[1];
		L[0]  = Temp;
		UT[0] = 0.0;
		Temp  = c*UT[2] + s*UT[1];
		UT[1] = c*UT[1] - s*UT[2];
		UT[2] = Temp;
		// Actualize eigenvectors
		TM[0] = s*V1[0] + c*V0[0];
		TM[1] = s*V1[1] + c*V0[1];
		TM[2] = s*V1[2] + c*V0[2];
		V1[0] = c*V1[0] - s*V0[0];
		V1[1] = c*V1[1] - s*V0[1];
		V1[2] = c*V1[2] - s*V0[2];
		V0[0] = TM[0];
		V0[1] = TM[1];
		V0[2] = TM[2];
		
		// i=1, j=2 ,r=0 (p=4)
		h = L[1]-L[2];
		if (fabs(h)<zero) t=1.0;
		else
		{
			th=0.5*h/UT[1];
			t=1.0/(fabs(th)+sqrt(th*th+1.0));
			if (th<0.0) t=-t;
		}
		c  = 1.0/sqrt(1.0+t*t);
		s  = c*t;
		cc = c*c;
		ss = s*s;
		// Zeroes term UT[1]
		Temp  = cc*L[1] + 2.0*c*s*UT[1] + ss*L[2];
		L[2]  = ss*L[1] - 2.0*c*s*UT[1] + cc*L[2];
		L[1]  = Temp;
		UT[1] = 0.0;
		Temp  = c*UT[0] + s*UT[2];
		UT[2] = c*UT[2] - s*UT[0];
		UT[0] = Temp;
		// Actualize eigenvectors
		TM[1] = s*V2[1] + c*V1[1];
		TM[2] = s*V2[2] + c*V1[2];
		TM[0] = s*V2[0] + c*V1[0];
		V2[1] = c*V2[1] - s*V1[1];
		V2[2] = c*V2[2] - s*V1[2];
		V2[0] = c*V2[0] - s*V1[0];
		V1[1] = TM[1];
		V1[2] = TM[2];
		V1[0] = TM[0];

		// i=0, j=2 ,r=1 (p=5)
		h = L[0]-L[2];
		if (fabs(h)<zero) t=1.0;
		else
		{
			th=0.5*h/UT[2];
			t=1.0/(fabs(th)+sqrt(th*th+1.0));
			if (th<0.0) t=-t;
		}
		c  = 1.0/sqrt(1.0+t*t);
		s  = c*t;
		cc = c*c;
		ss = s*s;
		// Zeroes term UT[2]
		Temp  = cc*L[0] + 2.0*c*s*UT[2] + ss*L[2];
		L[2]  = ss*L[0] - 2.0*c*s*UT[2] + cc*L[2];
		L[0]  = Temp;
		UT[2] = 0.0;
		Temp  = c*UT[0] + s*UT[1];
		UT[1] = c*UT[1] - s*UT[0];
		UT[0] = Temp;
		// Actualize eigenvectors
		TM[0] = s*V2[0] + c*V0[0];
		TM[2] = s*V2[2] + c*V0[2];
		TM[1] = s*V2[1] + c*V0[1];
		V2[0] = c*V2[0] - s*V0[0];
		V2[2] = c*V2[2] - s*V0[2];
		V2[1] = c*V2[1] - s*V0[1];
		V0[0] = TM[0];
		V0[2] = TM[2];
		V0[1] = TM[1];
	}

	return -1; // Did not converge

}

inline int JacobiRot(Tensors::Tensor2 const & T, Tensors::Tensor1 & L)
{
	const double maxIt  = 30;      // Max number of iterations
	const double errTol = 1.0e-15; // Erro: Tolerance
	const double zero   = 1.0e-10; // Erro: Tolerance
	
	double UT[3];         // Values of the Upper Triangle part of symmetric matrix A
	double th;            // theta = (Aii-Ajj)/2Aij
	double c;             // Cossine
	double s;             // Sine
	double cc;            // Cossine squared
	double ss;            // Sine squared
	double t;             // Tangent
	double Temp;          // Auxiliar variable
	double sq2=sqrt(2.0); // Useful for the conversion: Tensor (Mandel) ==> matriz A
	double sumUT;         // Sum of upper triangle of abs(A) that measures the error
	int    it;            // Iteration number
	double h;             // Difference L(i)-L(j)
	
	// Initialize eigenvalues which correnspond to the diagonal part of A
	L(0)=T(0); L(1)=T(1); L(2)=T(2);

	// Initialize Upper Triangle part of A matrix (array[3])
	UT[0]=T(3)/sq2; UT[1]=T(4)/sq2; UT[2]=T(5)/sq2;

	// Iterations
	for (it=1; it<=maxIt; ++it)
	{
		// Check error
		sumUT = fabs(UT[0])+fabs(UT[1])+fabs(UT[2]);
		if (sumUT<=errTol) return it;
		
		// i=0, j=1 ,r=2 (p=3)
		h=L(0)-L(1);
		if (fabs(h)<zero) t=1.0;
		else
		{
			th=0.5*h/UT[0];
			t=1.0/(fabs(th)+sqrt(th*th+1.0));
			if (th<0.0) t=-t;
		}
		c  = 1.0/sqrt(1.0+t*t);
		s  = c*t;
		cc = c*c;
		ss = s*s;
		// Zeroes term UT[0]
		Temp  = cc*L(0) + 2.0*c*s*UT[0] + ss*L(1);
		L(1)  = ss*L(0) - 2.0*c*s*UT[0] + cc*L(1);
		L(0)  = Temp;
		UT[0] = 0.0;
		Temp  = c*UT[2] + s*UT[1];
		UT[1] = c*UT[1] - s*UT[2];
		UT[2] = Temp;
		
		// i=1, j=2 ,r=0 (p=4)
		h = L(1)-L(2);
		if (fabs(h)<zero) t=1.0;
		else
		{
			th=0.5*h/UT[1];
			t=1.0/(fabs(th)+sqrt(th*th+1.0));
			if (th<0.0) t=-t;
		}
		c  = 1.0/sqrt(1.0+t*t);
		s  = c*t;
		cc = c*c;
		ss = s*s;
		// Zeroes term UT[1]
		Temp  = cc*L(1) + 2.0*c*s*UT[1] + ss*L(2);
		L(2)  = ss*L(1) - 2.0*c*s*UT[1] + cc*L(2);
		L(1)  = Temp;
		UT[1] = 0.0;
		Temp  = c*UT[0] + s*UT[2];
		UT[2] = c*UT[2] - s*UT[0];
		UT[0] = Temp;

		// i=0, j=2 ,r=1 (p=5)
		h = L(0)-L(2);
		if (fabs(h)<zero) t=1.0;
		else
		{
			th=0.5*h/UT[2];
			t=1.0/(fabs(th)+sqrt(th*th+1.0));
			if (th<0.0) t=-t;
		}
		c  = 1.0/sqrt(1.0+t*t);
		s  = c*t;
		cc = c*c;
		ss = s*s;
		// Zeroes term UT[2]
		Temp  = cc*L(0) + 2.0*c*s*UT[2] + ss*L(2);
		L(2)  = ss*L(0) - 2.0*c*s*UT[2] + cc*L(2);
		L(0)  = Temp;
		UT[2] = 0.0;
		Temp  = c*UT[0] + s*UT[1];
		UT[1] = c*UT[1] - s*UT[0];
		UT[0] = Temp;
	}

	return -1; // Did not converge

}

}; // namespace Tensors

#endif // MECHSYS_TENSORS_JACOBIROT_H
