/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo, Fernando Alonso                   *
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

/* Numerical evaluation of eigenvalues and eigenvectors for NxN SYMMETRIC
 *  matrices using the Jacobi-rotation Method.
 *
 * 2005-09-27 (C): Dorival M. Pedroso - Completely rewritten (3x3 => NxN); based on C numerical recipes
 * 2005-05-10 (C): Dorival M. Pedroso - C++
 * 2004-08-15 (C): Dorival M. Pedroso - Fortran 95
 * 2002-04-15 (C): Dorival M. Pedroso - Fortran 77
 *
 */

#ifndef MECHSYS_LINALG_JACOBIROT_H
#define MECHSYS_LINALG_JACOBIROT_H

// STL
#include <cmath>

// MechSys
#include "linalg/matrix.h"
#include "linalg/vector.h"
#include "util/exception.h"

namespace LinAlg
{

/** \namespace LinAlg::Sym Linear Algebra - Symmetric matrices. */

namespace Sym
{

/** Jacobi Transformation of a Symmetric Matrix<double>.
 * The Jacobi method consists of a sequence of orthogonal similarity transformations.
 * Each transformation (a Jacobi rotation) is just a plane rotation designed to annihilate one of the off-diagonal matrix elements.
 * Successive transformations undo previously set zeros, but the off-diagonal elements nevertheless get smaller and smaller.
 * Accumulating the product of the transformations as you go gives the matrix of eigenvectors (Q), while the elements of the final
 * diagonal matrix (A) are the eigenvalues.
 * The Jacobi method is absolutely foolproof for all real symmetric matrices.
 * For matrices of order greater than about 10, say, the algorithm is slower, by a significant constant factor, than the QR method.
 *
 * \param A In/Out: is the matrix we seek for the eigenvalues (SYMMETRIC and square, i.e. Rows=Cols)
 * \param Q Out: is a matrix which columns are the eigenvectors
 * \param v Out: is a vector with the eigenvalues
 * \return The number of iterations
 */
inline int JacobiRot(Matrix<double> & A, Matrix<double> & Q, Vector<double> & v)
{
	int N = A.Rows();
	if (N<2)         throw new Fatal("Sym::JacobiRot: A(%d,%d) matrix must be at least 2x2 on size",A.Rows(),A.Cols());
	if (A.Cols()!=N) throw new Fatal("Sym::JacobiRot: A(%d,%d) matrix must be squared",A.Rows(),A.Cols());
	if (Q.Rows()!=N) throw new Fatal("Sym::JacobiRot: Q(%d,%d) must have the same dimensions than A(%d,%d)",Q.Rows(),Q.Cols(),A.Rows(),A.Cols());
	if (Q.Cols()!=N) throw new Fatal("Sym::JacobiRot: Q(%d,%d) must have the same dimensions than A(%d,%d)",Q.Rows(),Q.Cols(),A.Rows(),A.Cols());

	const int    maxIt  = 30;      // Max number of iterations
	const double errTol = 1.0e-15; // Tolerance
	const double Zero   = 1.0e-15; // Tolerance

	int    j,p,q;
	double theta,tau,t,sm,s,h,g,c;
	double * b = new double [N];
	double * z = new double [N];

	for (p=0; p<N; p++) // Initialize Q to the identity matrix.
	{
		for (q=0; q<N; q++) Q(p,q)=0.0;
		Q(p,p)=1.0;
	}
	for (p=0; p<N; p++)
	{
		b[p]=v(p)=A(p,p); // Initialize b and v to the diagonal of A
		z[p]=0.0;         // This vector will accumulate terms of the form tapq as in equation (11.1.14).
	}

	for (int it=0; it<maxIt; it++)
	{
		sm=0.0;
		for (p=0; p<N-1; p++) // Sum off-diagonal elements.
		{
			for (q=p+1; q<N; q++)
			sm += fabs(A(p,q));
		}
		if (sm<errTol) // The normal return
		{
			delete [] b;
			delete [] z;
			return it+1;
		}
		for (p=0; p<N-1; p++)
		{
			for (q=p+1; q<N; q++)
			{
				h=v(q)-v(p);
				if (fabs(h)<=Zero) t=1.0;
				else
				{
					theta=0.5*h/(A(p,q));
					t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
					if (theta < 0.0) t=-t;
				}
				c=1.0/sqrt(1.0+t*t);
				s=t*c;
				tau=s/(1.0+c);
				h=t*A(p,q);
				z[p] -= h;
				z[q] += h;
				v(p) -= h;
				v(q) += h;
				A(p,q)=0.0;
				for (j=0; j<p; j++)  // Case of rotations 0 <= j < p.
				{ 
					g = A(j,p);
					h = A(j,q);
					A(j,p) = g - s*(h+g*tau);
					A(j,q) = h + s*(g-h*tau);
				}
				for (j=p+1; j<q; j++) // Case of rotations p < j < q.
				{ 
					g = A(p,j);
					h = A(j,q);
					A(p,j) = g - s*(h+g*tau);
					A(j,q) = h + s*(g-h*tau);
				}
				for (j=q+1; j<N; j++) //Case of rotations q < j < N.
				{ 
					g = A(p,j);
					h = A(q,j);
					A(p,j) = g - s*(h+g*tau);
					A(q,j) = h + s*(g-h*tau);
				}
				for (j=0; j<N; j++) // Q matrix
				{
					g = Q(j,p);
					h = Q(j,q);
					Q(j,p) = g - s*(h+g*tau);
					Q(j,q) = h + s*(g-h*tau);
				}
			}
		}
		for (p=0; p<N; p++)
		{
			b[p] += z[p];
			z[p]  = 0.0;   // reinitialize z.
			v(p)  = b[p];  // update v with the sum of tapq,
		}
	}

	delete [] b;
	delete [] z;
	std::cout << "ERROR: Too many iterations in routine jacobi\n";
	return maxIt+1;

}

}; // namespace Sym

}; // namespace LinAlg

#endif // MECHSYS_LINALG_JACOBIROT_H
