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

/* JacobiRot - Copyright (C) 2007 Dorival de Moraes Pedroso */

#ifndef MPM_JACOBIROT_H
#define MPM_JACOBIROT_H

// STL
#include <cmath>  // for sqrt ...
#include <cfloat> // for DBL_EPSILON

// MechSys
#include <mechsys/util/fatal.h>

// Local
#include <mechsys/mpm/defs.h>

namespace MPM {

inline int JacobiRot(STensor2 const & T, double L[3])
{
	const double maxIt  = 30;              // Max number of iterations
	const double tol    = DBL_EPSILON*100; // Erro: Tolerance
	const double zero   = DBL_EPSILON;     // Zero
	
	double UT[3];  // Values of the Upper Triangle part of symmetric matrix A
	double th;     // theta = (Aii-Ajj)/2Aij
	double c;      // Cossine
	double s;      // Sine
	double cc;     // Cossine squared
	double ss;     // Sine squared
	double t;      // Tangent
	double tmp;    // Auxiliar variable
	double sumUT;  // Sum of upper triangle of abs(A) that measures the error
	int    it;     // Iteration number
	double h;      // Difference L[i]-L[j]
	
	// Initialize eigenvalues which correnspond to the diagonal part of A
	L[0]=T(0); L[1]=T(1); L[2]=T(2);

	// Initialize Upper Triangle part of A matrix (array[3])
	UT[0]=T(3)/SQ2; UT[1]=T(4)/SQ2; UT[2]=T(5)/SQ2;

	// Iterations
	for (it=1; it<=maxIt; ++it)
	{
		// Check error
		sumUT = fabs(UT[0])+fabs(UT[1])+fabs(UT[2]);
		if (sumUT<=tol) return it;
		
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
		tmp   = cc*L[0] + 2.0*c*s*UT[0] + ss*L[1];
		L[1]  = ss*L[0] - 2.0*c*s*UT[0] + cc*L[1];
		L[0]  = tmp;
		UT[0] = 0.0;
		tmp   = c*UT[2] + s*UT[1];
		UT[1] = c*UT[1] - s*UT[2];
		UT[2] = tmp;
		
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
		tmp   = cc*L[1] + 2.0*c*s*UT[1] + ss*L[2];
		L[2]  = ss*L[1] - 2.0*c*s*UT[1] + cc*L[2];
		L[1]  = tmp;
		UT[1] = 0.0;
		tmp   = c*UT[0] + s*UT[2];
		UT[2] = c*UT[2] - s*UT[0];
		UT[0] = tmp;

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
		tmp   = cc*L[0] + 2.0*c*s*UT[2] + ss*L[2];
		L[2]  = ss*L[0] - 2.0*c*s*UT[2] + cc*L[2];
		L[0]  = tmp;
		UT[2] = 0.0;
		tmp   = c*UT[0] + s*UT[1];
		UT[1] = c*UT[1] - s*UT[0];
		UT[0] = tmp;
	}
	throw new Fatal("JacobiRot dit not converge for %d iterations",it);
}

inline int JacobiRot (STensor2 const & T, double V0[3], double V1[3], double V2[3], double L[3]) 
{
	const double maxIt  = 30;              // Max number of iterations
	const double tol    = DBL_EPSILON*100; // Erro: Tolerance
	const double zero   = DBL_EPSILON;     // Zero
	
	double UT[3];  // Values of the Upper Triangle part of symmetric matrix A
	double th;     // theta = (Aii-Ajj)/2Aij
	double c;      // Cossine
	double s;      // Sine
	double cc;     // Cossine squared
	double ss;     // Sine squared
	double t;      // Tangent
	double tmp;    // Auxiliar variable
	double TM[3];  // Auxiliar array
	double sumUT;  // Sum of upper triangle of abs(A) that measures the error
	int    it;     // Iteration number
	double h;      // Difference L[i]-L[j]
	
	// Initialize eigenvalues which correnspond to the diagonal part of A
	L[0]=T(0); L[1]=T(1); L[2]=T(2);

	// Initialize Upper Triangle part of A matrix (array[3])
	UT[0]=T(3)/SQ2; UT[1]=T(4)/SQ2; UT[2]=T(5)/SQ2;

	// Initialize eigenvectors
	V0[0]=1.0; V1[0]=0.0; V2[0]=0.0;
	V0[1]=0.0; V1[1]=1.0; V2[1]=0.0;
	V0[2]=0.0; V1[2]=0.0; V2[2]=1.0;

	// Iterations
	for (it=1; it<=maxIt; ++it)
	{
		// Check error
		sumUT = fabs(UT[0])+fabs(UT[1])+fabs(UT[2]);
		if (sumUT<=tol) return it;
		
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
		tmp   = cc*L[0] + 2.0*c*s*UT[0] + ss*L[1];
		L[1]  = ss*L[0] - 2.0*c*s*UT[0] + cc*L[1];
		L[0]  = tmp;
		UT[0] = 0.0;
		tmp   = c*UT[2] + s*UT[1];
		UT[1] = c*UT[1] - s*UT[2];
		UT[2] = tmp;
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
		tmp   = cc*L[1] + 2.0*c*s*UT[1] + ss*L[2];
		L[2]  = ss*L[1] - 2.0*c*s*UT[1] + cc*L[2];
		L[1]  = tmp;
		UT[1] = 0.0;
		tmp   = c*UT[0] + s*UT[2];
		UT[2] = c*UT[2] - s*UT[0];
		UT[0] = tmp;
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
		tmp   = cc*L[0] + 2.0*c*s*UT[2] + ss*L[2];
		L[2]  = ss*L[0] - 2.0*c*s*UT[2] + cc*L[2];
		L[0]  = tmp;
		UT[2] = 0.0;
		tmp   = c*UT[0] + s*UT[1];
		UT[1] = c*UT[1] - s*UT[0];
		UT[0] = tmp;
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
	throw new Fatal("JacobiRot dit not converge for %d iterations",it);
}

}; // namespace MPM

#endif // MPM_JACOBIROT_H
