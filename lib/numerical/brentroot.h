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

#ifndef MECHSYS_BRENTROOT_H
#define MECHSYS_BRENTROOT_H

// STL
#include <cmath>
#include <cfloat> // for DBL_EPSILON

// MechSys
#include <mechsys/util/fatal.h>

namespace Numerical
{

/** Find a root for F(x)=0 inside [A,B] using Brent's method.
Examples:
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/numerical/tst/tbrentroot.cpp?view=markup"> ttbrentroot.cpp Test Brent's method</a>
*/
template<typename Instance>
class BrentRoot
{
public:
	// Constants
	static int BRENT_MAXIT; ///< Max iterations for Brent's method

	// Typedefs
	typedef double (Instance::*pFun) (double x, void const * P1, void const * P2, void const * P3) const; ///< Callback function (P# => extra parameters)

	/** Constructor. */
	BrentRoot (Instance const * p2Inst, pFun p2Fun);

	// Methods
	double Solve  (double A, double B, void const * P1, void const * P2, void const * P3) const; ///< Find root (P# => extra parameters)
	void   SetTol (double Tol) { _tol = Tol; } ///< Set tolerance

private:
	Instance const * _p2inst; ///< Pointer to object
	pFun             _p2fun;  ///< Pointer to method
	double           _tol;    ///< error tolerance

}; // class BrentRoot

template<typename Instance>
int BrentRoot<Instance>::BRENT_MAXIT = 100;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


template<typename Instance>
inline BrentRoot<Instance>::BrentRoot(Instance const * p2Inst, pFun p2Fun)
	: _p2inst (p2Inst),
	  _p2fun  (p2Fun),
	  _tol    (sqrt(DBL_EPSILON))
{
}

template<typename Instance>
inline double BrentRoot<Instance>::Solve(double A, double B, void const * P1, void const * P2, void const * P3) const
{
	/* Based on ZEROIN C math library: http://www.netlib.org/c/
	
	   Algorithm:
	   G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
	   computations. M., Mir, 1980, p.180 of the Russian edition
	
	   The function makes use of the bissection procedure combined with
	   the linear or quadric inverse interpolation.
	   At every step program operates on three abscissae: a, b, and c:
	    a: the last but one approximation
	    b: the last and the best approximation to the root
	    c: the last but one or even earlier approximation than a such that:
	        1) |f(b)| <= |f(c)|
	        2) f(b) and f(c) have opposite signs, i.e. b and c confine the root

	   At every step this algorithm selects one of the two new approximations, the former being
	   obtained by the bissection procedure and the latter resulting in the interpolation.
	   If a,b and c are all different, the quadric interpolation is utilized, otherwise the linear one.
	   If the latter (i.e. obtained by the interpolation) point is reasonable (i.e. lies within
	   the current interval [b,c] not being too close to the boundaries) it is accepted.
	   The bissection result is used in the other case. Therefore, the range of uncertainty is
	   ensured to be reduced at least by the factor 1.6.
	*/

	double a  = A; // the last but one approximation
	double b  = B; // the last and the best approximation to the root
	double c  = a; // the last but one or even earlier approximation than a that
	double fa = (_p2inst->*_p2fun)(a,P1,P2,P3);
	double fb = (_p2inst->*_p2fun)(b,P1,P2,P3);
	double fc = fa;

	// Check input
	if ((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0.0))
		throw new Fatal(_("BrentRoot::Solve: Root must be bracketed."));

	// Solve
	for(int i=0; i<BRENT_MAXIT; ++i)
	{
		double prev_step = b-a; // Distance from the last but one to the last approximation
		double tol_act;         // Actual tolerance
		double p,q;             // Interpolation step is calculated in the form p/q; division operations is delayed until the last moment
		double new_step;        // Step at this iteration

		// Swap data for b to be the best approximation
		if (fabs(fc)<fabs(fb))
		{
			 a =  b;     b =  c;     c =  a;
			fa = fb;    fb = fc;    fc = fa;
		}
		tol_act  = 2.0*DBL_EPSILON*fabs(b) + _tol/2.0;
		new_step = (c-b)/2.0;

		// Check for convergence
		if (fabs(new_step)<=tol_act || fb==0.0)
			return b; // Acceptable approx. is found

		// Decide if the interpolation can be tried
		if (fabs(prev_step)>=tol_act && fabs(fa)>fabs(fb)) // If prev_step was large enough and was in true direction
		{
			// Interpolation may be tried
			double t1,cb,t2;
			cb = c-b;
			if(a==c) // If we have only two distinct points, linear interpolation can only be applied
			{
				t1 = fb/fa;
				p  = cb*t1;
				q  = 1.0-t1;
			}
			else // Quadric inverse interpolation
			{
				q = fa/fc;     t1 = fb/fc;     t2 = fb/fa;
				p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
				q = (q-1.0) * (t1-1.0) * (t2-1.0);
			}

			// p was calculated with the oposite sign; make p positive and assign possible minus to q
			if (p>0.0) q = -q;
			else       p = -p;

			// If b+p/q falls in [b,c] and isn't too large
			if (p<(0.75*cb*q-fabs(tol_act*q)/2.0) && p<fabs(prev_step*q/2.0))
				new_step = p/q;// it is accepted

			// If p/q is too large then the bissection procedure can reduce [b,c] range to more extent
		}

		// Adjust the step to be not less than tolerance
		if (fabs(new_step)<tol_act)
		{
			if (new_step>0.0) new_step =  tol_act;
			else              new_step = -tol_act;
		}

		// Save the previous approx
		 a =  b;
		fa = fb;

		// Do step to a new approxim.
		b += new_step;
		fb = (_p2inst->*_p2fun)(b,P1,P2,P3);

		// Adjust c for it to have a sign opposite to that of b
		if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0))
		{
			 c =  a;
			fc = fa;
		}
	}

	// Did not converge
	throw new Fatal(_("BrentRoot::Solve: did not converge for %d iterations"),BRENT_MAXIT);
	return 0;

}


}; // namespace Numerical

#endif // MECHSYS_BRENTROOT_H
