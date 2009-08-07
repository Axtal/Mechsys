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

// Std Lib
#include <stdlib.h>

// GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

// MechSys
#include "mechsys.h"
#include "util/exception.h"

/* Computation of the integral,

	I = int (dx dy dz)/(2pi)^3  1/(1-cos(x)cos(y)cos(z))

	over (-pi,-pi,-pi) to (+pi, +pi, +pi).  The exact answer
	is Gamma(1/4)^4/(4 pi^3).  This example is taken from
	C.Itzykson, J.M.Drouffe, "Statistical Field Theory -
	Volume 1", Section 1.1, p21, which cites the original
	paper M.L.Glasser, I.J.Zucker, Proc.Natl.Acad.Sci.USA 74
	1800 (1977)

	For simplicity we compute the integral over the region 
	(0,0,0) -> (pi,pi,pi) and multiply by 8
*/

double exact = 1.3932039296856768591842462603255;

double g (double * k, size_t dim, void * params)
{
	double A = 1.0 / (M_PI * M_PI * M_PI);
	return A / (1.0 - cos (k[0]) * cos (k[1]) * cos (k[2]));
}

void display_results (char const * title, double result, double error)
{
	printf ("%s ==================\n", title);
	printf ("result = % .6f\n", result);
	printf ("sigma  = % .6f\n", error);
	printf ("exact  = % .6f\n", exact);
	printf ("error  = % .6f = %.1g sigma\n", result - exact, fabs (result - exact) / error);
}

int main(int argc, char **argv) try
{
	double res, err;
	double xl[3] = { 0, 0, 0 };
	double xu[3] = { M_PI, M_PI, M_PI };

	const gsl_rng_type *T;
	gsl_rng *r;

	gsl_monte_function G = { &g, 3, 0 };

	size_t calls = 500000;

	gsl_rng_env_setup ();

	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	// Plain
	{
		gsl_monte_plain_state *s = gsl_monte_plain_alloc (3);
		gsl_monte_plain_integrate (&G, xl, xu, 3, calls, r, s, &res, &err);
		gsl_monte_plain_free (s);
		display_results ("plain", res, err);
	}

	// Miser
	{
		gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
		gsl_monte_miser_integrate (&G, xl, xu, 3, calls, r, s, &res, &err);
		gsl_monte_miser_free (s);
		display_results ("miser", res, err);
	}

	// Vegas
	{
		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
		gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s, &res, &err);
		display_results ("vegas warm-up", res, err);
		printf ("converging...\n");
		do
		{
			gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
			&res, &err);
			printf ("result = % .6f sigma = % .6f "
			"chisq/dof = %.1f\n", res, err, s->chisq);
		}
		while (fabs (s->chisq - 1.0) > 0.5);
		display_results ("vegas final", res, err);
		gsl_monte_vegas_free (s);
	}

	// Clean up
	gsl_rng_free (r);
	return 0;
}
catch (Exception * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout<<"Fatal: "<<m<<std::endl;  exit(1); }
catch (...) { std::cout << "Some exception (...) ocurred\n"; } 
