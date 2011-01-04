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

// STL
#include <iostream>
#include <cmath>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/linalg/jacobirot.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/util/stopwatch.h>

using std::cout;
using std::endl;
using Util::SQ2;

int main(int argc, char **argv) try
{
    // jacobirot
    printf("\n%s . . . using jacobi rotation%s\n\n",TERM_BLACK_WHITE,TERM_RST);
    Mat_t A(3,3), Q;
    Vec_t L;
    A = 1., 2., 3.,
        2., 3., 2.,
        3., 2., 2.;
    int it = JacobiRot (A, Q, L);
    printf("number of iterations = %d\n", it);
    cout << "A =\n" << PrintMatrix (A);
    cout << "Q =\n" << PrintMatrix (Q);
    cout << "L ="   << PrintVector (L);

    // GSL
    printf("\n%s . . . using GSL%s\n\n",TERM_BLACK_WHITE,TERM_RST);
    A = 1., 2., 3.,
        2., 3., 2.,
        3., 2., 2.;
    Vec_t l;
    Eig (A, l);
    cout << "A =\n" << PrintMatrix (A);
    cout << "l ="   << PrintVector (l);

    // GSL + eigenvectors
    printf("\n%s . . . using GSL with eigenvectors%s\n\n",TERM_BLACK_WHITE,TERM_RST);
    A = 1., 2., 3.,
        2., 3., 2.,
        3., 2., 2.;
    Vec_t ll;
    Mat_t qq;
    Eig (A, ll, &qq);
    cout << "A  =\n" << PrintMatrix (A);
    cout << "qq =\n" << PrintMatrix (qq);
    cout << "ll ="   << PrintVector (ll);

    // error
    double err0 = fabs(L(0)-l(1));
    double err1 = fabs(L(1)-l(0));
    double err2 = fabs(L(2)-l(2));
    double max_error = 0.0;
    if (err0>max_error) max_error = err0;
    if (err1>max_error) max_error = err1;
    if (err2>max_error) max_error = err2;

    // tensor
    printf("\n%s . . . . . . . . . . . . tensor ..%s\n\n",TERM_BLACK_WHITE,TERM_RST);
    Vec_t sig(6);
    Mat3_t msig, r;
    Vec3_t v, n, t, tt;
    sig = -1.5, -2.0, -0.5,  -0.1*SQ2, -0.5*SQ2, -0.8*SQ2;
    n   = 1.0, 1.2, -0.8;
    Ten2Mat   (sig, msig);
    JacobiRot (msig, r, v);
    //t = product (msig, n); // this does not work, since the (msig) matrix was modified by JacobiRot
    Mat3_t msig_new;
    Ten2Mat (sig, msig_new);
    t = product (msig_new, n);
    Mult (sig, n, tt);
    cout << "sig =\n" << PrintMatrix (msig);
    cout << "r  =\n"  << PrintMatrix (r);
    cout << "v  ="    << PrintVector (v);
    cout << "t  ="    << PrintVector (t);
    cout << "tt ="    << PrintVector (tt);
    double err3 = Norm(t-tt);
    if (err3>max_error) max_error = err3;

    double tol = 1.0e-14;
    printf("\n Max error = %s%g%s\n", (max_error>tol ? TERM_RED : TERM_GREEN), max_error, TERM_RST);

    if (max_error>tol) return 1;
    return 0;
}
MECHSYS_CATCH
