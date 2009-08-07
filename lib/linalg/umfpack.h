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

#ifndef MECHSYS_UMFPACK_H
#define MECHSYS_UMFPACK_H

// UMFPACK
extern "C" {
#include <umfpack.h>
}

// MechSys
#include "linalg/matrix.h"
#include "linalg/vector.h"
#include "linalg/sparse_matrix.h"
#include "util/fatal.h"
#include "util/util.h"

#ifdef DO_DEBUG
  using std::cout;
  using std::endl;
#endif

/** \namespace UMFPACK Unsymmetric multifrontal sparse LU factorization package.
  %UMFPACK is a set of routines for solving unsymmetric sparse linear systems, Ax=b, using the Unsymmetric MultiFrontal method.
  See <a href="http://www.cise.ufl.edu/research/sparse/umfpack/">Tim Davis' web page.</a>
  
  Examples:
   - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tumfpack.cpp?view=markup">tumfpack.cpp Test UMFPACK solver</a>
 */
namespace UMFPACK
{

/** %UMFPACK internal error message.
 * \param info %UMFPACK internal error code
 */
String ErrorMsg(int info)
{
	switch (info)
	{
		case UMFPACK_ERROR_out_of_memory           : return String(_("UMFPACK_ERROR_out_of_memory (-1)"));
		case UMFPACK_ERROR_invalid_Numeric_object  : return String(_("UMFPACK_ERROR_invalid_Numeric_object (-3)"));
		case UMFPACK_ERROR_invalid_Symbolic_object : return String(_("UMFPACK_ERROR_invalid_Symbolic_object (-4)"));
		case UMFPACK_ERROR_argument_missing        : return String(_("UMFPACK_ERROR_argument_missing (-5)"));
		case UMFPACK_ERROR_n_nonpositive           : return String(_("UMFPACK_ERROR_n_nonpositive (-6)"));
		case UMFPACK_ERROR_invalid_matrix          : return String(_("UMFPACK_ERROR_invalid_matrix (-8)"));
		case UMFPACK_ERROR_different_pattern       : return String(_("UMFPACK_ERROR_different_pattern (-11)"));
		case UMFPACK_ERROR_invalid_system          : return String(_("UMFPACK_ERROR_invalid_system (-13)"));
		case UMFPACK_ERROR_invalid_permutation     : return String(_("UMFPACK_ERROR_invalid_permutation (-15)"));
		case UMFPACK_ERROR_internal_error          : return String(_("UMFPACK_ERROR_internal_error (-911)"));
		case UMFPACK_ERROR_file_IO                 : return String(_("UMFPACK_ERROR_file_IO (-17)"));
		default                                    : { String r; r.Printf(_("Unknown UMFPACK_ERROR (%d)"),info); return r; };
	}
}

/** Solves \f$ {X} \leftarrow [A]^{-1}{B} \f$.
 * \param A %Sparse matrix
 * \param B Right-hand side vector
 * \param X Result
 */
inline void Solve(Sparse::Matrix<double,int> const & A, LinAlg::Vector<double> const & B, LinAlg::Vector<double> & X)
{
#ifndef DNDEBUG
	if (A.Rows()!=A.Cols()) throw new Fatal(_("UMFPACK::Solve: The matrix A (%d x %d) must be squared."),A.Rows(),A.Cols());
	if (B.Size()!=A.Cols()) throw new Fatal(_("UMFPACK::Solve: The vector B (%d x 1) must have a size equal to the number of columns of matrix A (%d)"),B.Size(),A.Cols());
#endif
	double *null = (double *)NULL;
	void *symbolic, *numeric;
	int info = 0;
	info = umfpack_di_symbolic      (A.Rows(), A.Rows(), A.GetApPtr(), A.GetAiPtr(), A.GetAxPtr(), &symbolic, null, null);              if (info<0) throw new Fatal(_("UMFPACK::Solve: umfpack_dl_symbolic failed. %s"),ErrorMsg(info).CStr());
	info = umfpack_di_numeric       (A.GetApPtr(), A.GetAiPtr(), A.GetAxPtr(), symbolic, &numeric, null, null);                         if (info<0) throw new Fatal(_("UMFPACK::Solve: umfpack_dl_numeric failed. %s"),ErrorMsg(info).CStr());
	       umfpack_di_free_symbolic (&symbolic);
	info = umfpack_di_solve         (UMFPACK_A, A.GetApPtr(), A.GetAiPtr(), A.GetAxPtr(), X.GetPtr(), B.GetPtr(), numeric, null, null); if (info<0) throw new Fatal(_("UMFPACK::Solve: umfpack_dl_solve failed. %s"),ErrorMsg(info).CStr());
	       umfpack_di_free_numeric  (&numeric);
}

}; // namespace UMFPACK

#endif // MECHSYS_UMFPACK_H
