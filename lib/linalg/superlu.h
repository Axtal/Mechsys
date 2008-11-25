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

#ifndef MECHSYS_SUPERLU_H
#define MECHSYS_SUPERLU_H

// SuperLU
#include <slu_ddefs.h>

// MechSys
#include "linalg/matrix.h"
#include "linalg/vector.h"
#include "linalg/sparse_matrix.h"
#include "util/exception.h"
#include "util/util.h"

#ifdef DO_DEBUG
  using std::cout;
  using std::endl;
#endif

/** \namespace SuperLU Nonsymmetric sparse solver.
  %SuperLU is a general purpose library for the direct solution of large, sparse, nonsymmetric systems of linear equations on high performance machines.
  See <a href="http://crd.lbl.gov/~xiaoye/SuperLU/">Xiaoye's web page.</a>
  
  Examples:
   - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tsuperlu.cpp?view=markup">tsuperlu.cpp  Test SuperLU linear solver</a>
 */

namespace SuperLU
{

/** Solves \f$ {X} \leftarrow [A]^{-1}{X} \f$.
 * \param A %Sparse matrix
 * \param X Right-hand side vector and solution vector
 */
inline void Solve(Sparse::Matrix<double,int> const & A, LinAlg::Vector<double> & X)
{
#ifndef DNDEBUG
	if (A.Rows()!=A.Cols()) throw new Fatal(_("SuperLU::Solve: The matrix A (%d x %d) must be squared."),A.Rows(),A.Cols());
	if (X.Size()!=A.Cols()) throw new Fatal(_("SuperLU::Solve: The vector X (%d x 1) must have a size equal to the number of columns of matrix A (%d)"),X.Size(),A.Cols());
#endif

	// Size
	int m  = A.Rows();
	int n  = A.Cols();
	int nz = A.nZ();

    // Create matrix a in the format expected by SuperLU
	NCformat    ncf = { nz, const_cast<double*>(A.GetAxPtr()), const_cast<int*>(A.GetAiPtr()), const_cast<int*>(A.GetApPtr()) };
	SuperMatrix a   = { /*Stype=column-wise/no-supernode*/SLU_NC, /*Dtype=double*/SLU_D, /*Mtype=general*/SLU_GE, m, n, &ncf };

    // Create right-hand side vector b
	DNformat    dnf = { m, X.GetPtr() };
	SuperMatrix b   = { /*Stype=dense-storage*/SLU_DN, /*Dtype=double*/SLU_D, /*Mtype=general*/SLU_GE, m, 1, &dnf };

	// Workspaces
    int * perm_r = new int [m]; // row permutations from partial pivoting
    int * perm_c = new int [n]; // column permutation vector

    // Set the default input options
    superlu_options_t options;
    set_default_options(&options);
    options.ColPerm = COLAMD; //NATURAL;

    // Initialize the statistics variables
    SuperLUStat_t stat;
    StatInit(&stat);
    
	// Solve
	int info;
    SuperMatrix l, u;
    dgssv(&options, &a, perm_c, perm_r, &l, &u, &b, &stat, &info);
    
#ifdef DO_DEBUGx
    dPrint_CompCol_Matrix   ("a", &a);
    dPrint_CompCol_Matrix   ("u", &u);
    dPrint_SuperNode_Matrix ("l", &l);
#endif

	// Cleanup
    Destroy_SuperNode_Matrix (&l);
    Destroy_CompCol_Matrix   (&u);
    StatFree                 (&stat);
	delete [] perm_r;
	delete [] perm_c;

}

}; // namespace SuperLU

#endif // MECHSYS_SUPERLU_H
