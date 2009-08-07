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

#ifndef MECHSYS_LINALG_LAWRAP_H
#define MECHSYS_LINALG_LAWRAP_H

// MechSys
#include "linalg/matrix.h"
#include "linalg/vector.h"
#include "util/fatal.h"
#include "util/util.h"

extern "C"
{
#include "linalg/cblas.h"
	// BLAS
	void dscal_(int const *N, double const *alpha, double *X, int const *incX);

	void dcopy_(int const *N, double const *X, int const *incX, double *Y, int const *incY);

	void daxpy_(int const *N, double const *alpha, double const *X, int const *incX, double *Y, int const *incY);

	double ddot_(int const *N, double const *X, int const *incX, double const *Y, int const *incY);

	void dsymv_(char const *Uplo, int    const *N, double const *alpha, double const *A,
	            int  const *lda , double const *X, int    const *incX , double const *beta,
	            double     *Y   , int    const *incY);

	void dgemv_(char   const *TransA , int    const *M   , int    const *N  , double const *alpha ,
	            double const *A      , int    const *lda , double const *X  , int    const *incX  ,
	            double const *beta   , double       *Y   , int    const *incY);

	void dgemm_(char   const *TransA  , char   const *TransB , int    const *M    , int    const *N   ,
	            int    const *K       , double const *alpha  , double const *A    , int    const *lda ,
	            double const *B       , int    const *ldb    , double const *beta , double       *C   , int const *ldc);

	void dger_(int const *M   , int    const *N, double const *alpha, double const *X,
	           int const *incX, double const *Y, int    const *incY , double       *A, int const *lda);

	// DGESV
	void dgesv_(int *Np, int *NRHSp,
	            double *A, int *LDAp, int *IPIVp,
	            double *B, int *LDBp, int *INFOp);

	// DGETRF - compute an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges
	void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);

	// DGETRI - compute the inverse of a matrix using the LU factorization computed by DGETRF
	void dgetri_(const int* n, double* a, const int* lda, int* ipiv, double* work, const int* lwork, int* info);

	// DGTSV and EIGENPROBLEMS
	/* not working with Ubuntu's LAPACK package
	double dlamch_ (char *CMACHp);
	void   dgtsv_  (int *Np, int *NRHSp, double *DL, double *D, double *DU, double *B, int *LDBp, int *INFOp);
	void   dsyevr_ (char *JOBZp, char *RANGEp, char *UPLOp, int *Np, double *A, int *LDAp, double *VLp, double *VUp,
	                int *ILp, int *IUp, double *ABSTOLp, int *Mp, double *W, double *Z, int *LDZp, int *ISUPPZ,
	                double *WORK, int *LWORKp, int *IWORK, int *LIWORKp, int *INFOp);
	*/
}

/** \namespace LinAlg Linear Algebra.
  Examples
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/teigen1.cpp?view=markup">     teigen1.cpp Test eigenvalues # 1</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/teigen.cpp?view=markup">      teigen.cpp Test eigenvalues</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/texpressions.cpp?view=markup">texpressions.cpp Test expressions</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tjacob.cpp?view=markup">      tjacob.cpp Test Jacoby rotation</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tlinsolver.cpp?view=markup">  tlinsolver.cpp Test Linear solver</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tmv.cpp?view=markup">         tmv.cpp Test matrix-vector multiplication</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tsparse.cpp?view=markup">     tsparse.cpp Test sparse matrices</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tsuperlu.cpp?view=markup">    tsuperlu.cpp  Test SuperLU linear solver</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tsv.cpp?view=markup">         tsv.cpp Test LAPACK solver</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tsymmv.cpp?view=markup">      tsymmv.cpp Test symmetric matrix-vector multiplication</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tumfpack.cpp?view=markup">    tumfpack.cpp Test UMFPACK solver</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tvv.cpp?view=markup">         tvv.cpp Test vector-vector operations</a>
 */

namespace LinAlg
{


/** Linear Solver Type. */
enum LinSol_T
{
	LAPACK_T   = 1, ///< DENSE
	UMFPACK_T  = 2, ///< SPARSE
	SuperLU_T  = 3, ///< SPARSE
	SuperLUd_T = 4  ///< SPARSE (distributed/parallel)
};

/** Returns the name and description of a linear solver. */
void LinSolNameDesc(LinSol_T LinSol, String & Name, String & Desc)
{
	     if (LinSol==LinAlg::LAPACK_T  ) { Name = "LAPACK";   Desc = "Linear Algebra PACKage";                                }
	else if (LinSol==LinAlg::UMFPACK_T ) { Name = "UMFPACK";  Desc = "Unsymmetric Multifrontal PACKage";                      }
	else if (LinSol==LinAlg::SuperLU_T ) { Name = "SuperLU";  Desc = "Super LU decomposition";                                }
	else if (LinSol==LinAlg::SuperLUd_T) { Name = "SuperLUd"; Desc = "Super LU decomposition, parallel (distributed-memory)"; }
	else                                 { Name = "Unknown";  Desc = "Unknown linear solver";                                 }
}


// ############################################################################ BLAS


/* LEVEL 1 (vector-vector) */

/// scale:  \f$ \{X\} \gets \alpha \{X\} \f$
inline void Scal(double const a, Vector<double> & X)
{
	int N    = X.Size();
	int incX = 1;
	dscal_(&N,&a,X.GetPtr(),&incX);
}

/// copy:  \f$ \{Y\} \gets \{X\} \f$
inline void Copy(Vector<double> const & X, Vector<double> & Y)
{
#ifndef DNDEBUG
	if (Y.Size()!=X.Size()) throw new Fatal(_("LinAlg::Copy: The size (%d) of vector Y must be equal to the size (%d) of vector X"),Y.Size(),X.Size());
#endif
	int N    = X.Size();
	int incX = 1;
	int incY = 1;
	dcopy_(&N,X.GetPtr(),&incX,Y.GetPtr(),&incY);
}

/// a*X plus Y:  \f$ \{Y\} \gets \alpha\{X\} + \{Y\} \f$
inline void Axpy(double const a, Vector<double> const & X, Vector<double> & Y)
{
#ifndef DNDEBUG
	if (Y.Size()!=X.Size()) throw new Fatal(_("LinAlg::Axpy: The size (%d) of vector Y must be equal to the size (%d) of vector X"),Y.Size(),X.Size());
#endif
	int N    = X.Size();
	int incX = 1;
	int incY = 1;
	daxpy_(&N,&a,X.GetPtr(),&incX,Y.GetPtr(),&incY);
}

/// internal product:  \f$ s \gets \{X\} \bullet \{Y\} \f$
inline double Dot(Vector<double> const & X, Vector<double> const & Y)
{
#ifndef DNDEBUG
	if (Y.Size()!=X.Size()) throw new Fatal(_("LinAlg::Axpy: The size (%d) of vector Y must be equal to the size (%d) of vector X"),Y.Size(),X.Size());
#endif
	int N    = X.Size();
	int incX = 1;
	int incY = 1;
	return ddot_(&N,X.GetPtr(),&incX,Y.GetPtr(),&incY);
}

/// Norm:  \f$ Norm\{X\} \gets sqrt(\{X\}^T\{X\}) \f$
inline double Norm(Vector<double> const & X)
{
	return sqrt(LinAlg::Dot(X,X));
}

/// copy and scale:  \f$ \{Y\} \gets \alpha\{X\} \f$ _IMPORTANT_: this function is not efficient
inline void CopyScal(double const a, Vector<double> const & X, Vector<double> & Y)
{
#ifndef DNDEBUG
	if (Y.Size()!=X.Size()) throw new Fatal(_("LinAlg::CopyScal: The size (%d) of vector Y must be equal to the size (%d) of vector X"),Y.Size(),X.Size());
#endif
	Copy(X,Y); // Y <- X;
	Scal(a,Y); // Y <- aY = aX
}

/// Add scaled vectors:  \f$ \{Z\} \gets \alpha\{X\} + \beta\{Y\} \f$ _IMPORTANT_: this function is not efficient
inline void AddScaled(double const a, Vector<double> const & X, double const b, Vector<double> const & Y, Vector<double> & Z)
{
#ifndef DNDEBUG
	if (Y.Size()!=X.Size()) throw new Fatal(_("LinAlg::AddScaled: The size (%d) of vector Y must be equal to the size (%d) of vector X"),Y.Size(),X.Size());
#endif
	Z.SetValues(0.0); // Z <- 0.0
	Axpy(a,X,Z);      // Z <- a*X + 0.0
	Axpy(b,Y,Z);      // Z <- b*Y + a*X
}


/* LEVEL 2 (matrix-vector) */

/// SYmmetric matrix * vector:  \f$ \{Y\} \gets \alpha[A]\{X\} + \beta\{Y\} \f$
inline void Symv(double const a, Matrix<double> const & A, Vector<double> const & X, double const b, Vector<double> & Y)
{
	/*   {Y} <- a*[A]*{X} + b*{Y} 
	                _                 _
	     / 0 \     |  00 01 02 ... 0N  | / 0 \     / 0 \
	     | 1 |     |  01 11 12 ... 1N  | | 1 |     | 1 |
	     | 2 | = a*|  02 12 22 ... 2N  |*| 2 | + b*| 2 |
	     |...|     |       ...         | |...|     |...|
	     \ N /     |_ 0N 1N 2N ... NN _| \ N /     \ N /
	
	    This function preserves [A] matrix after multplication
	*/
#ifndef DNDEBUG
	if (A.Cols()!=A.Rows()) throw new Fatal(_("LinAlg::Symv: The number of columns (%d) of matrix A must be equal to the number of rows (%d) of matrix A"),A.Cols(),A.Rows());
	if (X.Size()!=A.Rows()) throw new Fatal(_("LinAlg::Symv: The size (%d) of vector X must be equal to the number of rows (%d) of matrix A"),X.Size(),A.Rows());
	if (Y.Size()!=A.Rows()) throw new Fatal(_("LinAlg::Symv: The size (%d) of vector Y must be equal to the number of rows (%d) of matrix A"),Y.Size(),A.Rows());
#endif
	int  N    = A.Rows();
	char Uplo = 'U';
	int  incX = 1;
	int  incY = 1;
	dsymv_(&Uplo,         // CBLAS_UPLO
	       &N,            // N
	       &a,            // Alpha
	       A.GetPtr(),    // const double * A
	       &N,            // LDA
	       X.GetPtr(),    // const double * X
	       &incX,         // incX
	       &b,            // Beta
	       Y.GetPtr(),    // double * Y
	       &incY);        // incY
}

/// GEneral matrix * vector: \f$ \{Y\} \gets \alpha[A]\{X\} + \beta\{Y\}   \f$
inline void Gemv(double const a, Matrix<double> const & A, Vector<double> const & X, double const b, Vector<double> & Y)
{
	/*   {Y} <- a*[A]*{X} + b*{Y}
	                _                       _  / 0 \
	     / 0 \     |  00 01 02 03 04 ... 0N  | | 1 |     / 0 \
	     | 1 |     |  10 11 12 13 14 ... 1N  | | 2 |     | 1 |
	     | 2 | = a*|  20 21 22 23 24 ... 2N  |*| 3 | + b*| 2 |
	     |...|     |          ...            | | 4 |     |...|
	     \ M /     |_ M0 M1 M2 M3 M4 ... MN _| |...|     \ M /
	                                           \ N /
	
	    This function preserves [A] matrix after multplication
	*/
#ifndef DNDEBUG
	if (X.Size()!=A.Cols()) throw new Fatal(_("LinAlg::Gemv: The size (%d) of vector X must be equal to the number of columns (%d) of matrix A"),X.Size(),A.Cols());
	if (Y.Size()!=A.Rows()) throw new Fatal(_("LinAlg::Gemv: The size (%d) of vector Y must be equal to the number of rows (%d) of matrix A"),Y.Size(),A.Rows());
#endif
	int  M      = A.Rows();
	int  N      = A.Cols();
	char TransA = 'N';
	int  incX   = 1;
	int  incY   = 1;
	dgemv_(&TransA,        // CBLAS_TRANSPOSE A
	       &M,             // M
	       &N,             // N
	       &a,             // Alpha
	       A.GetPtr(),     // const double * A
	       &M,             // LDA
	       X.GetPtr(),     // const double * X
	       &incX,          // incX
	       &b,             // Beta
	       Y.GetPtr(),     // double * Y
	       &incY);         // incY
}

/// GEneral matrix^T * vector: \f$ \{Y\} \gets \alpha[A]^T\{X\} + \beta\{Y\}   \f$
inline void Gemtv(double const a, Matrix<double> const & A, Vector<double> const & X, double const b, Vector<double> & Y)
{
#ifndef DNDEBUG
	if (X.Size()!=A.Rows()) throw new Fatal(_("LinAlg::Gemtv: The size (%d) of vector X must be equal to the number of rows (%d) of matrix A"),X.Size(),A.Rows());
	if (Y.Size()!=A.Cols()) throw new Fatal(_("LinAlg::Gemtv: The size (%d) of vector Y must be equal to the number of columns (%d) of matrix A"),Y.Size(),A.Cols());
#endif
	int  M      = A.Rows();
	int  N      = A.Cols();
	char TransA = 'T';
	int  incX   = 1;
	int  incY   = 1;
	dgemv_(&TransA,        // CBLAS_TRANSPOSE A
	       &M,             // M
	       &N,             // N
	       &a,             // Alpha
	       A.GetPtr(),     // const double * A
	       &M,             // LDA
	       X.GetPtr(),     // const double * X
	       &incX,          // incX
	       &b,             // Beta
	       Y.GetPtr(),     // double * Y
	       &incY);         // incY
}

/// GEneral matrix makeR: \f$ [A] \gets \alpha\{X\}\{Y\}^T + [A]   \f$
inline void Ger(double const a, Vector<double> const & X, Vector<double> const & Y, Matrix<double> & A)
{
	/*  [A] <- a*{X}transp({Y}) + [A]
	     _                       _               _             _
	    |  00 01 02 03 04 ... 0N  |       / 0 \ |  0 1 2 ... N  |
	    |  10 11 12 13 14 ... 1N  |       | 1 |  -             -
	    |  20 21 22 23 24 ... 2N  | = a * | 2 |*                  + [A]
	    |          ...            |       |...|
	    |_ M0 M1 M2 M3 M4 ... MN _|       \ M /
	 */
#ifndef DNDEBUG
	if (X.Size()!=A.Rows()) throw new Fatal(_("LinAlg::Ger: The size (%d) of vector X must be equal to the number of rows (%d) of matrix A"),X.Size(),A.Rows());
	if (Y.Size()!=A.Cols()) throw new Fatal(_("LinAlg::Ger: The size (%d) of vector Y must be equal to the number of columns (%d) of matrix A"),Y.Size(),A.Cols());
#endif
	int M    = A.Rows();
	int N    = A.Cols();
	int incX = 1;
	int incY = 1;
	dger_(&M,             // M
	      &N,             // N
	      &a,             // Alpha
	      X.GetPtr(),     // const double * X
	      &incX,          // incX
	      Y.GetPtr(),     // const double * Y
	      &incY,          // incY
	      A.GetPtr(),     // double * A
	      &M);            // LDA
}

/// {Z} <- a*[A]*{X} + b*[B]*{Y}
inline void Gemvpmv(double const a, Matrix<double> const & A, Vector<double> const & X, double const b, Matrix<double> const & B, Vector<double> const & Y, Vector<double> & Z)
{
#ifndef NDEBUG
	if (B.Rows()!=A.Rows()) throw new Fatal(_("LinAlg::Gemvpmv: The number of rows (%d) of matrix B must be equal to the number of rows (%d) of matrix A"),B.Rows(),A.Rows());
	if (X.Size()!=A.Cols()) throw new Fatal(_("LinAlg::Gemvpmv: The size (%d) of vector X must be equal to the number of columns (%d) of matrix A"),X.Size(),A.Cols());
	if (Y.Size()!=B.Cols()) throw new Fatal(_("LinAlg::Gemvpmv: The size (%d) of vector Y must be equal to the number of columns (%d) of matrix B"),Y.Size(),B.Cols());
	if (Z.Size()!=A.Rows()) throw new Fatal(_("LinAlg::Gemvpmv: The size (%d) of vector Z must be equal to the number of rows (%d) of matrix A"),Z.Size(),A.Rows());
#endif
	Gemv(a,A,X, 0.0,Z);  // Z <- a*A*X + 0*Z
	Gemv(b,B,Y, 1.0,Z);  // Z <- b*B*Y + 1*Z
}


/* LEVEL 3 (matrix-matrix) */

/// GEneral matrix x matrix: \f$ [C] \gets \alpha[A][B] + \beta[C]    \f$
inline void Gemm(double const a, Matrix<double> const & A, Matrix<double> const & B, double const b, Matrix<double> & C)
{
/*  [C] <- a*[A]*[B] + b*[C] 
     _         _       _                       _  / 00 . 0N \      _         _
    |  00 . 0N  |     |  00 01 02 03 04 ... 0K  | | 10 . 1N |     |  00 . 0N  |
    |  10 . 1N  |     |  10 11 12 13 14 ... 1K  | | 20 . 2N |     |  10 . 1N  |
    |  20 . 2N  | = a*|  20 21 22 23 24 ... 2K  |*| 30 . 3N | + b*|  20 . 2N  |
    |    ...    |     |          ...            | | 40 . 4N |     |    ...    |
    |_ M0 . MN _|     |_ M0 M1 M2 M3 M4 ... MK _| |   ...   |     |_ M0 . MN _|
                                                  \ K0 . KN /
 */
#ifndef DNDEBUG
	if (B.Rows()!=A.Cols()) throw new Fatal(_("LinAlg::Gemm: The number of rows (%d) of matrix B must be equal to the number of columns (%d) of matrix A"),B.Rows(),A.Cols());
	if (C.Rows()!=A.Rows()) throw new Fatal(_("LinAlg::Gemm: The number of rows (%d) of matrix C must be equal to the number of rows (%d) of matrix A"),C.Rows(),A.Rows());
	if (C.Cols()!=B.Cols()) throw new Fatal(_("LinAlg::Gemm: The number of columns (%d) of matrix C must be equal to the number of columns (%d) of matrix B"),C.Cols(),B.Cols());
#endif
	int  M      = A.Rows();
	int  K      = A.Cols();
	int  N      = B.Cols();
	char TransA = 'N';
	char TransB = 'N';
	dgemm_(&TransA,       // CBLAS_TRANSPOSE A
	       &TransB,       // CBLAS_TRANSPOSE B
	       &M,            // M
	       &N,            // N
	       &K,            // K
	       &a,            // Alpha
	       A.GetPtr(),    // const double * A
	       &M,            // LDA
	       B.GetPtr(),    // const double * B
	       &K,            // LDB
	       &b,            // Beta
	       C.GetPtr(),    // double * C
	       &M);           // LDC
}

/// GEneral matrix^T x matrix: \f$ [C] \gets \alpha[A]^T[B] + \beta[C]    \f$
inline void Gemtm(double const a, Matrix<double> const & A, Matrix<double> const & B, double const b, Matrix<double> & C)
{
#ifndef DNDEBUG
	if (B.Rows()!=A.Rows()) throw new Fatal(_("LinAlg::Gemtm: The number of rows (%d) of matrix B must be equal to the number of rows (%d) of matrix A"),B.Rows(),A.Rows());
	if (C.Rows()!=A.Cols()) throw new Fatal(_("LinAlg::Gemtm: The number of rows (%d) of matrix C must be equal to the number of columns (%d) of matrix A"),C.Rows(),A.Cols());
	if (C.Cols()!=B.Cols()) throw new Fatal(_("LinAlg::Gemtm: The number of columns (%d) of matrix C must be equal to the number of columns (%d) of matrix B"),C.Cols(),B.Cols());
#endif
	int  K      = A.Rows();
	int  M      = A.Cols();
	int  N      = B.Cols();
	char TransA = 'T';
	char TransB = 'N';
	dgemm_(&TransA,       // CBLAS_TRANSPOSE A
	       &TransB,       // CBLAS_TRANSPOSE B
	       &M,            // M
	       &N,            // N
	       &K,            // K
	       &a,            // Alpha
	       A.GetPtr(),    // const double * A
	       &K,            // LDA
	       B.GetPtr(),    // const double * B
	       &K,            // LDB
	       &b,            // Beta
	       C.GetPtr(),    // double * C
	       &M);           // LDC
}

/// GEneral matrix x matrix^T: \f$ [C] \gets \alpha[A][B]^T+ \beta[C]    \f$
inline void Gemmt(double const a, Matrix<double> const & A, Matrix<double> const & B, double const b, Matrix<double> & C)
{
#ifndef DNDEBUG
	if (B.Cols()!=A.Cols()) throw new Fatal(_("LinAlg::Gemmt: The number of columns (%d) of matrix B must be equal to the number of columns (%d) of matrix A"),B.Cols(),A.Cols());
	if (C.Rows()!=A.Rows()) throw new Fatal(_("LinAlg::Gemmt: The number of rows (%d) of matrix C must be equal to the number of rows (%d) of matrix A"),C.Rows(),A.Rows());
	if (C.Cols()!=B.Rows()) throw new Fatal(_("LinAlg::Gemmt: The number of columns (%d) of matrix C must be equal to the number of rows (%d) of matrix B"),C.Cols(),B.Rows());
#endif
	int  M      = A.Rows();
	int  K      = A.Cols();
	int  N      = B.Rows();
	char TransA = 'N';
	char TransB = 'T';
	dgemm_(&TransA,       // CBLAS_TRANSPOSE A
	       &TransB,       // CBLAS_TRANSPOSE B
	       &M,            // M
	       &N,            // N
	       &K,            // K
	       &a,            // Alpha
	       A.GetPtr(),    // const double * A
	       &M,            // LDA
	       B.GetPtr(),    // const double * B
	       &N,            // LDB
	       &b,            // Beta
	       C.GetPtr(),    // double * C
	       &M);           // LDC
}

/// GEneral matrix^T x matrix^T: \f$ [C] \gets \alpha[A]^T[B]^T + \beta[C]    \f$
inline void Gemtmt(double const a, Matrix<double> const & A, Matrix<double> const & B, double const b, Matrix<double> & C)
{
#ifndef DNDEBUG
	if (B.Cols()!=A.Rows()) throw new Fatal(_("LinAlg::Gemmt: The number of columns (%d) of matrix B must be equal to the number of rows (%d) of matrix A"),B.Cols(),A.Rows());
	if (C.Rows()!=A.Cols()) throw new Fatal(_("LinAlg::Gemmt: The number of rows (%d) of matrix C must be equal to the number of columns (%d) of matrix A"),C.Rows(),A.Cols());
	if (C.Cols()!=B.Rows()) throw new Fatal(_("LinAlg::Gemmt: The number of columns (%d) of matrix C must be equal to the number of rows (%d) of matrix B"),C.Cols(),B.Rows());
#endif
	int  K      = A.Rows();
	int  M      = A.Cols();
	int  N      = B.Rows();
	char TransA = 'T';
	char TransB = 'T';
	dgemm_(&TransA,       // CBLAS_TRANSPOSE A
	       &TransB,       // CBLAS_TRANSPOSE B
	       &M,            // M
	       &N,            // N
	       &K,            // K
	       &a,            // Alpha
	       A.GetPtr(),    // const double * A
	       &K,            // LDA
	       B.GetPtr(),    // const double * B
	       &N,            // LDB
	       &b,            // Beta
	       C.GetPtr(),    // double * C
	       &M);           // LDC
}


// ########################################################################## Solver

/// GEneral SolVer: \f$ \{Y\} \gets [A]^{-1}\{Y\}   \f$   OBS.: This function does not preserve A
inline int Gesv(Matrix<double> & A, Vector<double> & Y)
{
#ifndef DNDEBUG
	if (A.Rows()!=A.Cols()) throw new Fatal(_("LinAlg::Gesv: The number of rows (%d) of matrix A must be equal to the number of columns (%d) of matrix A"),A.Rows(),A.Cols());
	if (Y.Size()!=A.Cols()) throw new Fatal(_("LinAlg::Gesv: The size (%d) of vector Y must be equal to the number of columns (%d) of matrix A"),Y.Size(),A.Cols());
#endif
	int   info = 0;
	int   N    = A.Rows(); // == A.Cols
	int   NRHS = 1;        //int NRHS = Y.Cols();
	int * ipiv = new int [N];
	dgesv_(&N,                 // A(N,N)
	       &NRHS,              // {Y}(N,NRHS) (RHS: Right Hand Side) 
	       A.GetPtr(),         // double * A
	       &N,                 // LDA
	       ipiv,               // Pivot Indices
	       Y.GetPtr(),         // double * Y
	       &N,                 // LDY
	       &info);
	delete [] ipiv;
	if (info!=0) throw new Fatal (_("LinAlg::Gesv: Linear solver could not find the solution (singular matrix?)"));
	return info;
}

/// GE Inverse
inline int Geinv(Matrix<double> & A)
{
#ifndef DNDEBUG
	if (A.Rows()!=A.Cols()) throw new Fatal(_("LinAlg::Geinv: The number of rows (%d) of matrix A must be equal to the number of columns (%d) of matrix A"),A.Rows(),A.Cols());
#endif
	int   info = 0;
	int   N    = A.Rows(); // == A.Cols
	int * ipiv = new int [N];

	// Factorization
	dgetrf_(&N,         // M
	        &N,         // N
	        A.GetPtr(), // double * A
	        &N,         // LDA
	        ipiv,       // Pivot indices
	        &info);
	if (info!=0) throw new Fatal ("LinAlg::Geinv: Matrix LU factorization did not work");

	int      NB    = 2;                  // Optimal blocksize: TODO: Use ILAENV to find best block size
	int      lwork = N*NB;               // Dimension of work >= max(1,N), optimal=N*NB
	double * work  = new double [lwork]; // Work

	// Inversion
	dgetri_(&N,         // N
	        A.GetPtr(), // double * A
	        &N,         // LDA
	        ipiv,       // Pivot indices
	        work,       // work
	        &lwork,     // dimension of work
	        &info);
	if (info!=0) throw new Fatal ("LinAlg::Geinv: Matrix inversion did not work");

	delete [] ipiv;
	delete [] work;
	return info;
}

/// General Tridiagonal SolVer: \f$ \{Y\} \gets [A]^{-1}\{Y\}   \f$
inline int Gtsv(Vector<double> & DL, Vector<double> & D, Vector<double> & DU, Vector<double> & Y)
{
	/*   {Y} <- inv([A]) * {Y}
	
	     DL = [l0, l1, ..., lN-1]
	     D  = [d0, d1, d2, ..., dN]
	     DU = [u0, u1, ..., uN-1]
	
	               [A]            {Y}       {Y}
	      / d0 u0  0  0  0  0 \  / x0 \   / y0 \
	      | l0 d1 u1  0  0  0 |  | x1 |   | y1 |
	      | 0  l1 d2 u2  0  0 |  | x2 | = | y2 |
	      | 0  0  l2 d3 u3  0 |  | x3 |   | y3 |
	      |          l3.d4 u4 |  | x4 |   | y4 |
	      \ 0  0  0  0  l4 d5 /  \ x5 /   \ y5 /
	*/
#ifndef DNDEBUG
	if (Y.Size()!=D.Size())    throw new Fatal(_("LinAlg::Gtsv: The size (%d) of vector Y must be equal to the size (%d) of vector D"),Y.Size(),D.Size());
	if (DL.Size()!=D.Size()-1) throw new Fatal(_("LinAlg::Gtsv: The size (%d) of vector DL must be equal to the size (%d) of vector D minus 1"),DL.Size(),D.Size()-1);
	if (DU.Size()!=D.Size()-1) throw new Fatal(_("LinAlg::Gtsv: The size (%d) of vector DU must be equal to the size (%d) of vector D minus 1"),DU.Size(),D.Size()-1);
#endif
	int info = 0;
	/* not working with Ubuntu's LAPACK package
	int N    = D.Size();
	int NRHS = 1;       //int NRHS = Y.Cols();
	dgtsv_(&N,          // A(N,N)
		   &NRHS,       // {Y}(N,NRHS) (RHS: Right Hand Side)
		   DL.GetPtr(), // double * DL
		   D.GetPtr(),  // double * D
		   DU.GetPtr(), // double * DU
		   Y.GetPtr(),  // double * Y
		   &N,          // LDY
		   &info);      // int * INFO
	*/
	return info;
}


// ################################################################### Eigenproblems

/// SYmmetric EigenvaValue/Vector RRR (Relatively Robust Representations): \f$ [Q,L] \gets eig([A])   \f$
inline int Syevr(Matrix<double> & A, Vector<double> & L, Matrix<double> & Q)
{
#ifndef DNDEBUG
	if (A.Cols()!=A.Rows()) throw new Fatal(_("LinAlg::Gesv: The number of columns (%d) of matrix A must be equal to the number of rows (%d) of matrix A"),A.Cols(),A.Rows());
	if (L.Size()!=A.Rows()) throw new Fatal(_("LinAlg::Gesv: The size (%d) of vector L must be equal to the number of rows (%d) of matrix A"),L.Size(),A.Rows());
	if (Q.Rows()!=A.Rows()) throw new Fatal(_("LinAlg::Gesv: The number of rows (%d) of matrix Q must be equal to the number of rows (%d) of matrix A"),Q.Rows(),A.Rows());
	if (Q.Cols()!=A.Rows()) throw new Fatal(_("LinAlg::Gesv: The number of columns (%d) of matrix Q must be equal to the number of rows (%d) of matrix A"),Q.Cols(),A.Rows());
#endif
	int      info   = 0;
	/* not working with Ubuntu's LAPACK package
	int      N      = A.Rows();
	char     jobz   = 'V';
	char     range  = 'A';
	char     uplo   = 'L';
	char     cmach  = 'S';
	double   abstol = dlamch_(&cmach);
	int    * isuppz = new int [2*N];
	int      lwork  = 26*N;
	double * work   = new double [lwork];
	int      liwork = 10*N;
	int    * iwork  = new int [liwork];
	double   vl,vu;
	int      il,iu;
	int      m;
	dsyevr_(&jobz,      // (IN) JOBZ:  'V'=compute eigenvals and eigenvecs
			&range,     // (IN) RANGE: 'A'=all eigenvals
			&uplo,      // (IN) UPLO:  'L'=lower triangle of A is stored
			&N,         // (IN) N
			A.GetPtr(), // (IN/OUT) double * A
			&N,         // (IN) LDA
			&vl,        // (IN) VL for RANGE='V' ** not used **
			&vu,        // (IN) VU for RANGE='V' ** not used **
			&il,        // (IN) IL for RANGE='I' ** not used **
			&iu,        // (IN) IU for RANGE='I' ** not used **
			&abstol,    // (IN) ABSTOL
			&m,         // (OUT) M = total number of eigenvals found (M=N)
			L.GetPtr(), // (OUT) sorted (ascending) eigenvalues
			Q.GetPtr(), // (OUT) columns contain the eigenvecs corresp. to eigenvs
			&N,         // (IN) LDQ
			isuppz,     // (OUT) int    * ISUPPZ(2*max(1,M))
			work,       // (OUT) double * WORK(LWORK)
			&lwork,     // (IN)  int    * LWORK >= max(1,26*N)
			iwork,      // (OUT) int    * IWORK(LIWORK)
			&liwork,    // (IN)  int    * LIWORK >= max(1,10*N)
			&info);     // (OUT) INFO
	delete [] isuppz;
	delete [] work;
	delete [] iwork;
	*/
	return info;
}

/// SYmmetric EigenProjectors: \f$ P_{i=1\cdots N} \gets eigenprojectors([A])   \f$
inline int Syep(Matrix<double> & A, Vector<double> & L, Matrix<double> P[]) ///< In/Out: A(N,N) => symmetric matrix Out: L(N) => sorted (ascending) eigenvalues Out: P(N,N)[N] => N eigenprojectors
{
#ifndef DNDEBUG
	if (A.Cols()!=A.Rows()) throw new Fatal(_("LinAlg::Gesv: The number of columns (%d) of matrix A must be equal to the number of rows (%d) of matrix A"),A.Cols(),A.Rows());
	if (L.Size()!=A.Rows()) throw new Fatal(_("LinAlg::Gesv: The size (%d) of vector L must be equal to the number of rows (%d) of matrix A"),L.Size(),A.Rows());
#endif
	int N = A.Rows();
	Matrix<double> Q(N,N); // columns are eigenVector<double>s
	int info = Syevr(A,L,Q);
	if (info==0)
	{
		for (int i=0; i<N; ++i)
		{
#ifndef DNDEBUG
	if (P[i].Rows()!=A.Rows()) throw new Fatal(_("LinAlg::Gesv: The number of rows (%d) of matrix P[%d] must be equal to the number of rows (%d) of matrix A"),P[i].Rows(),i,A.Rows());
	if (P[i].Cols()!=A.Rows()) throw new Fatal(_("LinAlg::Gesv: The number of columns (%d) of matrix P[%d] must be equal to the number of rows (%d) of matrix A"),P[i].Cols(),i,A.Rows());
#endif
			P[i].SetValues(0.0);
			cblas_dger(CblasColMajor, // CBLAS_ORDER
			           N,             // M
			           N,             // N
			           1,             // Alpha
			           &(Q(0,i)),     // const double * X
			           1,             // incX
			           &(Q(0,i)),     // const double * Y
			           1,             // incY
			           P[i].GetPtr(), // double * A
			           N);            // LDA
		}
	}
	return info;
}

/// SYmmetric Isotropic Function: \f$ [IF] \gets isotropicfunc([A])   \f$
inline int Syif(Matrix<double> & A, Matrix<double> & IF, double (*func) (double)) ///< In/Out: A(N,N). Out: IF(N,N). In: (*func) => Isotropic function, ex.: &sqrt, &exp, etc.
{
#ifndef DNDEBUG
	if (A.Cols()!=A.Rows())  throw new Fatal(_("LinAlg::Gesv: The number of columns (%d) of matrix A must be equal to the number of rows (%d) of matrix A"),A.Cols(),A.Rows());
	if (IF.Rows()!=A.Rows()) throw new Fatal(_("LinAlg::Gesv: The number of rows (%d) of matrix IF must be equal to the number of rows (%d) of matrix A"),IF.Rows(),A.Rows());
	if (IF.Cols()!=A.Rows()) throw new Fatal(_("LinAlg::Gesv: The number of columns (%d) of matrix IF must be equal to the number of rows (%d) of matrix A"),IF.Cols(),A.Rows());
#endif
	int N = A.Rows();
	Matrix<double> Q(N,N); // columns are eigenVector<double>s
	Vector<double> L(N);   // components are eigenvalues
	double alpha;
	int info = Syevr(A,L,Q);
	if (info==0)
	{
		IF.SetValues(0.0);
		for (int i=0; i<N; ++i)
		{
			alpha = (*func) (L(i));
			cblas_dger(CblasColMajor, // CBLAS_ORDER
			           N,             // M
			           N,             // N
			           alpha,         // Alpha
			           &(Q(0,i)),     // const double * X
			           1,             // incX
			           &(Q(0,i)),     // const double * Y
			           1,             // incY
			           &(IF(0,0)),    // double * A
			           N);            // LDA
		}
	}
	return info;
}


}; // namespace LinAlg

#endif //#define MECHSYS_LINALG_LAWRAP_H
