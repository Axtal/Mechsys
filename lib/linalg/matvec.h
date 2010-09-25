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

#ifndef MECHSYS_MATVEC_H
#define MECHSYS_MATVEC_H

// Std Lib
#include <iostream>
#include <cstring>   // for strcmp
#include <cmath>     // for sqrt, pow
#include <algorithm> // for min, max

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

// MechSys
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/util/array.h>

// LAPACK
extern "C"
{
    // DGETRF - compute an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges
    void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);

    // DGETRI - compute the inverse of a matrix using the LU factorization computed by DGETRF
    void dgetri_(const int* n, double* a, const int* lda, int* ipiv, double* work, const int* lwork, int* info);

    // DGESVD - computes the singular value decomposition of a real M-by-N matrix A, optionally computing the left and/or right singular vectors
    void dgesvd_(const char* jobu, const char* jobvt, const int* M, const int* N, double* A, const int* lda, double* S, double* U, const int* ldu, double* VT, const int* ldvt, double* work,const int* lwork, const int* info);

    // DGESV - double-general-solver
    void dgesv_(int *Np, int *NRHSp, double *A, int *LDAp, int *IPIVp, double *B, int *LDBp, int *INFOp);
}


/////////////////////////////////////////////////////////////////////////////////////////// General MatVec /////


#ifdef USE_MTL4

// Boost/MTL4
#include <boost/numeric/mtl/mtl.hpp>

/** Dense matrix (general). */
typedef mtl::dense2D<double, mtl::matrix::parameters<mtl::tag::col_major> > Mat_t;

/** Dense vector (general). */
typedef mtl::dense_vector<double> Vec_t;

#else

#include <mechsys/linalg/vector.h>
#include <mechsys/linalg/matrix.h>
#include <mechsys/linalg/laexpr.h>

typedef LinAlg::Vector<double> Vec_t;
typedef LinAlg::Matrix<double> Mat_t;

inline double dot         (Vec_t const & V, Vec_t const & W) { return LinAlg::Dot(V,W); }
inline size_t size        (Vec_t const & V) { return V.Size(); }
inline void   set_to_zero (Vec_t       & V) { V.SetValues(0.0); }
inline void   set_to_zero (Mat_t       & M) { M.SetValues(0.0); }
inline size_t num_rows    (Vec_t const & V) { return V.Size(); }
inline size_t num_rows    (Mat_t const & M) { return M.Rows(); }
inline size_t num_cols    (Mat_t const & M) { return M.Cols(); }

#endif

/** Print vector. */
inline String PrintVector (Vec_t const & V, char const * Fmt="%13g", Array<long> const * SkipR=NULL, double Tol=1.0e-13)
{
    int m = size(V);
    String lin;
    for (int i=0; i<m; ++i)
    {
        bool skip_row = false;
        if (SkipR!=NULL) skip_row = (SkipR->Find(i)<0 ? false : true);
        if (!skip_row)
        {
            double val = (fabs(V(i))<Tol ? 0.0 : V(i));
            String buf;  buf.Printf(Fmt,val);
            lin.append(buf);
        }
    }
    lin.append("\n");
    return lin;
}

/** Print matrix. */
inline String PrintMatrix (Mat_t const & M, char const * Fmt="%13g", Array<long> const * SkipRC=NULL, double Tol=1.0e-13, bool NumPy=false)
{
    int m = M.num_rows();
    int n = M.num_cols();
    String lin;
    if (NumPy) lin.append("matrix([[");
    for (int i=0; i<m; ++i)
    {
        bool skip_row = false;
        if (SkipRC!=NULL) skip_row = (SkipRC->Find(i)<0 ? false : true);
        if (!skip_row)
        {
            if (NumPy && i!=0) lin.append("        [");
            for (int j=0; j<n; ++j)
            {
                bool skip_col = false;
                if (SkipRC!=NULL) skip_col = (SkipRC->Find(j)<0 ? false : true);
                if (!skip_col)
                {
                    double val = (fabs(M(i,j))<Tol ? 0.0 : M(i,j));
                    String buf;  buf.Printf(Fmt,val);
                    lin.append(buf);
                    if (NumPy && j!=n-1) lin.append(",");
                }
            }
            if (NumPy && i!=m-1) lin.append("],");
            if (NumPy && i==m-1) lin.append("]])");
            else lin.append("\n");
        }
    }
    return lin;
}

/** Compare two vectors. */
inline double CompareVectors (Vec_t const & A, Vec_t const & B)
{
    size_t m = size(A);
    if (m!=size(B)) throw new Fatal("matvec.h:CompareVectors: vectors A_%d and B_%d must have the same size",m,size(B));
    double error = 0.0;
    for (size_t i=0; i<m; ++i)
        error += fabs(A(i)-B(i));
    return error;
}

/** Compare two matrices. */
inline double CompareMatrices (Mat_t const & A, Mat_t const & B)
{
    size_t m = A.num_rows();
    size_t n = A.num_cols();
    if ((m!=B.num_rows()) || (n!=B.num_cols())) throw new Fatal("matvec.h:CompareMatrices: matrices A_%dx%d and B_%dx%d must have the same number of rows and columns",m,n,B.num_rows(),B.num_cols());
    double error = 0.0;
    for (size_t i=0; i<m; ++i)
    for (size_t j=0; j<n; ++j)
        error += fabs(A(i,j)-B(i,j));
    return error;
}

/** Check if matrix is diagonal. */
inline double CheckDiagonal (Mat_t const & M, bool CheckUnitDiag=false)
{
    size_t m = M.num_rows();
    size_t n = M.num_cols();
    double error = 0.0;
    for (size_t i=0; i<m; ++i)
    for (size_t j=0; j<n; ++j)
    {
        if ((i==j) && CheckUnitDiag) error += fabs(M(i,j)-1.0);
        if  (i!=j)                   error += fabs(M(i,j));
    }
    return error;
}

/** Determinant. */
inline double Det (Mat_t const & M)
{
    int m = M.num_rows();
    int n = M.num_cols();
    if (m==1)
    {
        double res = 0;
        for (int i=0; i<n; i++) res += M(0,i)*M(0,i);
        return sqrt(res);
    }
    else if (m==2 && n==2)
    {
        return M(0,0)*M(1,1) - M(1,0)*M(0,1);
    }
    /*
    else if (m==2 && n==3)
    {
        double d1 = M(0,0)*M(1,1) - M(0,1)*M(1,0);
        double d2 = M(0,1)*M(1,2) - M(0,2)*M(1,1);
        double d3 = M(0,2)*M(1,0) - M(0,0)*M(1,2);
        return sqrt(d1*d1 + d2*d2 + d3*d3);
    }
    */
    else if (m==3 && n==3)
    {
        return  M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
              - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
              + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));
    }
    else if (m==n)
    {
        // factorization
        int   info = 0;
        int * ipiv = new int [m];
        Mat_t Mcpy(M);
        dgetrf_(&m,         // M
                &m,         // N
                Mcpy.data,  // double * A
                &m,         // LDA
                ipiv,       // Pivot indices
                &info);     // INFO
        if (info!=0) throw new Fatal ("matvec.h::Det: LAPACK: LU factorization failed");

        // determinant
        double det = 1.0;
        for (int i=0; i<m; ++i)
        {
            if (ipiv[i]!=(i+1)) det = -det * Mcpy(i,i);
            else                det =  det * Mcpy(i,i);
        }

        // end
        delete [] ipiv;
        return det;
    }
    else throw new Fatal("matvec.h:Det: Method is not implemented for (%d x %d) matrices yet",m,n);
}

/** Identity. */
inline void Identity (size_t NRows, Mat_t & I)
{
    I.change_dim (NRows,NRows);
    for (size_t i=0; i<NRows; ++i)
    for (size_t j=0; j<NRows; ++j)
    {
        if (i==j) I(i,j) = 1.0;
        else      I(i,j) = 0.0;
    }
}

/** Singular value decomposition. M = U_mxm * D_mxn * Vt_nxn   */
inline void Svd (Mat_t const & M, Mat_t & U, Vec_t & S, Mat_t & Vt)
{
    int  info   = 0;
    char job    = 'A';
    int  m      = M.num_rows();
    int  n      = M.num_cols();
    int  min_mn = (m<n ? m : n);
    int  max_mn = (m>n ? m : n);
    int  lwork  = 2.0*std::max(3*min_mn+max_mn, 5*min_mn);

    U. change_dim (m, m);
    Vt.change_dim (n, n); // trans(V)
    S. change_dim (min_mn);

    double * work  = new double [lwork]; // Work

    // decomposition
    Mat_t tmp(M);
    dgesvd_(&job,      // JOBU
            &job,      // JOBVT
            &m,        // M
            &n,        // N
            tmp.data,  // A
            &m,        // LDA
            S.data,    // S
            U.data,    // U
            &m,        // LDU
            Vt.data,   // VT
            &n,        // LDVT
            work,      // WORK
            &lwork,    // LWORK
            &info);    // INFO
    if (info!=0) throw new Fatal ("matvec::Svd: LAPACK: Decomposition failed");

    delete [] work;
}

/** Inverse. */
inline void Inv (Mat_t const & M, Mat_t & Mi, double Tol=1.0e-10)
{
    int m = M.num_rows();
    int n = M.num_cols();
    Mi.change_dim(m,n);
    if (m==2 && n==2)
    {
        double det = Det(M);
        if (fabs(det)<Tol) throw new Fatal("matvec.h:Inv: Cannot calculate inverse due to zero determinant. det = %g",det);

        Mi(0,0) =  M(1,1) / det;
        Mi(0,1) = -M(0,1) / det;

        Mi(1,0) = -M(1,0) / det;
        Mi(1,1) =  M(0,0) / det;
    }
    else if (m==3 && n==3)
    {
        double det = Det(M);
        if (fabs(det)<Tol) throw new Fatal("matvec.h:Inv: Cannot calculate inverse due to zero determinant. det = %g",det);

        Mi(0,0) = (M(1,1)*M(2,2) - M(1,2)*M(2,1)) / det;
        Mi(0,1) = (M(0,2)*M(2,1) - M(0,1)*M(2,2)) / det;
        Mi(0,2) = (M(0,1)*M(1,2) - M(0,2)*M(1,1)) / det;

        Mi(1,0) = (M(1,2)*M(2,0) - M(1,0)*M(2,2)) / det;
        Mi(1,1) = (M(0,0)*M(2,2) - M(0,2)*M(2,0)) / det;
        Mi(1,2) = (M(0,2)*M(1,0) - M(0,0)*M(1,2)) / det;

        Mi(2,0) = (M(1,0)*M(2,1) - M(1,1)*M(2,0)) / det;
        Mi(2,1) = (M(0,1)*M(2,0) - M(0,0)*M(2,1)) / det;
        Mi(2,2) = (M(0,0)*M(1,1) - M(0,1)*M(1,0)) / det;
    }
    else if (m==n) // square
    {
        int   info = 0;
        int * ipiv = new int [m];

        // factorization
        Mi = M;
        dgetrf_(&m,       // M
                &m,       // N
                Mi.data,  // double * A
                &m,       // LDA
                ipiv,     // Pivot indices
                &info);   // INFO
        if (info!=0) throw new Fatal ("matvec.h::Inv: LAPACK: LU factorization failed");

        int      NB    = 4;                  // Optimal blocksize ?
        int      lwork = m*NB;               // Dimension of work >= max(1,m), optimal=m*NB
        double * work  = new double [lwork]; // Work

        // inversion
        dgetri_(&m,       // N
                Mi.data,  // double * A
                &m,       // LDA
                ipiv,     // Pivot indices
                work,     // work
                &lwork,   // dimension of work
                &info);   // INFO
        if (info!=0) throw new Fatal ("matvec::Inv: LAPACK: Inversion failed");

        delete [] ipiv;
        delete [] work;
    }
    else // generalized (pseudo) inverse
    {
        Mat_t U; Vec_t S; Mat_t Vt;
        Svd(M, U, S, Vt);
        Mat_t Di(m,n);
        set_to_zero(Di);
        for (size_t i=0; i<size(S); ++i)
        {
            if (S(i)>Tol) Di(i,i) = 1.0/S(i);
        }
        Mat_t tmp(U*Di*Vt);
        Mi.change_dim(n,m);
        Mi = trans(tmp);
    }
}

/** Linear Solver. {X} = [M]^{-1}{X} (M is lost) (X initially has the contents of the right-hand side) */
inline void Sol (Mat_t & M, Vec_t & X)
{
    int  m = M.num_rows();
    int  n = M.num_cols();
    int mv = size(X);
    if (m!=n)  throw new Fatal("Sol: Matrix must be square");
    if (m!=mv) throw new Fatal("Sol: Vector X must have the same number of rows of matrix M");

    int   info = 0;
    int   nrhs = 1; // vector X has 1 column
    int * ipiv = new int [m];
    dgesv_(&m,      // A(m,m)
           &nrhs,   // {X}(m,1) (RHS: Right Hand Side) 
           M.data,  // double * A
           &m,      // LDA
           ipiv,    // Pivot Indices
           X.data,  // double * Y
           &m,      // LDY
           &info);  // info
    delete [] ipiv;

    if (info!=0) 
    {
        throw new Fatal ("Sol: Linear solver (DGESV) failed (singular matrix?)");
    }
}

/** Linear Solver. {X} = [M]^{-1}{B}  */
inline void Sol (Mat_t const & M, Vec_t const & B, Vec_t & X)
{
    Mat_t m(M);
    X = B;
    Sol (m, X);
}

/** Norm. */
inline double Norm (Vec_t const & V)
{
    return sqrt(dot(V,V));
}

/** Dyadic product. */
inline void Dyad (Vec_t const & A, Vec_t const & B, Mat_t & M)
{
    M.change_dim(size(A),size(B));
    for (size_t i=0; i<size(A); ++i)
    for (size_t j=0; j<size(B); ++j)
        M(i,j) = A(i) * B(j);
}

/** Left multiplication. {B} = {A}*[M]. NOTE: this is not efficient for large matrices.  */
inline void Mult (Vec_t const & A, Mat_t const & M, Vec_t & B)
{
    B.change_dim (M.num_cols());
    set_to_zero  (B);
    for (size_t i=0; i<M.num_rows(); ++i)
    for (size_t j=0; j<M.num_cols(); ++j)
        B(j) += A(i)*M(i,j);
}

/** Create column matrix from vector. */
inline void Vec2ColMat (Vec_t const & V, Mat_t & M)
{
    M.change_dim (size(V),1);
    for (size_t i=0; i<size(V); ++i) M(i,0) = V(i);
}


////////////////////////////////////////////////////////////////////////////////////////////// Tiny MatVec /////


/** 3x3 Matrix. */
typedef blitz::TinyMatrix<double,3,3> Mat3_t;

/** 3x1 Vector. */
typedef blitz::TinyVector<double,3> Vec3_t;
typedef blitz::TinyVector<size_t,3> iVec3_t;
typedef blitz::TinyVector<bool,3>   bVec3_t;

/** Print vector. */
inline String PrintVector (Vec3_t const & V, char const * Fmt="%13g", double Tol=1.0e-13)
{
    int m = 3;
    String lin;
    for (int i=0; i<m; ++i)
    {
        double val = (fabs(V(i))<Tol ? 0.0 : V(i));
        String buf;  buf.Printf(Fmt,val);
        lin.append(buf);
    }
    lin.append("\n");
    return lin;
}

/** Print matrix. */
inline String PrintMatrix (Mat3_t const & M, char const * Fmt="%13g", double Tol=1.0e-13)
{
    int m = 3;
    int n = 3;
    String lin;
    for (int i=0; i<m; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            double val = (fabs(M(i,j))<Tol ? 0.0 : M(i,j));
            String buf;  buf.Printf(Fmt,val);
            lin.append(buf);
        }
        lin.append("\n");
    }
    return lin;
}

/** Compare two vectors. */
inline double CompareVectors (Vec3_t const & A, Vec3_t const & B)
{
    double error = 0.0;
    for (size_t i=0; i<3; ++i)
        error += fabs(A(i)-B(i));
    return error;
}

/** Compare two matrices. */
inline double CompareMatrices (Mat3_t const & A, Mat3_t const & B)
{
    double error = 0.0;
    for (size_t i=0; i<3; ++i)
    for (size_t j=0; j<3; ++j)
        error += fabs(A(i,j)-B(i,j));
    return error;
}

/** Transpose.*/
inline void Trans (Mat3_t const & M, Mat3_t & Mt)
{
    Mt(0,0)=M(0,0);   Mt(0,1)=M(1,0);   Mt(0,2)=M(2,0);
    Mt(1,0)=M(0,1);   Mt(1,1)=M(1,1);   Mt(1,2)=M(2,1);
    Mt(2,0)=M(0,2);   Mt(2,1)=M(1,2);   Mt(2,2)=M(2,2);
}

/** Determinant.*/
inline double Det (Mat3_t const & M)
{
    double det =   M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
                 - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
                 + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));
    return det;
}

/** Inverse.*/
inline void Inv (Mat3_t const & M, Mat3_t & Mi, double Tol=1.0e-10)
{
    double det =   M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
                 - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
                 + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));

    if (fabs(det)<Tol)
    {
        std::ostringstream oss;
        oss << PrintMatrix(M);
        throw new Fatal("matvec.h::Inv: 3x3 matrix inversion failed with null (%g) determinat. M =\n%s",Tol,oss.str().c_str());
    }

    Mi(0,0)=(M(1,1)*M(2,2)-M(1,2)*M(2,1))/det;  Mi(0,1)=(M(0,2)*M(2,1)-M(0,1)*M(2,2))/det;  Mi(0,2)=(M(0,1)*M(1,2)-M(0,2)*M(1,1))/det;
    Mi(1,0)=(M(1,2)*M(2,0)-M(1,0)*M(2,2))/det;  Mi(1,1)=(M(0,0)*M(2,2)-M(0,2)*M(2,0))/det;  Mi(1,2)=(M(0,2)*M(1,0)-M(0,0)*M(1,2))/det;
    Mi(2,0)=(M(1,0)*M(2,1)-M(1,1)*M(2,0))/det;  Mi(2,1)=(M(0,1)*M(2,0)-M(0,0)*M(2,1))/det;  Mi(2,2)=(M(0,0)*M(1,1)-M(0,1)*M(1,0))/det;
}

/** Linear Solver. {X} = [M]^{-1}{B}  */
inline void Sol (Mat3_t const & M, Vec3_t const & B, Vec3_t & X, double Tol=1.0e-10)
{
    // determinant
    double det =   M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
                 - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
                 + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));
    if (fabs(det)<Tol) throw new Fatal("matvec.h:Sol: Cannot calculate inverse due to zero determinant. det = %g",det);

    // inverse matrix
    Mat3_t Mi;
    Mi(0,0) = (M(1,1)*M(2,2) - M(1,2)*M(2,1)) / det;
    Mi(0,1) = (M(0,2)*M(2,1) - M(0,1)*M(2,2)) / det;
    Mi(0,2) = (M(0,1)*M(1,2) - M(0,2)*M(1,1)) / det;

    Mi(1,0) = (M(1,2)*M(2,0) - M(1,0)*M(2,2)) / det;
    Mi(1,1) = (M(0,0)*M(2,2) - M(0,2)*M(2,0)) / det;
    Mi(1,2) = (M(0,2)*M(1,0) - M(0,0)*M(1,2)) / det;

    Mi(2,0) = (M(1,0)*M(2,1) - M(1,1)*M(2,0)) / det;
    Mi(2,1) = (M(0,1)*M(2,0) - M(0,0)*M(2,1)) / det;
    Mi(2,2) = (M(0,0)*M(1,1) - M(0,1)*M(1,0)) / det;

    // solve system
    X(0) = Mi(0,0)*B(0) + Mi(0,1)*B(1) + Mi(0,2)*B(2);
    X(1) = Mi(1,0)*B(0) + Mi(1,1)*B(1) + Mi(1,2)*B(2);
    X(2) = Mi(2,0)*B(0) + Mi(2,1)*B(1) + Mi(2,2)*B(2);
}

/** Eigenvalues. NOTE: This function changes the matrix M. */
inline void Eig (Mat3_t & M, Vec3_t & L)
{
    // calculate
    gsl_matrix_view m = gsl_matrix_view_array (M.data(), 3, 3);
    gsl_vector * eval = gsl_vector_alloc      (3);
    gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (3);
    gsl_eigen_symm (&m.matrix, eval, w);

    // eigenvalues
    L = gsl_vector_get(eval,0), gsl_vector_get(eval,1), gsl_vector_get(eval,2);

    // clean up
    gsl_eigen_symm_free (w);
    gsl_vector_free     (eval);
}

/** Eigenvalues and eigenvectors. NOTE: This function changes the matrix M. */
inline void Eig (Mat3_t & M, Vec3_t & L, Vec3_t & V0, Vec3_t & V1, Vec3_t & V2, bool SortAsc=false, bool SortDesc=false)
{
    // calculate
    gsl_matrix_view m = gsl_matrix_view_array (M.data(), 3, 3);
    gsl_vector * eval = gsl_vector_alloc      (3);
    gsl_matrix * evec = gsl_matrix_alloc      (3, 3);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
    gsl_eigen_symmv (&m.matrix, eval, evec, w);

    // sort
    if (SortAsc)  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    if (SortDesc) gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);

    // eigenvalues
    L = gsl_vector_get(eval,0), gsl_vector_get(eval,1), gsl_vector_get(eval,2);

    // eigenvectors
    gsl_vector_view ev = gsl_matrix_column (evec,0);
    V0 = gsl_vector_get    (&ev.vector,0), gsl_vector_get(&ev.vector,1), gsl_vector_get(&ev.vector,2);
    ev = gsl_matrix_column (evec,1);
    V1 = gsl_vector_get    (&ev.vector,0), gsl_vector_get(&ev.vector,1), gsl_vector_get(&ev.vector,2);
    ev = gsl_matrix_column (evec,2);
    V2 = gsl_vector_get    (&ev.vector,0), gsl_vector_get(&ev.vector,1), gsl_vector_get(&ev.vector,2);

    // clean up
    gsl_eigen_symmv_free (w);
    gsl_vector_free      (eval);
    gsl_matrix_free      (evec);
}

/** Norm. */
inline double Norm (Vec3_t const & V)
{
    return sqrt(blitz::dot(V,V));
}

/** Dyadic product. */
inline void Dyad (Vec3_t const & A, Vec3_t const & B, Mat3_t & M)
{
    M(0,0)=A(0)*B(0);  M(0,1)=A(0)*B(1);  M(0,2)=A(0)*B(2);
    M(1,0)=A(1)*B(0);  M(1,1)=A(1)*B(1);  M(1,2)=A(1)*B(2);
    M(2,0)=A(2)*B(0);  M(2,1)=A(2)*B(1);  M(2,2)=A(2)*B(2);
}

/** Matrix multiplication. */
inline void Mult (Mat3_t const & A, Mat3_t const & B, Mat3_t & M)
{
    M(0,0)=A(0,2)*B(2,0)+A(0,1)*B(1,0)+A(0,0)*B(0,0);  M(0,1)=A(0,2)*B(2,1)+A(0,1)*B(1,1)+A(0,0)*B(0,1);  M(0,2)=A(0,2)*B(2,2)+A(0,1)*B(1,2)+A(0,0)*B(0,2);
    M(1,0)=A(1,2)*B(2,0)+B(1,0)*A(1,1)+B(0,0)*A(1,0);  M(1,1)=A(1,2)*B(2,1)+A(1,1)*B(1,1)+B(0,1)*A(1,0);  M(1,2)=A(1,2)*B(2,2)+A(1,1)*B(1,2)+B(0,2)*A(1,0);
    M(2,0)=A(2,2)*B(2,0)+B(1,0)*A(2,1)+B(0,0)*A(2,0);  M(2,1)=B(2,1)*A(2,2)+B(1,1)*A(2,1)+B(0,1)*A(2,0);  M(2,2)=A(2,2)*B(2,2)+B(1,2)*A(2,1)+B(0,2)*A(2,0);
}

/** Clear vector. */
inline void set_to_zero (Vec3_t & V)
{
    V = 0.0, 0.0, 0.0;
}

/** Clear matrix. */
inline void set_to_zero (Mat3_t & M)
{
    M = 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0;
}

// Constants
namespace OrthoSys
{
    Vec3_t O;        ///< Origin
    Vec3_t e0,e1,e2; ///< Basis
    Mat3_t I;        ///< Identity

    int __init_ortho_sys()
    {
        O  = 0.0, 0.0, 0.0;
        e0 = 1.0, 0.0, 0.0;
        e1 = 0.0, 1.0, 0.0;
        e2 = 0.0, 0.0, 1.0;
        I  = 1.0, 0.0, 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, 1.0;
        return 0.0;
    }

    int __dummy_init_ortho_sys = __init_ortho_sys();
}

#ifdef USE_BOOST_PYTHON
inline Vec3_t Tup2Vec3 (BPy::tuple const & T3)
{
    return Vec3_t(BPy::extract<double>(T3[0])(), BPy::extract<double>(T3[1])(), BPy::extract<double>(T3[2])());
}
#endif


///////////////////////////////////////////////////////////////////////////////////////////// Tensors ////////////


/** Deviator of 2nd order symmetric tensor Ten. */
inline void Dev (Vec_t const & Ten, Vec_t & DevTen)
{
    double coef = (Ten(0)+Ten(1)+Ten(2))/3.0;
    DevTen     = Ten;
    DevTen(0) -= coef;
    DevTen(1) -= coef;
    DevTen(2) -= coef;
}

/** Trace of 2nd order symmetric tensor Ten. */
inline double Tra (Vec_t const & Ten)
{
    return Ten(0)+Ten(1)+Ten(2);
}

/** Creates the matrix representation of 2nd order symmetric tensor Ten (Mandel's representation). */
inline void Ten2Mat (Vec_t const & Ten, Mat3_t & Mat)
{
    // matrix of tensor
    size_t ncp = size(Ten);
    if (ncp==4)
    {
        Mat = Ten(0),            Ten(3)/Util::SQ2,     0.0,
              Ten(3)/Util::SQ2,  Ten(1),               0.0,
                           0.0,               0.0,  Ten(2);
    }
    else if (ncp==6)
    {
        Mat = Ten(0),            Ten(3)/Util::SQ2,  Ten(5)/Util::SQ2,
              Ten(3)/Util::SQ2,  Ten(1),            Ten(4)/Util::SQ2,
              Ten(5)/Util::SQ2,  Ten(4)/Util::SQ2,  Ten(2);
        
    }
    else throw new Fatal("matvec.h::Ten2Mat: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Creates a 2nd order symmetric tensor Ten (using Mandel's representation) based on its matrix components. NOTE: This method does not check if matrix is symmetric. */
inline void Mat2Ten (Mat3_t const & Mat, Vec_t & Ten, size_t NumComponents)
{
    // tensor from matrix
    Ten.change_dim (NumComponents);
    if (NumComponents==4)
    {
        Ten = Mat(0,0), Mat(1,1), Mat(2,2), Util::SQ2*Mat(0,1);
    }
    else if (NumComponents==6)
    {
        Ten = Mat(0,0), Mat(1,1), Mat(2,2), Util::SQ2*Mat(0,1), Util::SQ2*Mat(1,2), Util::SQ2*Mat(2,0);
    }
    else throw new Fatal("matvec.h::Mat2Ten: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** 2nd order symmetric tensor Ten raised to the power of 2. */
inline void Pow2 (Vec_t const & Ten, Vec_t & Ten2)
{
    size_t ncp = size(Ten);
    Ten2.change_dim (ncp);
    if (ncp==4)
    {
        Ten2 = Ten(0)*Ten(0) + Ten(3)*Ten(3)/2.0,
               Ten(1)*Ten(1) + Ten(3)*Ten(3)/2.0,
               Ten(2)*Ten(2),
               Ten(0)*Ten(3) + Ten(1)*Ten(3);
    }
    else if (ncp==6)
    {
        Ten2 = Ten(0)*Ten(0)           + Ten(3)*Ten(3)/2.0       + Ten(5)*Ten(5)/2.0,
               Ten(3)*Ten(3)/2.0       + Ten(1)*Ten(1)           + Ten(4)*Ten(4)/2.0,
               Ten(5)*Ten(5)/2.0       + Ten(4)*Ten(4)/2.0       + Ten(2)*Ten(2),
               Ten(0)*Ten(3)           + Ten(3)*Ten(1)           + Ten(5)*Ten(4)/Util::SQ2,
               Ten(3)*Ten(5)/Util::SQ2 + Ten(1)*Ten(4)           + Ten(4)*Ten(2),
               Ten(0)*Ten(5)           + Ten(3)*Ten(4)/Util::SQ2 + Ten(5)*Ten(2);
    }
    else throw new Fatal("matvec.h::Pow2: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Determinant of 2nd order symmetric tensor Ten. */
inline double Det (Vec_t const & Ten)
{
    size_t ncp = size(Ten);
    if (ncp==4)
    {
        return   Ten(0)*Ten(1)*Ten(2)
               - Ten(2)*Ten(3)*Ten(3)/2.0;
    }
    else if (ncp==6)
    {
        return   Ten(0)*Ten(1)*Ten(2) 
               + Ten(3)*Ten(4)*Ten(5)/Util::SQ2
               - Ten(0)*Ten(4)*Ten(4)/2.0
               - Ten(1)*Ten(5)*Ten(5)/2.0
               - Ten(2)*Ten(3)*Ten(3)/2.0;
    }
    else throw new Fatal("matvec.h::Det: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Inverse of 2nd order symmetric tensor T. */
inline void Inv (Vec_t const & T, Vec_t & Ti, double Tol=1.0e-10)
{
    size_t ncp = size(T);
    Ti.change_dim (ncp);
    if (ncp==4)
    {
        double det =  T(0)*T(1)*T(2)
                    - T(2)*T(3)*T(3)/2.0;

        if (fabs(det)<Tol)
        {
            std::ostringstream oss;  oss<<T;
            throw new Fatal("matvec.h::Inv: inverse of 2nd order symmetric tensor failed with null (%g) determinat. T =\n%s",Tol,oss.str().c_str());
        }

        Ti(0) = T(1)*T(2)/det;
        Ti(1) = T(0)*T(2)/det;
        Ti(2) = (T(0)*T(1)-2.0*T(3)*T(3))/det;
        Ti(3) = -2.0*T(2)*T(3)/det;
    }
    else if (ncp==6)
    {
        double det =  T(0)*T(1)*T(2) 
                    + T(3)*T(4)*T(5)/Util::SQ2
                    - T(0)*T(4)*T(4)/2.0
                    - T(1)*T(5)*T(5)/2.0
                    - T(2)*T(3)*T(3)/2.0;

        if (fabs(det)<Tol)
        {
            std::ostringstream oss;  oss<<T;
            throw new Fatal("matvec.h::Inv: inverse of 2nd order symmetric tensor failed with null (%g) determinat. T =\n%s",Tol,oss.str().c_str());
        }

        Ti(0) = (T(1)*T(2)-T(4)*T(4)/2.0)/det;
        Ti(1) = (T(0)*T(2)-T(5)*T(5)/2.0)/det;
        Ti(2) = (T(0)*T(1)-T(3)*T(3)/2.0)/det;
        Ti(3) = ((T(4)*T(5))/Util::SQ2-T(2)*T(3))/det;
        Ti(4) = ((T(3)*T(5))/Util::SQ2-T(0)*T(4))/det;
        Ti(5) = ((T(3)*T(4))/Util::SQ2-T(1)*T(5))/det;
    }
    else throw new Fatal("matvec.h::Inv: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Characteristic invariants of 2nd order symmetric tensor Ten. */
inline void CharInvs (Vec_t const & Ten, double & I1, double & I2, double & I3)
{
    I1 = Ten(0) + Ten(1) + Ten(2);
    I2 = Ten(0)*Ten(1) + Ten(1)*Ten(2) + Ten(2)*Ten(0) - Ten(3)*Ten(3)/2.0;
    I3 = Ten(0)*Ten(1)*Ten(2) - Ten(2)*Ten(3)*Ten(3)/2.0;
    size_t ncp = size(Ten);
    if (ncp>4)
    {
        I2 += (-Ten(4)*Ten(4)/2.0 - Ten(5)*Ten(5)/2.0);
        I3 += (Ten(3)*Ten(4)*Ten(5)/Util::SQ2 - Ten(0)*Ten(4)*Ten(4)/2.0 - Ten(1)*Ten(5)*Ten(5)/2.0);
    }
}

/** Characteristic invariants of 2nd order symmetric tensor Ten and their derivatives. */
inline void CharInvs (Vec_t const & Ten, double & I1, double & I2, double & I3, Vec_t & dI1dTen, Vec_t & dI2dTen, Vec_t & dI3dTen)
{
    I1 = Ten(0) + Ten(1) + Ten(2);
    I2 = Ten(0)*Ten(1) + Ten(1)*Ten(2) + Ten(2)*Ten(0) - Ten(3)*Ten(3)/2.0;
    I3 = Ten(0)*Ten(1)*Ten(2) - Ten(2)*Ten(3)*Ten(3)/2.0;
    size_t ncp = size(Ten);
    dI1dTen.change_dim (ncp);
    dI2dTen.change_dim (ncp);
    dI3dTen.change_dim (ncp);
    if (ncp>4)
    {
        I2 += (-Ten(4)*Ten(4)/2.0 - Ten(5)*Ten(5)/2.0);
        I3 += (Ten(3)*Ten(4)*Ten(5)/Util::SQ2 - Ten(0)*Ten(4)*Ten(4)/2.0 - Ten(1)*Ten(5)*Ten(5)/2.0);
        dI1dTen = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
    }
    else dI1dTen = 1.0, 1.0, 1.0, 0.0;
    Vec_t Ten2(ncp);
    Pow2 (Ten, Ten2);
    dI2dTen = I1*dI1dTen - Ten;
    dI3dTen = Ten2 - I1*Ten + I2*dI1dTen;
}

/** Eigenprojectors of 2nd order symmetric tensor Ten. */
inline void EigenProj (Vec_t const & Ten, Vec3_t & L, Vec_t & P0, Vec_t & P1, Vec_t & P2, bool SortAsc=false, bool SortDesc=false)
{
    // matrix of tensor
    Mat3_t ten;
    Ten2Mat (Ten, ten);

    // eigen-values and vectors
    Vec3_t v0,v1,v2;
    Eig (ten, L, v0, v1, v2, SortAsc, SortDesc);

    // eigen-projectors
    size_t ncp = size(Ten);
    P0.change_dim (ncp);
    P1.change_dim (ncp);
    P2.change_dim (ncp);
    if (ncp==4)
    {
        P0 = v0(0)*v0(0), v0(1)*v0(1), v0(2)*v0(2), v0(0)*v0(1)*Util::SQ2;
        P1 = v1(0)*v1(0), v1(1)*v1(1), v1(2)*v1(2), v1(0)*v1(1)*Util::SQ2;
        P2 = v2(0)*v2(0), v2(1)*v2(1), v2(2)*v2(2), v2(0)*v2(1)*Util::SQ2;
    }
    else
    {
        P0 = v0(0)*v0(0), v0(1)*v0(1), v0(2)*v0(2), v0(0)*v0(1)*Util::SQ2, v0(1)*v0(2)*Util::SQ2, v0(2)*v0(0)*Util::SQ2;
        P1 = v1(0)*v1(0), v1(1)*v1(1), v1(2)*v1(2), v1(0)*v1(1)*Util::SQ2, v1(1)*v1(2)*Util::SQ2, v1(2)*v1(0)*Util::SQ2;
        P2 = v2(0)*v2(0), v2(1)*v2(1), v2(2)*v2(2), v2(0)*v2(1)*Util::SQ2, v2(1)*v2(2)*Util::SQ2, v2(2)*v2(0)*Util::SQ2;
    }
}

/** Eigenprojectors of 2nd order symmetric tensor Ten. */
inline void EigenProjAnalytic (Vec_t const & Ten, Vec3_t & L, Vec_t & P0, Vec_t & P1, Vec_t & P2)
{
    // identity tensor
    size_t ncp = size(Ten);
    Vec_t I(ncp);
    set_to_zero(I);
    I(0)=1.0;  I(1)=1.0;  I(2)=1.0;

    // characteristics invariants
    double I1,I2,I3;
    CharInvs (Ten,I1,I2,I3);

    // inverse tensor
    Vec_t Ti;
    Inv (Ten, Ti);

    // eigen-values/projectors
    P0.change_dim (ncp);
    P1.change_dim (ncp);
    P2.change_dim (ncp);
    Vec_t * P[3] = {&P0, &P1, &P2}; // all three eigen projectors
    double alpha = 2.0*sqrt(I1*I1-3.0*I2);
    double theta = acos((2.0*I1*I1*I1-9.0*I1*I2+27.0*I3)/(2.0*pow(I1*I1-3.0*I2,3.0/2.0)));
    //std::cout << "I1 = " << I1 << ",  I2 = " << I2 << ",  I3 = " << I3 << ",  I1*I1-3*I2 = " << I1*I1-3.0*I2 << std::endl;
    //std::cout << "alpha = " << alpha << ",  theta = " << theta << std::endl;
    for (size_t k=0; k<3; ++k)
    {
        L(k) = (I1+alpha*cos((theta+2.0*Util::PI*(1.0+k))/3.0))/3.0;
        if (fabs(L(k)<1.0e-14)) throw new Fatal("matvec.h:EigenProjAnalytic: L(%d)=%g must be non-zero",k,L(k));
        double coef = L(k)/(2.0*L(k)*L(k)-L(k)*I1+I3/L(k));
        (*P[k]) = coef*(Ten+(L(k)-I1)*I+(I3/L(k))*Ti);
    }
}

/** Initialize second order identity tensor. */
inline void Calc_I (size_t NCp, Vec_t & I)
{
    I.change_dim(NCp);
    if      (NCp==4) I = 1.0, 1.0, 1.0, 0.0;
    else if (NCp==6) I = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
    else throw new Fatal("matvec.h::Calc_I: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Initialize fourth order identity tensor. */
inline void Calc_IIsym (size_t NCp, Mat_t & IIsym)
{
    IIsym.change_dim(NCp,NCp);
    if (NCp==4)
    {
        IIsym = 1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 1.0;
    }
    else if (NCp==6)
    {
        IIsym = 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    }
    else throw new Fatal("matvec.h::Calc_IIsym: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Initialize IdyI fourth order tensor. */
inline void Calc_IdyI (size_t NCp, Mat_t & IdyI)
{
    IdyI.change_dim(NCp,NCp);
    if (NCp==4)
    {
        IdyI = 1.0, 1.0, 1.0, 0.0,
               1.0, 1.0, 1.0, 0.0,
               1.0, 1.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 0.0;
    }
    else if (NCp==6)
    {
        IdyI = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
               1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
               1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    }
    else throw new Fatal("matvec.h::Calc_IdyI: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Initialize fourth order symmetric-deviatoric tensor. */
inline void Calc_Psd (size_t NCp, Mat_t & Psd)
{
    Psd.change_dim(NCp,NCp);
    if (NCp==4)
    {
        Psd =  2.0/3.0, -1.0/3.0, -1.0/3.0, 0.0,
              -1.0/3.0,  2.0/3.0, -1.0/3.0, 0.0,
              -1.0/3.0, -1.0/3.0,  2.0/3.0, 0.0,
                   0.0,      0.0,      0.0, 1.0;
    }
    else if (NCp==6)
    {
        Psd =  2.0/3.0, -1.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0,
              -1.0/3.0,  2.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0,
              -1.0/3.0, -1.0/3.0,  2.0/3.0, 0.0, 0.0, 0.0,
                   0.0,      0.0,      0.0, 1.0, 0.0, 0.0,
                   0.0,      0.0,      0.0, 0.0, 1.0, 0.0,
                   0.0,      0.0,      0.0, 0.0, 0.0, 1.0;
    }
    else throw new Fatal("matvec.h::Calc_Psd: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Initialize fourth order isotropic tensor. */
inline void Calc_Piso (size_t NCp, Mat_t & Piso)
{
    Piso.change_dim(NCp,NCp);
    if (NCp==4)
    {
        Piso = 1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0,
               1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0,
               1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0,
                   0.0,     0.0,     0.0, 0.0;
    }
    else if (NCp==6)
    {
        Piso = 1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0,
               1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0,
               1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0,
                   0.0,     0.0,     0.0, 0.0, 0.0, 0.0,
                   0.0,     0.0,     0.0, 0.0, 0.0, 0.0,
                   0.0,     0.0,     0.0, 0.0, 0.0, 0.0;
    }
    else throw new Fatal("matvec.h::Calc_Piso: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}


////////////////////////////////////////////////////////////////////////////////////////// Invariants ////////////


// Cambridge invariants
inline double Calc_pcam  (Vec_t const & Sig) { return -(Sig(0)+Sig(1)+Sig(2))/3.0; }
inline double Calc_ev    (Vec_t const & Eps) { return   Eps(0)+Eps(1)+Eps(2);      }
inline double Calc_qcam  (Vec_t const & Sig) { double m = (size(Sig)>4 ? pow(Sig(4),2.0)+pow(Sig(5),2.0) : 0.0); return sqrt(pow(Sig(0)-Sig(1),2.0) + pow(Sig(1)-Sig(2),2.0) + pow(Sig(2)-Sig(0),2.0) + 3.0*(pow(Sig(3),2.0)+m))/sqrt(2.0); }
inline double Calc_ed    (Vec_t const & Eps) { double m = (size(Eps)>4 ? pow(Eps(4),2.0)+pow(Eps(5),2.0) : 0.0); return sqrt(pow(Eps(0)-Eps(1),2.0) + pow(Eps(1)-Eps(2),2.0) + pow(Eps(2)-Eps(0),2.0) + 3.0*(pow(Eps(3),2.0)+m))*(sqrt(2.0)/3.0); }

// Octahedral invariants
inline double Calc_poct  (Vec_t const & Sig) { return -(Sig(0)+Sig(1)+Sig(2))/sqrt(3.0); }
inline double Calc_evoct (Vec_t const & Eps) { return  (Eps(0)+Eps(1)+Eps(2))/sqrt(3.0); }
inline double Calc_qoct  (Vec_t const & Sig) { double m = (size(Sig)>4 ? pow(Sig(4),2.0)+pow(Sig(5),2.0) : 0.0); return sqrt(pow(Sig(0)-Sig(1),2.0) + pow(Sig(1)-Sig(2),2.0) + pow(Sig(2)-Sig(0),2.0) + 3.0*(pow(Sig(3),2.0)+m))/sqrt(3.0); }
inline double Calc_edoct (Vec_t const & Eps) { double m = (size(Eps)>4 ? pow(Eps(4),2.0)+pow(Eps(5),2.0) : 0.0); return sqrt(pow(Eps(0)-Eps(1),2.0) + pow(Eps(1)-Eps(2),2.0) + pow(Eps(2)-Eps(0),2.0) + 3.0*(pow(Eps(3),2.0)+m))/sqrt(3.0); }


/** Octahedral invariants of Sig. */
inline void OctInvs (Vec_t const & Sig, double & p, double & q, double & t, double qTol=1.0e-8)
{
    size_t ncp = size(Sig);
    Vec_t s;
    Dev (Sig, s);
    double m = (ncp>4 ? pow(Sig(4),2.0)+pow(Sig(5),2.0) : 0.0);
    p = -(Sig(0)+Sig(1)+Sig(2))/Util::SQ3;
    q = sqrt(pow(Sig(0)-Sig(1),2.0) + pow(Sig(1)-Sig(2),2.0) + pow(Sig(2)-Sig(0),2.0) + 3.0*(pow(Sig(3),2.0)+m))/sqrt(3.0);
    t = 0.0;
    if (q>qTol)
    {
        double det_s = Det(s);
        t = -3.0*Util::SQ6*det_s/(q*q*q);
        if (t<=-1.0) t = -1.0;
        if (t>= 1.0) t =  1.0;
    }
}

/** Octahedral invariants of Sig (L = principal values). */
inline void OctInvs (Vec3_t const & L, double & p, double & q, double & t, double qTol=1.0e-8)
{
    p = -(L(0)+L(1)+L(2))/Util::SQ3;
    q = sqrt(pow(L(0)-L(1),2.0) + pow(L(1)-L(2),2.0) + pow(L(2)-L(0),2.0))/sqrt(3.0);
    t = 0.0;
    if (q>qTol)
    {
        Vec3_t S((2.*L(0)-L(1)-L(2))/3.,
                 (2.*L(1)-L(2)-L(0))/3.,
                 (2.*L(2)-L(0)-L(1))/3.);
        t = -3.0*Util::SQ6*S(0)*S(1)*S(2)/pow(q,3.0);
        if (t<=-1.0) t = -1.0;
        if (t>= 1.0) t =  1.0;
    }
}

/** Octahedral invariants of Sig (L = principal values). With derivatives */
inline void OctInvs (Vec3_t const & L, double & p, double & q, double & t, Vec3_t & dpdL, Vec3_t & dqdL, Vec3_t & dtdL, double qTol=1.0e-8)
{
    Vec3_t one(1.0,1.0,1.0), s;
    p    = -(L(0)+L(1)+L(2))/Util::SQ3;
    q    = sqrt(pow(L(0)-L(1),2.0) + pow(L(1)-L(2),2.0) + pow(L(2)-L(0),2.0))/sqrt(3.0);
    t    = 0.0;
    s    = L - ((L(0)+L(1)+L(2))/3.0)*one;
    dpdL = (-1.0/Util::SQ3)*one;
    dqdL = 0.0, 0.0, 0.0;
    if (q>qTol)
    {
        double q3 = q*q*q;
        double q5 = q3*q*q;
        double l  = (L(0)-L(1))*(L(1)-L(2))*(L(2)-L(0));
        Vec3_t B(L(2)-L(1), L(0)-L(2), L(1)-L(0));
        t    = -3.0*Util::SQ6*s(0)*s(1)*s(2)/q3;
        dqdL = (1.0/q)*s;
        dtdL = (-Util::SQ6*l/q5)*B;
        if (t<=-1.0) t = -1.0;
        if (t>= 1.0) t =  1.0;
    }
}

/** Calc principal values given octahedral invariants. (-pi <= th(rad) <= pi) */
inline void pqth2L (double p, double q, double th, Vec3_t & L, char const * Type="oct")
{
    if (strcmp(Type,"cam")==0)
    {
        L(0) = -p + 2.0*q*sin(th-2.0*Util::PI/3.0)/3.0;
        L(1) = -p + 2.0*q*sin(th)                 /3.0;
        L(2) = -p + 2.0*q*sin(th+2.0*Util::PI/3.0)/3.0;
    }
    else if (strcmp(Type,"oct")==0) // oct
    {
        L(0) = -p/Util::SQ3 + 2.0*q*sin(th-2.0*Util::PI/3.0)/Util::SQ6;
        L(1) = -p/Util::SQ3 + 2.0*q*sin(th)                 /Util::SQ6;
        L(2) = -p/Util::SQ3 + 2.0*q*sin(th+2.0*Util::PI/3.0)/Util::SQ6;
    }
    else throw new Fatal("pqTh2L: Method is not available for invariant Type==%s",Type);
}

inline void OctDerivs (Vec3_t const & L, double & p, double & q, double & t, Mat3_t & dpqthdL, double qTol=1.0e-8)
{
    OctInvs (L, p,q,t, qTol);
    if (q>qTol)
    {

        Vec3_t S((2.*L(0)-L(1)-L(2))/3., 
                 (2.*L(1)-L(2)-L(0))/3.,
                 (2.*L(2)-L(0)-L(1))/3.);
        Vec3_t devSS((2.*S(1)*S(2) - S(2)*S(0) - S(0)*S(1))/3.,
                     (2.*S(2)*S(0) - S(1)*S(2) - S(0)*S(1))/3.,
                     (2.*S(0)*S(1) - S(1)*S(2) - S(2)*S(0))/3.);
        double th = asin(t)/3.0;
        double c  = -1.0/(pow(q,3.0)*cos(3.0*th));
        dpqthdL = -1./Util::SQ3,                   -1./Util::SQ3,                   -1./Util::SQ3,
                  S(0)/q,                          S(1)/q,                          S(2)/q,
                  c*Util::SQ6*devSS(0)+c*q*t*S(0), c*Util::SQ6*devSS(1)+c*q*t*S(1), c*Util::SQ6*devSS(2)+c*q*t*S(2);
    }
    else
    {
        dpqthdL = -1./Util::SQ3, -1./Util::SQ3, -1./Util::SQ3,
                   0.,            0.,            0.,
                   0.,            0.,            0.;
    }
}

inline void InvOctDerivs (Vec3_t const & Lsorted, double & p, double & q, double & t, Mat3_t & dLdpqth, double qTol=1.0e-8)
{
    OctInvs (Lsorted, p,q,t, qTol);
    double th = asin(t)/3.0;
    dLdpqth = -1./Util::SQ3, (2.*sin(th-(2.*Util::PI)/3.))/Util::SQ6, (2.*q*cos(th-(2.*Util::PI)/3.))/Util::SQ6,
              -1./Util::SQ3, (2.*sin(th)                 )/Util::SQ6, (2.*q*cos(th)                 )/Util::SQ6,
              -1./Util::SQ3, (2.*sin(th+(2.*Util::PI)/3.))/Util::SQ6, (2.*q*cos(th+(2.*Util::PI)/3.))/Util::SQ6;
}

#ifdef USE_BOOST_PYTHON
inline BPy::tuple Pypqth2L (double p, double q, double th, BPy::str const & Type)
{
    Vec3_t l;
    pqth2L (p, q, th, l, BPy::extract<char const *>(Type)());
    return BPy::make_tuple (l(0), l(1), l(2));
}
#endif


///////////////////////////////////////////////////////////////////////// Failure criteria constants /////////////


/** Calculate M = max q/p at compression (phi: friction angle at compression (degrees)). */
inline double Phi2M (double Phi, char const * Type="oct")
{
    double sphi = sin(Phi*Util::PI/180.0);
    if      (strcmp(Type,"oct")==0) return 2.0*Util::SQ2*sphi/(3.0-sphi);
    else if (strcmp(Type,"cam")==0) return 6.0*sphi/(3.0-sphi);
    else if (strcmp(Type,"smp")==0)
    {
        double eta = 2.0*Util::SQ2*sphi/(3.0-sphi);
        double c   = sqrt((2.0+Util::SQ2*eta-2.0*eta*eta)/(3.0*Util::SQ3*(Util::SQ2*eta+2.0)));
        double a   = sqrt((2.0*eta+Util::SQ2)/Util::SQ6);
        double b   = sqrt((Util::SQ2-eta)/Util::SQ6);
        return sqrt((eta*eta+1.0)/(c*c*pow(a+2.0*b,2.0))-1.0);
    }
    else throw new Fatal("Phi2M: Method is not available for invariant Type==%s",Type);
}

/** Calculate phi (friction angle at compression (degrees)) given M (max q/p at compression). */
inline double M2Phi (double M, char const * Type="oct")
{
    double sphi;
    if      (strcmp(Type,"oct")==0) sphi = 3.0*M/(M+2.0*Util::SQ2);
    else if (strcmp(Type,"cam")==0) sphi = 3.0*M/(M+6.0);
    else throw new Fatal("M2Phi: Method is not available for invariant Type==%s",Type);
    return asin(sphi)*180.0/Util::PI;
}

#endif // MECHSYS_MATVEC_H
