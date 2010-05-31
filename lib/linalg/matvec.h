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
#include <cmath>     // for sqrt, pow
#include <algorithm> // for min, max

// Boost
#include <boost/numeric/mtl/mtl.hpp>

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


/** Dense matrix (general). */
typedef mtl::dense2D<double, mtl::matrix::parameters<mtl::tag::col_major> > Mat_t;

/** Dense vector (general). */
typedef mtl::dense_vector<double> Vec_t;

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

/** Left multiplication. {B} = {A}*[M]  */
inline void Mult (Vec_t const & A, Mat_t const & M, Vec_t & B)
{
    B.change_dim (M.num_cols());
    set_to_zero  (B);
    for (size_t i=0; i<M.num_rows(); ++i)
    for (size_t j=0; j<M.num_cols(); ++j)
        B(j) += A(i)*M(i,j);
}


////////////////////////////////////////////////////////////////////////////////////////////// Tiny MatVec /////


/** 3x3 Matrix. */
typedef blitz::TinyMatrix<double,3,3> Mat3_t;

/** 3x1 Vector. */
typedef blitz::TinyVector<double,3> Vec3_t;
typedef blitz::TinyVector<bool,3>  bVec3_t;

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

/** Determinant.*/
inline double Det (Mat3_t const & M)
{
    double det =   M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
                 - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
                 + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));
    return det;
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

/** Eigenvalues and eigenvectors. NOTE: This function changes the matrix M. */
inline void Eig (Mat3_t & M, Vec3_t & L, Vec3_t & V0, Vec3_t & V1, Vec3_t & V2)
{
    // calculate
    gsl_matrix_view m = gsl_matrix_view_array (M.data(), 3, 3);
    gsl_vector * eval = gsl_vector_alloc      (3);
    gsl_matrix * evec = gsl_matrix_alloc      (3, 3);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
    gsl_eigen_symmv (&m.matrix, eval, evec, w);

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
    for (size_t i=0; i<3; ++i)
    for (size_t j=0; j<3; ++j)
        M(i,j) = A(i) * B(j);
}

/** Compare two vectors. */
inline double CompareVectors (Vec3_t const & A, Vec3_t const & B)
{
    double error = 0.0;
    for (size_t i=0; i<3; ++i)
        error += fabs(A(i)-B(i));
    return error;
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


///////////////////////////////////////////////////////////////////////////////////////////// Tensors ////////////

// Cambridge invariants
inline double Calc_pcam  (Vec_t const & Sig) { return -(Sig(0)+Sig(1)+Sig(2))/3.0; }
inline double Calc_ev    (Vec_t const & Eps) { return   Eps(0)+Eps(1)+Eps(2);      }
inline double Calc_qcam  (Vec_t const & Sig) { double m = (size(Sig)>4 ? pow(Sig(4),2.0)+pow(Sig(5),2.0) : 0.0); return sqrt(pow(Sig(0)-Sig(1),2.0) + pow(Sig(1)-Sig(2),2.0) + pow(Sig(2)-Sig(0),2.0) + 3.0*(pow(Sig(3),2.0)+m))/sqrt(2.0); }
inline double Calc_ed    (Vec_t const & Eps) { double m = (size(Eps)>4 ? pow(Eps(4),2.0)+pow(Eps(5),2.0) : 0.0); return sqrt(pow(Eps(0)-Eps(1),2.0) + pow(Eps(1)-Eps(2),2.0) + pow(Eps(2)-Eps(0),2.0) + 3.0*(pow(Eps(3),2.0)+m))*(sqrt(2.0)/3.0); }

// Octahedral invariants
inline double Calc_poct  (Vec_t const & Sig) { return -(Sig(0)+Sig(1)+Sig(2))/sqrt(3.0); }
inline double Calc_qoct  (Vec_t const & Sig) { double m = (size(Sig)>4 ? pow(Sig(4),2.0)+pow(Sig(5),2.0) : 0.0); return sqrt(pow(Sig(0)-Sig(1),2.0) + pow(Sig(1)-Sig(2),2.0) + pow(Sig(2)-Sig(0),2.0) + 3.0*(pow(Sig(3),2.0)+m))/sqrt(3.0); }

/** Deviator of Sig. */
inline void Dev (Vec_t const & Sig, Vec_t & DevSig)
{
    double coef = (Sig(0)+Sig(1)+Sig(2))/3.0;
    DevSig     = Sig;
    DevSig(0) -= coef;
    DevSig(1) -= coef;
    DevSig(2) -= coef;
}

/** Trace of Sig. */
inline double Tra (Vec_t const & Sig)
{
    return Sig(0)+Sig(1)+Sig(2);
}

/** Eigenprojectors of Sig. */
inline void EigenProj (Vec_t const & Sig, Vec3_t & L, Vec_t & P0, Vec_t & P1, Vec_t & P2)
{
    size_t ncp = size(Sig);
    Mat3_t sig;
    if (ncp==4)
    {
        sig = Sig(0),            Sig(3)/Util::SQ2,     0.0,
              Sig(3)/Util::SQ2,  Sig(1),               0.0,
                           0.0,               0.0,  Sig(2);
    }
    else
    {
        sig = Sig(0),            Sig(3)/Util::SQ2,  Sig(5)/Util::SQ2,
              Sig(3)/Util::SQ2,  Sig(1),            Sig(4)/Util::SQ2,
              Sig(5)/Util::SQ2,  Sig(4)/Util::SQ2,  Sig(2);
        
    }
    Vec3_t v0,v1,v2;
    Eig (sig, L, v0, v1, v2);
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
        double det_s = s(0)*s(1)*s(2) - s(2)*s(3)*s(3)/2.0; // TODO: extend to 3D
        t = -3.0*Util::SQ6*det_s/(q*q*q);
        if (t<=-1.0) t = -1.0;
        if (t>= 1.0) t =  1.0;
    }
}

/** Octahedral invariants of Sig. */
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

/** Sig raised to the power of 2. */
inline void Pow2 (Vec_t const & Sig, Vec_t & Sig2)
{
    size_t ncp = size(Sig);
    Sig2.change_dim (ncp);
    if (ncp==4)
    {
        Sig2 = Sig(0)*Sig(0) + Sig(3)*Sig(3)/2.0,
               Sig(1)*Sig(1) + Sig(3)*Sig(3)/2.0,
               Sig(2)*Sig(2),
               Sig(0)*Sig(3) + Sig(1)*Sig(3);
    }
    else
    {
        Sig2 = Sig(0)*Sig(0)           + Sig(3)*Sig(3)/2.0       + Sig(5)*Sig(5)/2.0,
               Sig(3)*Sig(3)/2.0       + Sig(1)*Sig(1)           + Sig(4)*Sig(4)/2.0,
               Sig(5)*Sig(5)/2.0       + Sig(4)*Sig(4)/2.0       + Sig(2)*Sig(2),
               Sig(0)*Sig(3)           + Sig(3)*Sig(1)           + Sig(5)*Sig(4)/Util::SQ2,
               Sig(3)*Sig(5)/Util::SQ2 + Sig(1)*Sig(4)           + Sig(4)*Sig(2),
               Sig(0)*Sig(5)           + Sig(3)*Sig(4)/Util::SQ2 + Sig(5)*Sig(2);
    }
}

/** Characteristic invariants of Sig. */
inline void CharInvs (Vec_t const & Sig, double & I1, double & I2, double & I3, Vec_t & dI1dSig, Vec_t & dI2dSig, Vec_t & dI3dSig)
{
    I1 = Sig(0) + Sig(1) + Sig(2);
    I2 = Sig(0)*Sig(1) + Sig(1)*Sig(2) + Sig(2)*Sig(0) - Sig(3)*Sig(3)/2.0;
    I3 = Sig(0)*Sig(1)*Sig(2) - Sig(2)*Sig(3)*Sig(3)/2.0;
    size_t ncp = size(Sig);
    dI1dSig.change_dim (ncp);
    dI2dSig.change_dim (ncp);
    dI3dSig.change_dim (ncp);
    if (ncp>4)
    {
        I2 += (-Sig(4)*Sig(4)/2.0 - Sig(5)*Sig(5)/2.0);
        I3 += (Sig(3)*Sig(4)*Sig(5)/Util::SQ2 - Sig(0)*Sig(4)*Sig(4)/2.0 - Sig(1)*Sig(5)*Sig(5)/2.0);
        dI1dSig = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
    }
    else dI1dSig = 1.0, 1.0, 1.0, 0.0;
    Vec_t Sig2(ncp);
    Pow2 (Sig, Sig2);
    dI2dSig = I1*dI1dSig - Sig;
    dI3dSig = Sig2 - I1*Sig + I2*dI1dSig;
}

/** Calc principal values given octahedral invariants. */
inline void pqt2L (double p, double q, double t, Vec3_t & L, char const * Type="oct")
{
    double ttemp = (t<=-1.0 ? -1.0 : (t>=1.0 ? 1.0 : t));
    double th    = asin(ttemp)/3.0;
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
    else throw new Fatal("pqt2L: Method is not available for invariant Type==%s",Type);
}

/** Calc principal values given octahedral invariants. */
inline void pqTh2L (double p, double q, double Th, Vec3_t & L, char const * Type="oct")
{
    double t = sin(3.0*Th*M_PI/180.0);
    pqt2L (p, q, t, L, Type);
}

#ifdef USE_BOOST_PYTHON
inline Vec3_t Tup2Vec3 (BPy::tuple const & T3)
{
    return Vec3_t(BPy::extract<double>(T3[0])(), BPy::extract<double>(T3[1])(), BPy::extract<double>(T3[2])());
}

inline BPy::tuple Pypqt2L (double p, double q, double t, BPy::str const & Type)
{
    Vec3_t l;
    pqt2L (p, q, t, l, BPy::extract<char const *>(Type)());
    return BPy::make_tuple (l(0), l(1), l(2));
}

inline BPy::tuple PypqTh2L (double p, double q, double Th, BPy::str const & Type)
{
    double t = sin(3.0*Th*M_PI/180.0);
    Vec3_t l;
    pqt2L (p, q, t, l, BPy::extract<char const *>(Type)());
    return BPy::make_tuple (l(0), l(1), l(2));
}
#endif

/** Calc octahedral coordinates given principal values. */
inline void OctCoords (Vec3_t const & L, double & sa, double & sb, double & sc)
{
    Vec3_t l(L);
    Util::Sort (l(0),l(1),l(2));
    sa = Util::SQ2*(l(1)-l(2))/2.0;
    sb = (l(2)+l(1)-2.0*l(0))/Util::SQ6;
    sc = -(l(0)+l(1)+l(2))/Util::SQ3;
}

/** Calc principal values given octahedral coordinates. */
inline void PrincVals (double sa, double sb, double sc, Vec3_t & L)
{
    L(0) =               - 2.0*sb/Util::SQ6 - sc/Util::SQ3;
    L(1) =  sa/Util::SQ2 +     sb/Util::SQ6 - sc/Util::SQ3;
    L(2) = -sa/Util::SQ2 +     sb/Util::SQ6 - sc/Util::SQ3;
}

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

/** Calculate q_failure for given undrained cohesion (cu) */
inline double cu2qf (double cu, char const * qType="oct", bool psa=false)
{
    double coef = (psa ? sqrt(3.0) : 2.0);
    if      (strcmp(qType,"oct")==0) return sqrt(2.0/3.0)*coef*cu;
    else if (strcmp(qType,"cam")==0) return coef*cu;
    else throw new Fatal("cu2qf: Method is not available for invariant qType==%s",qType);
}

/** Calculate undrained cohesion (cu) for given q_failure */
inline double qf2cu (double qf, char const * qType="oct", bool psa=false)
{
    double coef = (psa ? sqrt(3.0) : 2.0);
    if      (strcmp(qType,"oct")==0) return qf/(sqrt(2.0/3.0)*coef);
    else if (strcmp(qType,"cam")==0) return qf/coef;
    else throw new Fatal("qf2cu: Method is not available for invariant qType==%s",qType);
}


Mat_t Isy_2d (4,4); ///< Idendity of 4th order
Mat_t Psd_2d (4,4); ///< Projector: sym-dev
Mat_t IdI_2d (4,4); ///< 2nd order I dyadic 2nd order I

int __init_tensors()
{
    Isy_2d = 1.0;

    Psd_2d = 1.0;
    Psd_2d(0,0)= 2.0/3.0;  Psd_2d(0,1)=-1.0/3.0;  Psd_2d(0,2)=-1.0/3.0;
    Psd_2d(1,0)=-1.0/3.0;  Psd_2d(1,1)= 2.0/3.0;  Psd_2d(1,2)=-1.0/3.0;
    Psd_2d(2,0)=-1.0/3.0;  Psd_2d(2,1)=-1.0/3.0;  Psd_2d(2,2)= 2.0/3.0;

    IdI_2d = 0.0;
    IdI_2d(0,0)=1.0;  IdI_2d(0,1)=1.0;  IdI_2d(0,2)=1.0;
    IdI_2d(1,0)=1.0;  IdI_2d(1,1)=1.0;  IdI_2d(1,2)=1.0;
    IdI_2d(2,0)=1.0;  IdI_2d(2,1)=1.0;  IdI_2d(2,2)=1.0;

    return 0;
}

int __dummy_init_tensors = __init_tensors();

#endif // MECHSYS_MATVEC_H
