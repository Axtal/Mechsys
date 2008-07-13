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

#ifndef MECHSYS_SCALAPACK_H
#define MECHSYS_SCALAPACK_H

// STL
#ifdef DO_DEBUG
  #include <map>
#endif

// MechSys
#include "linalg/matrix.h"
#include "linalg/vector.h"
#include "util/exception.h"
#include "util/util.h"

extern "C"
{
	void   Cblacs_pinfo    (int* mypnum, int* nprocs);
	void   Cblacs_get      (int context, int request, int* value);
	int    Cblacs_gridinit (int* context, char * order, int np_row, int np_col);
	void   Cblacs_gridinfo (int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
	void   Cblacs_gridexit (int context);
	void   Cblacs_exit     (int error_code);
	void   Cblacs_barrier  (int context, char *scope);
	void   Cigebs2d        (int context, char *scope, char *top, int m, int n, int    *A, int lda);
	void   Cigebr2d        (int context, char *scope, char *top, int m, int n, int    *A, int lda, int rsrc, int csrc);
	void   Cdgebs2d        (int context, char *scope, char *top, int m, int n, double *A, int lda);
	void   Cdgebr2d        (int context, char *scope, char *top, int m, int n, double *A, int lda, int rsrc, int csrc);
	void   Cpdgemr2d       (int M, int N,
                            double *A, int IA, int JA, int *ADESC,
                            double *B, int IB, int JB, int *BDESC,
                            int CTXT);
	int    numroc_         (int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
	void   descinit_       (int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info);
	int    descset_        (int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld);
	double pdlamch_        (int *ictxt , char *cmach);
	double pdlange_        (char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);
	void   pdlacpy_        (char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);
	void   pdgesv_         (int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv, double *B, int *ib, int *jb, int *descb, int *info);
	void   pdgemm_         (char *TRANSA, char *TRANSB, int * M, int * N, int * K, double * ALPHA,
	                        double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB,
	                        double * BETA, double * C, int * IC, int * JC, int * DESCC );
	int    indxg2p_        (int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
	int    indxl2g_        (int *indxloc , int *nb, int *iproc, int *isrcproc, int *nprocs);
	void   pdgeadd_        (char *TRANS, int * M, int * N,
                            double * ALPHA,
                            double * A, int * IA, int * JA, int * DESCA,
                            double * BETA,
                            double * C, int * IC, int * JC, int * DESCC);
}

namespace ScaLAPACK
{

inline int IndxL2G(int idxloc, int nb, int iproc, int isrcproc, int nprocs)
{

	// (FORTRAN) indxl2g_ function uses indexes which starts at 1 (one)
	int fortran_idxloc = idxloc + 1;
	return indxl2g_(&fortran_idxloc, &nb, &iproc, &isrcproc, &nprocs) - 1;

}

inline void PDLAFill(LinAlg::Matrix<double> const & AA, double *A, int *DescA, int iRread, int iCread, double *Work)
{

	/*
	This function was based in a function inside file pdlaread.f
	That function (PDLAREAD) has the following caption:
	
	 -- ScaLAPACK auxiliary routine (version 2.0) --
	    University of Tennessee, Knoxville, Oak Ridge National Laboratory,
	    and University of California, Berkeley.
	    August 12, 2001 
	
	Purpose
	=======
	
	 PDLAREAD reads from a file named FILNAM a matrix and distribute
	 it to the process grid.
	
	 Only the process of coordinates {IRREAD, ICREAD} read the file.
	
	 WORK must be of size >= MB_ = DESCA( MB_ ).
	
	Further Details
	===============
	
	Contributed by Song Jin, University of Tennessee, 1996.
	*/
	// Parameters
	//int block_cyclic_2d = 1;
	//int dtype_          = 1;
	int ctxt_           = 2;
	//int m_              = 3;
	//int n_              = 4;
	int mb_             = 5;
	//int nb_             = 6;
	//int rsrc_           = 7;
	//int csrc_           = 8;
	//int lld_            = 9;
	int one = 1;

	// Local scalars
	bool isioprocessor = false;
	int  ictxt         = DescA[ctxt_-1];
	int  lwork         = DescA[mb_-1];
	int  mycol         = 0;
	int  myrow         = 0;
	int  npcol         = 0;
	int  nprow         = 0;

	// Check if this is the IO processor
	Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);
	isioprocessor = ((myrow==iRread)&&(mycol==iCread));

	// Get number of rows and columns
	int iwork[2];
	if (isioprocessor)
	{
		iwork[0] = AA.Rows();
		iwork[1] = AA.Cols();
		//igebs2d_(&ictxt, "All", " ", 2, 1, &iwork, 2);
		Cigebs2d(ictxt, "All", " ", 2, 1, iwork, 2);
	}
	else
	{
		//igebr2d_(&ictxt, "All", " ", 2, 1, &iwork, 2, &iRread, &iCread);
		Cigebr2d(ictxt, "All", " ", 2, 1, iwork, 2, iRread, iCread);
	}
	int m = iwork[0];
	int n = iwork[1];

	// DESCSET initializes a descriptor vector 
	int descwork[9];
	int mm   = Util::Max(1, Util::Min(m, lwork));
	int nn   = Util::Max(1, static_cast<int>(lwork/mm));
	int mb   = mm;
	int nb   = nn;
	int rsrc = iRread;
	int csrc = iCread;
	int ldd  = Util::Max(1, mm);
	descset_(descwork, &mm, &nn, &mb, &nb, &rsrc, &csrc, &ictxt, &ldd);

	// Fill matrix
	for (int jstart=0; jstart<n; jstart+=nn)
	{
		int jend  = Util::Min(n, jstart+nn);
		int jsize = jend - jstart;
		for (int istart=0; istart<m; istart+=mm)
		{
			int    iend  = Util::Min(m, istart+mm);
			int    isize = iend - istart;
			double alpha = 1.0;
			double beta  = 0.0;
			if (isioprocessor)
			{
				for (int j=0; j<jsize; j++)
				for (int i=0; i<isize; i++)
					Work[i+j*ldd] = AA(i,j);
			}
			pdgeadd_("N", &isize, &jsize,
			         &alpha,
			         Work, &one, &one, descwork,
			         &beta,
			         A, &istart, &jstart, DescA);
		}
	}

	// Flag (?)
	Work[0] = DescA[mb_];

}

/// GEneral SolVer: \f$ \{Y\} \gets [A]^{-1}\{Y\}   \f$
/** \verbatim
 * 
 * {Y} <- inv([A]) * {Y}
 *
 * OBS.: This function preserves A
 *
 * \endverbatim */
inline int Gesv(LinAlg::Matrix<double> const & A, LinAlg::Vector<double> & Y)
{

	// Information
	int size = MPI::COMM_WORLD.Get_size(); // Number of process
	if (size<2) throw new Fatal(_("ScaLAPACK::Gesv: There must be at least 2 processors for this solver."));

	// Factor the number of processors into a rectangle quadratic shape preferred
	int nprow = -1;
	int npcol = -1;
	Util::FindBestSquare(size, nprow,npcol);
	if ((nprow<1) || (npcol<1)) throw new Fatal(_("ScaLAPACK::Gesv: Could not find nprow=%d and npcol=%d"),nprow,npcol);

	// Determine a good block size
	int N    = A.Rows();
	int NRHS = 1;
	int nb   = 2;
	if ((N/nprow)<(N/npcol)) nb=(N/nprow);
	else                     nb=(N/npcol);
	if (nb<1)  nb=1;
	if (nb>64) nb=64;

	// Debug
#ifdef DO_DEBUG
	nprow=2;
	npcol=3;
	nb   =2;
	if (size!=6) throw new Fatal(_("ScaLAPACK::Gesv: The number of process for DEBUG must be equal to 6"));
#endif

	// Check
	if (nb>N) nb=N;
	if (nprow*npcol>size) throw new Fatal(_("ScaLAPACK::Gesv: There is not enough processes (%d<%d) to make a %d-by-%d process grid"),size,nprow*npcol,nprow,npcol);

	// Initialize BLACS process grid
	int iam,   nprocs, ictxt;
	int myrow, mycol;
	Cblacs_pinfo    (&iam, &nprocs);
	Cblacs_get      (-1, 0, &ictxt);
	Cblacs_gridinit (&ictxt, /*order*/"Col-major", nprow, npcol);
	Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

#ifdef DO_DEBUG
	// Debug
	String msg, buf;
	msg.Printf("Hi, I'am {%d,%d}: nprow=%d, npcol=%d, nb=%d",myrow,mycol,nprow,npcol,nb);
	std::map<int,String> letters;
	letters[19]="S";  letters[-19]="-S";
	letters[ 3]="C";  letters[ -3]="-C";
	letters[ 1]="A";  letters[ -1]="-A";
	letters[12]="L";  letters[-12]="-L";
	letters[16]="P";  letters[-16]="-P";
	letters[11]="K";  letters[-11]="-K";
	msg.append(" A =\n");
	int rank = MPI::COMM_WORLD.Get_rank();
	if (rank==0)
	{
		for (int i=0; i<A.Rows(); ++i)
		{
			for (int j=0; j<A.Rows(); ++j)
			{
				buf.Printf("   %2s",letters[static_cast<int>(A(i,j))].CStr());
				msg.append(buf);
			}
			msg.append("\n");
		}
	}
#endif
	
	// Work only the process in the process grid
	int info = 0;
	if ((myrow>-1)&&(mycol>-1)&&(myrow<nprow)&&(mycol<npcol))
	{
		// Compute the size of the local matrices
		int isrc    = 0;
		int m_loc_A = numroc_(&N   , &nb, &myrow, &isrc, &nprow); // No of rows in local data block
		int n_loc_A = numroc_(&N   , &nb, &mycol, &isrc, &npcol); // No of columns in local data block
		int n_loc_b = numroc_(&NRHS, &nb, &mycol, &isrc, &npcol); // No of columns in local data block

		// m_loc_b must be equal to m_loc_A (TODO: check this ?)
		
#ifdef DO_DEBUG
		// Debug
		buf.Printf(" m_loc_A=%d, n_loc_A=%d, n_loc_b=%d\n", m_loc_A,n_loc_A,n_loc_b);
		msg.append(buf);
#endif

		// Allocate and fill the matrices [loc_A] and [loc_b]
		double * loc_A = new double [m_loc_A*n_loc_A]; // local A (sub) array
		double * loc_b = NULL;
		if (n_loc_b>0)
			loc_b = new double [m_loc_A*n_loc_b]; // local b (sub) array

		// Initialize the array descriptor for A and b arrays
		int loc_ld_A = Util::Max(1,m_loc_A); // local leading dimension
		int loc_ld_b = loc_ld_A;
		int desc_A[9];
		int desc_b[9];
		descinit_(desc_A, &N, &N   , &nb, &nb, &isrc, &isrc, &ictxt, &loc_ld_A, &info);
		descinit_(desc_b, &N, &NRHS, &nb, &nb, &isrc, &isrc, &ictxt, &loc_ld_b, &info);

		// Fill [loc_A] and [loc_b] (sub) arrays
		for (int i_loc=0; i_loc<m_loc_A; i_loc++)
		{
			int i_glob = IndxL2G(i_loc, nb, myrow, isrc, nprow);
			for (int j_loc=0; j_loc<n_loc_A; j_loc++)
			{
				int j_glob = IndxL2G(j_loc, nb, mycol, isrc, npcol);
#ifdef DO_DEBUG
				// Debug
				buf.Printf("  [%d,%d]=A(%d,%d)=%2s",i_loc,j_loc,i_glob,j_glob,letters[static_cast<int>(A(i_glob,j_glob))].CStr());
				msg.append(buf);
#endif
				// Fill [loc_A]
				loc_A[j_loc*m_loc_A+i_loc] = A(i_glob,j_glob);
			}
#ifdef DO_DEBUG
			// Debug
			msg.append("\n");
#endif
			// Fill [loc_b]
			if (n_loc_b>0) loc_b[i_loc] = Y(i_glob);
		}

		// Allocate a 'pivot' array
		int * ippiv = new int [m_loc_A+nb];

		// Call ScaLAPACK PDGESV routine
		int ione = 1;
#ifdef DO_DEBUG
		double tini = MPI::Wtime();
#endif

		// Solve
		pdgesv_(&N, &NRHS, loc_A, &ione, &ione, desc_A, ippiv, loc_b, &ione, &ione, desc_b, &info);

#ifdef DO_DEBUG
		// Debug
		double tfin = MPI::Wtime();
		if (n_loc_b>0)
		{
			msg.append("\n loc_b = ");
			for (int i_loc=0; i_loc<m_loc_A; i_loc++)
			{
				buf.Printf("%12f ",loc_b[i_loc]);
				msg.append(buf);
			}
			msg.append("\n\n");
		}
		buf.Printf(" Elapsed time = %f\n",tfin-tini);
		msg.append(buf);
#endif

		// Fill global solution vector (all-to-gather)
		int lld_Y = N;
		int desc_Y[9]; 
		int my = N*nprow;
		int ny = NRHS*npcol;
		descinit_(desc_Y, &my, &ny, &N, &NRHS, &isrc, &isrc, &ictxt, &lld_Y, &info);
		for (int i=0; i<nprow; ++i)
		{
			int istart = i*N+1;
			for (int j=0; j<npcol; ++j)
			{
				int jstart = j*NRHS+1;
				Cpdgemr2d(N,NRHS,loc_b,1,1,desc_b, Y.GetPtr(),istart,jstart,desc_Y,ictxt);
			}
		}

		// Wait until Y is filled
		Cblacs_barrier(ictxt, "All");

#ifdef DO_DEBUG
		// Debug
		msg.append("\n\n Solution: Y = ");
		for (int i=0; i<Y.Size(); ++i)
		{
			buf.Printf("%e ",Y(i));
			msg.append(buf);
		}
		msg.append("\n\n");
		std::cout<<msg<<std::endl;
#endif

		// Clean up
		delete [] loc_A;
		if (n_loc_b>0)
			delete [] loc_b;
		delete [] ippiv;

		// Grid exit
		Cblacs_gridexit(0);
	}

	// End
	return info;

}

}; // namespace ScaLAPACK

#endif // MECHSYS_SCALAPACK_H

/* MANUAL
 

   These notes were 'borrowed' from ScaLAPACK manual:

  -- ScaLAPACK routine (version 1.7) --
     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
     and University of California, Berkeley.
     Jan 30, 2006

 
  ======================================================================

  INTEGER FUNCTION INDXL2G( INDXLOC, NB, IPROC, ISRCPROC, NPROCS )

  Purpose
  =======

  INDXL2G computes the global index of a distributed matrix entry
  pointed to by the local index INDXLOC of the process indicated by
  IPROC.

  Arguments
  =========

  INDXLOC   (global input) INTEGER
            The local index of the distributed matrix entry.

  NB        (global input) INTEGER
            Block size, size of the blocks the distributed matrix is
            split into.

  IPROC     (local input) INTEGER
            The coordinate of the process whose local array row or
            column is to be determined.

  ISRCPROC  (global input) INTEGER
            The coordinate of the process that possesses the first
            row/column of the distributed matrix.

  NPROCS    (global input) INTEGER
            The total number processes over which the distributed
            matrix is distributed.

  ======================================================================


  ======================================================================

  INTEGER FUNCTION NUMROC(N, NB, IPROC, ISRCPROC, NPROCS)

  NUMROC computes the NUMber of Rows Or Columns of a distributed
  matrix owned by the process indicated by IPROC.

  Arguments
  =========

  N         (global input) INTEGER
            The number of rows/columns in distributed matrix.

  NB        (global input) INTEGER
            Block size, size of the blocks the distributed matrix is
            split into.

  IPROC     (local input) INTEGER
            The coordinate of the process whose local array row or
            column is to be determined.

  ISRCPROC  (global input) INTEGER
            The coordinate of the process that possesses the first
            row or column of the distributed matrix.

  NPROCS    (global input) INTEGER
            The total number processes over which the matrix is
            distributed.

  ======================================================================


  ======================================================================

  SUBROUTINE
  DESCINIT(DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO)

  DESCINIT initializes the descriptor vector with the 8 input arguments
  M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD.

  Notes
  =====

  Each global data object is described by an associated description
  vector.  This vector stores the information required to establish
  the mapping between an object element and its corresponding process
  and memory location.

  Let A be a generic term for any 2D block cyclicly distributed array.
  Such a global array has an associated description vector DESCA.
  In the following comments, the character _ should be read as
  "of the global array".

  NOTATION        STORED IN      EXPLANATION
  --------------- -------------- --------------------------------------
  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
                                 DTYPE_A = 1.
  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
                                 the BLACS process grid A is distribu-
                                 ted over. The context itself is glo-
                                 bal, but the handle (the integer
                                 value) may vary.
  M_A    (global) DESCA( M_ )    The number of rows in the global
                                 array A.
  N_A    (global) DESCA( N_ )    The number of columns in the global
                                 array A.
  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
                                 the rows of the array.
  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
                                 the columns of the array.
  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
                                 row of the array A is distributed.
  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
                                 first column of the array A is
                                 distributed.
  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
                                 array.  LLD_A >= MAX(1,LOCr(M_A)).

  Let K be the number of rows or columns of a distributed matrix,
  and assume that its process grid has dimension p x q.
  LOCr( K ) denotes the number of elements of K that a process
  would receive if K were distributed over the p processes of its
  process column.
  Similarly, LOCc( K ) denotes the number of elements of K that a
  process would receive if K were distributed over the q processes of
  its process row.
  The values of LOCr() and LOCc() may be determined via a call to the
  ScaLAPACK tool function, NUMROC:
          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
  An upper bound for these quantities may be computed by:
          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A

  Arguments
  =========

  DESC    (output) INTEGER array of dimension DLEN_.
          The array descriptor of a distributed matrix to be set.

  M       (global input) INTEGER
          The number of rows in the distributed matrix. M >= 0.

  N       (global input) INTEGER
          The number of columns in the distributed matrix. N >= 0.

  MB      (global input) INTEGER
          The blocking factor used to distribute the rows of the
          matrix. MB >= 1.

  NB      (global input) INTEGER
          The blocking factor used to distribute the columns of the
          matrix. NB >= 1.

  IRSRC   (global input) INTEGER
          The process row over which the first row of the matrix is
          distributed. 0 <= IRSRC < NPROW.

  ICSRC   (global input) INTEGER
          The process column over which the first column of the
          matrix is distributed. 0 <= ICSRC < NPCOL.

  ICTXT   (global input) INTEGER
          The BLACS context handle, indicating the global context of
          the operation on the matrix. The context itself is global.

  LLD     (local input)  INTEGER
          The leading dimension of the local array storing the local
          blocks of the distributed matrix. LLD >= MAX(1,LOCr(M)).

  INFO    (output) INTEGER
          = 0: successful exit
          < 0: if INFO = -i, the i-th argument had an illegal value

  Note
  ====

  If the routine can recover from an erroneous input argument, it will
  return an acceptable descriptor vector.  For example, if LLD = 0 on
  input, DESC(LLD_) will contain the smallest leading dimension
  required to store the specified M-by-N distributed matrix, INFO
  will be set  -9 in that case.

  ======================================================================


  ======================================================================

  SUBROUTINE
  PDGEADD (TRANS, M, N, ALPHA, A, IA, JA, DESCA, BETA, C, IC, JC, DESCC)

  Purpose
  =======

  PDGEADD  adds a matrix to another

     sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )

  where

     sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),  and, op( X )  is one  of

     op( X ) = X   or   op( X ) = X'.

  Thus, op( sub( A ) ) denotes A(IA:IA+M-1,JA:JA+N-1)   if TRANS = 'N',
                               A(IA:IA+N-1,JA:JA+M-1)'  if TRANS = 'T',
                               A(IA:IA+N-1,JA:JA+M-1)'  if TRANS = 'C'.

  Alpha  and  beta  are scalars, sub( C ) and op( sub( A ) ) are m by n
  submatrices.

  Notes
  =====

  A description  vector  is associated with each 2D block-cyclicly dis-
  tributed matrix.  This  vector  stores  the  information  required to
  establish the  mapping  between a  matrix entry and its corresponding
  process and memory location.

  In  the  following  comments,   the character _  should  be  read  as
  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
  block cyclicly distributed matrix.  Its description vector is DESC_A:

  NOTATION         STORED IN       EXPLANATION
  ---------------- --------------- ------------------------------------
  DTYPE_A (global) DESCA[ DTYPE_ ] The descriptor type.
  CTXT_A  (global) DESCA[ CTXT_  ] The BLACS context handle, indicating
                                   the NPROW x NPCOL BLACS process grid
                                   A  is  distributed over. The context
                                   itself  is  global,  but  the handle
                                   (the integer value) may vary.
  M_A     (global) DESCA[ M_     ] The  number of rows in the distribu-
                                   ted matrix A, M_A >= 0.
  N_A     (global) DESCA[ N_     ] The number of columns in the distri-
                                   buted matrix A, N_A >= 0.
  IMB_A   (global) DESCA[ IMB_   ] The number of rows of the upper left
                                   block of the matrix A, IMB_A > 0.
  INB_A   (global) DESCA[ INB_   ] The  number  of columns of the upper
                                   left   block   of   the  matrix   A,
                                   INB_A > 0.
  MB_A    (global) DESCA[ MB_    ] The blocking factor used to  distri-
                                   bute the last  M_A-IMB_A  rows of A,
                                   MB_A > 0.
  NB_A    (global) DESCA[ NB_    ] The blocking factor used to  distri-
                                   bute the last  N_A-INB_A  columns of
                                   A, NB_A > 0.
  RSRC_A  (global) DESCA[ RSRC_  ] The process row over which the first
                                   row of the matrix  A is distributed,
                                   NPROW > RSRC_A >= 0.
  CSRC_A  (global) DESCA[ CSRC_  ] The  process column  over  which the
                                   first column of  A  is  distributed.
                                   NPCOL > CSRC_A >= 0.
  LLD_A   (local)  DESCA[ LLD_   ] The  leading dimension  of the local
                                   array  storing  the  local blocks of
                                   the distributed matrix A,
                                   IF( Lc( 1, N_A ) > 0 )
                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
                                   ELSE
                                      LLD_A >= 1.

  Let K be the number of  rows of a matrix A starting at the global in-
  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
  receive if these K rows were distributed over NPROW processes.  If  K
  is the number of columns of a matrix  A  starting at the global index
  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
  these K columns were distributed over NPCOL processes.

  The values of Lr() and Lc() may be determined via a call to the func-
  tion PB_Cnumroc:
  Lr( IA, K ) = PB_Cnumroc( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
  Lc( JA, K ) = PB_Cnumroc( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )

  Arguments
  =========

  TRANS   (global input) CHARACTER*1
          On entry,  TRANS   specifies the form of op( sub( A ) ) to be
          used in the matrix addition as follows:

             TRANS = 'N' or 'n'   op( sub( A ) ) = sub( A ),

             TRANS = 'T' or 't'   op( sub( A ) ) = sub( A )',

             TRANS = 'C' or 'c'   op( sub( A ) ) = sub( A )'.

  M       (global input) INTEGER
          On entry,  M  specifies the number of rows of  the  submatrix
          sub( C ) and the number of columns of the submatrix sub( A ).
          M  must be at least zero.

  N       (global input) INTEGER
          On entry, N  specifies the number of columns of the submatrix
          sub( C ) and the number of rows of the submatrix sub( A ).  N
          must be at least zero.

  ALPHA   (global input) DOUBLE PRECISION
          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
          supplied  as  zero  then  the  local entries of  the array  A
          corresponding to the entries of the submatrix  sub( A )  need
          not be set on input.

  A       (local input) DOUBLE PRECISION array
          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
          at least Lc( 1, JA+M-1 ).  Before  entry, this array contains
          the local entries of the matrix A.

  IA      (global input) INTEGER
          On entry, IA  specifies A's global row index, which points to
          the beginning of the submatrix sub( A ).

  JA      (global input) INTEGER
          On entry, JA  specifies A's global column index, which points
          to the beginning of the submatrix sub( A ).

  DESCA   (global and local input) INTEGER array
          On entry, DESCA  is an integer array of dimension DLEN_. This
          is the array descriptor for the matrix A.

  BETA    (global input) DOUBLE PRECISION
          On entry,  BETA  specifies the scalar  beta.   When  BETA  is
          supplied  as  zero  then  the  local entries of  the array  C
          corresponding to the entries of the submatrix  sub( C )  need
          not be set on input.

  C       (local input/local output) DOUBLE PRECISION array
          On entry, C is an array of dimension (LLD_C, Kc), where Kc is
          at least Lc( 1, JC+N-1 ).  Before  entry, this array contains
          the local entries of the matrix C.
          On exit, the entries of this array corresponding to the local
          entries of the submatrix  sub( C )  are  overwritten  by  the
          local entries of the m by n updated submatrix.

  IC      (global input) INTEGER
          On entry, IC  specifies C's global row index, which points to
          the beginning of the submatrix sub( C ).

  JC      (global input) INTEGER
          On entry, JC  specifies C's global column index, which points
          to the beginning of the submatrix sub( C ).

  DESCC   (global and local input) INTEGER array
          On entry, DESCC  is an integer array of dimension DLEN_. This
          is the array descriptor for the matrix C.

  -- Written on April 1, 1998 by
     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.

  ======================================================================


  ======================================================================

  SUBROUTINE PDGESV( N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB, INFO )

  Purpose
  =======

  PDGESV computes the solution to a real system of linear equations

                        sub( A ) * X = sub( B ),

  where sub( A ) = A(IA:IA+N-1,JA:JA+N-1) is an N-by-N distributed
  matrix and X and sub( B ) = B(IB:IB+N-1,JB:JB+NRHS-1) are N-by-NRHS
  distributed matrices.

  The LU decomposition with partial pivoting and row interchanges is
  used to factor sub( A ) as sub( A ) = P * L * U, where P is a permu-
  tation matrix, L is unit lower triangular, and U is upper triangular.
  L and U are stored in sub( A ). The factored form of sub( A ) is then
  used to solve the system of equations sub( A ) * X = sub( B ).

  Notes
  =====

  Each global data object is described by an associated description
  vector.  This vector stores the information required to establish
  the mapping between an object element and its corresponding process
  and memory location.

  Let A be a generic term for any 2D block cyclicly distributed array.
  Such a global array has an associated description vector DESCA.
  In the following comments, the character _ should be read as
  "of the global array".

  NOTATION        STORED IN      EXPLANATION
  --------------- -------------- --------------------------------------
  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
                                 DTYPE_A = 1.
  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
                                 the BLACS process grid A is distribu-
                                 ted over. The context itself is glo-
                                 bal, but the handle (the integer
                                 value) may vary.
  M_A    (global) DESCA( M_ )    The number of rows in the global
                                 array A.
  N_A    (global) DESCA( N_ )    The number of columns in the global
                                 array A.
  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
                                 the rows of the array.
  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
                                 the columns of the array.
  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
                                 row of the array A is distributed.
  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
                                 first column of the array A is
                                 distributed.
  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
                                 array.  LLD_A >= MAX(1,LOCr(M_A)).

  Let K be the number of rows or columns of a distributed matrix,
  and assume that its process grid has dimension p x q.
  LOCr( K ) denotes the number of elements of K that a process
  would receive if K were distributed over the p processes of its
  process column.
  Similarly, LOCc( K ) denotes the number of elements of K that a
  process would receive if K were distributed over the q processes of
  its process row.
  The values of LOCr() and LOCc() may be determined via a call to the
  ScaLAPACK tool function, NUMROC:
          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
  An upper bound for these quantities may be computed by:
          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A

  This routine requires square block decomposition ( MB_A = NB_A ).

  Arguments
  =========

  N       (global input) INTEGER
          The number of rows and columns to be operated on, i.e. the
          order of the distributed submatrix sub( A ). N >= 0.

  NRHS    (global input) INTEGER
          The number of right hand sides, i.e., the number of columns
          of the distributed submatrix sub( B ). NRHS >= 0.

  A       (local input/local output) DOUBLE PRECISION pointer into the
          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
          On entry, the local pieces of the N-by-N distributed matrix
          sub( A ) to be factored. On exit, this array contains the
          local pieces of the factors L and U from the factorization
          sub( A ) = P*L*U; the unit diagonal elements of L are not
          stored.

  IA      (global input) INTEGER
          The row index in the global array A indicating the first
          row of sub( A ).

  JA      (global input) INTEGER
          The column index in the global array A indicating the
          first column of sub( A ).

  DESCA   (global and local input) INTEGER array of dimension DLEN_.
          The array descriptor for the distributed matrix A.

  IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
          This array contains the pivoting information.
          IPIV(i) -> The global row local row i was swapped with.
          This array is tied to the distributed matrix A.

  B       (local input/local output) DOUBLE PRECISION pointer into the
          local memory to an array of dimension
          (LLD_B,LOCc(JB+NRHS-1)).  On entry, the right hand side
          distributed matrix sub( B ). On exit, if INFO = 0, sub( B )
          is overwritten by the solution distributed matrix X.

  IB      (global input) INTEGER
          The row index in the global array B indicating the first
          row of sub( B ).

  JB      (global input) INTEGER
          The column index in the global array B indicating the
          first column of sub( B ).

  DESCB   (global and local input) INTEGER array of dimension DLEN_.
          The array descriptor for the distributed matrix B.

  INFO    (global output) INTEGER
          = 0:  successful exit
          < 0:  If the i-th argument is an array and the j-entry had
                an illegal value, then INFO = -(i*100+j), if the i-th
                argument is a scalar and had an illegal value, then
                INFO = -i.
          > 0:  If INFO = K, U(IA+K-1,JA+K-1) is exactly zero.
                The factorization has been completed, but the factor U
                is exactly singular, so the solution could not be
                computed.

  =====================================================================


  =====================================================================

     SUBROUTINE PDGEMR2D( M, N,
                         A, IA, JA, ADESC,
                         B, IB, JB, BDESC,
                         CTXT)
   Purpose
   =======

   PDGEMR2D copies a submatrix of A on a submatrix of B.
   A and B can have different distributions: they can be on different
   processor grids, they can have different blocksizes, the beginning
   of the area to be copied can be at a different places on A and B.

   The parameters can be confusing when the grids of A and B are
   partially or completly disjoint, in the case a processor calls
   this routines but is either not in the A context or B context, the
   ADESC[CTXT] or BDESC[CTXT] must be equal to -1, to ensure the
   routine recognise this situation.
   To summarize the rule:
   - If a processor is in A context, all parameters related to A must be valid.
   - If a processor is in B context, all parameters related to B must be valid.
   -  ADESC[CTXT] and BDESC[CTXT] must be either valid contexts or equal to -1.
   - M and N must be valid for everyone.
   - other parameters are not examined.

   Notes
   =====

   A description vector is associated with each 2D block-cyclicly dis-
   tributed matrix.  This vector stores the information required to
   establish the mapping between a matrix entry and its corresponding
   process and memory location.

   In the following comments, the character _ should be read as
   "of the distributed matrix".  Let A be a generic term for any 2D
   block cyclicly distributed matrix.  Its description vector is DESC_A:

  NOTATION        STORED IN      EXPLANATION
  --------------- -------------- --------------------------------------
  DT_A   (global) DESCA( DT_ )   The descriptor type.
  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
                                 the BLACS process grid A is distribu-
                                 ted over. The context itself is glo-
                                 bal, but the handle (the integer
                                 value) may vary.
  M_A    (global) DESCA( M_ )    The number of rows in the distributed
                                 matrix A.
  N_A    (global) DESCA( N_ )    The number of columns in the distri-
                                 buted matrix A.
  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
                                 the rows of A.
  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
                                 the columns of A.
  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
                                 row of the matrix A is distributed.
  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
                                 first column of A is distributed.
  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
                                 array storing the local blocks of the
                                 distributed matrix A.
                                 LLD_A >= MAX(1,LOCp(M_A)).
   Important notice
   ================
    The parameters of the routine have changed in April 1996
    There is a new last argument. It must be a context englobing
    all processors involved in the initial and final distribution.

    Be aware that all processors  included in this
     context must call the redistribution routine.

   Parameters
   ==========


   M        (input) INTEGER.
            On entry, M specifies the number of rows of the
            submatrix to be copied.  M must be at least zero.
            Unchanged on exit.

   N        (input) INTEGER.
            On entry, N specifies the number of cols of the submatrix
            to be redistributed.rows of B.  M must be at least zero.
            Unchanged on exit.

   A        (input) DOUBLE PRECISION
            On entry, the source matrix.
            Unchanged on exit.

   IA,JA    (input) INTEGER
            On entry,the coordinates of the beginning of the submatrix
            of A to copy.
            1 <= IA <= M_A - M + 1,1 <= JA <= N_A - N + 1,
            Unchanged on exit.

   ADESC    (input) A description vector (see Notes above)
            If the current processor is not part of the context of A
            the ADESC[CTXT] must be equal to -1.


   B        (output) DOUBLE PRECISION
            On entry, the destination matrix.
            The portion corresponding to the defined submatrix are updated.

   IB,JB    (input) INTEGER
            On entry,the coordinates of the beginning of the submatrix
            of B that will be updated.
            1 <= IB <= M_B - M + 1,1 <= JB <= N_B - N + 1,
            Unchanged on exit.

   BDESC    (input) B description vector (see Notes above)
            For processors not part of the context of B
            BDESC[CTXT] must be equal to -1.

   CTXT     (input) a context englobing at least all processors included
               in either A context or B context

  Memory requirement :
  ====================

  for the processors belonging to grid 0, one buffer of size block 0
  and for the processors belonging to grid 1, also one buffer of size
  block 1.

  =====================================================================

*/
