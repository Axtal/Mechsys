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


/*  __ Solver for (un)coupled problems. Serial/Parallel with Dense/Sparse matrices __
 *
 *  This solver can solve the following kind of problems:
 *    _           _                 _          _
 *   |   K    L1   | / dU/dt \     |   0    0   | / U \     / dF/dt \
 *   |             | |       |  +  |            | |   |  =  |       |
 *   |_  L2   M   _| \ dP/dt /     |_  0    H  _| \ P /     \ dQ/dt /
 *
 *               [C] { dU/dt }  +             [K] { U }  =  { dF/dt }
 *
 *  For example:
 *
 *    U = Displacement
 *    P = Porepressure
 *    F = Force
 *    Q = Water discharge
 *
 *    K = Stiffness
 *   L1 = CouplingMatrix1
 *   L2 = CouplingMatrix2
 *    M = MassMatrix
 *    H = Permeability
 *
 *    C = Order1Matrix (Applied to the first order derivative dU)
 *    K = Order0Matrix (Applied to the zeroth order derivative U)
 *
 *    num(Order1Matrix) = 4
 *    num(Order0Matrix) = 1
 *
 *  Solve:     ([C] + alpha*h*[K]) * {dU} = {dF} - h*[K]*{U}
 */


#ifndef MECHSYS_FEM_SOLVER_H
#define MECHSYS_FEM_SOLVER_H

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "fem/geometry.h"
#include "fem/parallel.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "linalg/sparse_triplet.h"
#include "linalg/sparse_matrix.h"
#include "linalg/sparse_crmatrix.h"
#ifdef HAVE_UMFPACK
  #include "linalg/umfpack.h"
#endif
#ifdef HAVE_SUPERLU
  #include "linalg/superlu.h"
#endif
#ifdef HAVE_SUPERLUD
  #include "linalg/superlud.h"
#endif

namespace FEM
{

/// Class which solves the global DAS (Differential Algebraic System)
class Solver
{
public:
	// Constructor
	Solver() : _num_div(1), _delta_time(0.0), _do_output(false), _g(NULL), _pd(NULL), _ndofs(0) { SetLinSol("UM"); }

	// Destructor
	virtual ~Solver() {}

	// Methods
	Solver * SetGeom      (FEM::Geom * G)    { _g  = G;  return this; }              ///< Set the geometry (Nodes/Elements) to be used during solution
	Solver * SetPData     (FEM::PData * PD)  { _pd = PD; return this; }              ///< Set structure with data for parallel computation
	Solver * SetLinSol    (char const * Key);                                        ///< LinSolName: LA=LAPACK, UM=UMFPACK, SLU=SuperLU, SLUd=SuperLUd
	Solver * SetNumDiv    (int    NumDiv)    { _num_div   =NumDiv;    return this; } ///< Set the number of divisions to be used during each call of Solve
	Solver * SetDeltaTime (double DeltaTime) { _delta_time=DeltaTime; return this; } ///< Set the time stepsize to be used during each call of Solve
	void     Solve        ();                                                        ///< Solve ([C]+alpha*h*[K])*{dU} = {dF}-h*[K]*{U} for the boundary conditions defined inside the nodes array

	// Access methods
	LinAlg::Vector<double> const & Resid () const { return _resid;  } ///< Resid = dFext - dFint
	size_t                         nDOF  () const { return _ndofs;  } ///< Number of DOFs

	// Methods to be overloaded by derived classes
	virtual Solver * SetCte (char const * Key, double Value) =0; ///< Set solver constant such as number of subincrements, DTOL, etc.
	virtual double   GetVar (char const * Key) const         =0; ///< Get solver variable such as Residuals or Relative error

protected:
	// Data
	int                      _num_div;    ///< The number of divisions to be used during each call of Solve
	double                   _delta_time; ///< The time stepsize to be used during each call of Solve
	bool                     _do_output;  ///< Do output?
	LinAlg::Vector<double>   _dF_ext;     ///< Increment of natural values, divided by the number of increments, which update the current state towards the condition at the end the stage being solved.
	LinAlg::Vector<double>   _dU_ext;     ///< Increment of essential values, divided by the number of increments, which update the current state towards the condition at the end the stage being solved.
	LinAlg::Vector<double>   _dF_int;     ///< Increment of internal natural (forces) values, divided by the number of increments, correspondent to the increment of external forces.
	LinAlg::Vector<double>   _resid;      ///< Residual: resid = dFext - dFint
	LinAlg::Vector<double>   _hKU;        ///< Linearized independent term of the differential equation.
	bool                     _has_hKU;    ///< Flag which says if any element has to contribute to the hKU vector. If _has_hKU==false, there is no need for the hKU vector, because there are no Order0Matrices in this stage of the simulation.

	// Methods
	void   _inv_G_times_dF_minus_hKU   (double h, LinAlg::Vector<double> & dF, LinAlg::Vector<double> & dU);            ///< Compute (linear solver) inv(G)*(dF-hKU), where G may be assembled by Order1 and Order0 matrices
	void   _update_nodes_and_elements  (double h, LinAlg::Vector<double> const & dF,LinAlg::Vector<double> const & dU); ///< Update nodes essential/natural values and elements internal values (such as stresses/strains)
	void   _backup_nodes_and_elements  ();                                                                              ///< Backup nodes essential/natural values and elements internal values (such as stresses/strains)
	void   _restore_nodes_and_elements ();                                                                              ///< Restore nodes essential/natural values and elements internal values (such as stresses/strains)
	double _norm_essential_vector      ();                                                                              ///< Compute the Euclidian norm of the essential (displacements) vector
	double _norm_natural_vector        ();                                                                              ///< Compute the Euclidian norm of the natural (forces) vector

	// To be overriden by derived classes
	virtual void _do_solve_for_an_increment(double dTime) =0;

private:
	// Data
	FEM::Geom            * _g;      ///< Geometry: Nodes and Elements to be used during solution
	FEM::PData           * _pd;     ///< Structure with data for parallel computation
	LinAlg::LinSol_T       _linsol; ///< The type of the linear solver to be used in the inversion of G11. Solves: X<-inv(A)*B
	int                    _nudofs; ///< Total number of unknown DOFs of a current stage
	int                    _npdofs; ///< Total number of prescribed DOFs of a current stage
	int                    _ndofs;  ///< Total number of DOFs
	Array<Node::DOF*>      _udofs;  ///< Unknown DOFs
	Array<Node::DOF*>      _pdofs;  ///< Prescribed DOFs
	LinAlg::Vector<double> _U_bkp;  ///< Place to store a temporary state of essential (displacements) values (backup).
	LinAlg::Vector<double> _F_bkp;  ///< Place to store a temporary state of natural (forces) values (backup).

	// Used only for "DENSE" solvers
	LinAlg::Matrix<double> _G;      ///< "Global-stiffness" dense matrix

	// Used only for "SPARSE" solvers
	Sparse::Triplet<double,int> _T11;      ///< Sparse G11 piece of the global matrix _G in the Triplet format
	Sparse::Triplet<double,int> _T12;      ///< Sparse G12 piece of the global matrix _G in the Triplet format
	Sparse::Triplet<double,int> _T21;      ///< Sparse G21 piece of the global matrix _G in the Triplet format
	Sparse::Triplet<double,int> _T22;      ///< Sparse G22 piece of the global matrix _G in the Triplet format
	int                         _T11_size; ///< Size of T11
	int                         _T12_size; ///< Size of T12
	int                         _T21_size; ///< Size of T21
	int                         _T22_size; ///< Size of T22

	// Private methods
	void _assemb_G_and_hKU     (double h);                                                                                                                                                   ///< Assemble matrix G for the actual state of elements/nodes. Note h == TimeInc
	void _mount_into_hKU       (LinAlg::Vector<double> const & V, Array<size_t> const & RowsMap);                                                                                            ///< Add a local vector, such as [U P]^T into the hKU vector
	void _mount_into_global    (LinAlg::Matrix<double> const & M, Array<size_t> const & RowsMap, Array<size_t> const & ColsMap, Array<bool> const & RowsPre, Array<bool> const & ColsPre);   ///< Add a local matrix, such as K,L1,L2 or M into the global DENSE matrix G or into the global SPARSE matrix G (T11,T12,T21,T22)
	void _compute_global_size  ();                                                                                                                                                           ///< Compute the sizes of T11,T12,T21 and T22
	void _increase_global_size (Array<size_t> const & RowsMap, Array<size_t> const & ColsMap, Array<bool> const & RowsPre, Array<bool> const & ColsPre);                                     ///< Auxiliar method to compute the sizes of T11,T12,T21 and T22
	void _copy_partial_matrix  (LinAlg::Matrix<double> const & Source, LinAlg::Vector<int> const & ValidCols, LinAlg::Vector<int> const & ValidRows, LinAlg::Matrix<double> & Target) const; ///< Auxiliar method to copy a piece of a matrix inside other matrix
	void _do_scatter           (LinAlg::Matrix<double> const & K, LinAlg::Vector<double> const & X, LinAlg::Vector<double> const & Y, LinAlg::Matrix<double> & K11, LinAlg::Matrix<double> & K12, LinAlg::Matrix<double> & K21, LinAlg::Matrix<double> & K22, LinAlg::Vector<double> & Y2, LinAlg::Vector<double> & X1); ///< Scatter "global" matrices to small pieces
	void _do_gather            (LinAlg::Vector<double> const & Y1, LinAlg::Vector<double> const & Y2, LinAlg::Vector<double> const & X1, LinAlg::Vector<double> const & X2, LinAlg::Vector<double> & Y, LinAlg::Vector<double> & X);                                                                                     ///< Gather small pieces into "global" matrices

}; // class Solver


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline Solver * Solver::SetLinSol(char const * Key)
{
	     if (strcmp(Key,"LA"  )==0) _linsol=LinAlg::LAPACK_T;
	else if (strcmp(Key,"UM"  )==0) _linsol=LinAlg::UMFPACK_T;
	else if (strcmp(Key,"SLU" )==0) _linsol=LinAlg::SuperLU_T;
	else if (strcmp(Key,"SLUd")==0) _linsol=LinAlg::SuperLU_T;
	else throw new Fatal("Solver::SetLinSol: Solver key==%s is invalid",Key);
	return this;
}

inline void Solver::Solve()
{
	// Solve:    ([C] + alpha*h*[K]) * {dU} = {dF} - h*[K]*{U}
	
	// Check geometry
	if (_g==NULL) throw new Fatal(_("Solver::Solve: Solver::SetGeom(FEM::Geom * G) must be called before calling this method (Solver::Solve())"));
	if (_g->Check()==false) throw new Fatal("Solver::Solve: FEM Geometry was not properly set (pointers may be NULL)");

	// Loop over all active elements
	double has_fvol = false;
	for (size_t i=0; i<_g->NElems(); ++i)
	{
		if (_g->Ele(i)->IsActive())
		{
			// Check if active elements are properly initialized
			String msg;
			if (_g->Ele(i)->Check(msg)==false)
				throw new Fatal("Solver::Solve Element # %d was not properly initialized. Error: %s",i,msg.CStr());

			// Check if any element has volumetric forces
			has_fvol = _g->Ele(i)->HasVolForces();
		}
	}

	// Number of divisions for each increment
	double dTime = _delta_time / _num_div; // Time increment

	// Number of DOFs
	_ndofs  = 0; // total
	_nudofs = 0; // unknown
	_npdofs = 0; // prescribed

	// Set the ID (equation number) of all DOFs and setup _dF_ext and _dU_ext vectors
	_udofs.Resize (0);
	_pdofs.Resize (0);
	for (size_t i=0; i<_g->NNodes(); ++i)
	{
		if (_g->Nod(i)->nSharedBy()>0)
		{
			// Loop along DOFs inside Node [i]
			for (size_t j=0; j<_g->Nod(i)->nDOF(); ++j)
			{
				// Save pointers to the DOFs
				if (_g->Nod(i)->DOFVar(j).IsEssenPresc) _pdofs.Push(&_g->Nod(i)->DOFVar(j));
				else                                    _udofs.Push(&_g->Nod(i)->DOFVar(j));
				
				// Assing an equation ID to DOF [j]
				_g->Nod(i)->DOFVar(j).EqID = _ndofs;

				// Next EqID
				_ndofs++;
			}
		}
	}
	_nudofs = static_cast<int>(_udofs.Size());
	_npdofs = static_cast<int>(_pdofs.Size());
	
	// Check
	if (_nudofs==0) throw new Fatal("Solver::Solve: The number of Unknowns DOFs must be greater than zero.\n  There might be an error with the boundary conditions' setting up.");
	if (_npdofs==0) throw new Fatal("Solver::Solve: The number of Prescribed DOFs must be greater than zero.\n  There might be an error with the boundary conditions' setting up.");

#ifdef HAVE_SUPERLUD
	if (_pd==NULL)     throw new Fatal(_("Solve::Solve: For parallel computation, PData (ParallelData) must be set before calling this method"));
	if (_ndofs!=nDOFs) throw new Fatal(_("Solver::Solve: Just computed _ndofs (%d) must be equal to nDOFs (%d)"),_ndofs,DOFs);
#endif

	// Allocate vectors
	_dF_ext.Resize (_ndofs);
	_dU_ext.Resize (_ndofs);
	_U_bkp .Resize (_ndofs);
	_F_bkp .Resize (_ndofs);
	_hKU   .Resize (_ndofs);
	_dF_int.Resize (_ndofs);
	_resid .Resize (_ndofs);

	// Set the vectors with increments of boundary conditions
	_dF_ext.SetValues (0.0);
	_dU_ext.SetValues (0.0);
	for (int i=0; i<_nudofs; ++i)
	{
		_dF_ext(_udofs[i]->EqID) = _udofs[i]->NaturalBry   / _num_div;
		_dU_ext(_udofs[i]->EqID) = _udofs[i]->EssentialBry / _num_div;
	}
	for (int i=0; i<_npdofs; ++i)
	{
		_dF_ext(_pdofs[i]->EqID) = _pdofs[i]->NaturalBry   / _num_div;
		_dU_ext(_pdofs[i]->EqID) = _pdofs[i]->EssentialBry / _num_div;
	}

	// Set component of external forces due to volume forces such as body forces, heat source, etc.
	if (has_fvol)
	{
		// Calculate fvol
		LinAlg::Vector<double> fvol;
		fvol.Resize    (_ndofs);
		fvol.SetValues (0.0);
		for (size_t i=0; i<_g->NElems(); ++i)
			if (_g->Ele(i)->IsActive()) _g->Ele(i)->AddVolForces(fvol); // add results to fvol

		// Add to dFext vector
		LinAlg::Axpy(1.0/_num_div,fvol, _dF_ext); // dFext <- dFext + fvol/numdiv
	}

	// Allocate stifness DENSE matrix G or Allocate stifness SPARSE matrix G
	if (_linsol==LinAlg::LAPACK_T) _G.Resize(_ndofs,_ndofs);
	else
	{
		// Compute the sizes of the triplets
		_compute_global_size();
		// Augment G11 according to the prescribed DOFs
		_T11_size += _npdofs;
		// Allocate sparse triplets of the G matrix
		_T11.AllocSpace(_ndofs, _ndofs, _T11_size);
		_T12.AllocSpace(_ndofs, _ndofs, _T12_size);
		_T21.AllocSpace(_ndofs, _ndofs, _T21_size);
		_T22.AllocSpace(_ndofs, _ndofs, _T22_size);
		// Put ones in the position of prescribed DOFs
		for (int i=0; i<_npdofs; ++i) _T11.PushEntry(_pdofs[i]->EqID, _pdofs[i]->EqID, 1.0);
	}

	// Output initial state
	//if (iStage==0 && _do_output) OutputIncrement(iStage, -1);

	// Loop over increments (0<increment<_num_div)
	for (int increment=0; increment<_num_div; ++increment)
	{
		// Solve for increment
		_do_solve_for_an_increment(dTime);

		// Output results for increment
		//if (_do_output) OutputIncrement(iStage, increment);
	}

	// Output end of Stage
	//if (_do_output) OutputStage(iStage);

}


/* protected */

inline void Solver::_inv_G_times_dF_minus_hKU(double h, LinAlg::Vector<double> & dF, LinAlg::Vector<double> & dU)
{
	/*   _               _
	    |  [G11]   [G12]  | / {dU1}=? \   / {dF1}   \  <== unknowns DOFs    size = _udofs
	    |                 | |         | = |         |
	    |_ [G21]   [G22] _| \ {dU2}   /   \ {dF2}=? /  <== prescribed DOFs  size = _pdofs

	 G may be:            and dF may be:

	  [G] = [C]+h[K]{U}    {dF} = {dF} - h[K]{U}

	 where:
	    _           _                 _          _
	   |   K    L1   | / dU/dt \     |   0    0   | / U \     / dF/dt \
	   |             | |       |  +  |            | |   |  =  |       |
	   |_  L2   M   _| \ dP/dt /     |_  0    H  _| \ P /     \ dQ/dt /

	               [C] { dU/dt }  +             [K] { U }  =  { dF/dt }
	 */

	// 0) Assemble the global stiffness matrix G and the hKU vector
	_assemb_G_and_hKU(h);                      // G  <- C + h*alpha*K
	if (_has_hKU) LinAlg::Axpy(-1.0,_hKU, dF); // dF <- dF - hKU

	if (_linsol==LinAlg::LAPACK_T)
	{
		// 1) Modify G for essential boundary conditions
		//        G11*dU1 + G12*dU2 = dF1
		//        G21*dU1 + G22*dU2 = dF2
		//        dU1 = inv(G11)*(dF1 - G12*dU2)
		//        dF2 = G21*dU1 + G22*dU2

		// 2) Allocate temporary matrices and vectors
		LinAlg::Matrix<double> G11,G12;
		LinAlg::Matrix<double> G21,G22;
		LinAlg::Vector<double> dU2;
		LinAlg::Vector<double> dF1;
		
		if (_nudofs > 0)
		{

			// 3) Scatter global G, dF and dU values into smaller pieces
			_do_scatter(_G,dF,dU, G11,G12,G21,G22, dU2, dF1);

			// 4) Solve for {dU1} and {dF2}
			LinAlg::Vector<double> W(dF1);                    // W <- dF1 (workspace)
			LinAlg::Gemv(-1.0,G12,dU2,1.0,W);               // W <- -1.0*G12*dU2 + 1.0*dF1
			LinAlg::Gesv(G11, W);                           // W <- inv(G11)*W  (G11 is lost)
			LinAlg::Vector<double> & dU1 = W;                 // dU1 <- W (reference)
			LinAlg::Vector<double> dF2(dU2.Size());           // allocate dF2
			LinAlg::Gemvpmv(1.0,G21,dU1, 1.0,G22,dU2, dF2); // dF2 <- 1*G21*dU1 + 1*G22*dU2

			// 5) Gather dU1, dU2, dF1 and dF2 into dU and dF
			_do_gather(dU1,dU2, dF1,dF2, dU,dF);
		}
		else
		{
			LinAlg::Gemv(1.0, _G, dU, 0.0, dF);
		}

	}
	else
	{
		// 1) Allocate workspace
		LinAlg::Vector<double> W(_ndofs); // W <- 0

		// 2) Copy prescribed dF1 values into the workspace
		//                {W} = ( {dF1} )
		//                      (   0   )
		for (int i=0; i<_nudofs; ++i) W(_udofs[i]->EqID) = dF(_udofs[i]->EqID);

		// 3) Copy prescribed dU2 values into the workspace
		//                       _            _
		//    [G11_augmented] = |  [G11]   0   |     {W} = ( {dF1} ) 
		//                      |_   0    [1] _|           ( {dU2} )
		for (int i=0; i<_npdofs; ++i) W(_pdofs[i]->EqID) = dU(_pdofs[i]->EqID);

		// 4) Modify W so inv(G11_augmented) can be multiplied to get dU2: {W} <- {W} - [G12]*{W}
		//                       _          _    
		//    {W} = ( {dF1} ) - |  0  [G12]  | * (   0   ) = ( {dF1} - [G12]*{dU2} )
		//          ( {dU2} )   |_ 0    0   _|   ( {dU2} )   (        {dU2}        )
		for (int k=0; k<_T12_size; ++k) W(_T12.Ai(k)) -= _T12.Ax(k)*dU(_T12.Aj(k));

		// Clear dF at positions of prescribed displacements
		for (int k=0; k<_T21_size; ++k) dF(_T21.Ai(k)) = 0.0;
		for (int k=0; k<_T22_size; ++k) dF(_T22.Ai(k)) = 0.0;

		if (_linsol==LinAlg::UMFPACK_T)
		{
			// 5) Allocate G11_augmented sparse matrix in compressed-COLUMN format (convert from triplet and remove duplicates)
			Sparse::Matrix<double,int> G11(_T11);

			// 6) Solve for essential (displacements) values: dU <- inv(G11)*W
			#ifdef HAVE_UMFPACK
			UMFPACK::Solve(G11,W, dU); // dU <- inv(G11)*W
			#else
			throw new Fatal(_("Solve::_inv_G_times_dF_minus_hKU: UMFPACK is not available"));
			#endif

			// 7) Solve for natural (forces) values dF2
			for (int k=0; k<_T21_size; ++k) dF(_T21.Ai(k)) += _T21.Ax(k)*dU(_T21.Aj(k)); // dF <- dF + G21*dU
			for (int k=0; k<_T22_size; ++k) dF(_T22.Ai(k)) += _T22.Ax(k)*dU(_T22.Aj(k)); // dF <- dF + G22*dU
		}
		else if (_linsol==LinAlg::SuperLU_T)
		{
			// 5) Allocate G11_augmented sparse matrix in compressed-COLUMN format (convert from triplet and remove duplicates)
			Sparse::Matrix<double,int> G11(_T11);

			// 6) Solve for essential (displacements) values: dU <- inv(G11)*W
			#ifdef HAVE_SUPERLU
			LinAlg::Copy   (W,   dU); // dU <- W
			SuperLU::Solve (G11, dU); // dU <- inv(G11)*dU
			#else
			throw new Fatal(_("Solve::_inv_G_times_dF_minus_hKU: SuperLU is not available"));
			#endif

			// 7) Solve for natural (forces) values dF2
			for (int k=0; k<_T21_size; ++k) dF(_T21.Ai(k)) += _T21.Ax(k)*dU(_T21.Aj(k)); // dF <- dF + G21*dU
			for (int k=0; k<_T22_size; ++k) dF(_T22.Ai(k)) += _T22.Ax(k)*dU(_T22.Aj(k)); // dF <- dF + G22*dU
		}
		else if (_linsol==LinAlg::SuperLUd_T)
		{
			// 5) Allocate G11_augmented sparse matrix in compressed-ROW format (convert from triplet and remove duplicates)
			Sparse::CRMatrix<double,int> G11(_T11);

			#ifdef HAVE_SUPERLUD
			// 6) Solve for essential (displacements) values: dU <- inv(G11)*W
			LinAlg::Copy    (W, dU);                                                                       // dU <- W
			SuperLUd::Solve (_pd->MyMinEq, _pd->MyMaxEq, _pd->MyNumEqs, _pd->nProcs, _pd->nDOFs, G11, dU); // dU <- inv(G11)*dU

			// 7) Distribute all pieces of dU to all processors
			MPI::COMM_WORLD.Barrier();
			MPI::COMM_WORLD.Allgatherv(&dU.GetPtr()[_pd->MyMinEq], _pd->MyNumEqs, MPI::DOUBLE, dU.GetPtr(), &_pd->AllNumEqs[0], &_pd->AllMinEq[0], MPI::DOUBLE);

			// 8) Solve for natural (forces) values dF2
			for (int k=0; k<_T21_size; ++k) { int I=_T21.Ai(k); if (I>=_pd->MyMinEq && I<=MyMaxEq) dF(I) += _T21.Ax(k)*dU(_T21.Aj(k)); } // dF <- dF + G21*dU
			for (int k=0; k<_T22_size; ++k) { int I=_T22.Ai(k); if (I>=_pd->MyMinEq && I<=MyMaxEq) dF(I) += _T22.Ax(k)*dU(_T22.Aj(k)); } // dF <- dF + G22*dU

			// 9) Distribute all pieces of dF to all processors
			MPI::COMM_WORLD.Barrier();
			MPI::COMM_WORLD.Allgatherv(&dF.GetPtr()[_pd->MyMinEq], _pd->MyNumEqs, MPI::DOUBLE, dF.GetPtr(), &_pd->AllNumEqs[0], &_pd->AllMinEq[0], MPI::DOUBLE);
			#else
			throw new Fatal(_("Solve::_inv_G_times_dF_minus_hKU: SuperLUd is not available"));
			#endif
		}
		else throw new Fatal(_("Solve::_inv_G_times_dF_minus_hKU: Linear solver #%d is not available"),_linsol);
	}

#ifdef DO_DEBUGx
	std::ostringstream oss;
	if (_linsol==LinAlg::LAPACK_T) oss << "LAPACK ";
	else
	{
		if (_linsol==LinAlg::UMFPACK_T ) oss << "UMFPACK ";
		if (_linsol==LinAlg::SuperLU_T ) oss << "SuperLU ";
		//std::cout<<"_T11=\n"<<_T11<<std::endl;
		//std::cout<<"_T12=\n"<<_T12<<std::endl;
		//std::cout<<"_T21=\n"<<_T21<<std::endl;
		//std::cout<<"_T22=\n"<<_T22<<std::endl;
		if (_T11.Size()!=_T11.Top()) throw new Fatal(_("Solve::_inv_G_times_dF: _T11.Size=%d must be equal to _T11.Top=%d"),_T11.Size(),_T11.Top());
		if (_T12.Size()!=_T12.Top()) throw new Fatal(_("Solve::_inv_G_times_dF: _T12.Size=%d must be equal to _T12.Top=%d"),_T12.Size(),_T12.Top());
		if (_T21.Size()!=_T21.Top()) throw new Fatal(_("Solve::_inv_G_times_dF: _T21.Size=%d must be equal to _T21.Top=%d"),_T21.Size(),_T21.Top());
		if (_T22.Size()!=_T22.Top()) throw new Fatal(_("Solve::_inv_G_times_dF: _T22.Size=%d must be equal to _T22.Top=%d"),_T22.Size(),_T22.Top());
		//LinAlg::Matrix<double> D11;                                     G11.GetDense(D11); D11.SetNS(Util::_6_3)
		//LinAlg::Matrix<double> D21; Sparse::Matrix<double,int> G21(_T21); G21.GetDense(D21); D21.SetNS(Util::_6_3)
		//LinAlg::Matrix<double> D22; Sparse::Matrix<double,int> G22(_T22); G22.GetDense(D22); D22.SetNS(Util::_6_3)
		if (_linsol==LinAlg::SuperLUd_T)
		{
			oss<<"[1;33mProcessor #"<<MyID<<"[0m";
			oss<<" MyNumEqs="<<_pd->MyNumEqs<<" MyMinEq="<<_pd->MyMinEq<<" MyMaxEq="<<MyMaxEq<<" MyElements.size="<<MyElements.Size();
			oss<<" MyElements={"; for (size_t i=0; i<MyElements.Size(); ++i) {oss<<MyElements[i];if(i!=MyElements.Size()-1)oss<<",";}; oss<<"}";
		}
		//oss<<"A11 ([1;33m#"<<MyID<<"[0m) == G11 =\n"<<D11<<"\n";
		//oss<<"A21 ([1;33m#"<<MyID<<"[0m) == G21 =\n"<<D21<<"\n";
		//oss<<"A22 ([1;33m#"<<MyID<<"[0m) == G22 =\n"<<D22<<"\n";
		//oss<<"B ([1;33m#"<<MyID<<"[0m) == W   = " <<W<<"\n";
	}
	oss<<" dU={"; for (int i=0; i<dU.Size(); ++i) {oss<<dU(i);if(i!=dU.Size()-1)oss<<",";}; oss<<"}";
	oss<<" dF={"; for (int i=0; i<dF.Size(); ++i) {oss<<dF(i);if(i!=dF.Size()-1)oss<<",";}; oss<<"}"; oss<<"\n";
	std::cout<<oss.str();
#endif

}

inline void Solver::_update_nodes_and_elements(double h, LinAlg::Vector<double> const & dF, LinAlg::Vector<double> const & dU)
{
	// Update all essential values
	for (int i=0; i<_nudofs; ++i) _udofs[i]->EssentialVal += dU(_udofs[i]->EqID);
	for (int i=0; i<_npdofs; ++i) _pdofs[i]->EssentialVal += dU(_pdofs[i]->EqID);

	// Update all elements
	_dF_int.SetValues(0.0);
	for (size_t i=0; i<_g->NElems(); ++i)
	{
		if (_g->Ele(i)->IsActive())
			_g->Ele(i)->UpdateState(h,dU, _dF_int); // sum results into dF_int
	}

	// Calculate residual: resid = dFext-dFint
	_resid = dF - _dF_int;

	// Distribute all pieces of _dF_int to all processors
	#ifdef HAVE_SUPERLUD
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allgatherv(&_dF_int.GetPtr()[_pd->MyMinEq], _pd->MyNumEqs, MPI::DOUBLE, _dF_int.GetPtr(), &_pd->AllNumEqs[0], &_pd->AllMinEq[0], MPI::DOUBLE);
	#endif

	// Update all natural values
	for (int i=0; i<_nudofs; ++i) _udofs[i]->NaturalVal += _dF_int(_udofs[i]->EqID);
	for (int i=0; i<_npdofs; ++i) _pdofs[i]->NaturalVal += _dF_int(_pdofs[i]->EqID);

}

inline void Solver::_backup_nodes_and_elements()
{
	// Backup all nodes
	for (int i=0; i<_nudofs; ++i)
	{
		 _U_bkp(_udofs[i]->EqID) = _udofs[i]->EssentialVal;
		 _F_bkp(_udofs[i]->EqID) = _udofs[i]->NaturalVal;
	}
	for (int i=0; i<_npdofs; ++i)
	{	
		 _U_bkp(_pdofs[i]->EqID) = _pdofs[i]->EssentialVal;
		 _F_bkp(_pdofs[i]->EqID) = _pdofs[i]->NaturalVal;
	}

	// Backup all elements
	for (size_t i=0; i<_g->NElems(); ++i)
		_g->Ele(i)->BackupState();

}

inline void Solver::_restore_nodes_and_elements()
{
	// Restore all nodes
	for (int i=0; i<_nudofs; ++i)
	{
		_udofs[i]->EssentialVal = _U_bkp(_udofs[i]->EqID);
		_udofs[i]->NaturalVal   = _F_bkp(_udofs[i]->EqID);
	}
	for (int i=0; i<_npdofs; ++i)
	{
		_pdofs[i]->EssentialVal = _U_bkp(_pdofs[i]->EqID);
		_pdofs[i]->NaturalVal   = _F_bkp(_pdofs[i]->EqID);
	}

	// Restore all elements
	for (size_t i=0; i<_g->NElems(); ++i)
		_g->Ele(i)->RestoreState();

}

inline double Solver::_norm_essential_vector()
{
	double norm=0.0;
	for (int i=0; i<_nudofs; ++i) norm += pow(_udofs[i]->EssentialVal,2.0);
	for (int i=0; i<_npdofs; ++i) norm += pow(_pdofs[i]->EssentialVal,2.0);
	return sqrt(norm);
}

inline double Solver::_norm_natural_vector()
{
	double norm=0.0;
	for (int i=0; i<_nudofs; ++i) norm += pow(_udofs[i]->NaturalVal,2.0);
	for (int i=0; i<_npdofs; ++i) norm += pow(_pdofs[i]->NaturalVal,2.0);
	return sqrt(norm);
}


/* private */

inline void Solver::_assemb_G_and_hKU(double h)
{
	// Clear DENSE stifness matrix G or Clear SPARSE stifness matrix G (T11,T12,T21,T22)
	if (_linsol==LinAlg::LAPACK_T) _G.SetValues(0.0);
	else
	{
		// Reset Top (position to insert new values) => clear the triplets
		_T11.ResetTop(_npdofs); // leave the "ones" related to the prescribed DOFs unchanged
		_T12.ResetTop();
		_T21.ResetTop();
		_T22.ResetTop();
	}

	// Clear hKU and set alfa (step-controller)
	_hKU.SetValues(0.0);
	_has_hKU    = false;
	double alfa = 1.0;
    
	// Variables to be used during assemblage
	LinAlg::Matrix<double> M;        // Local matrix to be mounted
	LinAlg::Vector<double> V;        // Local vector to be mounted
	Array<size_t>          rows_map; // Map for row positions of components
	Array<size_t>          cols_map; // Map for col positions of components
	Array<size_t>          vect_map; // Map for row positions of components
	Array<bool>            rows_pre; // Map for row positions correspondent to prescribed essential DOFs
	Array<bool>            cols_pre; // Map for col positions correspondent to prescribed essential DOFs
	
	// Loop along elements
	for (size_t i_ele=0; i_ele<_g->NElems(); ++i_ele)
	{
		// Get a poniter to the element
		FEM::Element const * const elem = _g->Ele(i_ele);
		if (elem->IsActive())
		{
			// Assemble zero order matrices
			for (size_t i=0; i<elem->nOrder0Matrices(); i++)
			{
				// Get local matrix (ex.: permeability H) and maps
				elem->Order0MatMap (i, rows_map, cols_map, rows_pre, cols_pre);
				elem->Order0Matrix (i, M);

				// Get local vector related to the zeroth order term of the differential equation (ex.: [U P]^T)
				elem->Order0VecMap (i, vect_map);
				elem->Order0Vector (i, V);
				
				// Compute scaled vector and matrix
				V = h*M*V;
				M = alfa*h*M;

				// Assemble into global
				_has_hKU = true;
				_mount_into_hKU    (V, vect_map);
				_mount_into_global (M, rows_map, cols_map, rows_pre, cols_pre);
			}

			// Assemble first order matrices
			for (size_t i=0; i<elem->nOrder1Matrices(); i++)
			{
				// Get local matrix (ex.: stiffness Ke) and maps
				elem->Order1MatMap (i, rows_map, cols_map, rows_pre, cols_pre);
				elem->Order1Matrix (i, M);

				// Assemble into global
				_mount_into_global (M, rows_map, cols_map, rows_pre, cols_pre);
			}
		}
	}

	#ifdef HAVE_SUPERLUD
	if (_has_hKU)
	{
		// Distribute all pieces of _hKU to all processors
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Allgatherv(&_hKU.GetPtr()[_pd->MyMinEq], _pd->MyNumEqs, MPI::DOUBLE, _hKU.GetPtr(), &_pd->AllNumEqs[0], &_pd->AllMinEq[0], MPI::DOUBLE);
	}
	#endif

	// Check
	#ifndef DNDEBUG
	//_G.SetNS(Util::_6_3);
	//std::cout << "G(global stiffness) = \n" << _G << std::endl;
	for (int i=0; i<_G.Rows(); ++i)
	for (int j=0; j<_G.Cols(); ++j)
		if (_G(i,j)!=_G(i,j)) throw new Fatal (_("Solver::_assemb_G_and_hKU: DENSE stiffness matrix has NaNs"));
	#endif
}

inline void Solver::_mount_into_hKU(LinAlg::Vector<double> const & V, Array<size_t> const & RowsMap)
{
	// Assemble local (LocalVector) DOFs into hKU matrix
	for (int i=0; i<V.Size(); ++i)
		_hKU(RowsMap[i]) += V(i);
}

inline void Solver::_mount_into_global(LinAlg::Matrix<double> const & M, Array<size_t> const & RowsMap, Array<size_t> const & ColsMap, Array<bool> const & RowsPre, Array<bool> const & ColsPre)
{
	// Assemble local matrix (M) into global "global" matrix G
	int I,J; // global indexes
	if (_linsol==LinAlg::LAPACK_T) // DENSE
	{
		for (int i=0; i<M.Rows(); ++i)
		for (int j=0; j<M.Cols(); ++j)
		{
			I = RowsMap[i];
			J = ColsMap[j];
			_G(I,J) += M(i,j);
		}
	}
	else // SPARSE
	{
		for (int i=0; i<M.Rows(); i++)
		for (int j=0; j<M.Cols(); j++)
		{
			I = RowsMap[i];
			J = ColsMap[j];
				 if (!RowsPre[i] && !ColsPre[j]) { _T11.PushEntry(I,J,M(i,j)); }
			else if (!RowsPre[i] &&  ColsPre[j]) { _T12.PushEntry(I,J,M(i,j)); }
			else if ( RowsPre[i] && !ColsPre[j]) { _T21.PushEntry(I,J,M(i,j)); }
			else if ( RowsPre[i] &&  ColsPre[j]) { _T22.PushEntry(I,J,M(i,j)); }
		}
	}
#ifdef DO_DEBUG
	if (M.Rows()!=static_cast<int>(RowsMap.Size())) throw new Fatal(_("Solver::_mount_into_global: M.Rows must be equal to RowsMap.size"));
	if (M.Cols()!=static_cast<int>(ColsMap.Size())) throw new Fatal(_("Solver::_mount_into_global: M.Cols must be equal to ColsMap.size"));
	if (M.Rows()!=static_cast<int>(RowsPre.Size())) throw new Fatal(_("Solver::_mount_into_global: M.Rows must be equal to RowsPre.size"));
	if (M.Cols()!=static_cast<int>(ColsPre.Size())) throw new Fatal(_("Solver::_mount_into_global: M.Cols must be equal to ColsPre.size"));
#endif
}

inline void Solver::_compute_global_size()
{
	// Clear
	_T11_size = 0;
	_T12_size = 0;
	_T21_size = 0;
	_T22_size = 0;

	// Variables to be used during computation of G size
	Array<size_t> rows_map; // Map for row positions of components
	Array<size_t> cols_map; // Map for col positions of components
	Array<size_t> vect_map; // Map for row positions of components
	Array<bool>   rows_pre; // Map for row positions correspondent to prescribed essential DOFs
	Array<bool>   cols_pre; // Map for col positions correspondent to prescribed essential DOFs
	
	// Loop along elements
	for (size_t i_ele=0; i_ele<_g->NElems(); ++i_ele)
	{
		// Get a poniter to an element
		FEM::Element const * const elem = _g->Ele(i_ele);
		if (elem->IsActive())
		{
			// Compute zero order matrices size
			for (size_t i=0; i<elem->nOrder0Matrices(); i++)
			{
				// Get map (ex.: permeability H)
				elem->Order0MatMap(i, rows_map, cols_map, rows_pre, cols_pre);

				// Increase global stiffness pieces size
				_increase_global_size(rows_map, cols_map, rows_pre, cols_pre);
			}

			// Compute first order matrices size
			for (size_t i=0; i<elem->nOrder1Matrices(); i++)
			{
				// Get map (ex.: stiffness Ke)
				elem->Order1MatMap(i, rows_map, cols_map, rows_pre, cols_pre);

				// Increase global stiffness pieces size
				_increase_global_size(rows_map, cols_map, rows_pre, cols_pre);
			}
		}
	}
}

inline void Solver::_increase_global_size(Array<size_t> const & RowsMap, Array<size_t> const & ColsMap, Array<bool> const & RowsPre, Array<bool> const & ColsPre)
{
	// Compute the size of global stiffness pieces T11,T12,T21 an dT22
	for (size_t i=0; i<RowsMap.Size(); i++)
	for (size_t j=0; j<ColsMap.Size(); j++)
	{
			 if (!RowsPre[i] && !ColsPre[j]) { _T11_size++; }
		else if (!RowsPre[i] &&  ColsPre[j]) { _T12_size++; }
		else if ( RowsPre[i] && !ColsPre[j]) { _T21_size++; }
		else if ( RowsPre[i] &&  ColsPre[j]) { _T22_size++; }
	}
}

inline void Solver::_copy_partial_matrix(LinAlg::Matrix<double> const & Source, LinAlg::Vector<int> const & ValidCols, LinAlg::Vector<int> const & ValidRows, LinAlg::Matrix<double> & Target) const
{
	int m = ValidRows.Size();
	int n = ValidCols.Size();
	int s = Source.Rows();
	assert(Target.Rows()==m);
	assert(Target.Cols()==n);
	assert(Source.Cols()==s);
	double const * ptrSource = Source.GetPtr();
	double       * ptrTarget = Target.GetPtr();
	int    const * ptrVR     = ValidRows.GetPtr();
	int    const * ptrVC     = ValidCols.GetPtr();
	//  
	for (int j =0; j<n; j++) 
		for (int i =0; i<m; i++)
			ptrTarget[i+m*j] = ptrSource[ ptrVR[i] + s*ptrVC[j]];

}

inline void Solver::_do_scatter(LinAlg::Matrix<double> const & K,  
                                LinAlg::Vector<double> const & X,
                                LinAlg::Vector<double> const & Y,
                                LinAlg::Matrix<double>       & K11,
                                LinAlg::Matrix<double>       & K12,
                                LinAlg::Matrix<double>       & K21,
                                LinAlg::Matrix<double>       & K22,
                                LinAlg::Vector<double>       & Y2,
                                LinAlg::Vector<double>       & X1)
{
	LinAlg::Vector<int> ueqs(_nudofs);
	LinAlg::Vector<int> peqs(_npdofs);

	for (int i=0; i<_nudofs; ++i) ueqs(i) = _udofs[i]->EqID;
	for (int i=0; i<_npdofs; ++i) peqs(i) = _pdofs[i]->EqID;

	// K11: 
	K11.Resize(_nudofs,_nudofs);
	_copy_partial_matrix(K, ueqs, ueqs, K11);

	// K12:
	K12.Resize(_nudofs,_npdofs);
	_copy_partial_matrix(K, peqs, ueqs, K12);

	// X1
	X1.Resize(_nudofs);
	for (int i=0; i<_nudofs; ++i)
		X1(i) = X(_udofs[i]->EqID);

	// K21:
	K21.Resize(_npdofs,_nudofs);
	_copy_partial_matrix(K, ueqs, peqs, K21);

	// K22:
	K22.Resize(_npdofs,_npdofs);
	_copy_partial_matrix(K, peqs, peqs, K22);

	// Y2
	Y2.Resize(_npdofs);
	for (int i=0; i<_npdofs; ++i)
		Y2(i) = Y(_pdofs[i]->EqID);
}

inline void Solver::_do_gather(LinAlg::Vector<double> const & Y1,
                               LinAlg::Vector<double> const & Y2,
                               LinAlg::Vector<double> const & X1,
                               LinAlg::Vector<double> const & X2,
                               LinAlg::Vector<double>       & Y,
                               LinAlg::Vector<double>       & X) 
{
	for (int i=0; i<_nudofs; ++i)
	{
		Y(_udofs[i]->EqID) = Y1(i); // Y1
		X(_udofs[i]->EqID) = X1(i); // X1
	}
	for (int i=0; i<_npdofs; ++i)
	{
		Y(_pdofs[i]->EqID) = Y2(i); // Y2
		X(_pdofs[i]->EqID) = X2(i); // X2
	}
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Define a pointer to a function that makes (allocate) a new solver
typedef Solver * (*SolverMakerPtr)();

// Typdef of the array map that contains all the pointers to the functions that makes solvers
typedef std::map<String, SolverMakerPtr, std::less<String> > SolverFactory_t;

// Instantiate the array map that contains all the pointers to the functions that makes solvers
SolverFactory_t SolverFactory;

// Allocate a new solver according to a string giving the name of the solver
Solver * AllocSolver(char const * SolverName)
{
	SolverMakerPtr ptr=NULL;
	ptr = SolverFactory[SolverName];
	if (ptr==NULL)
		throw new Fatal(_("FEM::AllocSolver: There is no < %s > solver implemented in this library"), SolverName);
	return (*ptr)();
}

}; // namespace FEM

// Boost::Python
#ifdef USE_BOOST_PYTHON
//{

namespace BPy = boost::python;
class PySolver
{
public:
	           PySolver     (BPy::str const & Name) { _sol = FEM::AllocSolver (BPy::extract<char const *>(Name)()); }
	          ~PySolver     ()                      { if (_sol!=NULL) delete _sol; }
	PySolver & SetGeom      (FEM::Geom & G)         { _sol->SetGeom      (&G);                                return (*this); }
	PySolver & SetLinSol    (BPy::str const & Key)  { _sol->SetLinSol    (BPy::extract<char const *>(Key)()); return (*this); }
	PySolver & SetNumDiv    (int Numdiv)            { _sol->SetNumDiv    (Numdiv);                            return (*this); }
	PySolver & SetDeltaTime (double DeltaTime)      { _sol->SetDeltaTime (DeltaTime);                         return (*this); }
	void       Solve        ()                      { _sol->Solve        (); }
private:
	FEM::Solver * _sol;
}; // class PySolver

//}
#endif // USE_BOOST_PYTHON

#endif // MECHSYS_FEM_SOLVER_H
