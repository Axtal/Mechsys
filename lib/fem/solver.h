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


/*  __ Solver for (un)coupled problems. Serial/Parallel with Dense/Sparse matrices __
 *
 *  This solver can solve the following kind of problems:
 *    _           _                 _          _
 *   |   K    L1   | / dU/dt \     |   0    0   | / U \     / dF/dt \
 *   |             | |       |  +  |            | |   |  =  |       |
 *   |_  L2   M   _| \ dP/dt /     |_  0    H  _| \ P /     \ dQ/dt /
 *
 *               [C] { dU/dt }  +             [H] { U }  =  { dF/dt }
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
 *    num(CMatrix) = 4
 *    num(HMatrix) = 1
 *
 *  Solve:     ([C] + alpha*h*[H]) * {dU} = {dF} - h*[K]*{U}
 */


#ifndef MECHSYS_FEM_SOLVER_H
#define MECHSYS_FEM_SOLVER_H

// Std Lib
#include <cstring>  // for strcmp
#include <ctime>    // for clock
#include <iostream> // for cout

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "fem/data.h"
#include "fem/output.h"
#include "fem/parallel.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "util/fatal.h"
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
	// enums
	enum Solver_T { FE_T, AME_T }; ///< FE:Forward-Euler, AME:Automatic-Modified-Euler

	// Constructor
	Solver (FEM::Data & D)                       { _initialize (&D, NULL);    }
	Solver (FEM::Data & D, char const * FileKey) { _initialize (&D, FileKey); }

	// Destructor
	~Solver() { if (_out!=NULL) delete _out; }

	// Methods
	Solver * SetPData      (FEM::PData * PD)  { _pd = PD; return this; }                             ///< Set structure with data for parallel computation
	Solver * SetType       (char const * Type);                                                      ///< Set solver type: FE:Forward-Euler, AME:Automatic-Modified-Euler
	Solver * SetCte        (char const * Key, double Value);                                         ///< Set scheme constant such as number of subincrements, DTOL, etc.
	Solver * SetLinSol     (char const * Key);                                                       ///< LinSolName: LA=LAPACK, UM=UMFPACK, SLU=SuperLU, SLUd=SuperLUd
	void     Solve         (int NDiv=1, double DTime=0.0, bool ClearDisp=false);                                           ///< Solve ([C]+alpha*h*[K])*{dU} = {dF}-h*[K]*{U} for the boundary conditions defined inside the nodes array
	void     SolveWithInfo (int NDiv=1, double DTime=0.0, int iStage=0, char const * MoreInfo=NULL, bool ClearDisp=false); ///< Solve and show time/resid information
	void     DynSolve      (double tIni, double tFin, double h, double dtOut, size_t MaxIt=10);      ///< TODO
	double   F             (double t) const;                                                         ///< TODO

	// Access methods
	double Time () const { return _data->Time(); } ///< Return current time
	size_t nDOF () const { return _ndofs;        } ///< Number of DOFs

	// Access output
	Output const & Out () const { return (*_out); }

#ifdef USE_BOOST_PYTHON
//{
	         Solver          (FEM::Data & D, BPy::str const & FileKey)         { _initialize   (&D, BPy::extract<char const *>(FileKey)()); }
	Solver & PySetType       (BPy::str const & Type)                           { SetType       (BPy::extract<char const *>(Type)());     return (*this); }
	Solver & PySetCte        (BPy::str const & Key, double Val)                { SetCte        (BPy::extract<char const *>(Key)(), Val); return (*this); }
	Solver & PySetLinSol     (BPy::str const & Key)                            { SetLinSol     (BPy::extract<char const *>(Key)());      return (*this); }
	void     PySolveWithInfo (int NDiv, double DTime, int iStage,
	                          BPy::str const & MoreInfo, bool ClearDisp=false) { SolveWithInfo (NDiv, DTime, iStage, BPy::extract<char const *>(MoreInfo)(), ClearDisp); }
//}
#endif // USE_BOOST_PYTHON

private:
	// Data
	FEM::Data              * _data;    ///< FEM Data structure
	FEM::PData             * _pd;      ///< Structure with data for parallel computation
	Solver_T                 _type;    ///< Solver type: FE:Forward-Euler, AME:Automatic-Modified-Euler
	int                      _inc;     ///< Current increment
	int                      _ndofs;   ///< Total number of DOFs
	LinAlg::LinSol_T         _linsol;  ///< The type of the linear solver to be used in the inversion of G11. Solves: X<-inv(A)*B
	bool                     _has_hHU; ///< Flag which says if any element has to contribute to the hHU vector. If _has_hHU==false, there is no need for the hHU vector, because there are no HMatrices in this stage of the simulation.
	bool                     _has_dFs; ///< Has delta Force (s)tar (to be added to dFext)
	int                      _nudofs;  ///< Total number of unknown DOFs of a current stage
	int                      _npdofs;  ///< Total number of prescribed DOFs of a current stage
	Array<Node::DOF*>        _udofs;   ///< Unknown DOFs
	Array<Node::DOF*>        _pdofs;   ///< Prescribed DOFs
	LinAlg::Vector<double>   _U_bkp;   ///< Place to store a temporary state of essential (displacements) values (backup).
	LinAlg::Vector<double>   _F_bkp;   ///< Place to store a temporary state of natural (forces) values (backup).
	LinAlg::Vector<double>   _dF_ext;  ///< Increment of natural values, divided by the number of increments, which update the current state towards the condition at the end the stage being solved.
	LinAlg::Vector<double>   _dU_ext;  ///< Increment of essential values, divided by the number of increments, which update the current state towards the condition at the end the stage being solved.
	LinAlg::Vector<double>   _dF_int;  ///< Increment of internal natural (forces) values, divided by the number of increments, correspondent to the increment of external forces.
	LinAlg::Vector<double>   _resid;   ///< Residual: resid = dFext - dFint
	LinAlg::Vector<double>   _hHU;     ///< Linearized independent term of the differential equation.
	Output                 * _out;     ///< Write output file

	// Data for Forward-Euler
	int    _nSI;        ///< Number of sub-increments
	double _norm_resid; ///< Residual Norm(dFext - dFint)

	// Data for Automatic-Modified-Euler
	int                    _maxSI; ///< Max number of subincrements
	double                 _DTOL;  ///< Global solver local tolerance
	double                 _dTini; ///< Delta T initial
	double                 _mMin;  ///< m (multiplier) mininum
	double                 _mMax;  ///< m (multiplier) maximum
	double                 _mCoef; ///< m coefficient
	double                 _ZTOL;  ///< Zero tolerance
	bool                   _Cconv; ///< Check convergence ?
	double                 _Rerr;  ///< Relative error
	LinAlg::Vector<double> _dF_1;
	LinAlg::Vector<double> _dU_1;
	LinAlg::Vector<double> _dF_2;
	LinAlg::Vector<double> _dU_2;
	LinAlg::Vector<double> _dU_ME;
	LinAlg::Vector<double> _dF_ME;
	LinAlg::Vector<double> _Err_U;
	LinAlg::Vector<double> _Err_F;

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

	// Methods
	void   _initialize                 (FEM::Data * D, char const * FileKey);                                ///< Initialize the solver
	void   _inv_G_times_dF_minus_hHU   (double h, LinAlg::Vector<double> & dF, LinAlg::Vector<double> & dU); ///< Compute (linear solver) inv(G)*(dF-hHU), where G may be assembled by [C] and [H] matrices
	void   _update_nodes_and_elements  (double h,LinAlg::Vector<double> const & dU);                         ///< Update nodes essential/natural values and elements internal values (such as stresses/strains)
	void   _backup_nodes_and_elements  ();                                                                   ///< Backup nodes essential/natural values and elements internal values (such as stresses/strains)
	void   _restore_nodes_and_elements ();                                                                   ///< Restore nodes essential/natural values and elements internal values (such as stresses/strains)
	double _norm_essential_vector      ();                                                                   ///< Compute the Euclidian norm of the essential (displacements) vector
	double _norm_natural_vector        ();                                                                   ///< Compute the Euclidian norm of the natural (forces) vector

	// Private methods
	void _assemb_G_and_hHU     (double h);                                                                                                                                                   ///< Assemble matrix G for the actual state of elements/nodes. Note h == TimeInc
	void _mount_into_global    (LinAlg::Matrix<double> const & M, Array<size_t> const & RowsMap, Array<size_t> const & ColsMap, Array<bool> const & RowsPre, Array<bool> const & ColsPre);   ///< Add a local matrix, such as K,L1,L2 or M into the global DENSE matrix G or into the global SPARSE matrix G (T11,T12,T21,T22)
	void _compute_global_size  ();                                                                                                                                                           ///< Compute the sizes of T11,T12,T21 and T22
	void _increase_global_size (Array<size_t> const & RowsMap, Array<size_t> const & ColsMap, Array<bool> const & RowsPre, Array<bool> const & ColsPre);                                     ///< Auxiliar method to compute the sizes of T11,T12,T21 and T22
	void _copy_partial_matrix  (LinAlg::Matrix<double> const & Source, LinAlg::Vector<int> const & ValidCols, LinAlg::Vector<int> const & ValidRows, LinAlg::Matrix<double> & Target) const; ///< Auxiliar method to copy a piece of a matrix inside other matrix
	void _do_scatter           (LinAlg::Matrix<double> const & K, LinAlg::Vector<double> const & X, LinAlg::Vector<double> const & Y, LinAlg::Matrix<double> & K11, LinAlg::Matrix<double> & K12, LinAlg::Matrix<double> & K21, LinAlg::Matrix<double> & K22, LinAlg::Vector<double> & Y2, LinAlg::Vector<double> & X1); ///< Scatter "global" matrices to small pieces
	void _do_gather            (LinAlg::Vector<double> const & Y1, LinAlg::Vector<double> const & Y2, LinAlg::Vector<double> const & X1, LinAlg::Vector<double> const & X2, LinAlg::Vector<double> & Y, LinAlg::Vector<double> & X);                                                                                     ///< Gather small pieces into "global" matrices

	// FE and AME schemes
	void _fe_solve_for_an_increment  (double dTime);
	void _ame_solve_for_an_increment (double dTime);
	void _calc_resid                 (LinAlg::Vector<double> & dFext);

}; // class Solver


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline Solver * Solver::SetType(char const * Type)
{
	     if (strcmp(Type,"FE" )==0) _type = FE_T;
	else if (strcmp(Type,"AME")==0) _type = AME_T;
	else throw new Fatal("Solver::SetType: Solver type==%s is invalid",Type);
	return this;
}

inline Solver * Solver::SetCte(char const * Key, double Value)
{
	if (_type==FE_T)
	{
		if (strcmp(Key,"nSI")==0) _nSI = Value;
		else throw new Fatal("Solver::SetCte: ForwardEuler: This solver does not have a constant named %s",Key);
	}
	else
	{
		     if (strcmp(Key,"maxSI")==0) _maxSI = Value;
		else if (strcmp(Key,"DTOL" )==0) _DTOL  = Value;
		else if (strcmp(Key,"dTini")==0) _dTini = Value;
		else if (strcmp(Key,"mMin" )==0) _mMin  = Value;
		else if (strcmp(Key,"mMax" )==0) _mMax  = Value;
		else if (strcmp(Key,"mCoef")==0) _mCoef = Value;
		else if (strcmp(Key,"ZTOL" )==0) _ZTOL  = Value;
		else if (strcmp(Key,"Cconv")==0) _Cconv = Value;
		else throw new Fatal("Solver::SetCte: Automatic-Modified-Euler: This solver does not have a constant named %s",Key);
	}
	return this;
}

inline Solver * Solver::SetLinSol(char const * Key)
{
	     if (strcmp(Key,"LA"  )==0) _linsol=LinAlg::LAPACK_T;
	else if (strcmp(Key,"UM"  )==0) _linsol=LinAlg::UMFPACK_T;
	else if (strcmp(Key,"SLU" )==0) _linsol=LinAlg::SuperLU_T;
	else if (strcmp(Key,"SLUd")==0) _linsol=LinAlg::SuperLU_T;
	else throw new Fatal("Solver::SetLinSol: Solver key==%s is invalid",Key);
	return this;
}

inline void Solver::Solve(int NDiv, double DTime, bool ClearDisp)
{
	// Solve:    ([C] + alpha*h*[K]) * {dU} = {dF} - h*[K]*{U}
	
	// Check geometry
	if (_data->Check()==false) throw new Fatal("Solver::Solve: FEM Geometry was not properly set (pointers may be NULL)");

	// Loop over all active elements
	_has_dFs = false;
	for (size_t i=0; i<_data->NElems(); ++i)
	{
		if (_data->Ele(i)->IsActive())
		{
			// Check if active elements are properly initialized
			String msg;
			if (_data->Ele(i)->Check(msg)==false)
				throw new Fatal("Solver::Solve Element # %d was not properly initialized. Error: %s",i,msg.CStr());

			// Check if element has component to be added to F_ext
			if (_data->Ele(i)->HasDFs()) _has_dFs = true;
		}
	}

	// Number of DOFs
	_ndofs  = 0; // total
	_nudofs = 0; // unknown
	_npdofs = 0; // prescribed

	// Set the ID (equation number) of all DOFs and setup _dF_ext and _dU_ext vectors
	_udofs.Resize (0);
	_pdofs.Resize (0);
	for (size_t i=0; i<_data->NNodes(); ++i)
	{
		if (_data->Nod(i)->nSharedBy()>0)
		{
			// Loop along DOFs inside Node [i]
			for (size_t j=0; j<_data->Nod(i)->nDOF(); ++j)
			{
				// Save pointers to the DOFs
				if (_data->Nod(i)->DOFVar(j).IsEssenPresc) _pdofs.Push(&_data->Nod(i)->DOFVar(j));
				else                                       _udofs.Push(&_data->Nod(i)->DOFVar(j));
				
				// Assing an equation ID to DOF [j]
				_data->Nod(i)->DOFVar(j).EqID = _ndofs;

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
	_hHU   .Resize (_ndofs);
	_dF_int.Resize (_ndofs);
	_resid .Resize (_ndofs);

	// Set the vectors with increments of boundary conditions
	_dF_ext.SetValues (0.0);
	_dU_ext.SetValues (0.0);
	for (int i=0; i<_nudofs; ++i)
	{
		_dF_ext(_udofs[i]->EqID) = _udofs[i]->NaturalBry   / NDiv;
		_dU_ext(_udofs[i]->EqID) = _udofs[i]->EssentialBry / NDiv;
	}
	for (int i=0; i<_npdofs; ++i)
	{
		_dF_ext(_pdofs[i]->EqID) = _pdofs[i]->NaturalBry   / NDiv;
		_dU_ext(_pdofs[i]->EqID) = _pdofs[i]->EssentialBry / NDiv;
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

	// Number of divisions for each increment
	double dTime = DTime / NDiv; // Time increment

	// Solve
	if (_type==FE_T) { for (_inc=0; _inc<NDiv; ++_inc) _fe_solve_for_an_increment (dTime); }
	else             { for (_inc=0; _inc<NDiv; ++_inc) _ame_solve_for_an_increment(dTime); }

	// Clear boundary conditions for a next stage
	for (size_t i=0; i<_data->NNodes(); ++i) _data->Nod(i)->ClearBryValues();

	// Clear displacements
	if (ClearDisp) _data->ClearDisp();

	// Write output
	if (_out!=NULL) _out->Write();
}

inline void Solver::SolveWithInfo(int NDiv, double DTime, int iStage, char const * MoreInfo, bool ClearDisp)
{
	double start = std::clock();
	Solve (NDiv, DTime, ClearDisp); // <<<<<<< Solve
	double total = std::clock() - start;
	std::cout << "[1;33m\n--- Stage # " << iStage << " --------------------------------------------------[0m\n";
	if (MoreInfo!=NULL) std::cout << MoreInfo;
	std::cout << "[1;36m    Time elapsed (FE solution) = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
	std::cout << "[1;35m    Norm(Resid=DFext-DFint)    = " << _norm_resid << "[0m\n";
	std::cout << "[1;32m    Number of DOFs             = " << nDOF()      << "[0m" << std::endl;
}

inline double Solver::F(double t) const
{
	if (t<0.1) return t*10000.0/0.1;
	else       return   10000.0;
}

inline void Solver::DynSolve(double tIni, double tFin, double h, double dtOut, size_t MaxIt)
{
	/*
	double t    = tIni;
	double tout = t + dtOut;
	while (t<tFin)
	{
		M        = calc_M();
		C_n      = calc_C(u_n, v_n);
		K_n      = calc_K(u_n, v_n);
		Fint_n   = C_n * v_n + K_n * u_n;                   // Eq. (2.363)
		v_np1    = v_n + (1.0-gam)*h*a_n;                   // Eq. (2.360)
		u_np1    = u_n + h*v_n + 0.5*h*h*(1.0-2.0*bet)*a_n; // Eq. (2.361)
		a_np1    = 0.0;                                     // Eq. (2.381)
		Fext_n   = F(t);
		Fext_np1 = F(t+h);
		for (size_t i=0; i<MaxIt; ++i)
		{
			C_np1   = calc_C(u_np1, v_np1);                                             // Eq. (2.365)
			K_np1   = calc_K(u_np1, v_np1);                                             // Eq. (2.366)
			Fint_st = C_np1 * v_np1 + K_np1 * u_np1;                                    // Eq. (2.384)
			F_st    = (1.0+alp)*Fext_np1 - alp*Fext_n - (1.0+alp)*Fint_st + alp*Fint_n; // Eq. (2.383)
			M_st    = M + (1.0+alp)*gam*h*C_np1 + (1.0+alp)*bet*h*h*K_np1;              // Eq. (2.373)
			R       = F_st - M * a_np1;                                                 // Eq. (2.382)
			normR   = Norm(R);
			if (normR<Tol) break;
			linsol(M_st,R,da); // Eq. (2.385)
			//if (i==1)
			//{
				//for (j in 1:length(u_given))
					//if (u_given[j]) Da[j] = (u_n[j]-u_np1[j])/(bet*hh)
			//}
			v_np1 += gam*h*da;   // Eq. (2.386)
			u_np1 += bet*h*h*da; // Eq. (2.387)
			a_np1 += da;         // Eq. (2.380)
		}
		t += h;
		if (t>=tOut)
		{
			std::cout << t << v_n << std::endl;
			tout = tout + dtOut;
		}
	}
	*/
}


/* private */

inline void Solver::_initialize(FEM::Data * D, char const * FileKey)
{
	_data    = D;
	_pd      = NULL;
	_inc     = 0;
	_ndofs   = 0;
	_has_hHU = false;
	_has_dFs = false;
	SetType   ("FE");
	SetLinSol ("UM");
	if (FileKey==NULL) _out = NULL;
	else               _out = new Output (_data, FileKey);

	// Constants
	_nSI    = 1;
	_maxSI  = 400;
	_DTOL   = 1.0e-3;
	_dTini  = 0.001;
	_mMin   = 0.01;
	_mMax   = 10;
	_mCoef  = 0.7;
	_ZTOL   = 1.0e-7;
	_Cconv  = true;
	_Rerr   = 0.0;
}

inline void Solver::_inv_G_times_dF_minus_hHU(double h, LinAlg::Vector<double> & dF, LinAlg::Vector<double> & dU)
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

	               [C] { dU/dt }  +             [H] { U }  =  { dF/dt }
	 */

	// 0) Assemble the global stiffness matrix G and the hHU vector
	_assemb_G_and_hHU(h);                      // G  <- C + h*alpha*K
	if (_has_hHU) LinAlg::Axpy(-1.0,_hHU, dF); // dF <- dF - hHU

	// Add extra force components to dFext
	if (_has_dFs)
	{
		// Add to Fext
		for (size_t i=0; i<_data->NElems(); ++i)
			if (_data->Ele(i)->IsActive()) _data->Ele(i)->AddToDFext (_data->Time(), h, dF);
	}

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

			//std::cout << "G11 = \n" << G11 << std::endl;

			// 4) Solve for {dU1} and {dF2}
			LinAlg::Vector<double> W(dF1);                    // W <- dF1 (workspace)
			LinAlg::Gemv(-1.0,G12,dU2,1.0,W);               // W <- -1.0*G12*dU2 + 1.0*dF1
			LinAlg::Gesv(G11, W);                           // W <- inv(G11)*W  (G11 is lost)
			LinAlg::Vector<double> & dU1 = W;                 // dU1 <- W (reference)
			LinAlg::Vector<double> dF2(dU2.Size());           // allocate dF2
			LinAlg::Gemvpmv(1.0,G21,dU1, 1.0,G22,dU2, dF2); // dF2 <- 1*G21*dU1 + 1*G22*dU2

			// 5) Gather dU1, dU2, dF1 and dF2 into dU and dF
			_do_gather(dU1,dU2, dF1,dF2, dU,dF);

			//std::cout << "dF  = \n" << dF << std::endl;
			//std::cout << "dU  = \n" << dU << std::endl;

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
			throw new Fatal(_("Solve::_inv_G_times_dF_minus_hHU: UMFPACK is not available"));
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
			throw new Fatal(_("Solve::_inv_G_times_dF_minus_hHU: SuperLU is not available"));
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
			throw new Fatal(_("Solve::_inv_G_times_dF_minus_hHU: SuperLUd is not available"));
			#endif
		}
		else throw new Fatal(_("Solve::_inv_G_times_dF_minus_hHU: Linear solver #%d is not available"),_linsol);
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

inline void Solver::_update_nodes_and_elements(double h, LinAlg::Vector<double> const & dU)
{
	// Update all essential values
	for (int i=0; i<_nudofs; ++i) _udofs[i]->EssentialVal += dU(_udofs[i]->EqID);
	for (int i=0; i<_npdofs; ++i) _pdofs[i]->EssentialVal += dU(_pdofs[i]->EqID);

	// Update all elements
	_dF_int.SetValues(0.0);
	for (size_t i=0; i<_data->NElems(); ++i)
	{
		if (_data->Ele(i)->IsActive())
			_data->Ele(i)->Update (_data->Time(), h, dU, _dF_int); // sum results into dF_int
	}

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
	for (size_t i=0; i<_data->NElems(); ++i)
		_data->Ele(i)->Backup ();

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
	for (size_t i=0; i<_data->NElems(); ++i)
		_data->Ele(i)->Restore ();

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

inline void Solver::_assemb_G_and_hHU(double h)
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

	// Clear hHU and set alfa (step-controller)
	_hHU.SetValues(0.0);
	_has_hHU    = false;
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
	for (size_t i_ele=0; i_ele<_data->NElems(); ++i_ele)
	{
		FEM::Element const * const elem = _data->Ele(i_ele);
		if (elem->IsActive())
		{
			// Assemble [H] matrices
			for (size_t i=0; i<elem->NHMats(); i++)
			{
				elem->HMatMap (i, rows_map, cols_map, rows_pre, cols_pre);
				elem->HMatrix (i, M);
				M = alfa*h*M;
				_mount_into_global (M, rows_map, cols_map, rows_pre, cols_pre);
			}

			// Assemble [C] matrices
			for (size_t i=0; i<elem->NCMats(); i++)
			{
				elem->CMatMap (i, rows_map, cols_map, rows_pre, cols_pre);
				elem->CMatrix (i, M);
				_mount_into_global (M, rows_map, cols_map, rows_pre, cols_pre);
			}

			// Assemble [U] Vectors
			size_t n_RHS = elem->NUVecs();
			for (size_t i=0; i<n_RHS; i++)
			{
				elem->UVecMap (i, vect_map);
				elem->UVector (i, V);
				_has_hHU = true;
				for (int j=0; j<V.Size(); ++j) _hHU(vect_map[j]) += h*V(j);
			}
		}
	}

	#ifdef HAVE_SUPERLUD
	if (_has_hHU)
	{
		// Distribute all pieces of _hHU to all processors
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Allgatherv(&_hHU.GetPtr()[_pd->MyMinEq], _pd->MyNumEqs, MPI::DOUBLE, _hHU.GetPtr(), &_pd->AllNumEqs[0], &_pd->AllMinEq[0], MPI::DOUBLE);
	}
	#endif

	// Check
	#ifndef DNDEBUG
	//_G.SetNS(Util::_6_3);
	//std::cout << "G(global stiffness) = \n" << _G << std::endl;
	/*
	for (int i=0; i<_G.Rows(); ++i)
	for (int j=0; j<_G.Cols(); ++j)
		if (_G(i,j)!=_G(i,j)) throw new Fatal (_("Solver::_assemb_G_and_hHU: DENSE stiffness matrix has NaNs"));
	*/
	#endif
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
	for (size_t i_ele=0; i_ele<_data->NElems(); ++i_ele)
	{
		// Get a poniter to an element
		FEM::Element const * const elem = _data->Ele(i_ele);
		if (elem->IsActive())
		{
			// Compute [H] matrices size
			for (size_t i=0; i<elem->NHMats(); i++)
			{
				elem->HMatMap (i, rows_map, cols_map, rows_pre, cols_pre);
				_increase_global_size (rows_map, cols_map, rows_pre, cols_pre);
			}

			// Compute [C] matrices size
			for (size_t i=0; i<elem->NCMats(); i++)
			{
				elem->CMatMap (i, rows_map, cols_map, rows_pre, cols_pre);
				_increase_global_size (rows_map, cols_map, rows_pre, cols_pre);
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
	if (Target.Rows()!=m) throw new Fatal("Solver::_copy_partial_matrix: Target.Rows()==%d must be equal to %d",Target.Rows(),m);
	if (Target.Cols()!=n) throw new Fatal("Solver::_copy_partial_matrix: Target.Cols()==%d must be equal to %d",Target.Cols(),n);
	if (Source.Cols()!=s) throw new Fatal("Solver::_copy_partial_matrix: Source.Cols()==%d must be equal to %d",Source.Cols(),s);
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

inline void Solver::_fe_solve_for_an_increment(double dTime)
{
	// Allocate and fill dF_ext and dU_ext
	int ndofs = _dF_ext.Size();
	LinAlg::Vector<double> dF_ext(ndofs);
	LinAlg::Vector<double> dU_ext(ndofs);
	LinAlg::CopyScal(1.0/_nSI,_dF_ext, dF_ext); // dF_ext <- _dF_ext/_nSI
	LinAlg::CopyScal(1.0/_nSI,_dU_ext, dU_ext); // dU_ext <- _dU_ext/_nSI
	double h = dTime/_nSI;

	// Start
	for (int i=0; i<_nSI; ++i)
	{
		// Assemble G matrix and calculate dU_ext
		_inv_G_times_dF_minus_hHU(h, dF_ext, dU_ext); // dU_ext <- inv(G)*(dF_ext - hHU)

		// Update nodes and elements state
		_update_nodes_and_elements(h, dU_ext); // AND calculate _resid

		// Update time
		 _data->UpdateTime (h);
	}

	// Calculate residual
	_calc_resid (dF_ext);
}

inline void Solver::_ame_solve_for_an_increment(double dTime)
{
	// Allocate and fill local dF_ext and dU_ext
	LinAlg::Vector<double> dF_ext(_dF_ext);
	LinAlg::Vector<double> dU_ext(_dU_ext);

	// Allocate auxiliar vectors
	if (_inc==0)
	{
		int ndofs = dF_ext.Size();
		_dF_1 .Resize (ndofs);
		_dU_1 .Resize (ndofs);
		_dF_2 .Resize (ndofs);
		_dU_2 .Resize (ndofs);
		_dU_ME.Resize (ndofs);
		_dF_ME.Resize (ndofs);
		_Err_U.Resize (ndofs);
		_Err_F.Resize (ndofs);
	}

	// Start substeps
	double  T = 0.0;
	double dT = _dTini;
	for (int k=0; k<_maxSI; ++k)
	{
		if (T>=1.0)
		{
			_calc_resid (dF_ext);
			return;
		}

		// Sub-divide timestep
		double h = dT*dTime;

		// Calculate scaled increment vectors for this sub-step
		LinAlg::CopyScal(dT,dF_ext, _dF_1); // _dF_1 <- dT*dF_ext
		LinAlg::CopyScal(dT,dU_ext, _dU_1); // _dU_1 <- dT*dU_ext;
		LinAlg::CopyScal(dT,dF_ext, _dF_2); // _dF_2 <- dT*dF_ext
		LinAlg::CopyScal(dT,dU_ext, _dU_2); // _dU_2 <- dT*dU_ext;

		// Backup element state (needed when updating disp. state for a ME increment)
		_backup_nodes_and_elements();

		// Forward-Euler: Assemble G matrix and calculate _dU_1
		_inv_G_times_dF_minus_hHU(h, _dF_1, _dU_1); // _dU_1 <- inv(G)*(dF_ext - hHU)
	
		// Forward-Euler: update nodes and elements state
		_update_nodes_and_elements(h, _dU_1); // AND calculate _resid

		// Update time
		_data->UpdateTime (h);

		// Modified-Euler: Assemble G matrix and calculate dU_2
		_inv_G_times_dF_minus_hHU(h, _dF_2, _dU_2); // _dU_2 <- inv(G)*(dF_ext - hKun)
	
		// Save the norm of essential and natural vectors
		double normU = _norm_essential_vector();
		double normF = _norm_natural_vector();

		// Calculate local error
		if (normF<=_ZTOL) throw new Fatal(_("AutoME::_do_solve_for_an_increment: k=%d: normF=%e cannot be equal to ZTOL (%e)"),k,normF,_ZTOL);
		LinAlg::AddScaled (0.5,_dU_2, -0.5,_dU_1, _Err_U); // Error on U (%)
		LinAlg::AddScaled (0.5,_dF_2, -0.5,_dF_1, _Err_F); // Error on F
		double Rerr_U = LinAlg::Norm(_Err_U)/(1.0+normU);  // (R)elative error on U
		double Rerr_F = LinAlg::Norm(_Err_F)/normF;        // (R)elative error on F
		_Rerr = (Rerr_U>Rerr_F ? Rerr_U : Rerr_F);      // (R)elative error
		double m = _mCoef*sqrt(_DTOL/_Rerr);

		// Restore nodes and element for initial disp. state given at the start of the increment
		_restore_nodes_and_elements();

		// Restore time
		_data->UpdateTime (-h);

		if (_Rerr<=_DTOL)
		{
			// Calculate Modified-Euler force and displacement increment vectors
			LinAlg::AddScaled(0.5,_dU_1, 0.5,_dU_2, _dU_ME);
			LinAlg::AddScaled(0.5,_dF_1, 0.5,_dF_2, _dF_ME);

			// Update nodes and elements state for a Modified-Euler evaluation of displacements
			_update_nodes_and_elements(h, _dU_ME); // AND calculate _resid

			// Update time
			 _data->UpdateTime (h);

			// Next pseudo time
			T = T + dT;
			if (m>_mMax) m=_mMax;
		}
		else
			if (m<_mMin) m=_mMin;

		// Next substep size
		dT = m*dT;
		if (dT>1.0-T) dT=1.0-T;
	}
	if (_Cconv) throw new Fatal(_("AutoME::_do_solve_for_an_increment: did not converge for %d substeps"), _maxSI);
}

void Solver::_calc_resid(LinAlg::Vector<double> & dFext)
{
	/*
	LinAlg::Axpy      (+1.0,_hHU, dFext);                 //  dFext <- dFext + hHU
	LinAlg::AddScaled (1.0,dFext, -1.0,_dF_int, _resid);  // _resid <- dFext - dFint
	*/

	if (_has_hHU) dFext += _hHU;

	_resid      = dFext - _dF_int;
	_norm_resid = LinAlg::Norm(_resid);
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

#endif // MECHSYS_FEM_SOLVER_H
