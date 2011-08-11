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

#ifndef MECHSYS_FEM_SOLVER_H
#define MECHSYS_FEM_SOLVER_H

// Std Lib
#include <cstring>  // for strcmp
#include <iostream> // for cout

// Blitz++
#include <blitz/tinyvec-et.h>

// MPI
#ifdef HAS_MPI
  #include <mpi.h>
#endif

// MechSys
#include <mechsys/fem/node.h>
#include <mechsys/fem/element.h>
#include <mechsys/fem/domain.h>
#include <mechsys/linalg/sparse_triplet.h>
#include <mechsys/linalg/sparse_matrix.h>
#include <mechsys/linalg/umfpack.h>
#include <mechsys/linalg/mumps.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/numerical/odesolver.h>

namespace FEM
{

class Solver
{
public:
    // typedefs
    typedef void (*pOutFun) (Solver * Sol, void * OutDat); ///< Pointer to output function

    // Constructor
    Solver (Domain & Dom, SDPair const & Flags, pOutFun OutFun=NULL, void * OutDat=NULL,
                                                pOutFun DbgFun=NULL, void * DbgDat=NULL); ///< Allocate solver object

    // Destructor
    virtual ~Solver () {}

    // Methods
    virtual String Name () const =0;

    // Auxiliary methods
    double Timestep (int i, int a=7, double L=100.0, int sch=0, double m=2.0) const; ///< Calculate nonlinear steps

    // Data
    Domain  & Dom;      ///< Domain
    pOutFun   OutFun;   ///< Output function (called during output)
    void    * OutDat;   ///< Debug data (to be used with either OutFun or DbgFun)
    pOutFun   DbgFun;   ///< Debug function (called everytime for some internal (update) methods)
    void    * DbgDat;   ///< Debug data (to be used with either OutFun or DbgFun)
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Solver::Solver (Domain & TheDom, SDPair const & Flags, pOutFun TheOutFun, void * TheOutDat, pOutFun TheDbgFun, void * TheDbgDat)
    : Dom      (TheDom),
      OutFun   (TheOutFun),
      OutDat   (TheOutDat),
      DbgFun   (TheDbgFun),
      DbgDat   (TheDbgDat)
{
}

inline double Solver::Timestep (int i, int a, double L, int sch, double m) const
{
    return (i<a ? pow(2.0,i)/L : (sch==0 ? pow(2.0,i-a) : pow(2.0,a)/L+m*(i-a)));
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


typedef Solver * (*SolverMakerPtr)(Domain & Dom, SDPair const & Flags, Solver::pOutFun OutFun, void * OutDat, Solver::pOutFun DbgFun, void * DbgDat);

typedef std::map<String, SolverMakerPtr> SolverFactory_t;

SolverFactory_t SolverFactory;

Solver * AllocSolver(String const & Name, Domain & Dom, SDPair const & Flags, Solver::pOutFun OutFun=NULL, void * OutDat=NULL, Solver::pOutFun DbgFun=NULL, void * DbgDat=NULL)
{
    SolverFactory_t::iterator it = SolverFactory.find(Name);
    if (it==SolverFactory.end()) throw new Fatal("AllocSolver: '%s' is not available", Name.CStr());
    Solver * ptr = (*it->second)(Dom,Flags,OutFun,OutDat,DbgFun,DbgDat);
    return ptr;
}


}; // namespace FEM

#endif // MECHSYS_FEM_SOLVER_H
