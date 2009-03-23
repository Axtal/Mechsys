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

// Std Lib
#include <iostream>
#include <stdlib.h>

// MechSys
#include "lbm/lattice.h"
#include "lbm/mixture.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	// Allocate lattice
	LBM::Mixture m( /*FileKey*/ "multi", /*Is3D*/false, /*NComp*/2,/*Nx*/50, /*Ny*/50);

	// Set walls (top and bottom)
	// Lattice 0
	for (size_t i=0; i<m.GetLattice(0)->Top()   .Size(); ++i) m.GetLattice(0)->Top()[i]->SetSolid();
	for (size_t i=0; i<m.GetLattice(0)->Bottom().Size(); ++i) m.GetLattice(0)->Bottom()[i]->SetSolid();
	// Lattice 1
	for (size_t i=0; i<m.GetLattice(1)->Top()   .Size(); ++i) m.GetLattice(1)->Top()[i]->SetSolid();
	for (size_t i=0; i<m.GetLattice(1)->Bottom().Size(); ++i) m.GetLattice(1)->Bottom()[i]->SetSolid();

	//// Initialize
	//for (size_t i=0; i<m.Nx(); ++i)
	//for (size_t j=0; j<m.Ny(); ++j)
	//{
	//	//double rho0 = 1.4 +(.02*rand())/RAND_MAX;
	//	Vec3_t V; V = 0.0, 0.0, 0.0;
	//	m.GetLattice(0)->GetCell(i,j)->Initialize (0.1, V);
	//	m.GetLattice(1)->GetCell(i,j)->Initialize (1.0, V);
	//}

	// lattice[0] : Water
	// lattice[1] : Oil

	// Properties

	// Set inner drop
	int obsX = m.Nx()/2;   // x position
	int obsY = m.Ny()/2;   // y position
	int radius = m.Nx()/5; // Inital drop radius

	for (size_t i=0; i<m.Nx(); ++i)
	for (size_t j=0; j<m.Ny(); ++j)
	{
		Vec3_t V;  V = 0.0, 0.0, 0.0;
		if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(radius,2.0)) // circle equation
		{
			m.GetLattice(0)->GetCell(i,j)->Initialize (0.01, V);
			m.GetLattice(1)->GetCell(i,j)->Initialize (1.0, V);
		}
		else
		{
			m.GetLattice(0)->GetCell(i,j)->Initialize (1.0, V);
			m.GetLattice(1)->GetCell(i,j)->Initialize (0.01, V);
		}
	}

	//// Properties
	m.GetLattice(0)->SetTau(1.0)->SetG(-3.0)->SetGSolid(-0.0);
	m.GetLattice(1)->SetTau(1.0)->SetG(-1.0)->SetGSolid(-0.0);
	
	m.SetGravity(0.0, -0.0005, 0.0);
	m.SetMixG(0.5);

	// Solve
	m.Solve(/*tIni*/0.0, /*tFin*/5000.0, /*dt*/1.0, /*dtOut*/25.0);

}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
