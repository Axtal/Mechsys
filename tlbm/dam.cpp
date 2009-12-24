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
#include <mechsys/lbm/lattice.h>
#include <mechsys/lbm/mixture.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	//allocate viscosity
	double nu[2];
	nu[0]=0.5*1./6.;
	nu[1]=1./6.;

	// Allocate lattice
	LBM::Mixture m( /*FileKey*/ "dam", /*Is3D*/false, /*NComp*/2,nu,/*Nx*/100, /*Ny*/150,1,1,1.);

	// Set walls (top and bottom)
	// Lattice 0
	for (size_t i=0; i<m.GetLattice(0)->Top()   .Size(); ++i) m.GetLattice(0)->Top()[i]->SetSolid();
	for (size_t i=0; i<m.GetLattice(0)->Bottom().Size(); ++i) m.GetLattice(0)->Bottom()[i]->SetSolid();		
	for (size_t i=0; i<m.GetLattice(0)->Left()  .Size(); ++i) m.GetLattice(0)->Left()[i]->SetSolid();
	for (size_t i=0; i<m.GetLattice(0)->Right() .Size(); ++i) m.GetLattice(0)->Right()[i]->SetSolid();	

	// Lattice 1
	for (size_t i=0; i<m.GetLattice(1)->Top()   .Size(); ++i) m.GetLattice(1)->Top()[i]->SetSolid();
	for (size_t i=0; i<m.GetLattice(1)->Bottom().Size(); ++i) m.GetLattice(1)->Bottom()[i]->SetSolid();
	for (size_t i=0; i<m.GetLattice(1)->Left()  .Size(); ++i) m.GetLattice(1)->Left()[i]->SetSolid();
	for (size_t i=0; i<m.GetLattice(1)->Right() .Size(); ++i) m.GetLattice(1)->Right()[i]->SetSolid();	
		
	for (size_t i=0; i<m.Nx(); ++i)
	for (size_t j=0; j<m.Ny(); ++j)
	{
		Vec3_t v0;  v0 = 0.0, 0.0, 0.0;
		if ((i>m.Nx()/2-10)&&(i<(m.Nx()/2+10))&&(j<m.Ny()/5)) {
			m.GetLattice(0)->GetCell(i,j)->SetSolid();
			m.GetLattice(0)->GetCell(i,j)->Initialize(0.0001,v0,m.GetLattice(0)->Cs());
			m.GetLattice(1)->GetCell(i,j)->SetSolid();
			m.GetLattice(1)->GetCell(i,j)->Initialize(0.01,v0,m.GetLattice(1)->Cs());

		}
		else if ((i<m.Nx()/3)&&(i>=1)&&(j<3*m.Ny()/4)&&(j>=1)) {
			m.GetLattice(0)->GetCell(i,j)->Initialize(2.22,v0,m.GetLattice(0)->Cs());
			m.GetLattice(1)->GetCell(i,j)->Initialize(0.01,v0,m.GetLattice(1)->Cs());
		}
		else {
			m.GetLattice(0)->GetCell(i,j)->Initialize(0.0001,v0,m.GetLattice(0)->Cs());
			m.GetLattice(1)->GetCell(i,j)->Initialize(0.01,v0,m.GetLattice(1)->Cs());
		}

	}

	m.GetLattice(0)->SetG(-3.5)->SetGSolid(-0.0)->Homogenize();
	m.GetLattice(1)->SetG(-0.0)->SetGSolid(-0.0)->Homogenize();
	m.SetMixG(0.0);
	
	
	

	// Solve
	m.WriteState(0);
	//l.Solve(/*tIni*/0.0, /*tFin*/200.0, /*dt*/1.0, /*dtOut*/10.0);

	m.SetGravity(0.0,-0.001,0.0);
	m.Solve(/*tIni*/0.0, /*tFin*/1500.0, /*dtOut*/10.);

}
MECHSYS_CATCH
