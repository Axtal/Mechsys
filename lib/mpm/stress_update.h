/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso                                *
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

namespace MPM {

inline void MPoints2D::_stress_update(double Dt, double & DsE) // 3)
{
	// Stress/strain update
	DsE = 0.0;
	for (size_t p=0; p<_p.Size(); ++p)
	{
		if (_mdl[p]!=NULL)
		{
			// Variables
			STensor2 sp0; sp0 = _mdl[p]->Sig(); // initial point stress
			ATensor2 L;   L   = 0.0;            // mat point velocity gradient

			// Compute the velocity gradient
			int k = 0;
			for (int i=0; i<4; ++i)
			for (int j=0; j<4; ++j)
			{
				int n = _lbn[p]+((i*_g->nCol())+j); // grid node number
				DyadUp (_g->v(n), _sg[p].G[k],  L); // L += v dyad G
				k++;
			}
			
			// Compute the strain increment == rate of deformation D*Dt
			STensor2 de; Sym (Dt, L,  de); // de = Dt*sym(L)

			// Update the stress state
			_mdl[p]->Update (de);

			// Strain energy
			double Vp = 4.0*_l[p](0)*_l[p](1); // particle volume
			DsE += 0.5*blitz::dot((sp0+_mdl[p]->Sig()),de)*Vp; // Strain Energy on points

			// Update deformation gradient
			ATensor2 DF; DF = 1.0,1.0,1.0, 0.0,0.0,0.0, 0.0,0.0,0.0; // Delta F
			DF += L*Dt;             // DF = I + L*Dt
			Dot (DF, _F[p], _F[p]); // F = DF*F
		}
	}
}

}; // namespace MPM
