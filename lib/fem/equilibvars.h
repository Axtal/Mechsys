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

#ifndef MECHSYS_FEM_EQUILIBVARS_H
#define MECHSYS_FEM_EQUILIBVARS_H

namespace FEM
{

// Biot
const size_t ND_BIOT_3D        = 4;
const char   UD_BIOT_3D[ 4][4] = {"ux", "uy", "uz", "pwp"};
const char   FD_BIOT_3D[ 4][4] = {"fx", "fy", "fz", "vol"};

const size_t ND_BIOT_2D        = 3;
const char   UD_BIOT_2D[ 3][4] = {"ux", "uy", "pwp"};
const char   FD_BIOT_2D[ 3][4] = {"fx", "fy", "vol"};

const size_t NL_BIOT_3D        = 26;
const char   LB_BIOT_3D[26][4] = {"Ex", "Ey", "Ez", "Exy", "Eyz", "Ezx", "Sx", "Sy", "Sz", "Sxy", "Syz", "Szx", "E1", "E2", "E3", "S1", "S2", "S3", "p", "q", "Ev", "Ed", "Vx", "Vy", "Vz", "H"};

const size_t NL_BIOT_2D        = 19;
const char   LB_BIOT_2D[19][4] = {"Ex", "Ey", "Ez", "Exy", "Sx", "Sy", "Sz", "Sxy", "E1", "E2", "S1", "S2", "p", "q", "Ev", "Ed", "Vx", "Vy", "H"};

// Equilib
const size_t ND_EQUILIB_3D        = 3;
const char   UD_EQUILIB_3D[ 3][4] = {"ux", "uy", "uz"};
const char   FD_EQUILIB_3D[ 3][4] = {"fx", "fy", "fz"};

const size_t ND_EQUILIB_2D        = 2;
const char   UD_EQUILIB_2D[ 2][4] = {"ux", "uy"};
const char   FD_EQUILIB_2D[ 2][4] = {"fx", "fy"};

const size_t NL_EQUILIB_3D        = 22;
const char   LB_EQUILIB_3D[22][4] = {"Ex", "Ey", "Ez", "Exy", "Eyz", "Ezx", "Sx", "Sy", "Sz", "Sxy", "Syz", "Szx", "E1", "E2", "E3", "S1", "S2", "S3", "p", "q", "Ev", "Ed"};

const size_t NL_PSTRAIN        = 16;
const char   LB_PSTRAIN[16][4] = {"Ex", "Ey", "Ez", "Exy", "Sx", "Sy", "Sz", "Sxy", "E1", "E2", "S1", "S2", "p", "q", "Ev", "Ed"};

const size_t NL_PSTRESS        = 10;
const char   LB_PSTRESS[10][4] = {"Ex", "Ey", "Exy", "Sx", "Sy", "Sxy", "E1", "E2", "S1", "S2" };

// Beam
const size_t ND_BEAM_3D       = 6;
const char   UD_BEAM_3D[6][4] = {"ux", "uy", "uz", "wx", "wy", "wz"};
const char   FD_BEAM_3D[6][4] = {"fx", "fy", "fz", "mx", "my", "mz"};

const size_t ND_BEAM_2D       = 3;
const char   UD_BEAM_2D[3][4] = {"ux", "uy", "wz"};
const char   FD_BEAM_2D[3][4] = {"fx", "fy", "mz"};

const size_t NL_BEAM_3D       = 5;
const char   LB_BEAM_3D[5][4] = {"Ea", "Sa", "N", "V", "M"};

const size_t NL_BEAM_2D       = 5;
const char   LB_BEAM_2D[5][4] = {"Ea", "Sa", "N", "V", "M"};

// Rod
const size_t NL_ROD       = 3;
const char   LB_ROD[3][4] = {"Ea", "Sa", "N"};

}; // namespace FEM

#endif // MECHSYS_FEM_EQUILIBVARS_H
