/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2013 William Oquendo                                   *
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


#include <vector>

#include <mechsys/dem/domain.h>
using std::cout;
using std::endl;
using DEM::Domain;

int main( void )
{
    Domain dom;
    dom.AddVoroPack (-1, 0.05, 10.0, 10.0, 10.0, 10, 10, 20, 1.0, true, 1200, 1.0,1.0);
    dom.Initialize();
    dom.Save("domainwrite");
    dom.WriteXDMF("domainwrite");
    
    return 0;
}
