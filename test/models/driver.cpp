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

// STL
#include <iostream>
#include <sstream>
#include <fstream>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/linalg/umfpack.h>
#include <mechsys/models/linelastic.h>
#include <mechsys/models/elastoplastic.h>
#include <mechsys/models/camclay.h>
#include <mechsys/models/unconv01.h>
#include <mechsys/models/unconv02.h>
#include <mechsys/models/unconv03.h>
#include <mechsys/models/unconv04.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/numstreams.h>

using std::cout;
using std::endl;
using Util::_6_3;
using Util::_8s;
using Util::_10_6;
using Util::SQ2;
using Util::SQ3;
using Util::SQ6;
using Util::PI;

void zTgIncs (Model const * Mdl, EquilibState const * Sta, double LodeDeg, double dp, double dez, Vec_t & deps, Vec_t & dsig, Vec_t & divs, double dexy=0., double deyz=0., double dezx=0.)
{
    Mat_t D;
    Mdl->Stiffness (Sta, D);

    double c = -SQ3*dp;
    double dsx, dsy, dsz, dsxy, dsyz, dszx, dex, dey; 
    if (fabs(fabs(LodeDeg)-90.0)<1.0e-14)
    {
        //printf("--------------------------------------\n");
        dsx  = -((D(0,0)*(D(1,5)*D(2,1)-D(1,1)*D(2,5))+D(0,1)*(D(1,0)*D(2,5)-D(1,5)*D(2,0))+D(0,5)*(D(1,1)*D(2,0)-D(1,0)*D(2,1)))*dezx+(D(0,0)*(D(1,2)*D(2,1)-D(1,1)*D(2,2))+D(0,1)*(D(1,0)*D(2,2)-D(1,2)*D(2,0))+D(0,2)*(D(1,1)*D(2,0)-D(1,0)*D(2,1)))*dez+ (D(0,0)*(D(1,4)*D(2,1)-D(1,1)*D(2,4))+D(0,1)*(D(1,0)*D(2,4)-D(1,4)*D(2,0))+D(0,4)*(D(1,1)*D(2,0)-D(1,0)*D(2,1)))*deyz+(D(0,0)*(D(1,3)*D(2,1)-D(1,1)*D(2,3))+D(0,1)*(D(1,0)*D(2,3)-D(1,3)*D(2,0))+D(0,3)*(D(1,1)*D(2,0)-D(1,0)*D(2,1)))*dexy+D(0,0)*D(1,1)*c-D(0,1)*D(1,0)*c)/(D(1,0)*D(2,1)+D(0,0)*(-D(2,1)-2*D(1,1))+D(0,1)*(D(2,0)+2*D(1,0))-D(1,1)*D(2,0));
        dsy  = -((D(0,0)*(D(1,5)*D(2,1)-D(1,1)*D(2,5))+D(0,1)*(D(1,0)*D(2,5)-D(1,5)*D(2,0))+D(0,5)*(D(1,1)*D(2,0)-D(1,0)*D(2,1)))*dezx+(D(0,0)*(D(1,2)*D(2,1)-D(1,1)*D(2,2))+D(0,1)*(D(1,0)*D(2,2)-D(1,2)*D(2,0))+D(0,2)*(D(1,1)*D(2,0)-D(1,0)*D(2,1)))*dez+ (D(0,0)*(D(1,4)*D(2,1)-D(1,1)*D(2,4))+D(0,1)*(D(1,0)*D(2,4)-D(1,4)*D(2,0))+D(0,4)*(D(1,1)*D(2,0)-D(1,0)*D(2,1)))*deyz+(D(0,0)*(D(1,3)*D(2,1)-D(1,1)*D(2,3))+D(0,1)*(D(1,0)*D(2,3)-D(1,3)*D(2,0))+D(0,3)*(D(1,1)*D(2,0)-D(1,0)*D(2,1)))*dexy+D(0,0)*D(1,1)*c-D(0,1)*D(1,0)*c)/(D(1,0)*D(2,1)+D(0,0)*(-D(2,1)-2*D(1,1))+D(0,1)*(D(2,0)+2*D(1,0))-D(1,1)*D(2,0));
        dsz  =  ((D(0,0)*(2*D(1,5)*D(2,1)-2*D(1,1)*D(2,5))+D(0,1)*(2*D(1,0)*D(2,5)-2*D(1,5)*D(2,0))+D(0,5)*(2*D(1,1)*D(2,0)-2*D(1,0)*D(2,1)))*dezx+ (D(0,0)*(2*D(1,2)*D(2,1)-2*D(1,1)*D(2,2))+D(0,1)*(2*D(1,0)*D(2,2)-2*D(1,2)*D(2,0))+D(0,2)*(2*D(1,1)*D(2,0)-2*D(1,0)*D(2,1)))*dez+(D(0,0)*(2*D(1,4)*D(2,1)-2*D(1,1)*D(2,4))+D(0,1)*(2*D(1,0)*D(2,4)-2*D(1,4)*D(2,0))+D(0,4)*(2*D(1,1)*D(2,0)-2*D(1,0)*D(2,1)))*deyz+ (D(0,0)*(2*D(1,3)*D(2,1)-2*D(1,1)*D(2,3))+D(0,1)*(2*D(1,0)*D(2,3)-2*D(1,3)*D(2,0))+D(0,3)*(2*D(1,1)*D(2,0)-2*D(1,0)*D(2,1)))*dexy+D(1,0)*D(2,1)*c-D(0,0)*D(2,1)*c-D(1,1)*D(2,0)*c+D(0,1)*D(2,0)*c)/(D(1,0)*D(2,1)+D(0,0)*(-D(2,1)-2*D(1,1))+D(0,1)*(D(2,0)+2*D(1,0))- D(1,1)*D(2,0));
        dsxy =  ((D(1,0)*(D(2,1)*D(3,5)-D(2,5)*D(3,1))+D(0,0)*(-D(2,1)*D(3,5)-2*D(1,1)*D(3,5)+D(2,5)*D(3,1)+2*D(1,5)*D(3,1))+D(0,1)*(D(2,0)*D(3,5)+2*D(1,0)*D(3,5)-D(2,5)*D(3,0)-2*D(1,5)*D(3,0))+D(1,1)*(D(2,5)*D(3,0)-D(2,0)*D(3,5))+D(1,5)*(D(2,0)*D(3,1)-D(2,1)*D(3,0)) +D(0,5)*(-D(2,0)*D(3,1)-2*D(1,0)*D(3,1)+D(2,1)*D(3,0)+2*D(1,1)*D(3,0)))*dezx+(D(1,0)*(D(2,1)*D(3,2)-D(2,2)*D(3,1))+D(0,0)*(-D(2,1)*D(3,2)-2*D(1,1)*D(3,2)+D(2,2)*D(3,1)+2*D(1,2)*D(3,1))+D(0,1)*(D(2,0)*D(3,2)+2*D(1,0)*D(3,2)-D(2,2)*D(3,0)-2*D(1,2)*D(3,0))+D(1,1)* (D(2,2)*D(3,0)-D(2,0)*D(3,2))+D(1,2)*(D(2,0)*D(3,1)-D(2,1)*D(3,0))+D(0,2)*(-D(2,0)*D(3,1)-2*D(1,0)*D(3,1)+D(2,1)*D(3,0)+2*D(1,1)*D(3,0)))*dez+(D(1,0)*(D(2,1)*D(3,4)-D(2,4)*D(3,1))+D(0,0)*(-D(2,1)*D(3,4)-2*D(1,1)*D(3,4)+D(2,4)*D(3,1)+2*D(1,4)*D(3,1))+D(0,1)* (D(2,0)*D(3,4)+2*D(1,0)*D(3,4)-D(2,4)*D(3,0)-2*D(1,4)*D(3,0))+D(1,1)*(D(2,4)*D(3,0)-D(2,0)*D(3,4))+D(1,4)*(D(2,0)*D(3,1)-D(2,1)*D(3,0))+D(0,4)*(-D(2,0)*D(3,1)-2*D(1,0)*D(3,1)+D(2,1)*D(3,0)+2*D(1,1)*D(3,0)))*deyz+(D(1,0)*(D(2,1)*D(3,3)-D(2,3)*D(3,1))+D(0,0)* (-D(2,1)*D(3,3)-2*D(1,1)*D(3,3)+D(2,3)*D(3,1)+2*D(1,3)*D(3,1))+D(0,1)*(D(2,0)*D(3,3)+2*D(1,0)*D(3,3)-D(2,3)*D(3,0)-2*D(1,3)*D(3,0))+D(1,1)*(D(2,3)*D(3,0)-D(2,0)*D(3,3))+D(1,3)*(D(2,0)*D(3,1)-D(2,1)*D(3,0))+D(0,3)* (-D(2,0)*D(3,1)-2*D(1,0)*D(3,1)+D(2,1)*D(3,0)+2*D(1,1)*D(3,0)))*dexy+D(1,0)*D(3,1)*c-D(0,0)*D(3,1)*c-D(1,1)*D(3,0)*c+D(0,1)*D(3,0)*c)/(D(1,0)*D(2,1)+D(0,0)*(-D(2,1)-2*D(1,1))+D(0,1)*(D(2,0)+2*D(1,0))-D(1,1)*D(2,0));
        dsyz =  ((D(1,0)*(D(2,1)*D(4,5)-D(2,5)*D(4,1))+D(0,0)*(-D(2,1)*D(4,5)-2*D(1,1)*D(4,5)+D(2,5)*D(4,1)+2*D(1,5)*D(4,1))+D(0,1)*(D(2,0)*D(4,5)+2*D(1,0)*D(4,5)-D(2,5)*D(4,0)-2*D(1,5)*D(4,0))+D(1,1)*(D(2,5)*D(4,0)-D(2,0)*D(4,5))+D(1,5)*(D(2,0)*D(4,1)-D(2,1)*D(4,0)) +D(0,5)*(-D(2,0)*D(4,1)-2*D(1,0)*D(4,1)+D(2,1)*D(4,0)+2*D(1,1)*D(4,0)))*dezx+(D(1,0)*(D(2,1)*D(4,2)-D(2,2)*D(4,1))+D(0,0)*(-D(2,1)*D(4,2)-2*D(1,1)*D(4,2)+D(2,2)*D(4,1)+2*D(1,2)*D(4,1))+D(0,1)*(D(2,0)*D(4,2)+2*D(1,0)*D(4,2)-D(2,2)*D(4,0)-2*D(1,2)*D(4,0))+D(1,1)* (D(2,2)*D(4,0)-D(2,0)*D(4,2))+D(1,2)*(D(2,0)*D(4,1)-D(2,1)*D(4,0))+D(0,2)*(-D(2,0)*D(4,1)-2*D(1,0)*D(4,1)+D(2,1)*D(4,0)+2*D(1,1)*D(4,0)))*dez+(D(1,0)*(D(2,1)*D(4,4)-D(2,4)*D(4,1))+D(0,0)*(-D(2,1)*D(4,4)-2*D(1,1)*D(4,4)+D(2,4)*D(4,1)+2*D(1,4)*D(4,1))+D(0,1)* (D(2,0)*D(4,4)+2*D(1,0)*D(4,4)-D(2,4)*D(4,0)-2*D(1,4)*D(4,0))+D(1,1)*(D(2,4)*D(4,0)-D(2,0)*D(4,4))+D(1,4)*(D(2,0)*D(4,1)-D(2,1)*D(4,0))+D(0,4)*(-D(2,0)*D(4,1)-2*D(1,0)*D(4,1)+D(2,1)*D(4,0)+2*D(1,1)*D(4,0)))*deyz+(D(1,0)*(D(2,1)*D(4,3)-D(2,3)*D(4,1))+D(0,0)* (-D(2,1)*D(4,3)-2*D(1,1)*D(4,3)+D(2,3)*D(4,1)+2*D(1,3)*D(4,1))+D(0,1)*(D(2,0)*D(4,3)+2*D(1,0)*D(4,3)-D(2,3)*D(4,0)-2*D(1,3)*D(4,0))+D(1,1)*(D(2,3)*D(4,0)-D(2,0)*D(4,3))+D(1,3)*(D(2,0)*D(4,1)-D(2,1)*D(4,0))+D(0,3)* (-D(2,0)*D(4,1)-2*D(1,0)*D(4,1)+D(2,1)*D(4,0)+2*D(1,1)*D(4,0)))*dexy+D(1,0)*D(4,1)*c-D(0,0)*D(4,1)*c-D(1,1)*D(4,0)*c+D(0,1)*D(4,0)*c)/(D(1,0)*D(2,1)+D(0,0)*(-D(2,1)-2*D(1,1))+D(0,1)*(D(2,0)+2*D(1,0))-D(1,1)*D(2,0));
        dszx =  ((D(1,0)*(D(2,1)*D(5,5)-D(2,5)*D(5,1))+D(0,0)*(-D(2,1)*D(5,5)-2*D(1,1)*D(5,5)+D(2,5)*D(5,1)+2*D(1,5)*D(5,1))+D(0,1)*(D(2,0)*D(5,5)+2*D(1,0)*D(5,5)-D(2,5)*D(5,0)-2*D(1,5)*D(5,0))+D(1,1)*(D(2,5)*D(5,0)-D(2,0)*D(5,5))+D(1,5)*(D(2,0)*D(5,1)-D(2,1)*D(5,0)) +D(0,5)*(-D(2,0)*D(5,1)-2*D(1,0)*D(5,1)+D(2,1)*D(5,0)+2*D(1,1)*D(5,0)))*dezx+(D(1,0)*(D(2,1)*D(5,2)-D(2,2)*D(5,1))+D(0,0)*(-D(2,1)*D(5,2)-2*D(1,1)*D(5,2)+D(2,2)*D(5,1)+2*D(1,2)*D(5,1))+D(0,1)*(D(2,0)*D(5,2)+2*D(1,0)*D(5,2)-D(2,2)*D(5,0)-2*D(1,2)*D(5,0))+D(1,1)* (D(2,2)*D(5,0)-D(2,0)*D(5,2))+D(1,2)*(D(2,0)*D(5,1)-D(2,1)*D(5,0))+D(0,2)*(-D(2,0)*D(5,1)-2*D(1,0)*D(5,1)+D(2,1)*D(5,0)+2*D(1,1)*D(5,0)))*dez+(D(1,0)*(D(2,1)*D(5,4)-D(2,4)*D(5,1))+D(0,0)*(-D(2,1)*D(5,4)-2*D(1,1)*D(5,4)+D(2,4)*D(5,1)+2*D(1,4)*D(5,1))+D(0,1)* (D(2,0)*D(5,4)+2*D(1,0)*D(5,4)-D(2,4)*D(5,0)-2*D(1,4)*D(5,0))+D(1,1)*(D(2,4)*D(5,0)-D(2,0)*D(5,4))+D(1,4)*(D(2,0)*D(5,1)-D(2,1)*D(5,0))+D(0,4)*(-D(2,0)*D(5,1)-2*D(1,0)*D(5,1)+D(2,1)*D(5,0)+2*D(1,1)*D(5,0)))*deyz+(D(1,0)*(D(2,1)*D(5,3)-D(2,3)*D(5,1))+D(0,0)* (-D(2,1)*D(5,3)-2*D(1,1)*D(5,3)+D(2,3)*D(5,1)+2*D(1,3)*D(5,1))+D(0,1)*(D(2,0)*D(5,3)+2*D(1,0)*D(5,3)-D(2,3)*D(5,0)-2*D(1,3)*D(5,0))+D(1,1)*(D(2,3)*D(5,0)-D(2,0)*D(5,3))+D(1,3)*(D(2,0)*D(5,1)-D(2,1)*D(5,0))+D(0,3)* (-D(2,0)*D(5,1)-2*D(1,0)*D(5,1)+D(2,1)*D(5,0)+2*D(1,1)*D(5,0)))*dexy+D(1,0)*D(5,1)*c-D(0,0)*D(5,1)*c-D(1,1)*D(5,0)*c+D(0,1)*D(5,0)*c)/(D(1,0)*D(2,1)+D(0,0)*(-D(2,1)-2*D(1,1))+D(0,1)*(D(2,0)+2*D(1,0))-D(1,1)*D(2,0));
        dex  =  ((D(1,1)*D(2,5)+D(0,1)*(-D(2,5)-2*D(1,5))+D(0,5)*(D(2,1)+2*D(1,1))-D(1,5)*D(2,1))*dezx+(D(1,1)*D(2,2)+D(0,1)*(-D(2,2)-2*D(1,2))+D(0,2)*(D(2,1)+2*D(1,1))-D(1,2)*D(2,1))*dez+ (D(1,1)*D(2,4)+D(0,1)*(-D(2,4)-2*D(1,4))+D(0,4)*(D(2,1)+2*D(1,1))-D(1,4)*D(2,1))*deyz+(D(1,1)*D(2,3)+D(0,1)*(-D(2,3)-2*D(1,3))+D(0,3)*(D(2,1)+2*D(1,1))-D(1,3)*D(2,1))*dexy-D(1,1)*c+D(0,1)*c)/(D(1,0)*D(2,1)+D(0,0)*(-D(2,1)-2*D(1,1))+D(0,1)* (D(2,0)+2*D(1,0))-D(1,1)*D(2,0));
        dey  = -((D(1,0)*D(2,5)+D(0,0)*(-D(2,5)-2*D(1,5))+D(0,5)*(D(2,0)+2*D(1,0))-D(1,5)*D(2,0))*dezx+(D(1,0)*D(2,2)+D(0,0)*(-D(2,2)-2*D(1,2))+D(0,2)*(D(2,0)+2*D(1,0))-D(1,2)*D(2,0))*dez+ (D(1,0)*D(2,4)+D(0,0)*(-D(2,4)-2*D(1,4))+D(0,4)*(D(2,0)+2*D(1,0))-D(1,4)*D(2,0))*deyz+(D(1,0)*D(2,3)+D(0,0)*(-D(2,3)-2*D(1,3))+D(0,3)*(D(2,0)+2*D(1,0))-D(1,3)*D(2,0))*dexy-D(1,0)*c+D(0,0)*c)/(D(1,0)*D(2,1)+D(0,0)*(-D(2,1)-2*D(1,1))+D(0,1)* (D(2,0)+2*D(1,0))-D(1,1)*D(2,0));
    }
    else
    {
        double lode = LodeDeg*PI/180.;
        double m    = tan(lode);
        double a    = (1.0-SQ3*m)/2.;
        double b    = (1.0+SQ3*m)/2.;
        dsx  =  ((D(0,1)*(D(1,0)*D(2,5)*(a+1)+D(1,5)*D(2,0)*(-a-1))+D(0,0)*(D(1,5)*D(2,1)*(a+1)+D(1,1)*D(2,5)*(-a-1))+D(0,5)*(D(1,1)*D(2,0)*(a+1)+D(1,0)*D(2,1)*(-a-1)))*dezx+ (D(0,1)*(D(1,0)*D(2,2)*(a+1)+D(1,2)*D(2,0)*(-a-1))+D(0,0)*(D(1,2)*D(2,1)*(a+1)+D(1,1)*D(2,2)*(-a-1))+D(0,2)*(D(1,1)*D(2,0)*(a+1)+D(1,0)*D(2,1)*(-a-1)))*dez+ (D(0,1)*(D(1,0)*D(2,4)*(a+1)+D(1,4)*D(2,0)*(-a-1))+D(0,0)*(D(1,4)*D(2,1)*(a+1)+D(1,1)*D(2,4)*(-a-1))+D(0,4)*(D(1,1)*D(2,0)*(a+1)+D(1,0)*D(2,1)*(-a-1)))*deyz+ (D(0,1)*(D(1,0)*D(2,3)*(a+1)+D(1,3)*D(2,0)*(-a-1))+D(0,0)*(D(1,3)*D(2,1)*(a+1)+D(1,1)*D(2,3)*(-a-1))+D(0,3)*(D(1,1)*D(2,0)*(a+1)+D(1,0)*D(2,1)*(-a-1)))*dexy+D(0,0)*(D(1,1)*a*c-D(2,1)*c)+D(0,1)*(D(2,0)*c-D(1,0)*a*c))/(D(0,1)* (D(1,0)*(b-a)+D(2,0)*(b+1))+D(0,0)*(D(1,1)*(a-b)+D(2,1)*(-b-1))+D(1,1)*D(2,0)*(a+1)+D(1,0)*D(2,1)*(-a-1));
        dsy  = -((D(0,1)*(D(1,0)*D(2,5)*(b+1)+D(1,5)*D(2,0)*(-b-1))+D(0,0)*(D(1,5)*D(2,1)*(b+1)+D(1,1)*D(2,5)*(-b-1))+D(0,5)*(D(1,1)*D(2,0)*(b+1)+D(1,0)*D(2,1)*(-b-1)))*dezx+ (D(0,1)*(D(1,0)*D(2,2)*(b+1)+D(1,2)*D(2,0)*(-b-1))+D(0,0)*(D(1,2)*D(2,1)*(b+1)+D(1,1)*D(2,2)*(-b-1))+D(0,2)*(D(1,1)*D(2,0)*(b+1)+D(1,0)*D(2,1)*(-b-1)))*dez+ (D(0,1)*(D(1,0)*D(2,4)*(b+1)+D(1,4)*D(2,0)*(-b-1))+D(0,0)*(D(1,4)*D(2,1)*(b+1)+D(1,1)*D(2,4)*(-b-1))+D(0,4)*(D(1,1)*D(2,0)*(b+1)+D(1,0)*D(2,1)*(-b-1)))*deyz+ (D(0,1)*(D(1,0)*D(2,3)*(b+1)+D(1,3)*D(2,0)*(-b-1))+D(0,0)*(D(1,3)*D(2,1)*(b+1)+D(1,1)*D(2,3)*(-b-1))+D(0,3)*(D(1,1)*D(2,0)*(b+1)+D(1,0)*D(2,1)*(-b-1)))*dexy+D(0,0)*D(1,1)*b*c-D(0,1)*D(1,0)*b*c+D(1,0)*D(2,1)*c-D(1,1)*D(2,0)*c)/(D(0,1)* (D(1,0)*(b-a)+D(2,0)*(b+1))+D(0,0)*(D(1,1)*(a-b)+D(2,1)*(-b-1))+D(1,1)*D(2,0)*(a+1)+D(1,0)*D(2,1)*(-a-1));
        dsz  =  ((D(0,1)*(D(1,0)*D(2,5)*(b-a)+D(1,5)*D(2,0)*(a-b))+D(0,0)*(D(1,5)*D(2,1)*(b-a)+D(1,1)*D(2,5)*(a-b))+D(0,5)*(D(1,1)*D(2,0)*(b-a)+D(1,0)*D(2,1)*(a-b)))*dezx+ (D(0,1)*(D(1,0)*D(2,2)*(b-a)+D(1,2)*D(2,0)*(a-b))+D(0,0)*(D(1,2)*D(2,1)*(b-a)+D(1,1)*D(2,2)*(a-b))+D(0,2)*(D(1,1)*D(2,0)*(b-a)+D(1,0)*D(2,1)*(a-b)))*dez+ (D(0,1)*(D(1,0)*D(2,4)*(b-a)+D(1,4)*D(2,0)*(a-b))+D(0,0)*(D(1,4)*D(2,1)*(b-a)+D(1,1)*D(2,4)*(a-b))+D(0,4)*(D(1,1)*D(2,0)*(b-a)+D(1,0)*D(2,1)*(a-b)))*deyz+ (D(0,1)*(D(1,0)*D(2,3)*(b-a)+D(1,3)*D(2,0)*(a-b))+D(0,0)*(D(1,3)*D(2,1)*(b-a)+D(1,1)*D(2,3)*(a-b))+D(0,3)*(D(1,1)*D(2,0)*(b-a)+D(1,0)*D(2,1)*(a-b)))*dexy-D(0,0)*D(2,1)*b*c+D(0,1)*D(2,0)*b*c-D(1,0)*D(2,1)*a*c+D(1,1)*D(2,0)*a*c)/(D(0,1)* (D(1,0)*(b-a)+D(2,0)*(b+1))+D(0,0)*(D(1,1)*(a-b)+D(2,1)*(-b-1))+D(1,1)*D(2,0)*(a+1)+D(1,0)*D(2,1)*(-a-1));
        dsxy =  ((D(0,1)*(D(1,0)*D(3,5)*(b-a)+D(2,0)*D(3,5)*(b+1)+D(1,5)*D(3,0)*(a-b)+D(2,5)*D(3,0)*(-b-1))+D(0,0)*(D(1,5)*D(3,1)*(b-a)+D(2,5)*D(3,1)*(b+1)+D(1,1)*D(3,5)*(a-b)+D(2,1)*D(3,5)*(-b-1))+D(0,5)* (D(1,1)*D(3,0)*(b-a)+D(2,1)*D(3,0)*(b+1)+D(1,0)*D(3,1)*(a-b)+D(2,0)*D(3,1)*(-b-1))+D(1,1)*(D(2,0)*D(3,5)*(a+1)+D(2,5)*D(3,0)*(-a-1))+D(1,0)*(D(2,5)*D(3,1)*(a+1)+D(2,1)*D(3,5)*(-a-1))+D(1,5)*(D(2,1)*D(3,0)*(a+1)+D(2,0)*D(3,1)*(-a-1)))*dezx+( D(0,1)*(D(1,0)*D(3,2)*(b-a)+D(2,0)*D(3,2)*(b+1)+D(1,2)*D(3,0)*(a-b)+D(2,2)*D(3,0)*(-b-1))+D(0,0)*(D(1,2)*D(3,1)*(b-a)+D(2,2)*D(3,1)*(b+1)+D(1,1)*D(3,2)*(a-b)+D(2,1)*D(3,2)*(-b-1))+D(0,2)* (D(1,1)*D(3,0)*(b-a)+D(2,1)*D(3,0)*(b+1)+D(1,0)*D(3,1)*(a-b)+D(2,0)*D(3,1)*(-b-1))+D(1,1)*(D(2,0)*D(3,2)*(a+1)+D(2,2)*D(3,0)*(-a-1))+D(1,0)*(D(2,2)*D(3,1)*(a+1)+D(2,1)*D(3,2)*(-a-1))+D(1,2)*(D(2,1)*D(3,0)*(a+1)+D(2,0)*D(3,1)*(-a-1)))*dez+( D(0,1)*(D(1,0)*D(3,4)*(b-a)+D(2,0)*D(3,4)*(b+1)+D(1,4)*D(3,0)*(a-b)+D(2,4)*D(3,0)*(-b-1))+D(0,0)*(D(1,4)*D(3,1)*(b-a)+D(2,4)*D(3,1)*(b+1)+D(1,1)*D(3,4)*(a-b)+D(2,1)*D(3,4)*(-b-1))+D(0,4)* (D(1,1)*D(3,0)*(b-a)+D(2,1)*D(3,0)*(b+1)+D(1,0)*D(3,1)*(a-b)+D(2,0)*D(3,1)*(-b-1))+D(1,1)*(D(2,0)*D(3,4)*(a+1)+D(2,4)*D(3,0)*(-a-1))+D(1,0)*(D(2,4)*D(3,1)*(a+1)+D(2,1)*D(3,4)*(-a-1))+D(1,4)*(D(2,1)*D(3,0)*(a+1)+D(2,0)*D(3,1)*(-a-1)))*deyz+( D(0,1)*(D(1,0)*D(3,3)*(b-a)+D(2,0)*D(3,3)*(b+1)+D(1,3)*D(3,0)*(a-b)+D(2,3)*D(3,0)*(-b-1))+D(0,0)*(D(1,3)*D(3,1)*(b-a)+D(2,3)*D(3,1)*(b+1)+D(1,1)*D(3,3)*(a-b)+D(2,1)*D(3,3)*(-b-1))+D(0,3)* (D(1,1)*D(3,0)*(b-a)+D(2,1)*D(3,0)*(b+1)+D(1,0)*D(3,1)*(a-b)+D(2,0)*D(3,1)*(-b-1))+D(1,1)*(D(2,0)*D(3,3)*(a+1)+D(2,3)*D(3,0)*(-a-1))+D(1,0)*(D(2,3)*D(3,1)*(a+1)+D(2,1)*D(3,3)*(-a-1))+D(1,3)*(D(2,1)*D(3,0)*(a+1)+D(2,0)*D(3,1)*(-a-1)))*dexy-D(0,0)* D(3,1)*b*c+D(0,1)*D(3,0)*b*c-D(1,0)*D(3,1)*a*c+D(1,1)*D(3,0)*a*c+D(2,0)*D(3,1)*c-D(2,1)*D(3,0)*c)/(D(0,1)*(D(1,0)*(b-a)+D(2,0)*(b+1))+D(0,0)*(D(1,1)*(a-b)+D(2,1)*(-b-1))+D(1,1)*D(2,0)*(a+1)+D(1,0)*D(2,1)*(-a-1));
        dsyz =  ((D(0,1)*(D(1,0)*D(4,5)*(b-a)+D(2,0)*D(4,5)*(b+1)+D(1,5)*D(4,0)*(a-b)+D(2,5)*D(4,0)*(-b-1))+D(0,0)*(D(1,5)*D(4,1)*(b-a)+D(2,5)*D(4,1)*(b+1)+D(1,1)*D(4,5)*(a-b)+D(2,1)*D(4,5)*(-b-1))+D(0,5)* (D(1,1)*D(4,0)*(b-a)+D(2,1)*D(4,0)*(b+1)+D(1,0)*D(4,1)*(a-b)+D(2,0)*D(4,1)*(-b-1))+D(1,1)*(D(2,0)*D(4,5)*(a+1)+D(2,5)*D(4,0)*(-a-1))+D(1,0)*(D(2,5)*D(4,1)*(a+1)+D(2,1)*D(4,5)*(-a-1))+D(1,5)*(D(2,1)*D(4,0)*(a+1)+D(2,0)*D(4,1)*(-a-1)))*dezx+( D(0,1)*(D(1,0)*D(4,2)*(b-a)+D(2,0)*D(4,2)*(b+1)+D(1,2)*D(4,0)*(a-b)+D(2,2)*D(4,0)*(-b-1))+D(0,0)*(D(1,2)*D(4,1)*(b-a)+D(2,2)*D(4,1)*(b+1)+D(1,1)*D(4,2)*(a-b)+D(2,1)*D(4,2)*(-b-1))+D(0,2)* (D(1,1)*D(4,0)*(b-a)+D(2,1)*D(4,0)*(b+1)+D(1,0)*D(4,1)*(a-b)+D(2,0)*D(4,1)*(-b-1))+D(1,1)*(D(2,0)*D(4,2)*(a+1)+D(2,2)*D(4,0)*(-a-1))+D(1,0)*(D(2,2)*D(4,1)*(a+1)+D(2,1)*D(4,2)*(-a-1))+D(1,2)*(D(2,1)*D(4,0)*(a+1)+D(2,0)*D(4,1)*(-a-1)))*dez+( D(0,1)*(D(1,0)*D(4,4)*(b-a)+D(2,0)*D(4,4)*(b+1)+D(1,4)*D(4,0)*(a-b)+D(2,4)*D(4,0)*(-b-1))+D(0,0)*(D(1,4)*D(4,1)*(b-a)+D(2,4)*D(4,1)*(b+1)+D(1,1)*D(4,4)*(a-b)+D(2,1)*D(4,4)*(-b-1))+D(0,4)* (D(1,1)*D(4,0)*(b-a)+D(2,1)*D(4,0)*(b+1)+D(1,0)*D(4,1)*(a-b)+D(2,0)*D(4,1)*(-b-1))+D(1,1)*(D(2,0)*D(4,4)*(a+1)+D(2,4)*D(4,0)*(-a-1))+D(1,0)*(D(2,4)*D(4,1)*(a+1)+D(2,1)*D(4,4)*(-a-1))+D(1,4)*(D(2,1)*D(4,0)*(a+1)+D(2,0)*D(4,1)*(-a-1)))*deyz+( D(0,1)*(D(1,0)*D(4,3)*(b-a)+D(2,0)*D(4,3)*(b+1)+D(1,3)*D(4,0)*(a-b)+D(2,3)*D(4,0)*(-b-1))+D(0,0)*(D(1,3)*D(4,1)*(b-a)+D(2,3)*D(4,1)*(b+1)+D(1,1)*D(4,3)*(a-b)+D(2,1)*D(4,3)*(-b-1))+D(0,3)* (D(1,1)*D(4,0)*(b-a)+D(2,1)*D(4,0)*(b+1)+D(1,0)*D(4,1)*(a-b)+D(2,0)*D(4,1)*(-b-1))+D(1,1)*(D(2,0)*D(4,3)*(a+1)+D(2,3)*D(4,0)*(-a-1))+D(1,0)*(D(2,3)*D(4,1)*(a+1)+D(2,1)*D(4,3)*(-a-1))+D(1,3)*(D(2,1)*D(4,0)*(a+1)+D(2,0)*D(4,1)*(-a-1)))*dexy-D(0,0)* D(4,1)*b*c+D(0,1)*D(4,0)*b*c-D(1,0)*D(4,1)*a*c+D(1,1)*D(4,0)*a*c+D(2,0)*D(4,1)*c-D(2,1)*D(4,0)*c)/(D(0,1)*(D(1,0)*(b-a)+D(2,0)*(b+1))+D(0,0)*(D(1,1)*(a-b)+D(2,1)*(-b-1))+D(1,1)*D(2,0)*(a+1)+D(1,0)*D(2,1)*(-a-1));
        dszx =  ((D(0,1)*(D(1,0)*D(5,5)*(b-a)+D(2,0)*D(5,5)*(b+1)+D(1,5)*D(5,0)*(a-b)+D(2,5)*D(5,0)*(-b-1))+D(0,0)*(D(1,5)*D(5,1)*(b-a)+D(2,5)*D(5,1)*(b+1)+D(1,1)*D(5,5)*(a-b)+D(2,1)*D(5,5)*(-b-1))+D(0,5)* (D(1,1)*D(5,0)*(b-a)+D(2,1)*D(5,0)*(b+1)+D(1,0)*D(5,1)*(a-b)+D(2,0)*D(5,1)*(-b-1))+D(1,1)*(D(2,0)*D(5,5)*(a+1)+D(2,5)*D(5,0)*(-a-1))+D(1,0)*(D(2,5)*D(5,1)*(a+1)+D(2,1)*D(5,5)*(-a-1))+D(1,5)*(D(2,1)*D(5,0)*(a+1)+D(2,0)*D(5,1)*(-a-1)))*dezx+( D(0,1)*(D(1,0)*D(5,2)*(b-a)+D(2,0)*D(5,2)*(b+1)+D(1,2)*D(5,0)*(a-b)+D(2,2)*D(5,0)*(-b-1))+D(0,0)*(D(1,2)*D(5,1)*(b-a)+D(2,2)*D(5,1)*(b+1)+D(1,1)*D(5,2)*(a-b)+D(2,1)*D(5,2)*(-b-1))+D(0,2)* (D(1,1)*D(5,0)*(b-a)+D(2,1)*D(5,0)*(b+1)+D(1,0)*D(5,1)*(a-b)+D(2,0)*D(5,1)*(-b-1))+D(1,1)*(D(2,0)*D(5,2)*(a+1)+D(2,2)*D(5,0)*(-a-1))+D(1,0)*(D(2,2)*D(5,1)*(a+1)+D(2,1)*D(5,2)*(-a-1))+D(1,2)*(D(2,1)*D(5,0)*(a+1)+D(2,0)*D(5,1)*(-a-1)))*dez+( D(0,1)*(D(1,0)*D(5,4)*(b-a)+D(2,0)*D(5,4)*(b+1)+D(1,4)*D(5,0)*(a-b)+D(2,4)*D(5,0)*(-b-1))+D(0,0)*(D(1,4)*D(5,1)*(b-a)+D(2,4)*D(5,1)*(b+1)+D(1,1)*D(5,4)*(a-b)+D(2,1)*D(5,4)*(-b-1))+D(0,4)* (D(1,1)*D(5,0)*(b-a)+D(2,1)*D(5,0)*(b+1)+D(1,0)*D(5,1)*(a-b)+D(2,0)*D(5,1)*(-b-1))+D(1,1)*(D(2,0)*D(5,4)*(a+1)+D(2,4)*D(5,0)*(-a-1))+D(1,0)*(D(2,4)*D(5,1)*(a+1)+D(2,1)*D(5,4)*(-a-1))+D(1,4)*(D(2,1)*D(5,0)*(a+1)+D(2,0)*D(5,1)*(-a-1)))*deyz+( D(0,1)*(D(1,0)*D(5,3)*(b-a)+D(2,0)*D(5,3)*(b+1)+D(1,3)*D(5,0)*(a-b)+D(2,3)*D(5,0)*(-b-1))+D(0,0)*(D(1,3)*D(5,1)*(b-a)+D(2,3)*D(5,1)*(b+1)+D(1,1)*D(5,3)*(a-b)+D(2,1)*D(5,3)*(-b-1))+D(0,3)* (D(1,1)*D(5,0)*(b-a)+D(2,1)*D(5,0)*(b+1)+D(1,0)*D(5,1)*(a-b)+D(2,0)*D(5,1)*(-b-1))+D(1,1)*(D(2,0)*D(5,3)*(a+1)+D(2,3)*D(5,0)*(-a-1))+D(1,0)*(D(2,3)*D(5,1)*(a+1)+D(2,1)*D(5,3)*(-a-1))+D(1,3)*(D(2,1)*D(5,0)*(a+1)+D(2,0)*D(5,1)*(-a-1)))*dexy-D(0,0)* D(5,1)*b*c+D(0,1)*D(5,0)*b*c-D(1,0)*D(5,1)*a*c+D(1,1)*D(5,0)*a*c+D(2,0)*D(5,1)*c-D(2,1)*D(5,0)*c)/(D(0,1)*(D(1,0)*(b-a)+D(2,0)*(b+1))+D(0,0)*(D(1,1)*(a-b)+D(2,1)*(-b-1))+D(1,1)*D(2,0)*(a+1)+D(1,0)*D(2,1)*(-a-1));
        dex  =  ((D(0,5)*(D(1,1)*(b-a)+D(2,1)*(b+1))+D(0,1)*(D(1,5)*(a-b)+D(2,5)*(-b-1))+D(1,5)*D(2,1)*(a+1)+D(1,1)*D(2,5)*(-a-1))*dezx+ (D(0,2)*(D(1,1)*(b-a)+D(2,1)*(b+1))+D(0,1)*(D(1,2)*(a-b)+D(2,2)*(-b-1))+D(1,2)*D(2,1)*(a+1)+D(1,1)*D(2,2)*(-a-1))*dez+(D(0,4)*(D(1,1)*(b-a)+D(2,1)*(b+1))+D(0,1)*(D(1,4)*(a-b)+D(2,4)*(-b-1))+D(1,4)*D(2,1)*(a+1)+D(1,1)*D(2,4)*(-a-1))* deyz+(D(0,3)*(D(1,1)*(b-a)+D(2,1)*(b+1))+D(0,1)*(D(1,3)*(a-b)+D(2,3)*(-b-1))+D(1,3)*D(2,1)*(a+1)+D(1,1)*D(2,3)*(-a-1))*dexy+D(0,1)*b*c+D(1,1)*a*c-D(2,1)*c)/(D(0,1)*(D(1,0)*(b-a)+D(2,0)*(b+1))+D(0,0)*(D(1,1)*(a-b)+D(2,1)*(-b-1))+ D(1,1)*D(2,0)*(a+1)+D(1,0)*D(2,1)*(-a-1));
        dey  = -((D(0,5)*(D(1,0)*(b-a)+D(2,0)*(b+1))+D(0,0)*(D(1,5)*(a-b)+D(2,5)*(-b-1))+D(1,5)*D(2,0)*(a+1)+D(1,0)*D(2,5)*(-a-1))*dezx+ (D(0,2)*(D(1,0)*(b-a)+D(2,0)*(b+1))+D(0,0)*(D(1,2)*(a-b)+D(2,2)*(-b-1))+D(1,2)*D(2,0)*(a+1)+D(1,0)*D(2,2)*(-a-1))*dez+(D(0,4)*(D(1,0)*(b-a)+D(2,0)*(b+1))+D(0,0)*(D(1,4)*(a-b)+D(2,4)*(-b-1))+D(1,4)*D(2,0)*(a+1)+D(1,0)*D(2,4)*(-a-1))* deyz+(D(0,3)*(D(1,0)*(b-a)+D(2,0)*(b+1))+D(0,0)*(D(1,3)*(a-b)+D(2,3)*(-b-1))+D(1,3)*D(2,0)*(a+1)+D(1,0)*D(2,3)*(-a-1))*dexy+D(0,0)*b*c+D(1,0)*a*c-D(2,0)*c)/(D(0,1)*(D(1,0)*(b-a)+D(2,0)*(b+1))+D(0,0)*(D(1,1)*(a-b)+D(2,1)*(-b-1))+ D(1,1)*D(2,0)*(a+1)+D(1,0)*D(2,1)*(-a-1));
    }

    dsig.change_dim(6);
    deps.change_dim(6);
    dsig = dsx, dsy, dsz, dsxy, dsyz, dszx;
    deps = dex, dey, dez, dexy, deyz, dezx;

    // calc divs
    Vec_t deps_tmp(deps), dsig_tmp(dsig);
    Mdl->TgIncs (Sta, deps_tmp, dsig_tmp, divs);

    //printf("deps=[%g, %g, %g],  dsig=[%g, %g, %g]\n",deps(0),deps(1),deps(2),dsig(0),dsig(1),dsig(2));
}

//size_t ITMAX = 0;

void kTgIncs (Model const * Mdl, EquilibState const * Sta, double LodeDeg, double k, double dez, Vec_t & deps, Vec_t & dsig, Vec_t & divs)
{
    Mat_t D;
    Mdl->Stiffness (Sta, D);

    double dsx,dsy,dsz,dex,dey;
    double m = tan(LodeDeg*PI/180.);
    double a = (1.+m*SQ3)/2.;
    double b = (1.-m*SQ3)/2.;
    double n = SQ3*sqrt(1.+m*m)/(SQ2*k);
    double c = n-1.;
    double d = n+1.;
    //printf("m = %g\n",m);
    //if (fabs(fabs(LodeDeg)-90.0)<1.0e-14)
    //{
        //dsx = -((D(0,0)*(D(1,2)*D(2,1)-D(1,1)*D(2,2))+D(0,1)*(D(1,0)*D(2,2)-D(1,2)*D(2,0))+D(0,2)*(D(1,1)*D(2,0)-D(1,0)*D(2,1)))*dez)/(D(0,1)*(D(1,0)*(d-c)+D(2,0))+D(0,0)*(D(1,1)*(c-d)-D(2,1))+D(1,0)*D(2,1)-D(1,1)*D(2,0));
        //dsy = -((D(0,0)*(D(1,2)*D(2,1)-D(1,1)*D(2,2))+D(0,1)*(D(1,0)*D(2,2)-D(1,2)*D(2,0))+D(0,2)*(D(1,1)*D(2,0)-D(1,0)*D(2,1)))*dez)/(D(0,1)*(D(1,0)*(d-c)+D(2,0))+D(0,0)*(D(1,1)*(c-d)-D(2,1))+D(1,0)*D(2,1)-D(1,1)*D(2,0));
        //dsz =  ((D(0,1)*(D(1,0)*D(2,2)*(d-c)+D(1,2)*D(2,0)*(c-d))+D(0,0)*(D(1,2)*D(2,1)*(d-c)+D(1,1)*D(2,2)*(c-d))+D(0,2)* (D(1,1)*D(2,0)*(d-c)+D(1,0)*D(2,1)*(c-d)))*dez)/(D(0,1)*(D(1,0)*(d-c)+D(2,0))+D(0,0)*(D(1,1)*(c-d)-D(2,1))+D(1,0)* D(2,1)-D(1,1)*D(2,0));
        //dex =  ((D(0,2)*(D(1,1)*(d-c)+D(2,1))+D(0,1)*(D(1,2)*(c-d)-D(2,2))+D(1,1)*D(2,2)-D(1,2)*D(2,1))*dez)/(D(0,1)*(D(1,0)*(d-c)+D(2,0))+D(0,0)*(D(1,1)*(c-d)-D(2,1))+D(1,0)*D(2,1)-D(1,1)*D(2,0));
        //dey = -((D(0,2)*(D(1,0)*(d-c)+D(2,0))+D(0,0)*(D(1,2)*(c-d)-D(2,2))+D(1,0)*D(2,2)-D(1,2)*D(2,0))*dez)/(D(0,1)*(D(1,0)*(d-c)+D(2,0))+D(0,0)*(D(1,1)*(c-d)-D(2,1))+D(1,0)*D(2,1)-D(1,1)*D(2,0));
    //}
    //else
    //{
        dsx =  ((D(0,1)*(D(1,0)*D(2,2)*(d+b)+D(1,2)*D(2,0)*(-d-b))+D(0,0)* (D(1,2)*D(2,1)*(d+b)+D(1,1)*D(2,2)*(-d-b))+D(0,2)*(D(1,1)*D(2,0)*(d+b)+D(1,0)*D(2,1)*(-d-b)))*dez)/(D(0,1)* (D(1,0)*(a*d+b*c)+D(2,0)*(a-c))+D(0,0)*(D(1,1)*(-a*d-b*c)+D(2,1)*(c-a))+D(1,1)*D(2,0)*(d+b)+D(1,0)*D(2,1)*(-d-b));
        dsy = -((D(0,0)*(D(1,1)*D(2,2)*(c-a)+D(1,2)*D(2,1)*(a-c))+D(0,2)*(D(1,0)*D(2,1)*(c-a)+D(1,1)*D(2,0)*(a-c))+D(0,1)* (D(1,2)*D(2,0)*(c-a)+D(1,0)*D(2,2)*(a-c)))*dez)/(D(0,1)*(D(1,0)*(a*d+b*c)+D(2,0)*(a-c))+D(0,0)* (D(1,1)*(-a*d-b*c)+D(2,1)*(c-a))+D(1,1)*D(2,0)*(d+b)+D(1,0)*D(2,1)*(-d-b));
        dsz =  ((D(0,1)*(D(1,0)*D(2,2)*(a*d+b*c)+D(1,2)*D(2,0)*(-a*d-b*c))+D(0,0)*(D(1,2)*D(2,1)*(a*d+b*c)+D(1,1)*D(2,2)*(-a*d-b*c))+D(0,2)* (D(1,1)*D(2,0)*(a*d+b*c)+D(1,0)*D(2,1)*(-a*d-b*c)))*dez)/(D(0,1)*(D(1,0)*(a*d+b*c)+D(2,0)*(a-c))+D(0,0)* (D(1,1)*(-a*d-b*c)+D(2,1)*(c-a))+D(1,1)*D(2,0)*(d+b)+D(1,0)*D(2,1)*(-d-b));
        dex =  ((D(0,2)*(D(1,1)*(a*d+b*c)+D(2,1)*(a-c))+D(0,1)*(D(1,2)*(-a*d-b*c)+D(2,2)*(c-a))+D(1,2)*D(2,1)*(d+b)+D(1,1)*D(2,2)*(-d-b))*dez)/(D(0,1)*(D(1,0)*(a*d+b*c)+D(2,0)*(a-c))+D(0,0)*(D(1,1)*(-a*d-b*c)+D(2,1)*(c-a))+D(1,1)*D(2,0)*(d+b)+D(1,0)*D(2,1)*(-d-b));
        dey = -((D(0,2)*(D(1,0)*(a*d+b*c)+D(2,0)*(a-c))+D(0,0)*(D(1,2)*(-a*d-b*c)+D(2,2)*(c-a))+D(1,2)*D(2,0)*(d+b)+D(1,0)*D(2,2)*(-d-b))*dez)/(D(0,1)*(D(1,0)*(a*d+b*c)+D(2,0)*(a-c))+D(0,0)*(D(1,1)*(-a*d-b*c)+D(2,1)*(c-a))+D(1,1)*D(2,0)*(d+b)+D(1,0)*D(2,1)*(-d-b));
    //}

    // trial
    /*
    deps = 0., 0., dez, 0., 0., 0.;
    dsig = D*deps;
    Vec_t sigf(Sta->Sig + dsig);

    // constants
    double sx0  = Sta->Sig(0);
    double sy0  = Sta->Sig(1);
    double sz0  = Sta->Sig(2);
    double A    = sx0+D(0,2)*dez;
    double B    = sy0+D(1,2)*dez;
    double C    = sz0+D(2,2)*dez;
    double sa0  = (sy0-sx0)/SQ2;
    double sb0  = (sy0+sx0-2.*sz0)/SQ6;
    double sqH0 = sqrt(pow(sx0-sy0,2.) + pow(sy0-sz0,2.) + pow(sz0-sx0,2.));
           c    = sqH0 + k*(sx0+sy0+sz0);
    double M;
    a = 1.-m*SQ3;
    b = 1.+m*SQ3;
    M = SQ6*(sb0-m*sa0);

    // initial values
    Vec_t x(5), R(5), dx(5);
    x = sigf(0), sigf(1), sigf(2), deps(0), deps(1);
    Mat_t J(5,5), Ji(5,5);

    // iterations
    double tolR  = 1.0e-11;
    size_t maxit = 10;
    size_t it;
    for (it=1; it<=maxit; ++it)
    {
        double sqHf = sqrt(pow(x(0)-x(1),2.) + pow(x(1)-x(2),2.) + pow(x(2)-x(0),2.));
        R = x(0) - D(0,0)*x(3) - D(0,1)*x(4) - A,
            x(1) - D(1,0)*x(3) - D(1,1)*x(4) - B,
            x(2) - D(2,0)*x(3) - D(2,1)*x(4) - C,
            b*x(0) + a*x(1) - 2.*x(2) - M,
            sqHf + k*(x(0)+x(1)+x(2)) - c;
        if (Norm(R)<tolR) break;
        J = 1., 0., 0., -D(0,0), -D(0,1),
            0., 1., 0., -D(1,0), -D(1,1),
            0., 0., 1., -D(2,0), -D(2,1),
            b,  a, -2., 0., 0.,
            (2.*x(0)-x(1)-x(2))/sqHf+k, (2.*x(1)-x(2)-x(0))/sqHf+k, (2.*x(2)-x(0)-x(1))/sqHf+k, 0., 0.;
        Inv (J, Ji);
        dx = Ji * R;
        x -= dx;
    }
    if (it>=maxit) throw new Fatal("kTgIncs failed after %d iterations",maxit);
    if (it>ITMAX) ITMAX = it;
    dsig = x(0)-Sta->Sig(0), x(1)-Sta->Sig(1), x(2)-Sta->Sig(2), 0., 0., 0.;
    deps = x(3), x(4), dez, 0., 0., 0.;
    */

    // results
    dsig = dsx, dsy, dsz, 0., 0., 0.;
    deps = dex, dey, dez, 0., 0., 0.;

    //printf("diff = %g, %g, %g, %g, %g\n",dsx-dsig(0), dsy-dsig(1), dsz-dsig(2), dex-deps(0), dey-deps(1));
    
    // calc divs
    Vec_t deps_tmp(deps), dsig_tmp(dsig);
    Mdl->TgIncs (Sta, deps_tmp, dsig_tmp, divs);
}

Sparse::Triplet<double,int> D11,D12,D21,D22;

void TgIncs (Model const * Mdl, EquilibState const * Sta, double dT, Vec_t const & DEps, Vec_t const & DSig, Array<bool> const & pDEps, Vec_t & deps, Vec_t & dsig, Vec_t & divs)
{
    // zero matrices
    D11.ResetTop ();
    D12.ResetTop ();
    D21.ResetTop ();
    D22.ResetTop ();
    set_to_zero  (deps);
    set_to_zero  (dsig);

    // get stiffness
    Mat_t D;
    Mdl->Stiffness (Sta, D);

    // assemble
    for (size_t i=0; i<6; ++i)
    for (size_t j=0; j<6; ++j)
    {
        if      (!pDEps[i] && !pDEps[j]) D11.PushEntry (i,j, D(i,j));
        else if (!pDEps[i] &&  pDEps[j]) D12.PushEntry (i,j, D(i,j));
        else if ( pDEps[i] && !pDEps[j]) D21.PushEntry (i,j, D(i,j));
        else if ( pDEps[i] &&  pDEps[j]) D22.PushEntry (i,j, D(i,j));
    }

    // modify (augment) D11 and W
    Vec_t W(6);
    for (size_t i=0; i<6; ++i) // for each row in [D] and {W}
    {
        if (!pDEps[i]) // deps not prescribed => {dsig1}
        {
            dsig(i) = dT*DSig(i); // set dsig1 equal to dT*DSig1
            W   (i) = dsig(i);    // set W1 equal to dsig1
        }
        else // deps prescribed => {dsig2}
        {
            dsig(i)   = 0.0;         // clear dsig2
            W   (i)   = dT*DEps(i);  // set W2 equal to dT*DEps2
            D11.PushEntry (i,i,1.0); // modify D11
        }
    }

    // calc deps and dsig
    Sparse::SubMult (D12,    W,    W); // W1    -= D12*deps2
    UMFPACK::Solve  (D11,    W, deps); // deps   = inv(D11)*W
    Sparse::AddMult (D21, deps, dsig); // dsig2 += D21*deps1
    Sparse::AddMult (D22, deps, dsig); // dsig2 += D22*deps2

    // calc divs
    Mdl->TgIncs (Sta, deps, dsig, divs);

    // check
    if (false)
    {
        double dsa = (dsig(1)-dsig(0))/SQ2;
        double dsb = (dsig(1)+dsig(0)-2.0*dsig(2))/SQ2/SQ3;
        double lod = (fabs(dsa)>1.0e-13 ? atan2(dsb,dsa)*180.0/PI : 90.0);
        double dp  = -(dsig(0)+dsig(1)+dsig(2))/SQ3;
        double dez = deps(2);
        Vec_t dsig2, deps2, divs2;
        zTgIncs (Mdl, Sta, lod, dp, dez, deps2, dsig2, divs2);
        Vec_t diff_dsig(dsig-dsig2), diff_deps(deps-deps2);
        double err_dsig = Norm(diff_dsig);
        double err_deps = Norm(diff_deps);
        if (err_dsig>1.0e-12) throw new Fatal("TgIncs:: Error with dsig = %g",err_dsig);
        if (err_deps>1.0e-14) throw new Fatal("TgIncs:: Error with deps = %g",err_deps);
    }
}

struct PathIncs
{
    double dsx, dsy, dsz, dsxy, dsyz, dszx; // stress increments
    double dex, dey, dez, dexy, deyz, dezx; // strain increments (percentage)
    double lode, dp;                        // path given Lode angle (deg), dpoct and dez (percentage)
    bool   zPath;                           // use lode, dp and dez ?
    int    ninc;                            // number of increments for this path. -1 => use general
    double k;                               // path given Lode, k=dqoct/dpoct, and dez
    bool   kPath;                           // with k=Dq/Dp
    PathIncs () : dsx(0.),dsy(0.),dsz(0.),dsxy(0.),dsyz(0.),dszx(0.), 
                  dex(0.),dey(0.),dez(0.),dexy(0.),deyz(0.),dezx(0.),
                  lode(0.),dp(0.),zPath(false),ninc(-1),k(0.),kPath(false) {}
};

std::ostream & operator<< (std::ostream & os, PathIncs const & P)
{
    os << "ninc="<<P.ninc << " ";
    if      (P.kPath) { os << "lode="<<P.lode << " kcam="<<P.k*3./SQ2 << " dez="<<P.dez << std::endl; }
    else if (P.zPath) { os << "lode="<<P.lode << " dpcam="<<P.dp/SQ3  << " dez="<<P.dez << std::endl; }
    else
    {
        os << "ds=["<<P.dsx << " "<<P.dsy << " "<<P.dsz << " "<<P.dsxy << " "<<P.dsyz << " "<<P.dszx << "] ";
        os << "de=["<<P.dex << " "<<P.dey << " "<<P.dez << " "<<P.dexy << " "<<P.deyz << " "<<P.dezx << "]\n";
    }
    return os;
}

std::ostream & operator<< (std::ostream & os, Array<PathIncs> const & A)
{
    for (size_t i=0; i<A.Size(); ++i) os << "  " << A[i];
    return os;
}

int main(int argc, char **argv) try
{
    // input filename
    String inp_fname("driver.inp");
    String mat_fname("materials.inp");
    if (argc>1) inp_fname = argv[1];
    if (argc>2) mat_fname = argv[2];

    // default values
    int    matID  = 0;     // material ID
    double pcam0  = 100.0; // pCam
    size_t gninc  = 10;    // general number of increments
    bool   cdrift = true;  // correct YS drift
    Array<PathIncs> path;  // path increments

    // parse input file
    std::fstream inp_file(inp_fname.CStr(), std::ios::in);
    if (!inp_file.is_open()) throw new Fatal("Could not open file <%s>",inp_fname.CStr());
    bool   reading_path = false;
    int    ndat         = -1;
    size_t line_num     = 1;
    int    idxdat       = 0;
    size_t idxpath      = 0;
    while (!inp_file.eof())
    {
        String line,key,equal;
        double val;
        std::getline (inp_file,line);
        std::istringstream iss(line);
        if (iss >> key >> equal >> val)
        {
            //cout << key << "  " << val << endl;
            if (key[0]=='#') { line_num++; continue; }
            if (reading_path)
            {
                if      (key=="ndat") ndat = (int)val;
                else if (ndat<0) throw new Fatal("Error in file <%s> @ line # %d: key 'ndat' must come before data. Key==%s is in the wrong position",inp_fname.CStr(),line_num,key.CStr());
                else if (key=="kcam")  { path[idxpath].k    = SQ2*val/3.; path[idxpath].kPath=true; path[idxpath].zPath=false; idxdat++; }
                else if (key=="dpcam") { path[idxpath].dp   = val*SQ3;    path[idxpath].zPath=true; path[idxpath].kPath=false; idxdat++; }
                else if (key=="lode")  { path[idxpath].lode = val;       idxdat++; if (val<30. || val>90.) throw new Fatal("Lode angle alpha must be inside [30,90]. Alpha==%g is invalid",val); }
                else if (key=="dex")   { path[idxpath].dex  = val/100.;  idxdat++; }
                else if (key=="dey")   { path[idxpath].dey  = val/100.;  idxdat++; }
                else if (key=="dez")   { path[idxpath].dez  = val/100.;  idxdat++; }
                else if (key=="dexy")  { path[idxpath].dexy = val/100.;  idxdat++; }
                else if (key=="deyz")  { path[idxpath].deyz = val/100.;  idxdat++; }
                else if (key=="dezx")  { path[idxpath].dezx = val/100.;  idxdat++; }
                else if (key=="dsx")   { path[idxpath].dsx  = val;       idxdat++; }
                else if (key=="dsy")   { path[idxpath].dsy  = val;       idxdat++; }
                else if (key=="dsz")   { path[idxpath].dsz  = val;       idxdat++; }
                else if (key=="dsxy")  { path[idxpath].dsxy = val;       idxdat++; }
                else if (key=="dsyz")  { path[idxpath].dsyz = val;       idxdat++; }
                else if (key=="dszx")  { path[idxpath].dszx = val;       idxdat++; }
                else if (key=="ninc")  { path[idxpath].ninc = (int)val;  idxdat++; }
                else throw new Fatal("Error in file<%s> @ line # %d: reading data of path # %d. Key==%s is invalid",inp_fname.CStr(),line_num,idxpath,key.CStr());
                if (idxdat==ndat)
                {
                    ndat   = -1;
                    idxdat = 0;
                    idxpath++;
                    if (idxpath==path.Size()) break;
                }
            }
            else
            {
                if      (key=="matID")  matID  = (int)val;
                else if (key=="pcam0")  pcam0  = val;
                else if (key=="ninc")   gninc  = (int)val;
                else if (key=="cdrift") cdrift = (bool)val;
                else if (key=="npath")
                {
                    path.Resize ((size_t)val);
                    reading_path = true;
                }
                else throw new Fatal("Error in file <%s> @ line #: Key==%s in invalid",inp_fname.CStr(),line_num,key.CStr());
            }
        }
        line_num++;
    }
    if (idxpath!=path.Size()) throw new Fatal("Error in file <%s>: Not all path data could be read for npath==%d",inp_fname.CStr(),path.Size());
    printf("Input data:\n");
    printf("  matID  = %d\n",matID);
    printf("  pcam0  = %g\n",pcam0);
    printf("  ninc   = %zd\n",gninc);
    printf("  cdrift = %d\n",cdrift);

    // default initial values
    Dict prms, inis;
    inis.Set (-1, "sx sy sz", -pcam0,-pcam0,-pcam0);

    // parse materials file
    String model_name;
    ReadMaterial (-1, matID, mat_fname.CStr(), model_name, prms, inis);
    printf("\nMaterial data:\n");
    printf("  name = %s\n",model_name.CStr());
    printf("  prms : "); cout << prms << endl;
    printf("  inis : "); cout << inis << endl;

    // show path
    printf("\nPath:\n");
    cout << path << endl;

    // allocate model
    Model * mdl = AllocModel (model_name, /*NDim*/3, prms(-1));

    // set initial state
    EquilibState sta(/*NDim*/3);
    mdl->InitIvs (inis(-1), &sta);

    // number of internal values
    size_t niv = size(sta.Ivs);

    // allocate memory
    D11.AllocSpace (6,6, 50);
    D12.AllocSpace (6,6, 50);
    D21.AllocSpace (6,6, 50);
    D22.AllocSpace (6,6, 50);

    // output initial state
    String res_fname = inp_fname.substr(0,inp_fname.size()-3) + "res";
    std::ofstream of(res_fname.CStr(), std::ios::out);
    sta.Output (of, /*header*/true, "%17.8e");

    // auxiliar variables
    Vec_t DEps   (6), DSig   (6);               // total increments for one stage
    Vec_t deps   (6), dsig   (6);               // subincrements
    Vec_t deps_tr(6), dsig_tr(6), divs_tr(niv); // trial increments
    Vec_t deps_1 (6), dsig_1 (6), divs_1 (niv); // intermediate increments
    Vec_t deps_2 (6), dsig_2 (6), divs_2 (niv); // ME increments
    Vec_t eps_dif(6), sig_dif(6);               // ME - FE strain/stress difference
    EquilibState sta_1 (/*NDim*/3);             // intermediate state
    EquilibState sta_ME(/*NDim*/3);             // Modified-Euler state
    Array<bool> pDEps(6);                       // prescribed strain components ?

    // constants
    double STOL   = 1.0e-5;
    double dTini  = 1.0;
    double mMin   = 0.1;
    double mMax   = 10.0;
    double MaxSS  = 2000;

    // for each state
    double lode, dp, dez, dexy, deyz, dezx, koct;
    for (size_t i=0; i<path.Size(); ++i)
    {
        // number of increments
        int ninc = (path[i].ninc<0 ? gninc : path[i].ninc);

        // set prescribed increments
        set_to_zero (DEps);
        set_to_zero (DSig);
        pDEps = false,false,false, false,false,false;
        bool zpath = path[i].zPath;
        bool kpath = path[i].kPath;
        if (kpath)
        {
            lode = path[i].lode;
            koct = path[i].k;
            dez  = path[i].dez/ninc;
        }
        else if (zpath)
        {
            lode = path[i].lode;
            dp   = path[i].dp  /ninc;
            dez  = path[i].dez /ninc;
            dexy = path[i].dexy/ninc;
            deyz = path[i].deyz/ninc;
            dezx = path[i].dezx/ninc;
        }
        else
        {
            DSig(0) = path[i].dsx;
            DSig(1) = path[i].dsy;
            DSig(2) = path[i].dsz;
            if (fabs(path[i].dex)>1.0e-10) { DEps(0) = path[i].dex;  pDEps[0] = true; }
            if (fabs(path[i].dey)>1.0e-10) { DEps(1) = path[i].dey;  pDEps[1] = true; }
            if (fabs(path[i].dez)>1.0e-10) { DEps(2) = path[i].dez;  pDEps[2] = true; }
            deps = DEps/ninc;
            dsig = DSig/ninc;
        }

        // for each increment
        cout << "Increment = ";
        for (int j=0; j<ninc; ++j)
        {
            // trial increments
            if      (kpath) kTgIncs (mdl, &sta, lode, koct, dez,              deps_tr, dsig_tr, divs_tr);
            else if (zpath) zTgIncs (mdl, &sta, lode, dp,   dez,              deps_tr, dsig_tr, divs_tr, dexy, deyz, dezx);
            else             TgIncs (mdl, &sta, /*dT*/1.0, deps, dsig, pDEps, deps_tr, dsig_tr, divs_tr);

            // loading-unloading ?
            double aint = -1.0; // no intersection
            bool   ldg  = mdl->LoadCond (&sta, deps_tr, aint); // returns true if there is loading (also when there is intersection)

            // with intersection ?
            if (aint>0.0 && aint<1.0)
            {
                // update to intersection
                if      (kpath) kTgIncs (mdl, &sta, lode,    koct, dez*aint,       deps_tr, dsig_tr, divs_tr);
                else if (zpath) zTgIncs (mdl, &sta, lode, dp*aint, dez*aint,       deps_tr, dsig_tr, divs_tr, dexy*aint, deyz*aint, dezx*aint);
                else             TgIncs (mdl, &sta, /*dT*/aint, deps, dsig, pDEps, deps_tr, dsig_tr, divs_tr);
                sta.Eps += deps_tr;
                sta.Sig += dsig_tr;
                sta.Ivs += divs_tr;
                deps = fabs(1.0-aint)*deps; // remaining of deps to be applied

                // drift correction
                if (cdrift) mdl->CorrectDrift (&sta);

                // output
                sta.Output (of, /*header*/false, "%17.8e");
            }

            // set loading flag (must be after intersection because the tgIncs during intersection must be calc with Ldg=false)
            sta   .Ldg = ldg;
            sta_1 .Ldg = ldg;
            sta_ME.Ldg = ldg;

            // update stress path in model
            mdl->UpdatePath (&sta, deps_tr, dsig_tr);

            // for each pseudo time T
            double T  = 0.0;
            double dT = dTini;
            size_t k  = 0;
            for (k=0; k<MaxSS; ++k)
            {
                // exit point
                if (T>=1.0) break;

                // FE and ME steps
                if      (kpath) kTgIncs (mdl, &sta, lode,  koct, dez*dT,   deps_1, dsig_1, divs_1);
                else if (zpath) zTgIncs (mdl, &sta, lode, dp*dT, dez*dT,   deps_1, dsig_1, divs_1, dexy*dT, deyz*dT, dezx*dT);
                else             TgIncs (mdl, &sta, dT, deps, dsig, pDEps, deps_1, dsig_1, divs_1);
                sta_1.Eps = sta.Eps + deps_1;
                sta_1.Sig = sta.Sig + dsig_1;
                sta_1.Ivs = sta.Ivs + divs_1;
                if      (kpath) kTgIncs (mdl, &sta_1, lode,  koct, dez*dT,   deps_2, dsig_2, divs_2);
                else if (zpath) zTgIncs (mdl, &sta_1, lode, dp*dT, dez*dT,   deps_2, dsig_2, divs_2, dexy*dT, deyz*dT, dezx*dT);
                else             TgIncs (mdl, &sta_1, dT, deps, dsig, pDEps, deps_2, dsig_2, divs_2);
                sta_ME.Eps = sta.Eps + 0.5*(deps_1+deps_2);
                sta_ME.Sig = sta.Sig + 0.5*(dsig_1+dsig_2);
                sta_ME.Ivs = sta.Ivs + 0.5*(divs_1+divs_2);

                // local error estimate
                eps_dif = sta_ME.Eps - sta_1.Eps;
                sig_dif = sta_ME.Sig - sta_1.Sig;
                double eps_err = Norm(eps_dif)/(1.0+Norm(sta_ME.Eps));
                double sig_err = Norm(sig_dif)/(1.0+Norm(sta_ME.Sig));
                double ivs_err = 0.0;
                for (size_t i=0; i<niv; ++i) ivs_err += fabs(sta_ME.Ivs(i)-sta_1.Ivs(i))/(1.0+fabs(sta_ME.Ivs(i)));
                double error = eps_err + sig_err + ivs_err;

                // step multiplier
                double m = (error>0.0 ? 0.9*sqrt(STOL/error) : mMax);

                // update
                if (error<STOL)
                {
                    // update state
                    T += dT;
                    sta.Eps = sta_ME.Eps;
                    sta.Sig = sta_ME.Sig;
                    sta.Ivs = sta_ME.Ivs;

                    // drift correction
                    if (cdrift) mdl->CorrectDrift (&sta);

                    // update stress path in model
                    mdl->UpdatePath (&sta, Vec_t(0.5*(deps_1+deps_2)), Vec_t(0.5*(dsig_1+dsig_2)));

                    // limit change on stepsize
                    if (m>mMax) m = mMax;

                    // output
                    sta.Output (of, /*header*/false, "%17.8e");
                }
                else if (m<mMin) m = mMin;

                // change next step size
                dT = m * dT;

                // check for last increment
                if (dT>1.0-T) dT = 1.0-T;
            }
            if (k>=MaxSS) throw new Fatal("main: Modified-Euler did not converge after %d substeps",k);
            cout << j << " ";
        }
        cout << "\n";
    }

    //cout << "ITMAX = " << ITMAX << endl;

    // end
    of.close();
    cout << "\nFile <" << TERM_CLR_BLUE << res_fname << TERM_RST << "> written\n\n";
    delete mdl;
    return 0;
}
MECHSYS_CATCH
