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

#ifndef MPM_PROBLEMS2D_H
#define MPM_PROBLEMS2D_H

// FLTK
#include <FL/Fl_Choice.H>

// MechSys
#include <mechsys/util/array.h>
#include<mechsys/models/linelastic.h>
#include<mechsys/models/neohookean.h>
#include<mechsys/models/camclay.h>
#include<mechsys/models/elastoplastic.h>

// Local
#include <mechsys/mpm/defs.h>
#include <mechsys/mpm/grid2d.h>
#include <mechsys/mpm/mpoints2d.h>

namespace MPM {

// Fixed nodes?
bool CbIsFixed1 (Vector3D const & N, FixType & FType)
{
	if (N(0)<1.01)              { FType = FIX_XY; return true; }
	if (N(1)<0.11 || N(1)>0.19) { FType = FIX_Y;  return true; }
	return false;
}
bool CbIsFixed2 (Vector3D const & N, FixType & FType)
{
	if (N(0)<0.01 && N(1)<0.01) { FType=FIX_XY; return true; }
	if (N(0)<0.01)              { FType=FIX_X;  return true; }
	if (N(1)<0.01)              { FType=FIX_Y;  return true; }
	return false;
}
bool CbIsFixed3 (Vector3D const & N, FixType & FType)
{
	/*
	if (N(0)>-0.051 && N(0)<-0.049 && N(1)>-0.051 && N(1)<-0.049) { FType=FIX_XY; return true; }
	if (N(0)> 0.99  && N(0)< 1.01  && N(1)>-0.051 && N(1)<-0.049) { FType=FIX_XY; return true; }
	if (N(0)> 0.99  && N(0)< 1.01  && N(1)> 0.99  && N(1)<1.01)   { FType=FIX_XY; return true; }
	if (N(0)>-0.051 && N(0)<-0.049 && N(1)> 0.99  && N(1)<1.01)   { FType=FIX_XY; return true; }
	if (N(0)>-0.051 && N(0)<-0.049 && N(1)>-0.049 && N(1)<1.01)   { FType=FIX_X;  return true; }
	if (N(0)< 1.01  && N(0)>-0.051 && N(1)>-0.051 && N(1)<-0.049) { FType=FIX_Y;  return true; }
	if (N(0)> 0.99  && N(0)< 1.01  && N(1)>-0.049 && N(1)<1.01)   { FType=FIX_X;  return true; }
	if (N(0)< 1.01  && N(0)>-0.051 && N(1)> 0.99  && N(1)<1.01)   { FType=FIX_Y;  return true; }
	*/
	if (N(0)>-0.051 && N(0)<-0.049 && N(1)>-0.051 && N(1)<-0.049) { FType=FIX_XY; return true; }
	if (N(0)> 0.69  && N(0)< 0.71  && N(1)>-0.051 && N(1)<-0.049) { FType=FIX_XY; return true; }
	if (N(0)> 0.69  && N(0)< 0.71  && N(1)> 0.69  && N(1)<0.71)   { FType=FIX_XY; return true; }
	if (N(0)>-0.051 && N(0)<-0.049 && N(1)> 0.69  && N(1)<0.71)   { FType=FIX_XY; return true; }
	if (N(0)>-0.051 && N(0)<-0.049 && N(1)>-0.049 && N(1)<0.71)   { FType=FIX_X;  return true; }
	if (N(0)< 0.71  && N(0)>-0.051 && N(1)>-0.051 && N(1)<-0.049) { FType=FIX_Y;  return true; }
	if (N(0)> 0.69  && N(0)< 0.71  && N(1)>-0.049 && N(1)<0.71)   { FType=FIX_X;  return true; }
	if (N(0)< 0.71  && N(0)>-0.051 && N(1)> 0.69  && N(1)<0.71)   { FType=FIX_Y;  return true; }
	return false;
}
bool CbIsFixed4 (Vector3D const & N, FixType & FType)
{
	if (N(0)<0.01 && N(1)<0.01) { FType=FIX_XY; return true; }
	if (N(0)<0.01)              { FType=FIX_X;  return true; }
	if (N(1)<0.01)              { FType=FIX_Y;  return true; }
	return false;
}
bool CbIsFixed5 (Vector3D const & N, FixType & FType) // micrometer
{
	if (N(0)>-0.1 && N(0)<0.1) { FType=FIX_X; return true; }
	return false;
}
bool CbIsFixed6 (Vector3D const & N, FixType & FType)
{
	/*
	if (N(0)>-0.01  && N(0)<0.01   && N(1)>-0.01  && N(1)<0.01)  { FType=FIX_XY; return true; }
	if (N(0)> 19.99 && N(0)< 20.01 && N(1)>-0.01  && N(1)<0.01)  { FType=FIX_XY; return true; }
	if (N(0)> 19.99 && N(0)< 20.01 && N(1)> 19.99 && N(1)<20.01) { FType=FIX_XY; return true; }
	if (N(0)>-0.01  && N(0)<0.01   && N(1)> 19.99 && N(1)<20.01) { FType=FIX_XY; return true; }
	if (N(0)>-0.01  && N(0)<0.01   && N(1)>0.01   && N(1)<20.01) { FType=FIX_X;  return true; }
	if (N(0)< 20.01 && N(0)>-0.01  && N(1)>-0.01  && N(1)<0.01)  { FType=FIX_Y;  return true; }
	if (N(0)> 19.99 && N(0)< 20.01 && N(1)>0.01   && N(1)<20.01) { FType=FIX_X;  return true; }
	if (N(0)< 20.01 && N(0)>-0.01  && N(1)> 19.99 && N(1)<20.01) { FType=FIX_Y;  return true; }
	*/
	return false;
}
bool CbIsFixed7 (Vector3D const & N, FixType & FType)
{
	if (N(0)<0.01) { FType=FIX_XY; return true; }
	else return false;
}
bool CbIsFixed8 (Vector3D const & N, FixType & FType)
{
	if (N(1)<0.01) { FType=FIX_XY; return true; }
	else return false;
}
bool CbIsFixed10(Vector3D const & N, FixType & FType)
{
	if (N(0)<0.01) { FType = FIX_XY; return true; }
	else           { FType = FIX_Y;  return true; }
}
bool CbIsFixed11(Vector3D const & N, FixType & FType)
{
	     if (N(0)<0.01 && N(1)<0.01) { FType = FIX_XY; return true; }
	else if (N(0)<0.01)              { FType = FIX_X;  return true; }
	else if (N(1)<0.01)              { FType = FIX_Y;  return true; }
	return false;
}
bool CbIsFixed12(Vector3D const & N, FixType & FType)
{
	if (N(1)<0.01) // bottom
	{
		if (fabs(N(0)-20.0)<0.01) { FType = FIX_XY; return true; }
		else                      { FType = FIX_Y;  return true; }
	}
	return false;
}
bool CbIsFixed13(Vector3D const & N, FixType & FType)
{
	if (N(1)<0.01) // bottom
	{
		if (N(0)<0.01 || N(0)>14.95) { FType = FIX_XY; return true; }
		else                         { FType = FIX_Y;  return true; }
	}
    else if (N(0)<0.01 || N(0)>14.95) // left or right
    {
        FType = FIX_X; return true;
    }
	return false;
}
bool CbIsFixed14(Vector3D const & N, FixType & FType)
{
    double L = 30.0;
	if (N(1)<0.01) // bottom
	{
		if (N(0)<0.01 || N(0)>L-0.5) { FType = FIX_XY; return true; }
		else                         { FType = FIX_Y;  return true; }
	}
    else if (N(0)<0.01 || N(0)>L-0.5) // left or right
    {
        FType = FIX_X; return true;
    }
	return false;
}
bool CbIsFixed15(Vector3D const & N, FixType & FType)
{
    if (fabs(N(0)-5.0)<0.01 && fabs(N(1)-5.0)<0.01)
    {
        FType = FIX_XY;
        return true;
    }
    return false;
}
bool CbIsFixed16(Vector3D const & N, FixType & FType)
{
    if (fabs(N(1)-9.5)<0.05)
    {
        FType = FIX_Y;
        return true;
    }
    return false;
}
bool CbIsFixed17(Vector3D const & N, FixType & FType)
{
    if (N(1)>9.4)
    {
        if (fabs(N(0)-1.0)<0.01) FType = FIX_XY;
        else                     FType = FIX_Y;
        return true;
    }
    return false;
}

// Is point inside geometry callbacks
bool CbIsPointInGeom1 (Vector3D const & P, int & ClrIdx)
{

	if (P(0)>=1.0 && P(0)<=2.0 && P(1)>=1.50 && P(1)<=2.0) return true;
	else                                                   return false;

}
bool CbIsPointInGeom2 (Vector3D const & P, int & ClrIdx)
{
	if (P(0)>0.0 && P(0)<1.0 && P(1)>0.0 && P(1)<1.0) return true;
	else                                              return false;
}
bool CbIsPointInGeom3 (Vector3D const & P, int & ClrIdx)
{
	// Cylinders
	/*
	double xc1=0.25; double yc1=0.25; double r1=0.14;
	double xc2=0.70; double yc2=0.70; double r2=0.14;
	*/
	double xc1=0.15; double yc1=0.15; double r1=0.14;
	double xc2=0.50; double yc2=0.50; double r2=0.14;

	// Bottom-left cylinder
	double d1 = sqrt(pow(P(0)-xc1,2.0)+pow(P(1)-yc1,2.0));

	// Upper-right cylinder
	double d2 = sqrt(pow(P(0)-xc2,2.0)+pow(P(1)-yc2,2.0));

	if (d1<=r1 || d2<=r2) return true;
	else return false;
	/*
	else
	{
		if (P(0)<0.05 && P(0)>0    && P(1)>0    && P(1)<0.95) return true; // left
		if (P(0)>0.9  && P(0)<0.95 && P(1)>0    && P(1)<0.95) return true; // right
		if (P(0)>0    && P(0)<0.95 && P(1)<0.05 && P(1)>0   ) return true; // bottom
		if (P(0)>0    && P(0)<0.95 && P(1)>0.9  && P(1)<0.95) return true; // top
		return false;
	}
	*/
}
bool CbIsPointInGeom4 (Vector3D const & P, int & ClrIdx)
{
	if (P(0)>0.0 && P(0)<1.0 && P(1)>0.0 && P(1)<1.0) return true;
	else                                              return false;
}
bool CbIsPointInGeom5 (Vector3D const & P, int & ClrIdx) // micrometer
{
	if (P(0)>0.0 && P(0)<60.0 && P(1)>0.0 && P(1)<40.0) return true;
	else                                                return false;
}
bool CbIsPointInGeom6 (Vector3D const & P, int & ClrIdx)
{
	// Cylinders
	double xc=10.0; double yc=10.0; double r=2.5;

	// Distance
	double d = sqrt(pow(P(0)-xc,2.0)+pow(P(1)-yc,2.0));

	if (d<=r) return true;
	else
	{
		if (P(0)> 2.0 && P(0)< 3.0 && P(1)> 2.0 && P(1)<18.0) return true;
		if (P(0)>17.0 && P(0)<18.0 && P(1)> 2.0 && P(1)<18.0) return true;
		if (P(0)> 2.0 && P(0)<18.0 && P(1)>17.0 && P(1)<18.0) return true;
		if (P(0)> 2.0 && P(0)<18.0 && P(1)> 2.0 && P(1)< 3.0) return true;
		return false;
	}
}
bool CbIsPointInGeom7 (Vector3D const & P, int & ClrIdx)
{
	if (P(0)>0.0 && P(0)<19.0 && P(1)>10.0 && P(1)<12.0) return true;
	else                                                return false;
}
bool CbIsPointInGeom8 (Vector3D const & P, int & ClrIdx)
{
	if (P(0)>5.0 && P(0)<7.0 && P(1)>0.0 && P(1)<12.0) return true;
	else                                               return false;
}
bool CbIsPointInGeom10(Vector3D const & P, int & ClrIdx)
{
	if (P(0)>=0.0 && P(0)<=25.0 && P(1)>=1.0 && P(1)<=2.0) return true;
	else                                                   return false;
}
bool CbIsPointInGeom11(Vector3D const & P, int & ClrIdx)
{
	if (P(0)>0.0 && P(0)<10.01 && P(1)>0.0 && P(1)<10.01)
	{
		double xc1=0.0; double yc1=0.0; double r1=1.0;
		double d1 = sqrt(pow(P(0)-xc1,2.0)+pow(P(1)-yc1,2.0));
		if (d1>=r1) return true;
	}
	return false;
}
bool CbIsPointInGeom12(Vector3D const & P, int & ClrIdx)
{
    if (P(1)<=30.0-P(0) && P(1)<=P(0)-10.0 && P(1)>=0.0) return true;
	return false;
}
bool CbIsPointInGeom13(Vector3D const & P, int & ClrIdx)
{
    if (P(0)>=0.0 && P(0)<=15.0 && P(1)>0.0 && P(1)<=10.0) return true;
	return false;
}
bool CbIsPointInGeom14(Vector3D const & P, int & ClrIdx)
{
    // ground
    double L = 30.0;
    double H = 10.0;
    if (P(0)>=0.0 && P(0)<=L && P(1)>0.0 && P(1)<=H)
    {
        ClrIdx = -((static_cast<int>(P(1)) % 2) + 1);
        return true;
    }

    // comet
    double xc = L/2.0;
    double yc = H+3.0;
    double r  = 2.0;
    if (sqrt(pow(P(0)-xc,2.)+pow(P(1)-yc,2.))<=r)
    {
        ClrIdx = -3;
        return true;
    }

    // default
	return false;
}
bool CbIsPointInGeom15(Vector3D const & P, int & ClrIdx)
{
    double xc = 5.0;
    double yc = 5.0;
    double r  = 4.0;
    if (sqrt(pow(P(0)-xc,2.)+pow(P(1)-yc,2.))<=r) return true;
    return false;
}
bool CbIsPointInGeom16(Vector3D const & P, int & ClrIdx)
{
    if ((P(0)>0.5 && P(0)<1.5) && (P(1)>8.5 && P(1)<9.5)) return true;
    return false;
}
bool CbIsPointInGeom17(Vector3D const & P, int & ClrIdx)
{
    if ((P(0)>0.5 && P(0)<1.5) && (P(1)>8.5 && P(1)<9.5))
    {
        double xmin = 0.5;
        double l    = 0.125;
        ClrIdx = -((static_cast<int>((P(0)-xmin)/l) % 2) + 1);
        return true;
    }
    return false;
}

// Density callbacks
double CbDensity1 (Vector3D const & P) { return 1.0; }
double CbDensity2 (Vector3D const & P) { return 1.0; }
double CbDensity3 (Vector3D const & P) { return 1.0; }
double CbDensity4 (Vector3D const & P) { return 1.0; }
double CbDensity5 (Vector3D const & P) { return 2.71; } // picogram/(micrometer*nanoseconds)
double CbDensity6 (Vector3D const & P) { return 1.0; }
double CbDensity7 (Vector3D const & P) { return 1.0; }
double CbDensity8 (Vector3D const & P) { return 1.0; }
double CbDensity10(Vector3D const & P) { return 1.0; }
double CbDensity11(Vector3D const & P) { return 1.0; }
double CbDensity12(Vector3D const & P) { return 1.0; }
double CbDensity13(Vector3D const & P) { return 1.0; }
double CbDensity14(Vector3D const & P) { return 1.0; }
double CbDensity15(Vector3D const & P) { return 1.0; }
double CbDensity16(Vector3D const & P) { return 1050.0; } // kg/m3
double CbDensity17(Vector3D const & P) { return 1050.0; } // kg/m3

const double DTRUE = 1.0;

// Allocate model callbacks
void CbAllocMdl1 (Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu psa", 4.0*PI*PI, 0.0, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name   = "NeoHookean";
}
void CbAllocMdl2 (Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu psa", 100.0, 0.25, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name   = "NeoHookean";
}
void CbAllocMdl3 (Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu psa", 1.0, 0.2, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name   = "NeoHookean";
}
void CbAllocMdl4 (Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu psa", 210.0, 0.3, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name   = "NeoHookean";
}
void CbAllocMdl5 (Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
	double E  = 175.8; // E  picogram/(micrometer*nanoseconds) == GPa
	double nu = 0.28;  // nu
    Prms.Set ("E nu psa", E, nu, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name   = "NeoHookean";
}
void CbAllocMdl6 (Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
	double xc=10.0; double yc=10.0;
	double d = sqrt(pow(P(0)-xc,2.0)+pow(P(1)-yc,2.0));

	if (d<=2.5)
	{
        Prms.Set ("E nu psa", 1.0, 0.49, DTRUE);
        Inis.Set ("zero", DTRUE);
        Name   = "NeoHookean";
	}
	else
	{
        Prms.Set ("E nu sY VM psa", 10000.0, 0.1, 2.0, DTRUE, DTRUE);
        Inis.Set ("zero", DTRUE);
        Name = "ElastoPlastic";
	}
}
void CbAllocMdl7 (Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu psa", 1000.0, 0.25, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name   = "NeoHookean";
}
void CbAllocMdl7b (Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu sY VM psa", 100.0, 0.25, 3.0, DTRUE, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name = "ElastoPlastic";
}
void CbAllocMdl8 (Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu psa", 100.0, 0.25, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name   = "NeoHookean";
}
void CbAllocMdl10(Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu psa", 100.0, 0.0, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name   = "NeoHookean";
}
void CbAllocMdl11(Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu psa", 6778.0, 0.2103, DTRUE);
    Inis.Set ("sx sy sz", -30.0, -30.0, -30.0);
    Name   = "NeoHookean";
}
void CbAllocMdl12(Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    //Prms.Set ("E nu psa", 10.0, 0.3, DTRUE);
    //Inis.Set ("zero", DTRUE);
    //Name   = "NeoHookean";

    Inis.Set ("zero", DTRUE);
    Prms.Set ("E nu sY VM psa", 100000.0, 0.3, 2.0, DTRUE, DTRUE);
    //Prms.Set ("E nu sY VM psa", 100.0, 0.3, 2.0, DTRUE, DTRUE);
    Name   = "ElastoPlastic";
}
void CbAllocMdl13(Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu psa", 10.0, 0.3, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name   = "NeoHookean";
}
void CbAllocMdl14(Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    // comet
    double L  = 30.0;
    double H  = 10.0;
    double xc = L/2.0;
    double yc = H+3.0;
    double r  = 2.0;
    if (sqrt(pow(P(0)-xc,2.)+pow(P(1)-yc,2.))<=r)
    {
        Prms.Set ("E nu psa", 100000.0, 0.1, DTRUE);
        Inis.Set ("zero", DTRUE);
        Name   = "NeoHookean";
        Tag    = -1;
    }

    // ground
    else
    {
        Inis.Set ("zero", DTRUE);
        Prms.Set ("E nu sY VM psa", 100.0, 0.3, 2.0, DTRUE, DTRUE);
        //Prms.Set ("E nu c phi  DP psa", 100.0, 0.3, 0.0, 20.0,  DTRUE, DTRUE);
        Name   = "ElastoPlastic";
        Tag    = -2;
    }
}
void CbAllocMdl15(Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu psa", 1000.0, 0.1, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name   = "NeoHookean";
}
void CbAllocMdl16(Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu", 1.0e+6, 0.3, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name   = "NeoHookean";
}
void CbAllocMdl17(Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme)
{
    Prms.Set ("E nu", 1.0e+6, 0.3, DTRUE);
    Inis.Set ("zero", DTRUE);
    Name   = "NeoHookean";
}

// Initial velocity
void CbIniVeloc1 (Vector3D const & P, Vector3D & v) { v = 0.1, 0.0, 0.0; }
void CbIniVeloc2 (Vector3D const & P, Vector3D & v) { v = 0.0; }
void CbIniVeloc3 (Vector3D const & P, Vector3D & v)
{
	/*
	     if (P(0)<0.5 && P(1)<0.5 && P(0)>0.1 && P(1)>0.1) v =  0.15,  0.15, 0.0;
	else if (P(0)>0.5 && P(1)>0.5 && P(0)<0.9 && P(1)<0.9) v = -0.15, -0.1, 0.0;
	else                                                   v =  0.0,  0.0, 0.0;
	*/
	/*
	     if (P(0)<0.5 && P(1)<0.5 && P(0)>0.1 && P(1)>0.1) v =  0.1,  0.1, 0.0;
	else if (P(0)>0.5 && P(1)>0.5 && P(0)<0.9 && P(1)<0.9) v = -0.1, -0.1, 0.0;
	else                                                   v =  0.0,  0.0, 0.0;
	*/
	if (P(0)<0.325 && P(1)<0.325) v =  0.1,  0.1, 0.0;
	else                          v = -0.1, -0.1, 0.0;
}
void CbIniVeloc4 (Vector3D const & P, Vector3D & v) { v = 0.0; }
void CbIniVeloc5 (Vector3D const & P, Vector3D & v) { v = 0.0; } // micrometer/nanoseconds
void CbIniVeloc6 (Vector3D const & P, Vector3D & v) { v = 0.0; }
void CbIniVeloc7 (Vector3D const & P, Vector3D & v) { v = 0.0; }
void CbIniVeloc8 (Vector3D const & P, Vector3D & v) { v = 0.0; }
void CbIniVeloc10(Vector3D const & P, Vector3D & v)
{
	double v0   = 0.1;                          // initial velocity
	//double E    = 100;                          // Young modulus
	//double rho0 = 1.0;                          // initial density
	double mode = 5.0;                          // Mode
	double bet  = ((2.0*mode-1.0)/2.0)*PI/25.0; // Eigenvalues
	//double c    = sqrt(E/rho0);                 // elastic wave speed
	//double w    = bet*c;                        // frequency
	double vel0 = v0*sin(bet*P(0));             // initial velocity
	v = vel0, 0.0, 0.0;
	//double xcm  = 12.5;                                // centre of mass
	//double velcm= function(t) v0*cos(w*t)*sin(bet*xcm) // correct velocity (center of mass)
}
void CbIniVeloc11 (Vector3D const & P, Vector3D & v) { v = 0.0; }
void CbIniVeloc12 (Vector3D const & P, Vector3D & v) { v = 0.0; }
void CbIniVeloc13 (Vector3D const & P, Vector3D & v) { v = 0.0; }
void CbIniVeloc14 (Vector3D const & P, Vector3D & v)
{
    // comet
    double L  = 30.0;
    double H  = 10.0;
    double xc = L/2.0;
    double yc = H+3.0;
    double r  = 2.0;
    if (sqrt(pow(P(0)-xc,2.)+pow(P(1)-yc,2.))<=r) v = 0.0, -10.0, 0.0;
    else v = 0.0;
}
void CbIniVeloc15 (Vector3D const & P, Vector3D & v)
{
    // comet
    double xc   = 5.0;
    double yc   = 5.0;
    double dfdx = 2.0*(P(0)-xc);
    double dfdy = 2.0*(P(1)-yc);
    double d    = sqrt(pow(P(0)-xc,2.)+pow(P(1)-yc,2.));
    v = -dfdy, dfdx, 0.0;
    v *= (0.1*d*NormV(v));
}
void CbIniVeloc16 (Vector3D const & P, Vector3D & v) { v = 0.0; }
void CbIniVeloc17 (Vector3D const & P, Vector3D & v) { v = 0.0; }

/*   x  sqrt(x)  2*sqrt(x)  log2(2*sqrt(x))  1/sqrt(x)  1-1/sqrt(x)
 *   1        1          2                1       1            0.0
 *   4        2          4                2       0.5          0.5
 *  16        4          8                3       0.25         0.75
 *
 *        1                    4                     16
 *  1_ _______           1_ _______             1_ _______ 
 *    |       |        0.5_|   | _ |_3L/4    0.75>|_|_|_|_|<7L/8
 *  0_|   _   |_L/2      0_|___|___|            0_|_|_|_|_|
 *    |       |            |   | _ |_L/4          |_|_|_|_|
 * -1_|_______|         -1_|___|___|           -1_|_|_|_|_|<L/8    */

// Applied tractions
bool CbHasTraction1 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasTraction2 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t)
{
	double xmax = 1.0;
	double ymax = 1.0;
	if (N[0](1)>ymax-1.01*L(1) && N[1](0)<xmax/1.99) // half/upper cells
	{
		double sqnpc = sqrt(static_cast<double>(nPCell));
		if (P(1)>N[3](1)-0.51*L(1)/sqnpc)
		{
			double q = 5.0;
			t = 0.0,  -q*L(0)/sqnpc,  0.0;
			return true;
		}
	}
	return false;
}
bool CbHasTraction3 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasTraction4 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t)
{
	if (N[0](1)>0.99-L(1))
	{
		if (P(1)>N[0](1)+L(1)/2.01)
		{
			double p = 1.0;
			if (nPCell==1) t = 0.0, p*L(0)    , 0.0;
			else           t = 0.0, p*L(0)/2.0, 0.0;
			return true;
		}
	}
	return false;
}
bool CbHasTraction5 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t)
{
	if (N[0](0)>59.9-L(0))
	{
		if (P(0)>59.9-L(0)/2.0)
		{
			double p = 400.0; // 40 picogram/(micrometer*nanoseconds) == GPa
            double val = p*L(1)/sqrt(nPCell);
			t = val, 0.0, 0.0;
			return true;
		}
	}
	return false;
}
bool CbHasTraction6 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t)
{
	double xc=10.0; double yc=10.0;
	double f = 1000.0;
	double d = sqrt(pow(P(0)-xc,2.0)+pow(P(1)-yc,2.0));

	if (d<=2.5)
	{
		Vector3D c, r;
		c = xc, yc, 0.0;
		r = P - c;
		double n = sqrt(dot(r,r));
		if (n>0.0) r /= n;
		else       r = 0.0;
		t = f*r;
		return true;
	}
	else return false;
}
bool CbHasTraction7 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t)
{
	double npc = sqrt(static_cast<double>(nPCell));
	if (P(0)>19.0-L(0)/npc)
	{
		t = 0.0, -10.0, 0.0;
		return true;
	}
	return false;
}
bool CbHasTraction8 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t)
{
	double npc = sqrt(static_cast<double>(nPCell));
	if (P(0)>5.0 && P(0)<5.0+L(0)/npc && P(1)>12.0-L(1)/npc)
	{
		t = 0.0, -10.0, 0.0;
		return true;
	}
	return false;
}
bool CbHasTraction10(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasTraction11(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t)
{
	double xmax = 10.0;
	double ymax = 10.0;
	double sqnpc = sqrt(static_cast<double>(nPCell));
	if (N[0](1)>ymax-1.01*L(1) && 0.99*N[1](0)<xmax) // upper cells
	{
		if (P(1)>N[3](1)-0.51*L(1)/sqnpc)
		{
			double q = 30.0;
			t = 0.0,  -q*L(0)/sqnpc,  0.0;
			return true;
		}
	}
	if (N[0](0)>xmax-1.01*L(0) && 0.99*N[1](1)<ymax) // right cells
	{
		if (P(0)>N[1](0)-0.51*L(0)/sqnpc)
		{
			double q = 30.0;
			t = -q*L(0)/sqnpc,  0.0,  0.0;
			return true;
		}
	}
	return false;
}
bool CbHasTraction12(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasTraction13(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasTraction14(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasTraction15(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasTraction16(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasTraction17(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }

// Applied displacements
bool CbHasAppDisp1 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp2 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp3 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp4 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp5 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp6 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp7 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp8 (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp10(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp11(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp12(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp13(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t)
{
    if (P(0)>=0.0 && P(0)<=0.5 && P(1)>9.7)
    {
        t = 0., -0.1, 0.; // downward velocity
        return true;
    }
    return false;
}
bool CbHasAppDisp14(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp15(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp16(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }
bool CbHasAppDisp17(Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t) { return false; }

// Body mass callbacks
static inline void CbB1 (double t, Vector3D & B) { B=0.0; }
static inline void CbB2 (double t, Vector3D & B) { B=0.0; }
static inline void CbB3 (double t, Vector3D & B) { B=0.0; }
static inline void CbB4 (double t, Vector3D & B) { B=0.0; }
static inline void CbB5 (double t, Vector3D & B) { B=0.0; }
static inline void CbB6 (double t, Vector3D & B) { B=0.0; }//, -1.0, 0.0; }
static inline void CbB7 (double t, Vector3D & B) { B=0.0; }//, -1.0, 0.0; }
static inline void CbB8 (double t, Vector3D & B) { B=0.0; }//, -1.0, 0.0; }
static inline void CbB10(double t, Vector3D & B) { B=0.0; }
static inline void CbB11(double t, Vector3D & B) { B=0.0; }
static inline void CbB12(double t, Vector3D & B)
{
    //double tsw = 5.0; // very funny!
    double tsw = 20.0;
    double mul = (t<tsw ? sin(PI*t/(2.0*tsw)) : 1.0);
    B = 0.0, -mul, 0.0;
}
static inline void CbB13(double t, Vector3D & B) { B=0.0; }
static inline void CbB14(double t, Vector3D & B) { B=0.0; }
static inline void CbB15(double t, Vector3D & B) { B=0.0; }
static inline void CbB16(double t, Vector3D & B) { B=0.0, -1000.0, 0.0; } // m/s2
static inline void CbB17(double t, Vector3D & B) { B=0.0, -3500.0, 0.0; } // m/s2

// Loading multiplier callbacks
static inline void CbLdM1 (double t, double & M) { M=0.0;                           } // t in seconds (s)
static inline void CbLdM2 (double t, double & M) { M=t;                             } // t in seconds (s)
static inline void CbLdM3 (double t, double & M) { M=0.0;                           } // t in seconds (s)
static inline void CbLdM4 (double t, double & M) { M=(t<0.5 ? t/0.5 : 1.0);         } // t in seconds (s)
static inline void CbLdM5 (double t, double & M) { M=(t<60.0 ? t/60.0 : 1.0);       } // t in nano-seconds (ns)
//static inline void CbLdM5 (double t, double & M) { M=(t<10.0 ? t/10.0 : 1.0);       } // t in nano-seconds (ns)
//static inline void CbLdM6 (double t, double & M) { M=(t>0.5 ? (t<1 ? 1.0 : 0.0) : 0.0); } // t in seconds (s)
static inline void CbLdM6 (double t, double & M) { M=1.0; }//(t>0.01 ? (t<0.11 ? 1.0 : 0.0) : 0.0); } // t in seconds (s)
static inline void CbLdM7 (double t, double & M) { M=(t<10.0 ? t/10.0 : 1.0);       } // t in seconds (s)
static inline void CbLdM8 (double t, double & M) { M=(t<5.0 ? t/5.0 : 0.0);       } // t in seconds (s)
static inline void CbLdM10(double t, double & M) { M=0.0;                           } // t in seconds (s)
static inline void CbLdM11(double t, double & M) { M=(t<100.0 ? t/100.0 : 1.0);     } // t in seconds (s)
static inline void CbLdM12(double t, double & M) { M=1.0;                           } // t in seconds (s)
static inline void CbLdM13(double t, double & M) { M=(t<10.0 ? t/10.0 : 1.0);       } // t in seconds (s)
static inline void CbLdM14(double t, double & M) { M=0.0;                           } // t in seconds (s)
static inline void CbLdM15(double t, double & M) { M=0.0;                           } // t in seconds (s)
static inline void CbLdM16(double t, double & M) { M=0.0;                           } // t in seconds (s)
static inline void CbLdM17(double t, double & M) { M=0.0;                           } // t in seconds (s)

// Correct velocities callbacks
static inline void CbVel1 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.1*cos(2.0*PI*t),0,0; }
static inline void CbVel2 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel3 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel4 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel5 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel6 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel7 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel8 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel10(double t, Vector3D const & XY, Vector3D & Vel)
{
	double v0    = 0.1;                          // initial velocity
	double E     = 100;                          // Young modulus
	double rho0  = 1.0;                          // initial density
	double mode  = 5.0;                          // Mode
	double bet   = ((2.0*mode-1.0)/2.0)*PI/25.0; // Eigenvalues
	double c     = sqrt(E/rho0);                 // elastic wave speed
	double w     = bet*c;                        // frequency
	double xcm   = 12.5;                         // centre of mass
	double velcm = v0*cos(w*t)*sin(bet*xcm);     // correct velocity (center of mass)
	Vel = velcm, 0.0, 0.0;
}
static inline void CbVel11 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel12 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel13 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel14 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel15 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel16 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }
static inline void CbVel17 (double t, Vector3D const & XY, Vector3D & Vel) { Vel=0.0; }

// Correct stresses callbacks
static inline void CbStress1 (double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress2 (double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress3 (double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress4 (double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress5 (double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress6 (double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress7 (double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress8 (double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress10(double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress11(double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress12(double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress13(double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress14(double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress15(double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress16(double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }
static inline void CbStress17(double t, Vector3D const & XY, STensor2 & Sig) { Sig=0.0; }

/** Problem related data. */
struct Problem2D
{
	// Data
	int               Id;      ///< Problem number
	int               NPcell;  ///< Number of points per cell
	Grid2D          * G;       ///< The grid
	MPoints2D       * P;       ///< The material points
	ptIsPointInGeom   pGeom;   ///< Pointer to a function which determines whether the point is inside the geometry
	ptB               pB;      ///< Pointer to a function which returns the body mass B
	ptLdM             pM;      ///< Pointer to a function which returns the loadings multiplier M
	ptVeloc           pVeloc;  ///< Pointer to a function which returns the correct velocity
	ptStress          pStress; ///< Pointer to a function which returns the correct stress
	bool              CPDI;    ///< use CPDI
	bool              MPM;     ///< use original MPM
	bool              USF;     ///< update stress first ?
	bool              SmallSt; ///< small strains ?
	bool              FwEuler; ///< Use Forward Euler ?
	bool              SaveVTU; ///< Save vtu at each output timestep?
	bool              SavePNG; ///< Save png at each output timestep?
    bool              SavePointData; ///< Save point data
	double            t;       ///< current time
	double            tf;      ///< final time
	double            Dt;      ///< time increment
	double            Dtout;   ///< time increment for output
	size_t            Outp;    ///< point/particle number to output data
	Vector3D          B;       ///< current b-body mass
	double            M;       ///< Current multiplier for external forces

	// Output arrays
	Array<double>         Out_t;       // time array
	Array<double>         Out_Wx_t;    // Body weights (in time)
	Array<double>         Out_sE_t;    // Out strain Energy (in time)
	Array<double>         Out_kE_t;    // Out kinetic Energy (in time)
	Array<double>         Out_tE_t;    // Out total Energy (in time)
	Array<Array<double> > Out_p_vx_t;  // Out particle x-velocity (in time) Size=np
	Array<Array<double> > Out_p_vy_t;  // Out particle y-velocity (in time) Size=np
	Array<Array<double> > Out_p_cvx_t; // Out correct particle x-velocity (in time) Size=np
	Array<Array<double> > Out_p_sxx_t; // Out particle xx-stress (in time) Size=np
	Array<Array<double> > Out_p_syy_t; // Out particle yy-stress (in time) Size=np
	Array<Array<double> > Out_p_exx_t; // Out particle xx-strain (in time) Size=np
	Array<Array<double> > Out_p_eyy_t; // Out particle yy-strain (in time) Size=np
	Array<Array<double> > Out_p_fx_t;  // Out external force x (in time) Size=np
	Array<Array<double> > Out_p_fy_t;  // Out external force y (in time) Size=np
	Array<Array<double> > Out_p_ux_t;  // Out displacement x (in time) Size=np
	Array<Array<double> > Out_p_uy_t;  // Out displacement y (in time) Size=np
};

/** (Re)set output arrays. */
void ReSetOutputArrays (Problem2D & Prob)
{
	// Erase output arrays
	Prob.Out_t      .Resize(0);
	Prob.Out_Wx_t   .Resize(0);
	Prob.Out_sE_t   .Resize(0);
	Prob.Out_kE_t   .Resize(0);
	Prob.Out_tE_t   .Resize(0);
	Prob.Out_p_vx_t .Resize(Prob.P->nPoints());
	Prob.Out_p_vy_t .Resize(Prob.P->nPoints());
	Prob.Out_p_cvx_t.Resize(Prob.P->nPoints());
	Prob.Out_p_sxx_t.Resize(Prob.P->nPoints());
	Prob.Out_p_syy_t.Resize(Prob.P->nPoints());
	Prob.Out_p_exx_t.Resize(Prob.P->nPoints());
	Prob.Out_p_eyy_t.Resize(Prob.P->nPoints());
	Prob.Out_p_fx_t .Resize(Prob.P->nPoints());
	Prob.Out_p_fy_t .Resize(Prob.P->nPoints());
	Prob.Out_p_ux_t .Resize(Prob.P->nPoints());
	Prob.Out_p_uy_t .Resize(Prob.P->nPoints());
}

Fl_Menu_Item Problems[] =
{
	{"1. Linear vibration"                , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"2. Indentation (WCCM8)"             , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"3. Bouncing disks"                  , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"4. 2D vibration - small"            , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"5. Particle separation problem"     , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"6. Explosion"                       , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"7. Cantilever"                      , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"8. Column (buckling)"               , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"9. Cantilever (plastic)"            , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"10.Linear vibration (continuum)"    , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"11.Cylindrical Hole"                , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"12.Spoil pile"                      , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"13.Footing penetration"             , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"14.Comet impact"                    , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"15.Rotating disk"                   , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"16. SBB-02a: gravity"               , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{"17. SBB-02b: gravity"               , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0}, 
	{0,0,0,0,0,0,0,0,0}
};
size_t DefaultProblem = 3;

/** Set problem. */
void SetProblem (Problem2D & Prob)
{
	// Deallocate previous grid and matpoints
	if (Prob.G!=NULL) delete Prob.G;
	if (Prob.P!=NULL) delete Prob.P;

	// Set problem data
	switch (Prob.Id)
	{
		case 1:
		{
			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/0, /*yMin*/0, /*nRow*/4, /*nCol*/4, /*Lx*/1.0, /*Ly*/1.0, &CbIsFixed1);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom1, &CbDensity1, &CbAllocMdl1, &CbIniVeloc1, &CbHasTraction1, &CbHasAppDisp1);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom1;
			Prob.pB      = &CbB1;
			Prob.pM      = &CbLdM1;
			Prob.pVeloc  = &CbVel1;
			Prob.pStress = &CbStress1;

			// Set constants
			Prob.t     = 0.0;          // current time
			Prob.tf    = 5.0;          // final time
			Prob.Dt    = 0.015;// 0.1/(2.0*PI); // time increment
			Prob.Dtout = 0.015;        // time increment for output
			Prob.Outp  = 0;            // point/particle number to output data
			Prob.B     = 0.0;          // current b-body mass
			Prob.M     = 0.0;          // multiplier for external forces

			break;
		}
		case 2:
		{
			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-0.25, /*yMin*/-0.25, /*nRow*/7, /*nCol*/7, /*Lx*/0.25, /*Ly*/0.25, &CbIsFixed2);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom2, &CbDensity2, &CbAllocMdl2, &CbIniVeloc2, &CbHasTraction2, &CbHasAppDisp2);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom2;
			Prob.pB      = &CbB2;
			Prob.pM      = &CbLdM2;
			Prob.pVeloc  = &CbVel2;
			Prob.pStress = &CbStress2;

			// Set constants
			Prob.t     = 0.0;    // current time
			Prob.tf    = 8.0;    // final time
			Prob.Dt    = 0.01;   // time increment
			Prob.Dtout = 0.1;    // time increment for output
			Prob.Outp  = 0;      // point/particle number to output data
			Prob.B     = 0.0;    // current b-body mass
			Prob.M     = 0.0;    // multiplier for external forces

			break;
		}
		case 3:
		{
			// Allocate grid
			//Prob.G = new Grid2D (/*xMin*/-0.10, /*yMin*/-0.10, /*nRow*/24, /*nCol*/24, /*Lx*/0.05, /*Ly*/0.05, &CbIsFixed3);
			Prob.G = new Grid2D (/*xMin*/-0.10, /*yMin*/-0.10, /*nRow*/18, /*nCol*/18, /*Lx*/0.05, /*Ly*/0.05, &CbIsFixed3);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom3, &CbDensity3, &CbAllocMdl3, &CbIniVeloc3, &CbHasTraction3, &CbHasAppDisp3);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom3;
			Prob.pB      = &CbB3;
			Prob.pM      = &CbLdM3;
			Prob.pVeloc  = &CbVel3;
			Prob.pStress = &CbStress3;

			// Set constants
			Prob.t     = 0.0;   // current time
			Prob.tf    = 10.0;  // final time
			Prob.Dt    = 0.001; // time increment
			Prob.Dtout = 0.1;   // time increment for output
			Prob.Outp  = 0;     // point/particle number to output data
			Prob.B     = 0.0;   // current b-body mass
			Prob.M     = 0.0;   // multiplier for external forces

			/*
			// Boundary conditions
			for (int i=0; i<Prob.G->nNodes(); ++i)
			{
				if (Prob.G->N(i)(0)<0.10) Prob.G->SetFixed(i,FIX_XY); // left
				if (Prob.G->N(i)(0)>0.89) Prob.G->SetFixed(i,FIX_XY); // right
				if (Prob.G->N(i)(1)<0.10) Prob.G->SetFixed(i,FIX_XY); // bottom
				if (Prob.G->N(i)(1)>0.89) Prob.G->SetFixed(i,FIX_XY); // top
			}
			*/

			break;
		}
		case 4:
		{
			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-0.25, /*yMin*/-0.25, /*nRow*/7, /*nCol*/7, /*Lx*/0.25, /*Ly*/0.25, &CbIsFixed4);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom4, &CbDensity4, &CbAllocMdl4, &CbIniVeloc4, &CbHasTraction4, &CbHasAppDisp4);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom4;
			Prob.pB      = &CbB4;
			Prob.pM      = &CbLdM4;
			Prob.pVeloc  = &CbVel4;
			Prob.pStress = &CbStress4;

			// Set constants
			Prob.t     = 0.0;    // current time
			Prob.tf    = 4.0;    // final time
			Prob.Dt    = 0.0005; // time increment
			Prob.Dtout = 0.01;   // time increment for output
			Prob.Outp  = 0;      // point/particle number to output data
			Prob.B     = 0.0;    // current b-body mass
			Prob.M     = 0.0;    // multiplier for external forces

			break;
		}
		case 5:
		{
			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-2.0, /*yMin*/-2.0, /*nRow*/23, /*nCol*/123, /*Lx*/2.0, /*Ly*/2.0, &CbIsFixed5);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom5, &CbDensity5, &CbAllocMdl5, &CbIniVeloc5, &CbHasTraction5, &CbHasAppDisp5);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom5;
			Prob.pB      = &CbB5;
			Prob.pM      = &CbLdM5;
			Prob.pVeloc  = &CbVel5;
			Prob.pStress = &CbStress5;

			// Set constants
			Prob.t     = 0.0;  // current time (nano-seconds)
			Prob.tf    = 60.0; // final time (nano-seconds)
			Prob.Dt    = 0.02; // time increment (nano-seconds)
			Prob.Dtout = 0.1;  // time increment for output (nano-seconds)
			Prob.Outp  = 0;    // point/particle number to output data
			Prob.B     = 0.0;  // current b-body mass
			Prob.M     = 0.0;  // multiplier for external forces

			break;
		}
		case 6:
		{
			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-1.0, /*yMin*/-1.0, /*nRow*/23, /*nCol*/23, /*Lx*/1.0, /*Ly*/1.0, &CbIsFixed6);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom6, &CbDensity6, &CbAllocMdl6, &CbIniVeloc6, &CbHasTraction6, &CbHasAppDisp6);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom6;
			Prob.pB      = &CbB6;
			Prob.pM      = &CbLdM6;
			Prob.pVeloc  = &CbVel6;
			Prob.pStress = &CbStress6;

			// Set constants
			Prob.t     = 0.0;    // current time
			Prob.tf    = 1.0;    // final time
			Prob.Dt    = 0.0001; // time increment
			Prob.Dtout = 0.01;   // time increment for output
			Prob.Outp  = 0;      // point/particle number to output data
			Prob.B     = 0.0;    // current b-body mass
			Prob.M     = 0.0;    // multiplier for external forces

			break;
		}
		case 7:
		{
			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-1.0, /*yMin*/-1.0, /*nRow*/16, /*nCol*/23, /*Lx*/1.0, /*Ly*/1.0, &CbIsFixed7);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom7, &CbDensity7, &CbAllocMdl7, &CbIniVeloc7, &CbHasTraction7, &CbHasAppDisp7);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom7;
			Prob.pB      = &CbB7;
			Prob.pM      = &CbLdM7;
			Prob.pVeloc  = &CbVel7;
			Prob.pStress = &CbStress7;

			// Set constants
			Prob.t     = 0.0;   // current time
			Prob.tf    = 10.0;  // final time
			Prob.Dt    = 0.001; // time increment
			Prob.Dtout = 0.1;   // time increment for output
			Prob.Outp  = 0;     // point/particle number to output data
			Prob.B     = 0.0;   // current b-body mass
			Prob.M     = 0.0;   // multiplier for external forces

			break;
		}
		case 8:
		{
			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-3.0, /*yMin*/-1.0, /*nRow*/16, /*nCol*/20, /*Lx*/1.0, /*Ly*/1.0, &CbIsFixed8);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom8, &CbDensity8, &CbAllocMdl8, &CbIniVeloc8, &CbHasTraction8, &CbHasAppDisp8);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom8;
			Prob.pB      = &CbB8;
			Prob.pM      = &CbLdM8;
			Prob.pVeloc  = &CbVel8;
			Prob.pStress = &CbStress8;

			// Set constants
			Prob.t     = 0.0;   // current time
			Prob.tf    = 20.0;  // final time
			Prob.Dt    = 0.001; // time increment
			Prob.Dtout = 0.1;   // time increment for output
			Prob.Outp  = 0;     // point/particle number to output data
			Prob.B     = 0.0;   // current b-body mass
			Prob.M     = 0.0;   // multiplier for external forces

			break;
		}
		case 9:
		{
			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-1.0, /*yMin*/-1.0, /*nRow*/16, /*nCol*/23, /*Lx*/1.0, /*Ly*/1.0, &CbIsFixed7);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom7, &CbDensity7, &CbAllocMdl7b, &CbIniVeloc7, &CbHasTraction7, &CbHasAppDisp7);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom7;
			Prob.pB      = &CbB7;
			Prob.pM      = &CbLdM7;
			Prob.pVeloc  = &CbVel7;
			Prob.pStress = &CbStress7;

			// Set constants
			Prob.t     = 0.0;   // current time
			Prob.tf    = 10.0;  // final time
			Prob.Dt    = 0.001; // time increment
			Prob.Dtout = 0.1;   // time increment for output
			Prob.Outp  = 0;     // point/particle number to output data
			Prob.B     = 0.0;   // current b-body mass
			Prob.M     = 0.0;   // multiplier for external forces

			break;
		}
		case 10:
		{
			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-1.0, /*yMin*/0, /*nRow*/4, /*nCol*/28, /*Lx*/1.0, /*Ly*/1.0, &CbIsFixed10);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom10, &CbDensity10, &CbAllocMdl10, &CbIniVeloc10, &CbHasTraction10, &CbHasAppDisp10);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom10;
			Prob.pB      = &CbB10;
			Prob.pM      = &CbLdM10;
			Prob.pVeloc  = &CbVel10;
			Prob.pStress = &CbStress10;

			// Set constants
			Prob.t     = 0.0;          // current time
			Prob.tf    = 100.0;        // final time
			Prob.Dt    = 0.01;         // time increment
			Prob.Dtout = 1.0;          // time increment for output
			Prob.Outp  = 12;           // point/particle number to output data
			Prob.B     = 0.0;          // current b-body mass
			Prob.M     = 0.0;          // multiplier for external forces

			break;
		}
		case 11:
		{
			/* Midas: GeoStrucAnal_01.pdf */

			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-1.0, /*yMin*/-1.0, /*nRow*/13, /*nCol*/13, /*Lx*/1.0, /*Ly*/1.0, &CbIsFixed11);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom11, &CbDensity11, &CbAllocMdl11, &CbIniVeloc11, &CbHasTraction11, &CbHasAppDisp11);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom11;
			Prob.pB      = &CbB11;
			Prob.pM      = &CbLdM11;
			Prob.pVeloc  = &CbVel11;
			Prob.pStress = &CbStress11;

			// Set constants
			Prob.t     = 0.0;          // current time
			Prob.tf    = 100.0;        // final time
			Prob.Dt    = 0.01;         // time increment
			Prob.Dtout = 1.0;          // time increment for output
			Prob.Outp  = 0;            // point/particle number to output data
			Prob.B     = 0.0;          // current b-body mass
			Prob.M     = 0.0;          // multiplier for external forces

			break;
		}
		case 12:
		{
			double H     = 10.0;   // height
			double L     = 40.0;   // length
			int    ndivx = 40;
			int    ndivy = 10;
			double Lx    = L/ndivx;
			double Ly    = H/ndivy;

			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-Lx, /*yMin*/-Ly, /*nRow*/ndivy+3, /*nCol*/ndivx+3, /*Lx*/Lx, /*Ly*/Ly, &CbIsFixed12);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom12, &CbDensity12, &CbAllocMdl12, &CbIniVeloc12, &CbHasTraction12, &CbHasAppDisp12);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom12;
			Prob.pB      = &CbB12;
			Prob.pM      = &CbLdM12;
			Prob.pVeloc  = &CbVel12;
			Prob.pStress = &CbStress12;

			// Set constants
			Prob.t     = 0.0;          // current time
			Prob.tf    = 30.0;         // final time
			Prob.Dt    = 0.01;         // time increment
			Prob.Dtout = 0.1;          // time increment for output
			Prob.Outp  = 0;            // point/particle number to output data
			Prob.B     = 0.0;          // current b-body mass
			Prob.M     = 0.0;          // multiplier for external forces

			break;
		}
		case 13:
		{
			double L     = 15.0;   // length
			double H     = 10.0;   // height
			int    ndivx = 60;
			int    ndivy = 40;
			double Lx    = L/ndivx;
			double Ly    = H/ndivy;

			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-Lx, /*yMin*/-Ly, /*nRow*/ndivy+3+5, /*nCol*/ndivx+3, /*Lx*/Lx, /*Ly*/Ly, &CbIsFixed13);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom13, &CbDensity13, &CbAllocMdl13, &CbIniVeloc13, &CbHasTraction13, &CbHasAppDisp13);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom13;
			Prob.pB      = &CbB13;
			Prob.pM      = &CbLdM13;
			Prob.pVeloc  = &CbVel13;
			Prob.pStress = &CbStress13;

			// Set constants
			Prob.t     = 0.0;          // current time
			Prob.tf    = 30.0;         // final time
			Prob.Dt    = 0.01;         // time increment
			Prob.Dtout = 0.1;          // time increment for output
			Prob.Outp  = 0;            // point/particle number to output data
			Prob.B     = 0.0;          // current b-body mass
			Prob.M     = 0.0;          // multiplier for external forces

			break;
		}
		case 14:
		{
			double L     = 30.0;   // length
			double H     = 15.0;   // height
			int    ndivx = 30;
			int    ndivy = 15;
			double Lx    = L/ndivx;
			double Ly    = H/ndivy;

			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-Lx, /*yMin*/-Ly, /*nRow*/ndivy+3, /*nCol*/ndivx+3, /*Lx*/Lx, /*Ly*/Ly, &CbIsFixed14);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom14, &CbDensity14, &CbAllocMdl14, &CbIniVeloc14, &CbHasTraction14, &CbHasAppDisp14);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom14;
			Prob.pB      = &CbB14;
			Prob.pM      = &CbLdM14;
			Prob.pVeloc  = &CbVel14;
			Prob.pStress = &CbStress14;

			// Set constants
			Prob.t     = 0.0;          // current time
			Prob.tf    = 0.7;          // final time
			Prob.Dt    = 0.00001;      // time increment
			//Prob.tf    = 30.0;         // final time
			//Prob.Dt    = 0.001;        // time increment
			Prob.Dtout = 0.01;         // time increment for output
			Prob.Outp  = 0;            // point/particle number to output data
			Prob.B     = 0.0;          // current b-body mass
			Prob.M     = 0.0;          // multiplier for external forces

			break;
		}
		case 15:
		{
			double L     = 10.0;   // length
			double H     = 10.0;   // height
			int    ndivx = 20;
			int    ndivy = 20;
			double Lx    = L/ndivx;
			double Ly    = H/ndivy;

			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-Lx, /*yMin*/-Ly, /*nRow*/ndivy+3, /*nCol*/ndivx+3, /*Lx*/Lx, /*Ly*/Ly, &CbIsFixed15);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom15, &CbDensity15, &CbAllocMdl15, &CbIniVeloc15, &CbHasTraction15, &CbHasAppDisp15);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom15;
			Prob.pB      = &CbB15;
			Prob.pM      = &CbLdM15;
			Prob.pVeloc  = &CbVel15;
			Prob.pStress = &CbStress15;

			// Set constants
			Prob.t     = 0.0;          // current time
			Prob.tf    = 5.0;          // final time
			Prob.Dt    = 0.0005;       // time increment
			Prob.Dtout = 0.01;         // time increment for output
			Prob.Outp  = 0;            // point/particle number to output data
			Prob.B     = 0.0;          // current b-body mass
			Prob.M     = 0.0;          // multiplier for external forces

			break;
		}
		case 16:
		{
			double L     = 2.0;   // length
			double H     = 10.0;   // height
			int    ndivx = 4;
			int    ndivy = 20;
			double Lx    = L/ndivx;
			double Ly    = H/ndivy;

			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-Lx, /*yMin*/-Ly, /*nRow*/ndivy+3, /*nCol*/ndivx+3, /*Lx*/Lx, /*Ly*/Ly, &CbIsFixed16);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom16, &CbDensity16, &CbAllocMdl16, &CbIniVeloc16, &CbHasTraction16, &CbHasAppDisp16);

			// Set callbacks
			Prob.pGeom   = &CbIsPointInGeom16;
			Prob.pB      = &CbB16;
			Prob.pM      = &CbLdM16;
			Prob.pVeloc  = &CbVel16;
			Prob.pStress = &CbStress16;

			// Set constants
			Prob.t     = 0.0;          // current time
			Prob.tf    = 0.25;         // final time
			Prob.Dt    = 0.00006;      // time increment
			Prob.Dtout = 0.001;        // time increment for output
			Prob.Outp  = 0;            // point/particle number to output data
			Prob.B     = 0.0;          // current b-body mass
			Prob.M     = 0.0;          // multiplier for external forces

			break;
		}
		case 17:
		{
			double L     = 2.0;   // length
			double H     = 10.0;   // height
			int    ndivx = 4;
			int    ndivy = 20;
			double Lx    = L/ndivx;
			double Ly    = H/ndivy;

			// Allocate grid
			Prob.G = new Grid2D (/*xMin*/-Lx, /*yMin*/-Ly, /*nRow*/ndivy+3, /*nCol*/ndivx+3, /*Lx*/Lx, /*Ly*/Ly, &CbIsFixed17);

			// Allocate material points
			Prob.P = new MPoints2D (Prob.NPcell, Prob.G, &CbIsPointInGeom17, &CbDensity17, &CbAllocMdl17, &CbIniVeloc17, &CbHasTraction17, &CbHasAppDisp17);

			// Set callbacks
            Prob.NPcell  = 9;
			Prob.pGeom   = &CbIsPointInGeom17;
			Prob.pB      = &CbB17;
			Prob.pM      = &CbLdM17;
			Prob.pVeloc  = &CbVel17;
			Prob.pStress = &CbStress17;

			// Set constants
			Prob.t     = 0.0;          // current time
			Prob.tf    = 0.25;         // final time
			Prob.Dt    = 0.00006;      // time increment
			Prob.Dtout = 0.001;        // time increment for output
			Prob.Outp  = 0;            // point/particle number to output data
			Prob.B     = 0.0;          // current b-body mass
			Prob.M     = 0.0;          // multiplier for external forces

			break;
		}
		default: { throw new Fatal("Prob # %d is not available.",Prob.Id); }
	}
	ReSetOutputArrays (Prob);
}

}; // namespace MPM

#endif // MPM_PROBLEMS2D_H
