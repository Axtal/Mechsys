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

#ifndef MECHSYS_TUNNEL_H
#define MECHSYS_TUNNEL_H

// MechSys
#include "util/array.h"
#include "util/util.h"
#include "util/exception.h"

class ElasticCtes
{
public:
	ElasticCtes() { Set_E_nu(1.0,0.25); }

	// Set methods
	ElasticCtes & Set_K_G  (double K, double G)  { _K = K;  _G  = G;   _E  = 9.0*_K*_G/(3.0*_K+_G);   _nu = (3.0*_K-2.0*_G)/(6.0*_K+2.0*_G);  return (*this); }
	ElasticCtes & Set_E_nu (double E, double nu) { _E = E;  _nu = nu;  _G  = _E/(2.0*(1.0+_nu));      _K  = _E/(3.0*(1.0-2.0*_nu));           return (*this); }
	ElasticCtes & Set_E_G  (double E, double G)  { _E = E;  _G  = G;   _nu = _E/(2.0*_G)-1.0;         _K  = _E/(3.0*(1.0-2.0*_nu));           return (*this); }
	ElasticCtes & Set_G_nu (double G, double nu) { _G = G;  _nu = nu;  _E  = 2.0*G*(1.0+_nu);         _K  = _E/(3.0*(1.0-2.0*_nu));           return (*this); }

	// Access
	double E  () const { return _E;  }
	double nu () const { return _nu; }
	double K  () const { return _K;  }
	double G  () const { return _G;  }

private:
	double _E;  // Young modulus
	double _nu; // Poission coefficient
	double _K;  // Bulk modulus
	double _G;  // Shear modulus

}; // class ElasticCtes


class Kirsch
{
public:
	Kirsch(double a, double px, double py, ElasticCtes const & C) : _a(a), _px(px), _py(py), _c(C) { _pm=(py+px)/2.0; _pd=(py-px)/2.0; }
	double Calc_srr(double r, double t)
	{
		if (r<_a) throw new Fatal("Kirsch::Srr: r==%g must be greater than/equal to a==%g",r,_a);
		double m  = _a*_a/(r*r);
		double mm = m*m;
		return _pm*(1.0-m)-_pd*(1.0-4.0*m+3.0*mm)*cos(2.0*t);
	}
	double Calc_stt(double r, double t)
	{
		if (r<_a) throw new Fatal("Kirsch::Stt: r==%g must be greater than/equal to a==%g",r,_a);
		double m  = _a*_a/(r*r);
		double mm = m*m;
		return _pm*(1.0+m)+_pd*(1.0+3.0*mm)*cos(2.0*t);
	}
	double Calc_srt(double r, double t)
	{
		if (r<_a) throw new Fatal("Kirsch::Srt: r==%g must be greater than/equal to a==%g",r,_a);
		double m  = _a*_a/(r*r);
		double mm = m*m;
		return _pd*(1.0+2.0*m-3.0*mm)*sin(2.0*t);
	}
	double Calc_sxx(double r, double t)
	{
		double srr = Calc_srr(r,t);
		double stt = Calc_stt(r,t);
		double srt = Calc_srt(r,t);
		double c   = cos(t);
		double s   = sin(t);
		return s*s*stt + c*c*srr - 2.0*c*s*srt;
	}
	double Calc_syy(double r, double t)
	{
		double srr = Calc_srr(r,t);
		double stt = Calc_stt(r,t);
		double srt = Calc_srt(r,t);
		double c   = cos(t);
		double s   = sin(t);
		return c*c*stt + s*s*srr + 2.0*c*s*srt;
	}
	double Calc_sxy(double r, double t)
	{
		double srr = Calc_srr(r,t);
		double stt = Calc_stt(r,t);
		double srt = Calc_srt(r,t);
		double c   = cos(t);
		double s   = sin(t);
		return -s*s*srt + c*c*srt - c*s*stt + c*s*srr;
	}
	double Calc_ur(double r, double t)
	{
		//return _a*(1.0-_c.nu()*_c.nu())*( _px+_py+2.0*(_py-_px)*cos(2.0*t) )/_c.E();
		//return (_a/_c.E())*(_px+_py-_c.nu()*_px+2.0*(1.0-_c.nu()*_c.nu())*cos(2.0*t));
		double k = _a*_a/r;
		double m = _a*_a/(r*r);
		return (0.5*_pm/_c.G())*k - (0.5*_pd/_c.G())*k*(4.0*(1.0-_c.nu())-m)*cos(2.0*t);
	}
	void ValsAlongDist(double t, Array<double> & R, Array<double> & Srr, Array<double> & Stt, Array<double> & Srt, Array<double> & Sxx, Array<double> & Syy, Array<double> & Sxy, Array<double> & ur, double Max_r=-1, double NPoints=100)
	{
		if (Max_r<0) Max_r = _a+10.0*_a;
		double dr = (Max_r-_a)/(NPoints+1.0);
		for (size_t i=0; i<NPoints; ++i)
		{
			double r = _a+i*dr;
			R.Push   (r);
			Srr.Push (Calc_srr(r,t));
			Stt.Push (Calc_stt(r,t));
			Srt.Push (Calc_srt(r,t));
			Sxx.Push (Calc_sxx(r,t));
			Syy.Push (Calc_syy(r,t));
			Sxy.Push (Calc_sxy(r,t));
			ur .Push (Calc_ur (r,t));
		}
	}
	void ValsAlongHole(Array<double> & T, Array<double> & Srr, Array<double> & Stt, Array<double> & Srt, Array<double> & Sxx, Array<double> & Syy, Array<double> & Sxy, Array<double> & ur, double NDiv=90)
	{
		double dt = (Util::PI/180.0)*90.0/NDiv;
		for (size_t i=0; i<NDiv+1; ++i)
		{
			double t = i*dt;
			T  .Push (t*180.0/Util::PI);
			Srr.Push (Calc_srr(_a,t));
			Stt.Push (Calc_stt(_a,t));
			Srt.Push (Calc_srt(_a,t));
			Sxx.Push (Calc_sxx(_a,t));
			Syy.Push (Calc_syy(_a,t));
			Sxy.Push (Calc_sxy(_a,t));
			ur .Push (Calc_ur (_a,t));
		}
	}

private:
	double              _a;  // Radius
	double              _px;
	double              _py;
	ElasticCtes const & _c; // Parameters
	double              _pm; // Mean pressure
	double              _pd; // Deviatoric pressure
}; // class Kirsch


class Tunnel
{
public:
	Tunnel() : _a(1.0), _zc(1.0) {}

	Tunnel & SetRadius(double a ) { _a  = a;  return (*this); }
	Tunnel & SetDepth (double zc) { _zc = zc; return (*this); }
private:
	double _a;  // radius
	double _zc; // depth
}; // class Tunnel

#endif // MECHSYS_TUNNEL_H
