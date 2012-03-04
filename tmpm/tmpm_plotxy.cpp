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

/* tplotxy - Copyright (C) 2007 Dorival de Moraes Pedroso */

// STL
#include <iostream>
#include <cmath>

// FLTK
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>

// Local
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mpm/plotxy.h>

using std::cout;
using std::endl;
using MPM::PlotXY;

/////////////////////////////////////////////////////////////////////////////////////////////// MainWindow /////

// MainWindow
class MainWindow : public Fl_Double_Window
{
public:
	// Constructor
	MainWindow(int width=1024, int height=824)
		: Fl_Double_Window (width,height,"Test PlotXY - (c) 2007 Dorival M. Pedroso") // {{{
	{
		// Add widgets to window
		begin();
			// Widgets
			_wp0 = new PlotXY (/*xmin*/50,  /*ymin*/50,  /*width*/400, /*height*/200);
			_wp1 = new PlotXY (/*xmin*/260, /*ymin*/260, /*width*/400, /*height*/400);
		end();

		// Set curve # 0
		for (int i=0; i<=20; ++i)
		{
			// Curve # 0
			double x = static_cast<double>(i)/10;
			double y = 2.0*x;
			_X0.Push(x);
			_Y0.Push(y);

			// Curve # 1
			x = static_cast<double>(i)/10;
			y = x*x;
			_X1.Push(x);
			_Y1.Push(y);

			// Curve # 2
			x = static_cast<double>(i);
			y = 20.0*x;
			_X2.Push(x);
			_Y2.Push(y);

			// Curve # 3
			x = static_cast<double>(i);
			y = x*x;
			_X3.Push(x);
			_Y3.Push(y);
		}

		// Set plotXY # 0
		_wp0->EqScales (false).RecalcSF(true);
		_wp0->AddCurve (&_X0, &_Y0, "0"); _wp0->SetCurve(0).Clr=FL_RED;  _wp0->SetCurve(0).Lwd=2; _wp0->SetCurve(0).Pch=16;
		_wp0->AddCurve (&_X1, &_Y1, "1"); _wp0->SetCurve(1).Clr=FL_BLUE; _wp0->SetCurve(1).Lwd=2; _wp0->SetCurve(1).Typ=MPM::CT_BOTH;

		// Set plotXY # 1
		_wp1->EqScales (false).RecalcSF(true);
		_wp1->AddCurve (&_X2, &_Y2, "2"); _wp1->SetCurve(0).Clr=FL_GREEN;   _wp1->SetCurve(0).Lwd=1; _wp1->SetCurve(0).Typ=MPM::CT_LINES;
		_wp1->AddCurve (&_X3, &_Y3, "3"); _wp1->SetCurve(1).Clr=FL_MAGENTA; _wp1->SetCurve(1).Lwd=1; _wp1->SetCurve(1).Typ=MPM::CT_BOTH;
		_wp1->AddCurve (&_X3, &_Y3, "3"); _wp1->SetCurve(2).Clr=FL_BLUE;    _wp1->SetCurve(2).Lwd=1; _wp1->SetCurve(2).Typ=MPM::CT_BOTH;
		_wp1->AddCurve (&_X3, &_Y3, "3"); _wp1->SetCurve(3).Clr=FL_CYAN;    _wp1->SetCurve(3).Lwd=2; _wp1->SetCurve(3).Typ=MPM::CT_BOTH;

		// Set resizable and show window
		resizable (this);
		show      ();

	} // }}}

private:
	// Widgets
	PlotXY * _wp0; 
	PlotXY * _wp1; 

	// Variables
	Array<double> _X0; // Curve # 0
	Array<double> _Y0; // 
	Array<double> _X1; // Curve # 1
	Array<double> _Y1; // 
	Array<double> _X2; // Curve # 2
	Array<double> _Y2; // 
	Array<double> _X3; // Curve # 3
	Array<double> _Y3; // 

}; // class MainWindow

/////////////////////////////////////////////////////////////////////////////////////////////// MainWindow /////

int main(int argc, char **argv) try
{
	MainWindow win;
	Fl::run();
	return 0;
}
catch (Fatal * f) //{{{ 
{
	f->Cout();
	delete f;
	exit (1);
}
catch (char const * m)
{
	std::cout << "Fatal: " << m << std::endl;
	exit (1);
}
catch (...)
{
	std::cout << "Some exception (...) ocurred\n";
} //}}} 

// vim:fdm=marker
