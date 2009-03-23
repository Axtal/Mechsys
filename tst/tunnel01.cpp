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
#include <fstream>
#include <cmath>

// FLTK
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>

// MechSys
#include "analytical/tunnel.h"
#include "util/util.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "gui/plotxy.h"

using std::cout;
using std::endl;
using Util::_4;
using Util::_8s;
using Util::PI;

class MainWindow : public Fl_Double_Window
{
public:
	// Constructor
	MainWindow(int width=1024, int height=824) : Fl_Double_Window (width,height,"ttunnel")
	{
		// Widget
		begin();
			_plt0 = new PlotXY (/*x*/0,  /*y*/0,    /*width*/600, /*height*/200, "theta=0");
			_plt1 = new PlotXY (/*x*/0,  /*y*/200,  /*width*/600, /*height*/200, "theta=45");
			_plt2 = new PlotXY (/*x*/0,  /*y*/400,  /*width*/600, /*height*/200, "theta=90");
			_plt3 = new PlotXY (/*x*/0,  /*y*/600,  /*width*/600, /*height*/200, "displacement");
		end();

		// Data
		ElasticCtes ec;  ec.Set_E_nu(1,0.25);
		Kirsch      ki(/*a*/1.0, /*px*/1.0, /*py*/1.0, ec);
		ki.ValsAlongDist (/*t*/   0.0, _r_t00, _srr_t00, _stt_t00, _srt_t00, _sxx_t00, _syy_t00, _sxy_t00, _ur_t00);
		ki.ValsAlongDist (/*t*/PI/4.0, _r_t45, _srr_t45, _stt_t45, _srt_t45, _sxx_t45, _syy_t45, _sxy_t45, _ur_t45);
		ki.ValsAlongDist (/*t*/PI/2.0, _r_t90, _srr_t90, _stt_t90, _srt_t90, _sxx_t90, _syy_t90, _sxy_t90, _ur_t90);
		ki.ValsAlongHole (             _t    , _srr    , _stt    , _srt    , _sxx    , _syy    , _sxy    , _ur    );

		// Plot # 0
		_plt0->EqScales (false).RecalcSF(true);
		_plt0->AddCurve (&_r_t00, &_srr_t00, "Srr");  _plt0->SetCurve(0).Clr=FL_RED;    _plt0->SetCurve(0).Lwd=2;  _plt0->SetCurve(0).Typ=PlotXY::CT_LINES;
		_plt0->AddCurve (&_r_t00, &_stt_t00, "Stt");  _plt0->SetCurve(1).Clr=FL_BLUE;   _plt0->SetCurve(1).Lwd=2;  _plt0->SetCurve(1).Typ=PlotXY::CT_LINES;
		_plt0->AddCurve (&_r_t00, &_srt_t00, "Srt");  _plt0->SetCurve(2).Clr=FL_GREEN;  _plt0->SetCurve(2).Lwd=2;  _plt0->SetCurve(2).Typ=PlotXY::CT_LINES;
		_plt0->AddCurve (&_r_t00, &_sxx_t00, "Sxx");  _plt0->SetCurve(3).Clr=FL_RED;    _plt0->SetCurve(3).Lwd=2;  _plt0->SetCurve(3).Typ=PlotXY::CT_POINTS;
		_plt0->AddCurve (&_r_t00, &_syy_t00, "Syy");  _plt0->SetCurve(4).Clr=FL_BLUE;   _plt0->SetCurve(4).Lwd=2;  _plt0->SetCurve(4).Typ=PlotXY::CT_POINTS;
		_plt0->AddCurve (&_r_t00, &_sxy_t00, "Sxy");  _plt0->SetCurve(5).Clr=FL_GREEN;  _plt0->SetCurve(5).Lwd=2;  _plt0->SetCurve(5).Typ=PlotXY::CT_POINTS;

		// Plot # 1
		_plt1->EqScales (false).RecalcSF(true);
		_plt1->AddCurve (&_r_t45, &_srr_t45, "Srr");  _plt1->SetCurve(0).Clr=FL_RED;    _plt1->SetCurve(0).Lwd=2;  _plt1->SetCurve(0).Typ=PlotXY::CT_LINES;
		_plt1->AddCurve (&_r_t45, &_stt_t45, "Stt");  _plt1->SetCurve(1).Clr=FL_BLUE;   _plt1->SetCurve(1).Lwd=2;  _plt1->SetCurve(1).Typ=PlotXY::CT_LINES;
		_plt1->AddCurve (&_r_t45, &_srt_t45, "Srt");  _plt1->SetCurve(2).Clr=FL_GREEN;  _plt1->SetCurve(2).Lwd=2;  _plt1->SetCurve(2).Typ=PlotXY::CT_LINES;
		_plt1->AddCurve (&_r_t45, &_sxx_t45, "Sxx");  _plt1->SetCurve(3).Clr=FL_RED;    _plt1->SetCurve(3).Lwd=2;  _plt1->SetCurve(3).Typ=PlotXY::CT_POINTS;
		_plt1->AddCurve (&_r_t45, &_syy_t45, "Syy");  _plt1->SetCurve(4).Clr=FL_BLUE;   _plt1->SetCurve(4).Lwd=2;  _plt1->SetCurve(4).Typ=PlotXY::CT_POINTS;
		_plt1->AddCurve (&_r_t45, &_sxy_t45, "Sxy");  _plt1->SetCurve(5).Clr=FL_GREEN;  _plt1->SetCurve(5).Lwd=2;  _plt1->SetCurve(5).Typ=PlotXY::CT_POINTS;

		// Plot # 2
		_plt2->EqScales (false).RecalcSF(true);
		_plt2->AddCurve (&_r_t90, &_srr_t90, "Srr");  _plt2->SetCurve(0).Clr=FL_RED;    _plt2->SetCurve(0).Lwd=2;  _plt2->SetCurve(0).Typ=PlotXY::CT_LINES;
		_plt2->AddCurve (&_r_t90, &_stt_t90, "Stt");  _plt2->SetCurve(1).Clr=FL_BLUE;   _plt2->SetCurve(1).Lwd=2;  _plt2->SetCurve(1).Typ=PlotXY::CT_LINES;
		_plt2->AddCurve (&_r_t90, &_srt_t90, "Srt");  _plt2->SetCurve(2).Clr=FL_GREEN;  _plt2->SetCurve(2).Lwd=2;  _plt2->SetCurve(2).Typ=PlotXY::CT_LINES;
		_plt2->AddCurve (&_r_t90, &_sxx_t90, "Sxx");  _plt2->SetCurve(3).Clr=FL_RED;    _plt2->SetCurve(3).Lwd=2;  _plt2->SetCurve(3).Typ=PlotXY::CT_POINTS;
		_plt2->AddCurve (&_r_t90, &_syy_t90, "Syy");  _plt2->SetCurve(4).Clr=FL_BLUE;   _plt2->SetCurve(4).Lwd=2;  _plt2->SetCurve(4).Typ=PlotXY::CT_POINTS;
		_plt2->AddCurve (&_r_t90, &_sxy_t90, "Sxy");  _plt2->SetCurve(5).Clr=FL_GREEN;  _plt2->SetCurve(5).Lwd=2;  _plt2->SetCurve(5).Typ=PlotXY::CT_POINTS;

		// Plot # 3
		_plt3->EqScales (false).RecalcSF(true);
		_plt3->AddCurve (&_t, &_ur, "ur");  _plt3->SetCurve(0).Clr=FL_RED;    _plt3->SetCurve(0).Lwd=2;  _plt3->SetCurve(0).Typ=PlotXY::CT_LINES;

		// Set resizable and show window
		resizable (this);
		show      ();
	}

private:
	PlotXY        * _plt0;
	PlotXY        * _plt1;
	PlotXY        * _plt2;
	PlotXY        * _plt3;
	Array<double>   _r_t00,   _r_t45,   _r_t90,   _t;
	Array<double>   _srr_t00, _srr_t45, _srr_t90, _srr;
	Array<double>   _stt_t00, _stt_t45, _stt_t90, _stt;
	Array<double>   _srt_t00, _srt_t45, _srt_t90, _srt;
	Array<double>   _sxx_t00, _sxx_t45, _sxx_t90, _sxx;
	Array<double>   _syy_t00, _syy_t45, _syy_t90, _syy;
	Array<double>   _sxy_t00, _sxy_t45, _sxy_t90, _sxy;
	Array<double>   _ur_t00,  _ur_t45,  _ur_t90,  _ur;

}; // class MainWindow


int main(int argc, char **argv) try
{
	MainWindow win;
	Fl::run();
	return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
