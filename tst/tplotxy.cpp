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
#include <cmath>

// FLTK
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>

// MechSys
#include "util/array.h"
#include "util/fatal.h"
#include "gui/plotxy.h"

using std::cout;
using std::endl;

class MainWindow : public Fl_Double_Window
{
public:
    // Constructor
    MainWindow(int width=1024, int height=824)
        : Fl_Double_Window (width,height,"Testing PlotXY")
    {
        // Add widgets to window
        begin();
            // Widgets
            P0 = new GUI::PlotXY (/*xmin*/50,  /*ymin*/50,  /*width*/400, /*height*/300, "First Graph");
            P1 = new GUI::PlotXY (/*xmin*/100, /*ymin*/400, /*width*/400, /*height*/400, "Second Graph");
        end();

        // Calc data
        const int np = 20;
        X .Resize(np);
        Y0.Resize(np);
        Y1.Resize(np);
        Y2.Resize(np);
        Y3.Resize(np);
        Y4.Resize(np);
        Y5.Resize(np);
        for (int i=0; i<np; ++i)
        {
            X [i] = static_cast<double>(i)/np;
            Y0[i] = 2.0*X[i];
            Y1[i] = X[i]*X[i];
            Y2[i] = X[i];
            Y3[i] = 0.5*X[i];
            Y4[i] = log(1.0+X[i]);
            Y5[i] = sqrt(X[i]);
        }

        // P0
        P0->EqScales (false).RecalcSF(true);
        P0->AddCurve (&X, &Y0, "2*x");  P0->_C[0].Clr=FL_RED;   P0->_C[0].Lwd=2;  P0->_C[0].Pch=16;
        P0->AddCurve (&X, &Y1, "x^2");  P0->_C[1].Clr=FL_BLUE;  P0->_C[1].Lwd=2;  P0->_C[1].Typ=GUI::CT_BOTH;

        // P1
        P1->EqScales (false).RecalcSF(true);
        P1->AddCurve (&X, &Y2, "x");         P1->_C[0].Clr=FL_GREEN;    P1->_C[0].Lwd=1;  P1->_C[0].Typ=GUI::CT_LINES;
        P1->AddCurve (&X, &Y3, "x/2");       P1->_C[1].Clr=FL_MAGENTA;  P1->_C[1].Lwd=1;  P1->_C[1].Typ=GUI::CT_BOTH;
        P1->AddCurve (&X, &Y4, "log(1+x)");  P1->_C[2].Clr=FL_BLUE;     P1->_C[2].Lwd=1;  P1->_C[2].Typ=GUI::CT_BOTH;
        P1->AddCurve (&X, &Y5, "sqrt(x)");   P1->_C[3].Clr=FL_CYAN;     P1->_C[3].Lwd=2;  P1->_C[3].Typ=GUI::CT_BOTH;

        // Set resizable and show window
        resizable (this);
        show      ();
    }

    // Data
    GUI::PlotXY * P0;
    GUI::PlotXY * P1;
    Array<double> X, Y0, Y1, Y2, Y3, Y4, Y5;
};

int main(int argc, char **argv) try
{
    MainWindow win;
    Fl::run();
    return 0;
}
MECHSYS_CATCH
