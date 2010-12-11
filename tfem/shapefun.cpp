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

// GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

// MechSys
#include <mechsys/fem/fem.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
using Util::_4_2;
using Util::_8s;
using Util::_2;

struct Prms
{
    FEM::GeomElem * ele;   // element
    int             i;     // shape function of node i
    int             j;     // node at which shape fun is returned
    bool            r_cte; // r-direction constant ?
    Array<double>   R;     // natural coordinates of nodes
    Array<double>   S;     // natural coordinates of nodes
};

double CalcShape (double x, void * prms)
{
    Prms * p = static_cast<Prms*>(prms);
    if (p->r_cte) p->ele->Shape (p->R[p->i], x, 0);
    else          p->ele->Shape (x, p->S[p->i], 0);
    return p->ele->N(p->j);
}

int main(int argc, char **argv) try
{
    Array<String> names(6);
    names = "Lin2", "Tri3", "Tri6", "Tri15", "Quad4", "Quad8";

    bool has_error = false;
    for (size_t iname=0; iname<names.Size(); ++iname)
    {
        cout << "\n=========================  " << names[iname] << "  ===========================" << endl;

        // parameters
        Prms p;

        // natural coordinates of nodes
        size_t ndim = 2;
        if (names[iname]=="Lin2")
        {
            p.R.Resize (2);
            p.S.Resize (2);
            //       0    1
            p.R = -1.0, 1.0;
            p.S =  0.0, 0.0;
        }
        else if (names[iname]=="Tri3")
        {
            p.R.Resize (3);
            p.S.Resize (3);
            //      0    1    2
            p.R = 0.0, 1.0, 0.0;
            p.S = 0.0, 0.0, 1.0;
        }
        else if (names[iname]=="Tri6")
        {
            p.R.Resize (6);
            p.S.Resize (6);
            //      0    1    2    3    4    5
            p.R = 0.0, 1.0, 0.0, 0.5, 0.5, 0.0;
            p.S = 0.0, 0.0, 1.0, 0.0, 0.5, 0.5;
        }
        else if (names[iname]=="Tri15")
        {
            p.R.Resize (15);
            p.S.Resize (15);
            //      0    1    2    3    4    5     6     7     8     9   10    11     12   13     14
            p.R = 0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 0.25, 0.75, 0.75, 0.25, 0.0,  0.0,  0.25, 0.5,  0.25;
            p.S = 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.0,  0.0,  0.25, 0.75, 0.75, 0.25, 0.25, 0.25, 0.50;
        }
        else if (names[iname]=="Quad4")
        {
            p.R.Resize (4);
            p.S.Resize (4);
            //       0     1     2     3
            p.R = -1.0,  1.0,  1.0, -1.0;
            p.S = -1.0, -1.0,  1.0,  1.0;
        }
        else if (names[iname]=="Quad8")
        {
            p.R.Resize (8);
            p.S.Resize (8);
            //       0     1     2     3     4     5     6     7
            p.R = -1.0,  1.0,  1.0, -1.0,  0.0,  1.0,  0.0, -1.0;
            p.S = -1.0, -1.0,  1.0,  1.0, -1.0,  0.0,  1.0,  0.0;
        }

        // element
        p.ele = FEM::AllocGeomElem (names[iname].CStr(), ndim);

        // shape
        double error_N = 0.0;
        double tol_N   = 1.0e-16;
        for (size_t i=0; i<p.ele->NN; ++i)
        {
            p.ele->Shape (p.R[i], p.S[i], 0);
            cout << "N"<<_2<<i<<"("<<_4_2<<p.R[i]<<","<<_4_2<<p.S[i]<<") = " << PrintVector(p.ele->N,"%1g");
            for (size_t j=0; j<p.ele->NN; ++j)
            {
                if (i==j) error_N += fabs(p.ele->N(j)-1.0);
                else      error_N += fabs(p.ele->N(j));
            }
        }
        cout << "error (N)    = " << (error_N>tol_N ? "[1;31m" : "[1;32m") << _8s<<error_N << "[0m\n";

        // derivatives
        double error_dNdr = 0.0;
        double error_dNds = 0.0;
        double tol_dNdr   = 1.0e-6;
        double tol_dNds   = 1.0e-6;
        double tol_dNdR   = 1.0e-5;
        double tol_dNdS   = 1.0e-5;
        double err_dNdr, err_dNds;
        double dNdr, dNds, dNdr_num, dNds_num, abserr;
        for (size_t i=0; i<p.ele->NN; ++i)
        {
            p.ele->Derivs (p.R[i], p.S[i], 0);
            p.i = i;
            for (size_t j=0; j<p.ele->NN; ++j)
            {
                gsl_function fun;
                fun.function = &CalcShape;
                fun.params   = &p;
                p.j          = j;

                p.r_cte = false; // varying r (s-constant)
                gsl_deriv_central (&fun, /*around r*/p.R[i], /*h*/1.0e-8, &dNdr_num, &abserr);

                if (abserr>tol_dNdr) cout << "abserr =" << _8s<<abserr << ", tol_dNdr =" << _8s<<tol_dNdr << endl;

                p.r_cte = true; // varying s (r-constant)
                gsl_deriv_central (&fun, /*around s*/p.S[i], /*h*/1.0e-8, &dNds_num, &abserr);

                if (abserr>tol_dNds) cout << "abserr =" << _8s<<abserr << ", tol_dNds =" << _8s<<tol_dNds << endl;

                dNdr        = p.ele->dNdR(0,j);
                dNds        = p.ele->dNdR(1,j);
                err_dNdr    = fabs(dNdr_num - dNdr);
                err_dNds    = fabs(dNds_num - dNds);
                error_dNdr += err_dNdr;
                error_dNds += err_dNds;

                if (err_dNdr>tol_dNdr) cout<<"dN("<<_2<<i<<")dr("<<_2<<j<<") ="<<_8s<<dNdr<<", dNdr_num ="<<_8s<<dNdr_num<<", error =[1;31m"<<_8s<<err_dNdr<<"[0m\n";
                if (err_dNds>tol_dNds) cout<<"dN("<<_2<<i<<")ds("<<_2<<j<<") ="<<_8s<<dNds<<", dNds_num ="<<_8s<<dNds_num<<", error =[1;31m"<<_8s<<err_dNds<<"[0m\n";
            }
        }
        cout << "error (dNdr) = " << (error_dNdr>tol_dNdR ? "[1;31m" : "[1;32m") << _8s<<error_dNdr << "[0m\n";
        cout << "error (dNds) = " << (error_dNdr>tol_dNdS ? "[1;31m" : "[1;32m") << _8s<<error_dNds << "[0m\n";

        if (error_N>tol_N)       has_error = true;
        if (error_dNdr>tol_dNdR) has_error = true;
        if (error_dNds>tol_dNdS) has_error = true;

        delete p.ele;
    }

    return has_error;
}
MECHSYS_CATCH
