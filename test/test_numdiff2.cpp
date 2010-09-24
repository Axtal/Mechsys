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

// Std Lib
#include <iostream>

// GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// MechSys
#include <mechsys/numerical/numdiff.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/numstreams.h>

using std::cout;
using std::endl;
using Util::_8s;
using Util::PI;
using Util::SQ2;

class Problem
{
public:
    Problem () { Sig0.change_dim(6);  Sig0 = 1.5,2.,0.5,0.1*SQ2,0.5*SQ2,0.8*SQ2; }
    void CalcSig (double t, Vec_t & Sig) const
    {
        Sig = (1.+t*t)*Sig0;
    }
    double q (double t) const
    {
        Vec_t sig;
        CalcSig (t, sig);
        return Calc_qoct (sig);
    }
    double dqdt (double t) const
    {
        Vec_t sig, S, dsigdt;
        CalcSig (t, sig);
        Dev     (sig, S);
        dsigdt = (2.*t)*Sig0;
        return dot(S,dsigdt)/Calc_qoct(sig);
        //Vec_t dSdt;
        //Dev (dsigdt, dSdt);
        //return Norm(dSdt);
    }
    Vec_t Sig0;
};

int main(int argc, char **argv) try
{
    Problem                  prob;
    Numerical::Diff<Problem> nd(&prob);

    size_t ndiv = 20;
    printf("%6s %12s %12s %12s %16s %16s\n","t","q","dqdt_num","dqdt","diff","error");
    for (size_t i=0; i<ndiv+1; ++i)
    {
        double t = (double)i/(double)ndiv;
        double dqdt = nd.DyDx (&Problem::q, t);
        printf("%6.3f %12.8f %12.8f %12.8f %16.8e %16.8e\n",t,prob.q(t),dqdt,prob.dqdt(t),fabs(dqdt-prob.dqdt(t)),nd.LastAbsErr);
    }

    // end
    return 0;
}
MECHSYS_CATCH
