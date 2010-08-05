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
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/numstreams.h>

using std::cout;
using std::endl;
using Util::_6_3;
using Util::_8s;
using Util::_10_6;
using Util::SQ2;

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
}

struct Increments
{
    double dsx, dsy, dsz; // stress increments
    double dex, dey, dez; // strain increments
    Increments () : dsx(0.),dsy(0.),dsz(0.), dex(0.),dey(0.),dez(0.) {}
};

int main(int argc, char **argv) try
{
    cout << Util::_8s << endl;

    // input
    int  tst       = 4;
    bool isotrop   = false;
    bool pcte      = false;
    bool pstrain   = false;
    int  ninc      = 10;
    bool cor_drift = true;
    if (argc>1) tst       = atoi(argv[1]);
    if (argc>2) isotrop   = atoi(argv[2]);
    if (argc>3) pcte      = atoi(argv[3]);
    if (argc>4) pstrain   = atoi(argv[4]);
    if (argc>5) ninc      = atoi(argv[5]);
    if (argc>6) cor_drift = atoi(argv[6]);

    // path
    Array<Increments> incs;

    // default constants (for tests 1-4)
    double E   = 6000.0;
    double nu  = 0.3;
    double p0  = 100.0;
    double M   = 1.0;
    double pf  = 3.0*p0/(3.0-M);
    double qf  = M*pf;
    double phi = M2Phi(M,"cam");
    double K   = E/(3.0*(1.0-2.0*nu));
    double G   = E/(2.0*(1.0+nu));
    incs.Resize (1);
    if (isotrop)
    {
        incs[0].dsx = -p0;
        incs[0].dsy = -p0;
        incs[0].dsz = -p0;
    }
    else
    {
        if (pcte)
        {
            incs[0].dsz = -80.0;
            incs[0].dsx = -incs[0].dsz/2.0;
            incs[0].dsy = -incs[0].dsz/2.0;
        }
        else if (pstrain)
        {
            cout << "pstrain" << endl;
            /* TODO: check this
            incs[0].dez = -0.001;
            incs[0].dex = 1.0e-9;
            incs[0].dey = 1.0e-9;
            */
            double pf  = p0;
            double qf  = 100.0;
            double thf = 0.0;
            Vec3_t L;
            pqth2L (pf,qf,thf,L,"cam");
            incs[0].dsx = L(0)-(-p0);
            incs[0].dsy = L(1)-(-p0);
            incs[0].dsz = L(2)-(-p0);
        }
        else incs[0].dey = -0.04;
    }

    // model, parameters, and initial values
    String name;
    SDPair prms, inis;
    switch (tst)
    {
        case 1: // Elastic
        {
            cout << "[1;33m================================ (1) Elastic ===================================[0m\n";
            name = "LinElastic";
            prms.Set ("E nu", E, nu);
            inis.Set ("sx sy sz", -p0,-p0,-p0);
            break;
        }
        case 2: // von Mises
        {
            cout << "[1;33m================================ (2) von Mises =================================[0m\n";
            name = "ElastoPlastic";
            prms.Set ("E nu fc sY Hp", E, nu, FAILCRIT("VM"), qf, 0.0);
            inis.Set ("sx sy sz", -p0,-p0,-p0);
            break;
        }
        case 3: // Mohr-Coulomb
        {
            cout << "[1;33m================================ (3) Mohr-Coulomb ==============================[0m\n";
            name = "ElastoPlastic";
            prms.Set ("E nu fc c phi psi NonAssoc", E, nu, FAILCRIT("MC"), 0.1, phi, 1.0, 1.);
            inis.Set ("sx sy sz", -p0,-p0,-p0);
            break;
        }
        case 4: // CamClay
        {
            cout << "[1;33m================================ (4) Cam-Clay ==================================[0m\n";
            name = "CamClay";
            prms.Set ("lam kap nu phi", 0.01, 0.001, nu, phi);
            inis.Set ("sx sy sz v0", -p0,-p0,-p0, 2.0);
            break;
        }
        case 5:
        {
            cout << "[1;33m================================ (5) Cam-Clay ==================================[0m\n";
            // model
            name = "CamClay";
            prms.Set ("lam kap nu phi", 0.0891, 0.0196, 0.2, 31.6);
            inis.Set ("sx sy sz v0", -196.0,-196.0,-196.0, 1.691);

            // increments
            incs.Resize (3);
            incs[0].dsz = -196.0;
            incs[1].dsz =  294.0;
            incs[2].dsz = -481.0;
            break;
        }
        case 6: // Unconv01
        {
            cout << "[1;33m================================ (6) Unconv 01 =================================[0m\n";
            name = "Unconv01";
            double l0    = 0.001;
            double l1    = 0.005;
            double l3    = 0.008;
            double betb  = 100.0;
            double betbb = 100.0;
            double v0    = 2.0;
            double xR10  = log(p0*sqrt(3.0));
            double xR30  = xR10+0.1;
            prms.Set ("l0 l1 l3 betb betbb phi K G", l0, l1, l3, betb, betbb, phi, K, G);
            inis.Set ("sx sy sz v0 xR10 xR30", -p0,-p0,-p0, v0, xR10, xR30);
            break;
        }
        case 7: // Unconv02
        {
            cout << "[1;33m================================ (7) Unconv 02 =================================[0m\n";
            name = "Unconv02";
            double k0    = 0.05;
            double k1    = -0.5;
            double betk  = 5.0;
            double ev1   = -0.01/100.;
            double l0    = 0.001;
            double l1    = 0.005;
            double l3    = 0.02;
            double betb  = 50.0;
            double betbb = 50.0;
            double v0    = 2.0;
            double xR10  = log(p0*sqrt(3.0));
            double xR30  = xR10;
            prms.Set ("l0 l1 l3 betb betbb phi nu k0 k1 betk ev1", l0, l1, l3, betb, betbb, phi, nu, k0, k1, betk, ev1);
            inis.Set ("sx sy sz v0 xR10 xR30", -p0,-p0,-p0, v0, xR10, xR30);
            break;
        }
        case 8: // Unconv03
        {
            cout << "[1;33m================================ (8) Unconv 03 =================================[0m\n";
            name = "Unconv03";
            double l0    = 0.001;
            double l1    = 0.005;
            double l3    = 0.02;
            double betb  = 40.0;
            double betbb = 40.0;
            double v0    = 2.0;
            double xR10  = log(p0*sqrt(3.0));
            double xR30  = xR10+0.5;
            prms.Set ("l0 l1 l3 betb betbb phi nu", l0, l1, l3, betb, betbb, phi, nu);
            inis.Set ("sx sy sz v0 xR10 xR30", -p0,-p0,-p0, v0, xR10, xR30);
            break;
        }
        default: throw new Fatal("main: Test = %d is not available",tst);
    }
    Model * mdl = AllocModel (name, /*NDim*/3, prms);

    // set initial state
    EquilibState sta(/*NDim*/3);
    mdl->InitIvs (inis, &sta);

    // number of internal values
    size_t niv = size(sta.Ivs);

    // allocate memory
    D11.AllocSpace (6,6, 50);
    D12.AllocSpace (6,6, 50);
    D21.AllocSpace (6,6, 50);
    D22.AllocSpace (6,6, 50);

    // output initial state
    std::ostringstream oss;
    oss << _8s<< "sx"   << _8s<< "sy"   << _8s<< "sz" << _8s<< "sxy" << _8s<< "syz" << _8s<< "szx";
    oss << _8s<< "ex"   << _8s<< "ey"   << _8s<< "ez" << _8s<< "exy" << _8s<< "eyz" << _8s<< "ezx";
    for (size_t i=0; i<niv; ++i) oss << _8s<< mdl->IvNames[i]; oss << "\n";
    for (size_t i=0; i<6;   ++i) oss << _8s<< sta.Sig(i);
    for (size_t i=0; i<6;   ++i) oss << _8s<< sta.Eps(i);
    for (size_t i=0; i<niv; ++i) oss << _8s<< sta.Ivs(i);
    oss << "\n";

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
    for (size_t i=0; i<incs.Size(); ++i)
    {
        // set prescribed increments
        set_to_zero (DEps);
        set_to_zero (DSig);
        pDEps = false,false,false, false,false,false;
        DSig(0) = incs[i].dsx;
        DSig(1) = incs[i].dsy;
        DSig(2) = incs[i].dsz;
        if (fabs(incs[i].dex)>1.0e-10) { DEps(0) = incs[i].dex;  pDEps[0] = true; }
        if (fabs(incs[i].dey)>1.0e-10) { DEps(1) = incs[i].dey;  pDEps[1] = true; }
        if (fabs(incs[i].dez)>1.0e-10) { DEps(2) = incs[i].dez;  pDEps[2] = true; }

        // for each increment
        deps = DEps/ninc;
        dsig = DSig/ninc;
        for (int j=0; j<ninc; ++j)
        {
            // trial increments
            TgIncs (mdl, &sta, /*dT*/1.0, deps, dsig, pDEps, deps_tr, dsig_tr, divs_tr);

            // loading-unloading ?
            double aint = -1.0; // no intersection
            bool   ldg  = mdl->LoadCond (&sta, deps_tr, aint); // returns true if there is loading (also when there is intersection)

            // with intersection ?
            if (aint>0.0 && aint<1.0)
            {
                // update to intersection
                TgIncs (mdl, &sta, /*dT*/aint, deps, dsig, pDEps, deps_tr, dsig_tr, divs_tr);
                sta.Eps += deps_tr;
                sta.Sig += dsig_tr;
                sta.Ivs += divs_tr;
                deps = fabs(1.0-aint)*deps; // remaining of deps to be applied

                // drift correction
                if (cor_drift) mdl->CorrectDrift (&sta);

                // output
                for (size_t i=0; i<6;   ++i) oss << _8s<< sta.Sig(i);
                for (size_t i=0; i<6;   ++i) oss << _8s<< sta.Eps(i);
                for (size_t i=0; i<niv; ++i) oss << _8s<< sta.Ivs(i);
                oss << "\n";
            }

            // set loading flag (must be after intersection because the TgIncs during intersection must be calc with Ldg=false)
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
                TgIncs (mdl, &sta, dT, deps, dsig, pDEps, deps_1, dsig_1, divs_1);
                sta_1.Eps = sta.Eps + deps_1;
                sta_1.Sig = sta.Sig + dsig_1;
                sta_1.Ivs = sta.Ivs + divs_1;
                TgIncs (mdl, &sta_1, dT, deps, dsig, pDEps, deps_2, dsig_2, divs_2);
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
                    if (cor_drift) mdl->CorrectDrift (&sta);

                    // update stress path in model
                    mdl->UpdatePath (&sta, Vec_t(0.5*(deps_1+deps_2)), Vec_t(0.5*(dsig_1+dsig_2)));

                    // limit change on stepsize
                    if (m>mMax) m = mMax;

                    // output
                    for (size_t i=0; i<6;   ++i) oss << _8s<< sta.Sig(i);
                    for (size_t i=0; i<6;   ++i) oss << _8s<< sta.Eps(i);
                    for (size_t i=0; i<niv; ++i) oss << _8s<< sta.Ivs(i);
                    oss << "\n";
                }
                else if (m<mMin) m = mMin;

                // change next step size
                dT = m * dT;

                // check for last increment
                if (dT>1.0-T) dT = 1.0-T;
            }
            if (k>=MaxSS) throw new Fatal("main: Modified-Euler did not converge after %d substeps",k);

            cout << "increment = " << j << endl;
        }
    }

    // output file
    String fn;
    fn.Printf("test_models_%d.res",tst);
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();

    // clean up
    delete mdl;
    return 0;
}
MECHSYS_CATCH