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

/* Output - Copyright (C) 2007 Dorival de Moraes Pedroso */

#ifndef MPM_OUTPUT_H
#define MPM_OUTPUT_H

// STL
#include <iostream>
#include <fstream>
#include <string>

// External
#include <sys/stat.h> // for mkdir()

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>

// Local
#include <mechsys/mpm/defs.h>

namespace MPM {

/** Break down a line using separators. */
inline void SeparatedLine (std::string const & Line, std::string const & Separator, Array<std::string> & Res)
{
    /* Ex.:    /home/dorival/teste/A File.txt 
     *          R[0]  R[1]   R[2]     R[3]
     *
     *           ->1234->456a->789as
     *             R[0]  R[1]  R[2]
     */

    // Clear resulting array
    Res.Resize (0);

    // Fill resulting array
    if (Separator.empty()) Res.Push (Line);
    else
    {
        // Loop along the line
        size_t start = 0;
        size_t end   = Line.find (Separator, start);
        while (end!=std::string::npos)
        {
            Res.Push (Line.substr(start, end-start));
            start = end + Separator.size();
            end   = Line.find (Separator, start);
        }
        Res.Push (Line.substr(start));
    }
}

/** Create directory (recursively). */
bool MkdirP (char const * Path)
{
    // Tree
    Array<std::string> tree;
    SeparatedLine (Path, "/", tree);

    int res = -1;
    if (tree.Size()>0)
    {
        std::string left = "";
        for (size_t i=0; i<tree.Size(); ++i)
        {
            left += tree[i] + "/";
            // Upon successful completion, mkdir() shall return 0.
            // Otherwise, -1 shall be returned, no directory shall be created, and errno shall be set to indicate the error.
            res = mkdir (left.c_str(),S_IRWXU); //if (res!=) // directory may be existent, so, it is ok...
        }
    }

    return res==0 ? true : false;
}

/** VTK XML file format - Unstructured. */
namespace VTU
{

namespace Out
{

/** Output vectors at each point. */
inline void PointDataVectors (Array<Vector3D> const & V, char const * Name, std::ostream & F)
{
    F << "        <DataArray type=\"Float64\" Name=\"" << Name << "\" format=\"ascii\" NumberOfComponents=\"3\">" << std::endl;
    int k = 0; F << "        ";
    for (size_t i=0; i<V.Size(); ++i)
    {
        F << "  " << V[i](0) << " ";
        F <<         V[i](1) << " ";
        F <<         V[i](2);
        k++;
        if (k>7) { F<<std::endl; if (i<V.Size()-1) F<<"        "; k=0; }
        else { if (i==V.Size()-1) F<<std::endl; }
    }
    F << "        </DataArray>" << std::endl;
}

/** Output stress component at each point. */
inline void PointDataStressComp (Array<EquilibState*> const & S, int iComp, std::ostream & F)
{
    char   name[32];
    double coef = 1.0;
    switch (iComp)
    {
        case 0: { snprintf(name, 32, "stress_xx");           break; }
        case 1: { snprintf(name, 32, "stress_yy");           break; }
        case 2: { snprintf(name, 32, "stress_zz");           break; }
        case 3: { snprintf(name, 32, "stress_xy"); coef=SQ2; break; }
        case 4: { snprintf(name, 32, "stress_yz"); coef=SQ2; break; }
        case 5: { snprintf(name, 32, "stress_xz"); coef=SQ2; break; }
        default: { throw new Fatal("VTU::Out::PointsStressComp: Invalid iComp=%d. Must be 0<=iComp<=5.",iComp); }
    }
    F << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">" << std::endl;
    int k = 0; F << "        ";
    for (size_t i=0; i<S.Size(); ++i)
    {
        if (iComp>(int)size(S[i]->Sig)-1) F << " " << 0.0;
        else                              F << " " << S[i]->Sig(iComp)*coef;
        k++;
        if (k>3*7) { F<<std::endl; if (i<S.Size()-1) F<<"        "; k=0; }
        else { if (i==S.Size()-1) F<<std::endl; }
    }
    F << "        </DataArray>" << std::endl;
}

/** Output strain component at each point. */
inline void PointDataStrainComp (Array<EquilibState*> const & S, int iComp, std::ostream & F)
{
    char   name[32];
    double coef = 1.0;
    switch (iComp)
    {
        case 0: { snprintf(name, 32, "strain_xx");           break; }
        case 1: { snprintf(name, 32, "strain_yy");           break; }
        case 2: { snprintf(name, 32, "strain_zz");           break; }
        case 3: { snprintf(name, 32, "strain_xy"); coef=SQ2; break; }
        case 4: { snprintf(name, 32, "strain_yz"); coef=SQ2; break; }
        case 5: { snprintf(name, 32, "strain_xz"); coef=SQ2; break; }
        default: { throw new Fatal("VTU::Out::PointsStrainComp: Invalid iComp=%d. Must be 0<=iComp<=5.",iComp); }
    }
    F << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">" << std::endl;
    int k = 0; F << "        ";
    for (size_t i=0; i<S.Size(); ++i)
    {
        if (iComp>(int)size(S[i]->Eps)-1) F << " " << 0.0;
        else                              F << " " << S[i]->Eps(iComp)*coef;
        k++;
        if (k>3*7) { F<<std::endl; if (i<S.Size()-1) F<<"        "; k=0; }
        else { if (i==S.Size()-1) F<<std::endl; }
    }
    F << "        </DataArray>" << std::endl;
}

/** Output assymmetric tensors at each point. */
inline void PointDataTensors (Array<ATensor2> const & T, char const * Name, std::ostream & F)
{
    F << "        <DataArray type=\"Float64\" Name=\"" << Name << "\" format=\"ascii\" NumberOfComponents=\"9\">" << std::endl;
    int k = 0; F << "        ";
    for (size_t i=0; i<T.Size(); ++i)
    {
        F << "  " << T[i](0) << " " << T[i](3) << " " << T[i](5) << " ";
        F << "  " << T[i](6) << " " << T[i](1) << " " << T[i](4) << " ";
        F << "  " << T[i](8) << " " << T[i](7) << " " << T[i](2);
        k++;
        if (k>7) { F<<std::endl; if (i<T.Size()-1) F<<"        "; k=0; }
        else { if (i==T.Size()-1) F<<std::endl; }
    }
    F << "        </DataArray>" << std::endl;
}

/** Output symmetric tensors at each point. */
inline void PointDataTensors (Array<STensor2> const & T, char const * Name, std::ostream & F)
{
    F << "        <DataArray type=\"Float64\" Name=\"" << Name << "\" format=\"ascii\" NumberOfComponents=\"9\">" << std::endl;
    int k = 0; F << "        ";
    for (size_t i=0; i<T.Size(); ++i)
    {
        F << "  " << T[i](0)     << " " << T[i](3)/SQ2 << " " << T[i](5)/SQ2 << " ";
        F << "  " << T[i](3)/SQ2 << " " << T[i](1)     << " " << T[i](4)/SQ2 << " ";
        F << "  " << T[i](5)/SQ2 << " " << T[i](4)/SQ2 << " " << T[i](2);
        k++;
        if (k>7) { F<<std::endl; if (i<T.Size()-1) F<<"        "; k=0; }
        else { if (i==T.Size()-1) F<<std::endl; }
    }
    F << "        </DataArray>" << std::endl;
}

/** Output tensors at each point. */
inline void PointDataSig (Array<EquilibState*> const & S, char const * Name, std::ostream & F)
{
    F << "        <DataArray type=\"Float64\" Name=\"" << Name << "\" format=\"ascii\" NumberOfComponents=\"9\">" << std::endl;
    int k = 0; F << "        ";
    for (size_t i=0; i<S.Size(); ++i)
    {
        if (size(S[i]->Sig)==4)
        {
            F << "  " << S[i]->Sig(0)     << " " << S[i]->Sig(3)/SQ2 << " " << 0.0              << " ";
            F << "  " << S[i]->Sig(3)/SQ2 << " " << S[i]->Sig(1)     << " " << 0.0              << " ";
            F << "  " << 0.0              << " " << 0.0              << " " << S[i]->Sig(2);
        }
        else
        {
            F << "  " << S[i]->Sig(0)     << " " << S[i]->Sig(3)/SQ2 << " " << S[i]->Sig(5)/SQ2 << " ";
            F << "  " << S[i]->Sig(3)/SQ2 << " " << S[i]->Sig(1)     << " " << S[i]->Sig(4)/SQ2 << " ";
            F << "  " << S[i]->Sig(5)/SQ2 << " " << S[i]->Sig(4)/SQ2 << " " << S[i]->Sig(2);
        }
        k++;
        if (k>7) { F<<std::endl; if (i<S.Size()-1) F<<"        "; k=0; }
        else { if (i==S.Size()-1) F<<std::endl; }
    }
    F << "        </DataArray>" << std::endl;
}

/** Output internal variables at each point. */
inline void PointDataIvs (Array<EquilibState*> const & S, char const * Name, std::ostream & F)
{
    size_t max_nivs = 0;
    for (size_t i=0; i<S.Size(); ++i) { size_t nivs = size(S[i]->Ivs); if (nivs>max_nivs) max_nivs = nivs; }
    F << "        <DataArray type=\"Float64\" Name=\"ivs\" format=\"ascii\" NumberOfComponents=\""<<max_nivs<<"\">" << std::endl;
    int k = 0; F << "        ";
    for (size_t i=0; i<S.Size(); ++i)
    {
        for (size_t iv=0; iv<max_nivs; ++iv)
        {
            if (size(S[i]->Ivs)>iv) F << " " << S[i]->Ivs(iv);
            else                    F << " " << 0.0;
        }
        k++;
        if (k>7) { F<<std::endl; if (i<S.Size()-1) F<<"        "; k=0; }
        else { if (i==S.Size()-1) F<<std::endl; }
    }
    F << "        </DataArray>" << std::endl;
}

}; // namespace Out

}; // namespace VTU

}; // namespace MPM

#endif // MPM_OUTPUT_H
