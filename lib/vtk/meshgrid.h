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

#ifndef MECHSYS_MESHGRID_H
#define MECHSYS_MESHGRID_H

// MechSys
#include <mechsys/util/fatal.h>

class MeshGrid
{
public:
    // Constructors
    MeshGrid (double StartX, double EndX, double StepX,
              double StartY, double EndY, double StepY);

    MeshGrid (double StartX, double EndX, int nX,
              double StartY, double EndY, int nY); // nX and nY must be equal to or greater than 1

    MeshGrid (double StartX, double EndX, int nX,
              double StartY, double EndY, int nY,
              double StartZ, double EndZ, int nZ); // nX and nY must be equal to or greater than 1

    // Destructor
    ~MeshGrid ();

    // Methods
    int nX     () const { return _nX;     }
    int nY     () const { return _nY;     }
    int nZ     () const { return _nZ;     }
    int Length () const { return _length; }

    // Data
    double * X;
    double * Y;
    double * Z;

private:
    int _nX;
    int _nY;
    int _nZ;
    int _length;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline MeshGrid::MeshGrid (double StartX, double EndX, double StepX, double StartY, double EndY, double StepY)
    : X(NULL), Y(NULL), Z(NULL), _nX(0), _nY(0), _nZ(0), _length(0)
{
    double NumX = (EndX-StartX+StepX)/StepX;
    double NumY = (EndY-StartY+StepY)/StepY;

    if (NumX<1) throw new Fatal("MeshGrid::MeshGrid: __internal_error__ NumX must be greater than 0");
    if (NumY<1) throw new Fatal("MeshGrid::MeshGrid: __internal_error__ NumY must be greater than 0");

    int nX = static_cast<int>(NumX);
    int nY = static_cast<int>(NumY);

    _nX     = nX;
    _nY     = nY;
    _length = nX*nY;

    X = new double [_length];
    Y = new double [_length];

    int    k=0;
    double x;
    double y = StartY;
    for (int j=0; j<nY; ++j)
    {
        x = StartX;
        for (int i=0; i<nX; ++i)
        {
            X[k] = x;
            Y[k] = y;
            x += StepX;
            k++;
        }
        y += StepY;
    }
}

inline MeshGrid::MeshGrid(double StartX, double EndX, int nX, double StartY, double EndY, int nY)
    : X(NULL), Y(NULL), Z(NULL), _nX(0), _nY(0), _nZ(0), _length(0)
{
    if (nX<1) throw new Fatal("MeshGrid::MeshGrid: nX must be greater than 0");
    if (nY<1) throw new Fatal("MeshGrid::MeshGrid: nY must be greater than 0");

    double StepX = ( nX>1 ? (EndX-StartX)/(nX-1.0) : 0.0 );
    double StepY = ( nY>1 ? (EndY-StartY)/(nY-1.0) : 0.0 );

    _nX     = nX;
    _nY     = nY;
    _length = nX*nY;

    X = new double [_length];
    Y = new double [_length];

    int k = 0;
    double x;
    double y = StartY;
    for (int j=0; j<nY; ++j)
    {
        x = StartX;
        for (int i=0; i<nX; ++i)
        {
            X[k] = x;
            Y[k] = y;
            x += StepX;
            k++;
        }
        y += StepY;
    }
}

inline MeshGrid::MeshGrid(double StartX, double EndX, int nX, double StartY, double EndY, int nY, double StartZ, double EndZ, int nZ)
    : X(NULL), Y(NULL), Z(NULL), _nX(0), _nY(0), _nZ(0), _length(0)
{
    if (nX<1) throw new Fatal("MeshGrid::MeshGrid: nX must be greater than 0");
    if (nY<1) throw new Fatal("MeshGrid::MeshGrid: nY must be greater than 0");
    if (nZ<1) throw new Fatal("MeshGrid::MeshGrid: nZ must be greater than 0");

    double StepX = ( nX>1 ? (EndX-StartX)/(nX-1.0) : 0.0 );
    double StepY = ( nY>1 ? (EndY-StartY)/(nY-1.0) : 0.0 );
    double StepZ = ( nZ>1 ? (EndZ-StartZ)/(nZ-1.0) : 0.0 );

    _nX     = nX;
    _nY     = nY;
    _nZ     = nZ;
    _length = nX*nY*nZ;

    X = new double [_length];
    Y = new double [_length];
    Z = new double [_length];

    int p = 0;
    double x;
    double y;
    double z = StartZ;
    for (int k=0; k<nZ; ++k)
    {
        y = StartY;
        for (int j=0; j<nY; ++j)
        {
            x = StartX;
            for (int i=0; i<nX; ++i)
            {
                X[p] = x;
                Y[p] = y;
                Z[p] = z;
                x += StepX;
                p++;
            }
            y += StepY;
        }
        z += StepZ;
    }
}

inline MeshGrid::~MeshGrid ()
{
    if (X!=NULL) delete [] X;
    if (Y!=NULL) delete [] Y;
    if (Z!=NULL) delete [] Z;
}

std::ostream & operator<< (std::ostream & os, MeshGrid & mg)
{
    for (int i=0; i<mg.Length(); ++i)
    {
        if (mg.nX()>0) os << mg.X[i] << "   ";
        if (mg.nY()>0) os << mg.Y[i] << "   ";
        if (mg.nZ()>0) os << mg.Z[i] << "   ";
        os << std::endl;
    }
    return os;
}

#endif // MECHSYS_MESHGRID_H
