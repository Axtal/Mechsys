#ifndef LATTICE2D_H
#define LATTICE2D_H


#include <iostream>
#include <sstream>
#include <fstream>

#include "lbm/matrix.h"
#include "lbm/node.h"

using namespace std;

class lattice2D
{

public:
	lattice2D(int nx, int ny)
	{
		// Allocate
		_nx = nx;
		_ny = ny;
		_size = nx*ny;
		_array.Resize(nx, ny);
	}

	node & get(int ix, int iy)
	{
		return _array(ix,iy);
	}

	void collide()
	{
		applyBryConds();
		// collide
		for (int i=0; i<_size; i++)
		{
			if (_array(i).getType() != node::SOLID)
				_array(i).collide();
		}

		bounceBack(); // ok.. good position for bounceBack
	}

	void stream()
	{
		// Local coordinates
		int c_local[9][2] = { 
		{ 0,  0},
		{ 1,  0},
		{ 0,  1},
		{-1,  0},
		{ 0, -1},
		{ 1,  1},
		{-1,  1},
		{-1, -1},
		{ 1, -1},
		};

		// stream
		for (int i=0; i<_nx; i++)
			for (int j=0; j<_ny; j++)
			{
				for (int k=0; k<9; k++)
				{
					int nextX = i + c_local[k][0];
					int nextY = j + c_local[k][1];
					if (nextX==-1)  nextX = _nx-1;
					if (nextX==_nx) nextX = 0;
					if (nextY==-1)  nextY = _ny-1;
					if (nextY==_ny) nextY = 0;
//					if (nextX<0 || nextX==_nx || nextY<0 || nextY==_ny) continue;
					_array(nextX,nextY).tmpF(k) = _array(i,j).F(k);
				}
			}
		
		// swap lattice2D
		for (int i=0; i<_size; i++)
			for (int k=0; k<9; k++)
			{
				node & N = _array(i);
				N.F(k) = N.tmpF(k);
			}
	}

	void applyBryConds()
	{
		// Set boundary conditions
		for (int i=0; i<_size; i++)
			_array(i).applyBryCond();
	}

	void bounceBack()
	{
		// Set bounce back boundary condition
		for (int i=0; i<_size; i++)
			if (_array(i).getType() == node::SOLID)
				_array(i).bounceBack();
	}

	void writeState(int it)
	{
		// Open/create file
		std::ofstream ofile;
		std::ostringstream oss_file;
		oss_file << "stream" << it << ".vtk";
		ofile.open(oss_file.str().c_str(), std::ios::out);

		// Temporary variable
		ostringstream oss;
		oss << "# vtk DataFile Version 2.0" << endl;
		
		oss << "Stream " << it << endl;
		oss << "ASCII" << endl;
		oss << "DATASET STRUCTURED_POINTS" << endl;
		oss << "DIMENSIONS " << _nx << " " << _ny << " " <<  1 << endl;
		oss << "ORIGIN " << 0 << " " << 0 << " " << 0 << endl;
		oss << "SPACING 1 1 1" << endl;
		oss << "POINT_DATA " << _nx*_ny << endl;

		//oss << "SCALARS Geom float 1" << endl;
		//oss << "LOOKUP_TABLE default" << endl;
		//for (int j=0; j<_ny; j++)
		//	for (int i=0; i<_nx; i++)
		//	{
		//		double value=0.0;
		//		if (_array(i,j).getType()==node::SOLID) value = 3.0;
		//		else value = 2.0;
		//		oss << value << endl;
		//	}

		oss << "SCALARS Density float 1" << endl;
		oss << "LOOKUP_TABLE default" << endl;
		for (int j=0; j<_ny; j++)
			for (int i=0; i<_nx; i++)
				oss << _array(i,j).getDensity() << endl;

		oss << "SCALARS SpeedX float 1" << endl; 
		oss << "LOOKUP_TABLE default" << endl;
		for (int j=0; j<_ny; j++)
			for (int i=0; i<_nx; i++)
			{
				double vx, vy;
				_array(i,j).getSpeed(vx, vy);
				if (_array(i,j).getType()==node::SOLID) vx = 0.0;
				oss << vx << endl;
			}

		oss << "SCALARS SpeedY float 1" << endl;
		oss << "LOOKUP_TABLE default" << endl;
		for (int j=0; j<_ny; j++)
			for (int i=0; i<_nx; i++)
			{
				double vx, vy;
				_array(i,j).getSpeed(vx, vy);
				if (_array(i,j).getType()==node::SOLID) vy = 0.0;
				oss << vy << endl;
			}
		
		//oss << "VECTORS speed float" << endl;
		//for (int j=0; j<_ny; j++)
		//	for (int i=0; i<_nx; i++)
		//	{
		//		double vx, vy;
		//		_array(i,j).getSpeed(vx, vy);
		//		if (_array(i,j).getType()==node::SOLID)
		//		{
		//			vx = 0.0;
		//			vy = 0.0;
		//		}
		//		oss << vx << " " << vy << " " << 0.0 << endl;
		//	}

		// Write to file
		ofile << oss.str();

		// Close file
		ofile.close();

	}

	double maxVel()
	{
		double vmax = 0.0;
		double vx, vy;
		for (int i=0; i<_size; i++)
		{
			_array(i).getSpeed(vx, vy);
			double vxy = vx*vx + vy*vy;
			if (vxy>vmax) vmax = vxy;
		}
		return sqrt(vmax);

	}

private:
	Matrix<node> _array;
	int _nx;
	int _ny;
	int _size;
	int _it;

};

#endif //LATTICE2D_H
