#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

#include "lbm/lattice2D.h"

using namespace std;

// Analysis constants
double u_max = 0.1;
double Re    = 100;
int    nx    = 200;
int    ny    = 50;
int    radius = ny/10 + 1;

double calcViscosity()
{
	double kin_visc = u_max*(2*radius)/Re; // nu
	return 3*kin_visc + 0.5;
}

void calcInitSpeed(int x, int y, double & vx, double & vy)
{
	double L = ny - 2;
	double yp = y - 1.5;
	vx = u_max*4/(L*L)*(L*yp - yp*yp);
	vy = 0.0;
}


int main()
{
	lattice2D lb(nx,ny);

	// Geometry
	//  Walls
	for (int i=0; i<nx; i++)
		for (int j=0; j<ny; j++)
		{
			if (j==0)     lb.get(i,j).setSolid();
			if (j==ny-1)  lb.get(i,j).setSolid();
		}

	//  Obstacle
	int obsX = nx/5;       // x position
	int obsY = ny/2+3;     // y position

	for (int i=0; i<nx; i++)
		for (int j=0; j<ny; j++)
			if ( pow(i-obsX,2.0) + pow(j-obsY,2.0) <= pow(radius,2.0) ) // circle equation
				lb.get(i,j).setSolid();

	// Define boundary conditions:
	for (int j=0; j<ny; j++)
	{
		double vx, vy;
		calcInitSpeed(0, j, vx, vy);
		
		lb.get(   0,j).setNeumannBC  (node::WEST, vx, vy);
		lb.get(nx-1,j).setDirichletBC(node::EAST, 1.0);
	}

	// Define Initial conditions: velocity speed and density
	for (int i=0; i<nx; i++)
		for (int j=0; j<ny; j++)
		{
			double vx, vy;
			double tau = calcViscosity();
			double rho = 1.0;
			calcInitSpeed(i, j, vx, vy);
			//lb.get(i,j).initiallize(tau, vx, vy, rho);
			lb.get(i,j).initiallize(tau, 0.08, 0, rho);
		}
	
	int max_it = 10000;
	lb.writeState(0);
	for (int i=1; i<max_it; i++)
	{

		if (i % 100 == 0)
		{
			cout << "  it: " << i;
			cout << "  maxVel: " << lb.maxVel();
			cout << endl;
			lb.writeState(i);
		}

		lb.collide();
		lb.stream();

	}

	return 0;

}
