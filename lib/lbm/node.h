#ifndef NODE_H
#define NODE_H

// Local velocities 
const double v_local[9][2] = { { 0,  0}, { 1,  0}, { 0,  1}, {-1,  0},
							   { 0, -1}, { 1,  1}, {-1,  1}, {-1, -1}, { 1, -1},
							   };
// Local weights
const double w_cent = 0.44444444444; // 4/9
const double w_side = 0.11111111111; // 1/9
const double w_diag = 0.02777777777; // 1/36
const double w[9] = { w_cent, w_side, w_side, w_side, w_side, w_diag, w_diag, w_diag, w_diag };

class node
{
public:
	enum     cellT { FLUID, SOLID};
	enum   bcSideT { NORTH, EAST, SOUTH, WEST};
	enum   bcTypeT { NEUMANN, DIRICHLET};


	node() { _cellType = FLUID; }

	void   setViscosity(double tau) { _tau = tau; }
	//double getX()              { return _x;};                     // x coordinate
	//double getY()              { return _y;};                     // y coordinate
	double getDensity()           
	{ 
		return _f[0]+_f[1]+_f[2]+_f[3]+_f[4]+_f[5]+_f[6]+_f[7]+_f[8];
	};                  // density
	double get_f(int idx)      { return _f[idx]; };               // distribution functions array
	cellT  getType()           { return _cellType; };             // get if node is boundary

	//void   updateEqDF(double vx, double vy, double rho)
	//{
	//	// Calc equilibrium distribution
	//	for (int k=0; k<9; k++)
	//	{
	//		double vxy  = vx*v_local[k][0] + vy*v_local[k][1];
	//		double vsqr = vx*vx + vy*vy;
	//		_f_eq[k] = w[k] * rho * ( 1.0 + 3.0*vxy + 4.5*vxy*vxy - 1.5*vsqr);
	//	}
	//}

	double calcEqDF(int index, double vx, double vy, double rho)
	{
		// Calc equilibrium distribution
		double vxy  = vx*v_local[index][0] + vy*v_local[index][1];
		double vsqr = vx*vx + vy*vy;
		return w[index] * rho * ( 1.0 + 3.0*vxy + 4.5*vxy*vxy - 1.5*vsqr);
	}

	//void   setDF (int index, double value) { _f[index] = value; }
	double & F   (int index) { return _f[index]; }
	double & tmpF(int index) { return _f_tmp[index]; }
	//double getEqDF(int index) { return _f_eq[index]; }
	void   initiallize(double tau, double vx, double vy, double rho)
	{
		setViscosity(tau);
		for (int k=0; k<9; k++)  // Set initial DF to equilibrium
			_f[k] = calcEqDF(k, vx, vy, rho);
	}

	void   getSpeed(double & vx, double & vy) 
	{
		double rho = getDensity();
		vx = 0.0; vy = 0.0;
		//if (_cellType != SOLID)
		//{
			vx = (_f[1] - _f[3] + _f[5] - _f[6] - _f[7] + _f[8]) / rho;
			vy = (_f[2] - _f[4] + _f[5] + _f[6] - _f[7] - _f[8]) / rho;
		//}
	}
	void setSolid() { _cellType = SOLID; };             // define the type of the node
	void setCell(cellT type) { _cellType = type; };             // define the type of the node
	void setNeumannBC(bcSideT bcSide, double vx0, double vy0)
	{
		_vx0 = vx0;
		_vy0 = vy0;
		_bcSide = bcSide;
		_bcType = NEUMANN;
	}
	void setDirichletBC(bcSideT bcSide, double rho)
	{
		_rho0 = rho;
		_bcSide = bcSide;
		_bcType = DIRICHLET;
	}
	void applyBryCond()
	{
		if (_bcType==NEUMANN)
		{
			if (_bcSide==WEST) neumannWest();
		}
		
		if (_bcType==DIRICHLET)
		{
			if (_bcSide==EAST) dirichletEast();
		}

	}
	void   neumannWest()
	{
		double vx = _vx0;
		double vy = _vy0;

		double rho = (_f[0]+_f[2]+_f[4] + 2.0*(_f[3]+_f[6]+_f[7]))/(1.0 - vx);

		_f[1] = _f[3] + 2.0/3.0*rho*vx;
		_f[5] = _f[7] + 1.0/6.0*rho*vx + 0.5*rho*vy - 0.5*(_f[2]-_f[4]);
		_f[8] = _f[6] + 1.0/6.0*rho*vx - 0.5*rho*vy + 0.5*(_f[2]-_f[4]);

	}
	void   dirichletEast()
	{
		double rho = _rho0;

		double vx = -1.0 + (_f[0]+_f[2]+_f[4] + 2.0*(_f[1]+_f[5]+_f[8]))/rho;
		double vy = 0.0;

		_f[3] = _f[1] - 2.0/3.0*rho*vx; 
		_f[7] = _f[5] - 1.0/6.0*rho*vx - 0.5*rho*vy + 0.5*(_f[2]-_f[4]);
		_f[6] = _f[8] - 1.0/6.0*rho*vx + 0.5*rho*vy - 0.5*(_f[2]-_f[4]);

	}

	void bounceBack()
	{
		double f_tmp[9];
		for (int i=1; i<9; i++) f_tmp[i] = _f[i];
		_f[1] = f_tmp[3]; _f[2] = f_tmp[4];
		_f[3] = f_tmp[1]; _f[4] = f_tmp[2];
		_f[5] = f_tmp[7]; _f[6] = f_tmp[8];
		_f[7] = f_tmp[5]; _f[8] = f_tmp[6];

	}

	void collide()
	{
		// Get information from lattice
		double vx, vy, rho;
		getSpeed(vx, vy);
		rho = getDensity();

		//Collision
		double omega = 1.0/_tau;
		//double omega = 1.86;
		for (int k=0; k<9; k++) 
		{
			double eqDF = calcEqDF(k, vx, vy, rho);
			_f[k] = (1.0-omega)*_f[k] + omega*eqDF;
		}
	}
private:
	//double  _x;    // x coordinate
	//double  _y;    // y coordinate
	double  _tau;
	double  _f[9];      // distribution functions array
	double  _f_tmp[9];  // temporary distribution functions array
	bcSideT _bcSide;    // define if node is boundary
	bcTypeT _bcType;    // define if node is boundary
	cellT   _cellType;  // define if node is boundary
	double  _vx0;
	double  _vy0;
	double  _rho0;
};

#endif // NODE_H
