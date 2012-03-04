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

/* Tiled - Copyright (C) 2007 Dorival de Moraes Pedroso */

#ifndef MPM_TILED_H
#define MPM_TILED_H

// STL
#include <cfloat> // for DBL_MAX

// FLTK
#include <FL/Fl.H>
#include <FL/Fl_Tile.H>
#include <FL/Fl_Group.H>
#include <FL/fl_draw.H>

namespace MPM {

typedef Fl_Group * (*pAllocGroup) (int X, int Y, int W, int H, void * Extra);

/////////////////////////////////////////////////////////////////////////////////////////// TwoTiledVert /////

class TwoTiledVert : public Fl_Tile
{
public:
	// Constructor
	TwoTiledVert (int X, int Y, int W, int H, void * Extra, pAllocGroup pAllocG1, pAllocGroup pAllocG2, double Hfactor=0.5) : Fl_Tile (X, Y, W, H)
	{
		int htop = static_cast<int>(Hfactor*H);
		int hbot = H-htop;
		_top = (*pAllocG1) (x(), y(),      W, htop, Extra);
		_bot = (*pAllocG2) (x(), y()+htop, W, hbot, Extra);
		end();
	}

	// Methods
	Fl_Group * T() { return _top; }
	Fl_Group * B() { return _bot; }

protected:
	// Resize proportionally
	void resize(int X, int Y, int W, int H)
	{
		// Get old proportions so we can preserve through resize
		double hfactor = static_cast<double> (_top->h())/h();
		int    htop    = static_cast<int>    (H*hfactor+0.5);

		// Resize our widget via Fl_Widget (to prevent children resizing)
		Fl_Widget::resize (X, Y, W, H);

		// Resize children with custom computations
		_top->resize (X, Y     , W, htop);
		_bot->resize (X, Y+htop, W, H-htop);
	}

private:
	// Data
	Fl_Group * _top;
	Fl_Group * _bot;
};

/////////////////////////////////////////////////////////////////////////////////////////// TwoTiledHorz /////

class TwoTiledHorz : public Fl_Tile
{
public:
	// Constructor
	TwoTiledHorz (int X, int Y, int W, int H, void * Extra, pAllocGroup pAllocG1, pAllocGroup pAllocG2, double Wfactor=0.5) : Fl_Tile (X, Y, W, H)
	{
		int wlef = static_cast<int>(Wfactor*W);
		int wrig = W-wlef;
		_lft = (*pAllocG1) (x(),      y(), wlef, H, Extra);
		_rig = (*pAllocG2) (x()+wlef, y(), wrig, H, Extra);
		end();
	}

	// Methods
	Fl_Group * L() { return _lft; }
	Fl_Group * R() { return _rig; }

protected:
	// Resize proportionally
	void resize(int X, int Y, int W, int H)
	{
		// Get old proportions so we can preserve through resize
		double wfactor = static_cast<double> (_lft->w())/w();
		int    wlef    = static_cast<int>    (W*wfactor+0.5);

		// Resize our widget via Fl_Widget (to prevent children resizing)
		Fl_Widget::resize (X, Y, W, H);

		// Resize children with custom computations
		_lft->resize (X      , Y , wlef   , H);
		_rig->resize (X+wlef , Y , W-wlef , H);
	}

private:
	// Data
	Fl_Group * _lft;
	Fl_Group * _rig;
};

////////////////////////////////////////////////////////////////////////////////////////////// FourTiled /////

class FourTiled : public Fl_Tile
{
public:
	// Constructor
	FourTiled (int X, int Y, int W, int H, void * Extra, pAllocGroup pAllocG1, pAllocGroup pAllocG2, pAllocGroup pAllocG3, pAllocGroup pAllocG4, double Wfactor=0.5, double Hfactor=0.5) : Fl_Tile (X, Y, W, H)
	{
		int htop = static_cast<int>(Hfactor*H);
		int hbot = H-htop;
		int wlef = static_cast<int>(Wfactor*W);
		int wrig = W-wlef;
		_toplft = (*pAllocG1) (x(),      y(),      wlef, htop, Extra);
		_toprig = (*pAllocG2) (x()+wlef, y(),      wrig, htop, Extra);
		_botlft = (*pAllocG3) (x(),      y()+htop, wlef, hbot, Extra);
		_botrig = (*pAllocG4) (x()+wlef, y()+htop, wrig, hbot, Extra);
		end();
	}

	// Methods
	Fl_Group * TL() { return _toplft; }
	Fl_Group * TR() { return _toprig; }
	Fl_Group * BL() { return _botlft; }
	Fl_Group * BR() { return _botrig; }

protected:
	// Resize proportionally
	void resize(int X, int Y, int W, int H)
	{
		// Get old proportions so we can preserve through resize
		double wfactor = static_cast<double> (_toplft->w())/w();
		double hfactor = static_cast<double> (_toplft->h())/h();
		int    wlef    = static_cast<int>    (W*wfactor+0.5);
		int    htop    = static_cast<int>    (H*hfactor+0.5);

		// Resize our widget via Fl_Widget (to prevent children resizing)
		Fl_Widget::resize (X, Y, W, H);

		// Resize children with custom computations
		_toplft->resize (X      , Y      , wlef   , htop);
		_toprig->resize (X+wlef , Y      , W-wlef , htop);
		_botlft->resize (X      , Y+htop , wlef   , H-htop);
		_botrig->resize (X+wlef , Y+htop , W-wlef , H-htop);
	}

private:
	// Data
	Fl_Group * _toplft;
	Fl_Group * _toprig;
	Fl_Group * _botlft;
	Fl_Group * _botrig;
};

}; // namespace MPM

#endif // MPM_TILED_H
