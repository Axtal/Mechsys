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

/* DrawArea2D - Copyright (C) 2007 Dorival de Moraes Pedroso */

#ifndef MPM_DRAWAREA2D_H
#define MPM_DRAWAREA2D_H

// STL
#include <cmath>  // for ceil and floor
#include <cfloat> // for DBL_EPSILON

// FLTK
#include <FL/Fl.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Group.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Pixmap.H>
#include <FL/Enumerations.H> // for Fl_Color

// libpng
#include <png.h>

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>

// Local
#include <mechsys/mpm/defs.h>
#include <mechsys/mpm/grid2d.h>
#include <mechsys/mpm/tensors.h>
#include <mechsys/mpm/mpoints2d.h>
#include <mechsys/mpm/infobox.h>
#include <mechsys/mpm/colormaps.h>

// Pixmaps
#include <mechsys/mpm/fixx.xpm>
#include <mechsys/mpm/fixy.xpm>
#include <mechsys/mpm/fixxy.xpm>

namespace MPM {

/* Structures. */
Fl_Menu_Item FldItems[] =
{
	{"Strain-xx" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Strain-yy" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Strain-zz" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Strain-xy" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Stress-xx" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Stress-yy" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Stress-zz" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Stress-xy" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Cam-p"     , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Cam-q"     , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Vel-norm"  , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"ClrIdx"    , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{0,0,0,0,0,0,0,0,0}
};
Fl_Menu_Item PtType[] =
{
	{"BoxBord" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Box"     , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{"Circles" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
	{0,0,0,0,0,0,0,0,0}
};

/* Drawing area 2D. */
class DrawArea2D : public Fl_Group
{
public:
	/* Constructor. */
	DrawArea2D (int xmin, int ymin, int width, int height); // Screen coordinates

	/* Destructor. */
	virtual ~DrawArea2D ();

	// Methods
	void   SetGrid   (Grid2D    const * G) { _g = G; } ///< Set grid
	void   SetPoints (MPoints2D const * P) { _p = P; } ///< Set material points
	void   SetTime   (double    const * t) { _t = t; } ///< Set time
	void   SetM      (double    const * M) { _M = M; } ///< Set multiplier for external forces
	void   CalcSF    ();                               ///< Calculate scale factors
	size_t SelP      () const        { return _selp; } ///< Return the selected mat point id
	size_t SelN      () const        { return _seln; } ///< Return the selected grid node id
	void   ResetSelP (size_t i);                       ///< Set the selected mat point for output (reset previous to -1 => must call redraw())

	// Set methods
	DrawArea2D & EqScales  (bool EQScales=true) { _eqsf=EQScales;  return (*this); } ///< Set equal scale factors
	DrawArea2D & RecalcSF  (bool REcalcSF=true) { _recsf=REcalcSF; return (*this); } ///< Recalculate scale factors during draw ?
	DrawArea2D & WithFrame (bool WFrame  =true) { _wframe=WFrame;  return (*this); } ///< Draw an all-around frame ?

    // Output method
    void SavePNG (char const * Filename)
    {
        FILE * fp = fopen (Filename, "wb");
        if (fp!=NULL)
        {
            int xmin = _x(_g->xMin());
            int xmax = _x(_g->xMax());
            int ymin = _y(_g->yMin());
            int ymax = _y(_g->yMax());
            int ww   = xmax - xmin;
            int hh   = ymin - ymax;
            if (ww>0 && hh>0)
            {
                // write png header information
                png_structp png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, (void*)NULL, NULL, NULL);
                png_infop   inf_ptr = png_create_info_struct  (png_ptr);
                png_init_io    (png_ptr, fp);
                png_set_IHDR   (png_ptr, inf_ptr, ww, hh, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
                png_write_info (png_ptr, inf_ptr);
                
                // write pixels
                uchar * pixels = fl_read_image (NULL, xmin, ymax, ww, hh, 0);
                for (int k=0; k<hh; k++) png_write_row (png_ptr, pixels + 3*k*ww);
                delete [] pixels;
                png_write_end (png_ptr, inf_ptr);
                png_destroy_write_struct (&png_ptr, &inf_ptr);
            }
        }
        fclose (fp);
    }

	/** Draw method (internal). */
	void draw ();

protected:
	/** Catch events. */
	int handle (int E);

private:
	// Widgets
	Fl_Input        * _w_psz; ///< Point size
	Fl_Choice       * _w_pty; ///< Point type
	Fl_Choice       * _w_cmt; ///< Color map type
	Fl_Choice       * _w_fld; ///< Field to be coloured
	Fl_Button       * _w_set; ///< Set all data
	Fl_Check_Button * _w_inv; ///< Invert colormap?
	Fl_Check_Button * _w_ldg; ///< Loads
	Fl_Check_Button * _w_vel; ///< Velocity
	Fl_Check_Button * _w_sel; ///< Velocity
	Fl_Check_Button * _w_cbr; ///< Colorbar
	
	// Data
	double     _sfX;    ///< Scale factor: x=xmin+(X-Xmin)*sf and y=ymin+ymax-(Y-Ymin)*sf, where x and y are screen coordinates
	double     _sfY;    ///< Scale factor: x=xmin+(X-Xmin)*sf and y=ymin+ymax-(Y-Ymin)*sf, where x and y are screen coordinates
	bool       _eqsf;   ///< Equal scale factors ?
	bool       _recsf;  ///< Recalculate scale factors during draw ?
	bool       _sfok;   ///< Is scale factor ok (set true after a first time calculation) ?
	bool       _wframe; ///< With frame? draw an all-around frame ?
	double     _Xmin;   ///< Minimum x value (real coordinates)
	double     _Ymin;   ///< Minimum y value (real coordinates)
	double     _Xmax;   ///< Maximum x value (real coordinates)
	double     _Ymax;   ///< Maximum y value (real coordinates)
	int        _hb;     ///< Increment for horizontal borders (to the inside) (screen coordinates)
	int        _vb;     ///< Increment for vertical borders (to the inside) (screen coordinates)
	int        _lfz;    ///< Labels font size
	int        _selp;   ///< Selected material point
	int        _seln;   ///< Selected grid node
	int        _pselp;  ///< Previous selected material point
	int        _pseln;  ///< Previous selected grid node
	int        _ctrh;   ///< Control height
	int        _cmbw;   ///< ColorBar width
	int        _cmbh;   ///< ColorBar height
	ClrMapType _cmtype; ///< ColorMap type
	FieldType  _field;  ///< Field
	bool       _cmrinv; ///< Invert colormap?
    bool       _loads;  ///< Show loads ?
    bool       _veloc;  ///< Show velocity ?
    bool       _selec;  ///< Show selection ?
    bool       _cbar;   ///< Show colorbar ?

	// Elements
	Grid2D    const * _g;  ///< Grid
	MPoints2D const * _p;  ///< Material points
	double    const * _t;  ///< Current time
	double    const * _M;  ///< Multiplier for external forces
	CurveProps        _cg; ///< Curve properties for drawing nodes (the grid)
	CurveProps        _cp; ///< Curve properties for drawing material points
	
	// XPMs
	Fl_Pixmap * _fixx;  ///< Pixmap representing a X fixing icon
	Fl_Pixmap * _fixy;  ///< Pixmap representing a Y fixing icon
	Fl_Pixmap * _fixxy; ///< Pixmap representing a XY fixing icon
	int         _hf;    ///< Fix pixmaps horizontal width
	int         _vf;    ///< Fix pixmaps vertical width

	// Private methods
	int    _x (double X) const { return static_cast<int>((x()+1+_hb)                          +_sfX*(X-_Xmin));     } ///< X in screen coordinates
	int    _y (double Y) const { return static_cast<int>((y()+1+_vb+_ctrh)+(h()-2-2*_vb-_ctrh)-_sfY*(Y-_Ymin));     } ///< Y in screen coordinates
	int    _l (double L) const { return static_cast<int>((_sfX>_sfY?_sfY:_sfX)*L);                      } ///< Length in screen coordinates
	double _X (int   xx) const { return _Xmin+(static_cast<double>(xx)-(x()+1+_hb))/_sfX;               } ///< x in real coordinates
	double _Y (int   yy) const { return _Ymin+((y()+1+_vb)+(h()-2-2*_vb)-static_cast<double>(yy))/_sfY; } ///< y in real coordinates

	// Private methods
	void _draw_vector_field (Array<Vector3D> const & P, Array<Vector3D> const & V, Fl_Color Clr); ///< Draw vector field
	void _draw_arrow        (int xi, int yi, int xf, int yf, Fl_Color Clr=FL_BLUE);               ///< Draw arrow
	void _draw_loaded_nodes ();                                                                   ///< Draw loaded nodes
	void _draw_coords       ();                                                                   ///< Draw coordinates in a small rectangle, without redrawing entire screen
	void _draw_selected     ();                                                                   ///< Draw coordinates in a small rectangle, without redrawing entire screen
	void _draw_loaded_points();
	void _draw_disp_points  ();
	void _draw_colorbar     ();
	void _draw_colorbar_labels(double min, double max);

	// Callbacks (must be static)
	static void _set_cb (Fl_Widget * o, void * v) { ((DrawArea2D*)v)->_set(); }

	// Control methods
	void _set();

}; // class DrawArea2D


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline DrawArea2D::DrawArea2D(int xmin, int ymin, int width, int height)
	: Fl_Group (xmin,ymin,width,height,0), // {
	  _sfX    (1.0),
	  _sfY    (1.0),
	  _eqsf   (true),
	  _recsf  (false),
	  _sfok   (false),
	  _wframe (true),
	  _Xmin   (0.0),
	  _Ymin   (0.0),
	  _Xmax   (1.0),
	  _Ymax   (1.0),
	  _lfz    (8),
	  _selp   (-1),
	  _seln   (-1),
	  _pselp  (-1),
	  _pseln  (-1),
	  _ctrh   (25),
	  _cmbw   (50),
	  _cmbh   (200),
	  _cmtype (CMT_JET),
	  _field  (FLD_TAG),//FLD_CAMQ),
	  _cmrinv (false),
	  _g      (NULL),
	  _p      (NULL),
	  _t      (NULL),
	  _M      (NULL),
	  _hf     (11),
	  _vf     (11)
{
	_w_psz = new Fl_Input        (x(), y(), 0, 0, "Sz");
	_w_pty = new Fl_Choice       (x(), y(), 0, 0, "Ty");
	_w_ldg = new Fl_Check_Button (x(), y(), 0, 0, "Loads");
	_w_vel = new Fl_Check_Button (x(), y(), 0, 0, "Veloc");
	_w_sel = new Fl_Check_Button (x(), y(), 0, 0, "Sel");
	_w_cbr = new Fl_Check_Button (x(), y(), 0, 0, "CBar");
	_w_cmt = new Fl_Choice       (x(), y(), 0, 0, "CMap");
	_w_fld = new Fl_Choice       (x(), y(), 0, 0, "Field");
	_w_set = new Fl_Button       (x(), y(), 0, 0, "Set");
	_w_inv = new Fl_Check_Button (x(), y(), 0, 0, "Inv");
	end();

	// Set curve properties for drawing nodes
	_cg.Typ = CT_BOTH;
	_cg.Clr = FL_GRAY;
	_cg.Lty = FL_SOLID;
	_cg.Lwd = 1;
	_cg.Pch = 1;
	_cg.Psz = 4;

	// Set curve properties for drawing material points
	_cp.Typ = CT_POINTS;
	_cp.Clr = FL_RED;
	_cp.Lty = FL_SOLID;
	_cp.Lwd = 1;
	_cp.Pch = 0;
	_cp.Psz = 8;

	// Borders
	_hb = 6+_hf+1; // 6 for border + hf for icons (xpms)
	_vb = 6+_vf+1; // 6 for border + vf for icons (xpms)

	// Set pixmaps
	_fixx  = new Fl_Pixmap (fixx_xpm);
	_fixy  = new Fl_Pixmap (fixy_xpm);
	_fixxy = new Fl_Pixmap (fixxy_xpm);

	// Initialize widgets
	_w_pty->menu  (PtType);
	_w_cmt->menu  (ClrMap::Items);
	_w_cmt->value (static_cast<int>(_cmtype));
	_w_fld->menu  (FldItems);
	_w_fld->value (static_cast<int>(_field));
	char buf[256];
	snprintf(buf, 256, "%d", _cp.Psz); _w_psz->value(buf);

	// Callbacks
	_w_set->callback (_set_cb, this);
} // }

inline DrawArea2D::~DrawArea2D()
{
	delete _fixx;
	delete _fixy;
	delete _fixxy;
}

inline void DrawArea2D::CalcSF()
{
	// Check input
	if (_g==NULL)       return;
	if (_g->nNodes()<2) return;

	// Bounding box
	_Xmin = _g->X(0)(0);
	_Ymin = _g->X(0)(1);
	_Xmax = _g->X(1)(0);
	_Ymax = _g->X(1)(1);
	for (size_t i=0; i<_g->nNodes(); ++i)
	{
		if (_g->X(i)(0)<_Xmin) _Xmin = _g->X(i)(0);
		if (_g->X(i)(0)>_Xmax) _Xmax = _g->X(i)(0);
		if (_g->X(i)(1)<_Ymin) _Ymin = _g->X(i)(1);
		if (_g->X(i)(1)>_Ymax) _Ymax = _g->X(i)(1);
	}

	// Scale factors
	_sfX = static_cast<double>((w()-2-2*_hb      )/(_Xmax-_Xmin));
	_sfY = static_cast<double>((h()-2-2*_vb-_ctrh)/(_Ymax-_Ymin));
	if (_eqsf)
	{
		double sf = (_sfX>_sfY ? _sfY : _sfX);
		_sfX = sf;
		_sfY = sf;
	}
}

inline void DrawArea2D::ResetSelP(size_t i)
{
	_pselp = -1;
	_selp  = i;
}

inline void DrawArea2D::draw()
{
	// User damage ONLY ?
	if (damage()==FL_DAMAGE_USER1)
	{
		// Just draw coords and done
		_draw_coords();
		return;
	}

	// Just draw selected points or nodes and done
	if (damage()==FL_DAMAGE_USER2)
	{
		_draw_selected();
		return;
	}

	// Clear control background
	fl_color (FL_BACKGROUND_COLOR);
	fl_rectf (x(), y(), w(), _ctrh);

	// Let group draw itself
	_w_psz->resize(x()+20,      y(), 40, 25);
	_w_pty->resize(x()+20+70,   y(), 90, 25);
    _w_ldg->resize(x()+20+170,  y(), 80, 25);
    _w_vel->resize(x()+20+235,  y(), 80, 25);
    _w_sel->resize(x()+20+300,  y(), 70, 25);
    _w_cbr->resize(x()+w()-460, y(), 80, 25);
	_w_cmt->resize(x()+w()-350, y(), 80, 25);
	_w_inv->resize(x()+w()-260, y(), 60, 25);
	_w_fld->resize(x()+w()-170, y(),100, 25);
	_w_set->resize(x()+w()-62,  y(), 60, 25);
	Fl_Group::draw();

	// Clear the background
	fl_color (FL_WHITE);
	fl_rectf (x(), y()+_ctrh, w(), h()-_ctrh);

	// Draw an all-around frame
	if (_wframe)
	{
		fl_color      (FL_BLACK);
		fl_line_style (FL_SOLID, 1);
		fl_rect       (x(), y()+_ctrh, w(), h()-_ctrh);
	}

	// Check input
	if (_g==NULL)       return;
	if (_g->nNodes()<2) return;

	// Calculate scale factors and draw rulers
	if (_recsf || !_sfok)
	{
		CalcSF ();
		_sfok = true;
	}

    // Set clipping in order to avoid drawing outside drawing area
    fl_push_clip (x()+1, y()+_ctrh+1, w()-2, h()-_ctrh-2);

	// Draw lines
	if (_cg.Typ==CT_LINES || _cg.Typ==CT_BOTH)
	{
		fl_color      (_cg.Clr);
		fl_line_style (_cg.Lty, _cg.Lwd);

		for (size_t i=0; i<_g->nRow(); ++i)
			fl_line (_x(_g->xMin()), _y(_g->yMin()+i*_g->L(1)), _x(_g->xMax()), _y(_g->yMin()+i*_g->L(1)));

		for (size_t j=0; j<_g->nCol(); ++j)
			fl_line (_x(_g->xMin()+j*_g->L(0)), _y(_g->yMin()), _x(_g->xMin()+j*_g->L(0)), _y(_g->yMax()));
	}

	// Draw fixed nodes
	for (size_t i=0; i<_g->nFix(); ++i)
	{
		int x = _x(_g->FixN(i)(0));
		int y = _y(_g->FixN(i)(1));
		switch (_g->FixT(i))
		{
			case FIX_X   : { _fixx->draw(x-_hf,y-_vf/2); break; }
			case FIX_Y   : { _fixy->draw(x-_hf/2,y);     break; }
			case FIX_Z   : {                             break; }
			case FIX_XY  : { _fixxy->draw(x-_hf/2,y);    break; }
			case FIX_YZ  : {                             break; }
			case FIX_ZX  : {                             break; }
			case FIX_XYZ : {                             break; }
		}
	}

	// Draw cell ids
	char buf[256];
	/*
	fl_font  (0, 8);
	fl_color (FL_BLACK);
	for (size_t i=0; i<_g->nCells(); ++i)
	{
		int p0 = _g->C(i)(0);
		int p1 = _g->C(i)(1);
		int p2 = _g->C(i)(2);
		int p3 = _g->C(i)(3);
		int xc = _x(0.5*(_g->X(p0)(0)+_g->X(p2)(0)));
		int yc = _y(0.5*(_g->X(p0)(1)+_g->X(p2)(1)));
		snprintf (buf, 256, "%d\n(%d,%d,%d,%d)", i,p0,p1,p2,p3); // format text
		fl_draw  (buf, xc, yc, 0, 0, FL_ALIGN_CENTER);           // draw text
	}
	*/

	// Draw node ids
	/*
	for (size_t i=0; i<_g->nNodes(); ++i)
	{
		int xi = _x(_g->X(i)(0));
		int yi = _y(_g->X(i)(1));
		snprintf (buf, 256, "%d", i);
		fl_draw  (buf, xi, yi, 0, 0, FL_ALIGN_CENTER); // draw text
	}
	*/

	// Check input
	if (_p==NULL) return;

	// Draw points
	double min = 0;
	double max = 0;
	if (_cp.Typ==CT_POINTS || _cp.Typ==CT_BOTH)
	{
		int r = _cp.Psz/2;
		switch (_field)
		{
			case FLD_EPSXX: { _p->MinMaxE(0,min,max);                     break; }
			case FLD_EPSYY: { _p->MinMaxE(1,min,max);                     break; }
			case FLD_EPSZZ: { _p->MinMaxE(2,min,max);                     break; }
			case FLD_EPSXY: { _p->MinMaxE(3,min,max); min*=SQ2; max*=SQ2; break; }
			case FLD_SIGXX: { _p->MinMaxS(0,min,max);                     break; }
			case FLD_SIGYY: { _p->MinMaxS(1,min,max);                     break; }
			case FLD_SIGZZ: { _p->MinMaxS(2,min,max);                     break; }
			case FLD_SIGXY: { _p->MinMaxS(3,min,max); min*=SQ2; max*=SQ2; break; }
			case FLD_CAMP:  { _p->MinMaxCamp(min,max);                    break; }
			case FLD_CAMQ:  { _p->MinMaxCamq(min,max);                    break; }
			case FLD_VELN:  { _p->MinMaxVelN(min,max);                    break; }
            case FLD_TAG:   { break; }
		}
		if (!(max-min>1.0e-9)) max=min;
		fl_line_style (FL_SOLID, 1);
		for (size_t i=0; i<_p->nPoints(); ++i)
		{
			/*
			int xi = _x(_p->P(i)(0)-_p->l(i)(0));
			int yi = _y(_p->P(i)(1)-_p->l(i)(1));
			int xf = _x(_p->P(i)(0)+_p->l(i)(0));
			int yf = _y(_p->P(i)(1)+_p->l(i)(1));
			fl_color (FL_RED);
			fl_rectf (xi, yf, xf-xi, yi-yf);
			fl_color (FL_BLACK);
			fl_rect  (xi, yf, xf-xi, yi-yf);
			*/
			Fl_Color clr = FL_WHITE;
			double   val = max;
			if (max-min>1.0e-9)
			{
				switch (_field)
				{
					case FLD_EPSXX: { val=(_p->e(i)(0)    -min)/(max-min); break; }
					case FLD_EPSYY: { val=(_p->e(i)(1)    -min)/(max-min); break; }
					case FLD_EPSZZ: { val=(_p->e(i)(2)    -min)/(max-min); break; }
					case FLD_EPSXY: { val=(_p->e(i)(3)*SQ2-min)/(max-min); break; }
					case FLD_SIGXX: { val=(_p->s(i)(0)    -min)/(max-min); break; }
					case FLD_SIGYY: { val=(_p->s(i)(1)    -min)/(max-min); break; }
					case FLD_SIGZZ: { val=(_p->s(i)(2)    -min)/(max-min); break; }
					case FLD_SIGXY: { val=(_p->s(i)(3)*SQ2-min)/(max-min); break; }
					case FLD_CAMP:  { val=(Calc_pcam(_p->s(i))-min)/(max-min); break; }
					case FLD_CAMQ:  { val=(Calc_qcam(_p->s(i))-min)/(max-min); break; }
					case FLD_VELN:  { val=(NormV(_p->v(i))-min)/(max-min); break; }
                    case FLD_TAG:   { break; }
				}
			}
            if (_field==FLD_TAG) clr = ClrMap::GetColor(_p->Clr(i));
            else
            {
                switch (_cmtype)
                {
                    case CMT_JET: { int idx=static_cast<int>(_cmrinv?ClrMap::JET_NCOLORS-1-(ClrMap::JET_NCOLORS-1.0)*val:(ClrMap::JET_NCOLORS-1.0)*val); clr=ClrMap::JET[idx]; break; }
                    case CMT_HOT: { int idx=static_cast<int>(_cmrinv?ClrMap::HOT_NCOLORS-1-(ClrMap::HOT_NCOLORS-1.0)*val:(ClrMap::HOT_NCOLORS-1.0)*val); clr=ClrMap::HOT[idx]; break; }
                    case CMT_BWP: { int idx=static_cast<int>(_cmrinv?ClrMap::BWP_NCOLORS-1-(ClrMap::BWP_NCOLORS-1.0)*val:(ClrMap::BWP_NCOLORS-1.0)*val); clr=ClrMap::BWP[idx]; break; }
                }
            }
			//if (_p->HasT(i)) clr=FL_BLACK;
			fl_color   (clr);
			if (_cp.Pch==0) // Real cube
			{
                _p->CalcC(i);
				fl_polygon ( _x(_p->C[0](0)),
				             _y(_p->C[0](1)),
				             _x(_p->C[1](0)),
				             _y(_p->C[1](1)),
				             _x(_p->C[2](0)),
				             _y(_p->C[2](1)),
				             _x(_p->C[3](0)),
				             _y(_p->C[3](1)) );
				fl_color (FL_BLACK);
				fl_loop  ( _x(_p->C[0](0)),
				           _y(_p->C[0](1)),
				           _x(_p->C[1](0)),
				           _y(_p->C[1](1)),
				           _x(_p->C[2](0)),
				           _y(_p->C[2](1)),
				           _x(_p->C[3](0)),
				           _y(_p->C[3](1)) );
			}
			else if (_cp.Pch==1) // Real cube with borders
			{
                _p->CalcC(i);
				fl_polygon ( _x(_p->C[0](0)),
				             _y(_p->C[0](1)),
				             _x(_p->C[1](0)),
				             _y(_p->C[1](1)),
				             _x(_p->C[2](0)),
				             _y(_p->C[2](1)),
				             _x(_p->C[3](0)),
				             _y(_p->C[3](1)) );
			}
			else if (_cp.Pch==2) // Filled circle
			{
				fl_pie (_x(_p->P(i)(0))-r, _y(_p->P(i)(1))-r, 2*r+1, 2*r+1, 0,360);
			}
		}
	}

	// Draw points lengths
	/*
	for (size_t i=0; i<_p->nPoints(); ++i)
	{
		_draw_arrow (_x(_p->P(i)(0)), _y(_p->P(i)(1)), _x(_p->P(i)(0)+_p->lx(i)(0)), _y(_p->P(i)(1)+_p->lx(i)(1)));
		_draw_arrow (_x(_p->P(i)(0)), _y(_p->P(i)(1)), _x(_p->P(i)(0)+_p->ly(i)(0)), _y(_p->P(i)(1)+_p->ly(i)(1)));
	}
	*/

	// Draw mat points ids
	/*
	fl_font  (0, 8);
	fl_color (FL_BLACK);
	for (size_t i=0; i<_p->nPoints(); ++i)
	{
		int xc = _x(_p->P(i)(0));
		int yc = _y(_p->P(i)(1));
		snprintf (buf, 256, "%d", i);                  // format text
		fl_draw  (buf, xc, yc, 0, 0, FL_ALIGN_CENTER); // draw text
	}
	*/
	
	// Draw velocity field
    if (_veloc) _draw_vector_field (_p->Ps(), _p->v(), FL_BLUE);

	// Draw boundary mat points
	/*
	if (_p->nPointsOnBry()>0)
	{
		int r = _cp.Psz/2+1;
		fl_color      (FL_BLACK);
		fl_line_style (FL_SOLID, 2);
		for (size_t i=0; i<_p->nPointsOnBry(); ++i)
			fl_pie (_x(_p->B(i)(0))-r, _y(_p->B(i)(1))-r, 2*r+1, 2*r+1, 0,360);
	}
	*/

	// Draw time
	if (_t==NULL) return;
	snprintf (buf, 256, "t=%g", (*_t)); // format text
	fl_font  (0, 16);                   // Set font for labels
	fl_color (FL_RED);                  // set color
	fl_draw  (buf, x()+w()-148, y()+h()-18-18, 0, 0, FL_ALIGN_LEFT_TOP); // draw text
	/*
	fl_font  (0, 22);                      // Set font for labels
	fl_color (FL_RED);                    // set color
	fl_draw  (buf, x()+w()-2, y(), 0, 0, FL_ALIGN_RIGHT_TOP); // draw text
	*/

	// Draw coords
	_draw_coords();

	// Draw loaded nodes and points
    if (_loads)
    {
        //_draw_loaded_nodes ();
        _draw_loaded_points();
    }

	// Draw selected
    if (_selec) _draw_selected();

    // Draw points with disp (velocity) prescribed
	_draw_disp_points();

	// Draw colorbar
    if (_cbar)
    {
        _draw_colorbar();
        _draw_colorbar_labels(min,max);
    }

    // Pop clipping
    fl_pop_clip ();
}

inline int DrawArea2D::handle(int E)
{
	int ret = Fl_Group::handle(E);
	switch (E)
	{
		case FL_MOVE:
		{
			damage (FL_DAMAGE_USER1);
			return 1;
		}
		case FL_PUSH:
		{
			if (_g==NULL) return ret;

			/*
			char   buf[256];
			double X = _X(Fl::event_x());
			double Y = _Y(Fl::event_y());
			sprintf(buf, "Selected: x=%.4g y=%.4g", X, Y);
			*/

			// Search mat points
			_pselp = _selp; // previous selected
			_selp  = -1;    // current selected
			for (size_t i=0; i<_p->nPoints(); ++i)
			{
				long dx = _x(_p->P(i)(0))-static_cast<long>(Fl::event_x());
				long dy = _y(_p->P(i)(1))-static_cast<long>(Fl::event_y());
				long d  = static_cast<long>(sqrt(dx*dx+dy*dy));
				if (d<10) { _selp=i; break; }
			}

			// Search nodes
			/*
			_pseln = _seln; // previous selected
			_seln  = -1;    // current selected
			if (_selp<0)
			{
				for (size_t i=0; i<_g->nNodes(); ++i)
				{
					long dx = _x(_g->X(i)(0))-static_cast<long>(Fl::event_x());
					long dy = _y(_g->X(i)(1))-static_cast<long>(Fl::event_y());
					long d  = static_cast<long>(sqrt(dx*dx+dy*dy));
					if (d<4) { _seln=i; break; }
				}
			}
			*/

			// Execute selection action
			if (_selp>=0) (*callback()) (this, user_data());

			// Redraw
			damage (FL_DAMAGE_USER2);
			return 1;
		}
	}
	return ret;
}


/* private */

inline void DrawArea2D::_draw_vector_field (Array<Vector3D> const & P, Array<Vector3D> const & V, Fl_Color Clr)
{
	// Check
	if (P.Size()<1) return;

	// Boundaries
    double maxvnorm = NormV(V[0]);
	for (size_t i=0; i<P.Size(); ++i)
	{
        double vnorm = NormV(V[i]);
        if (vnorm>maxvnorm) maxvnorm = vnorm;
	}

	// Draw
    if (maxvnorm>0.0)
    {
        for (size_t i=0; i<P.Size(); ++i)
        {
            int xi = _x(P[i](0));
            int yi = _y(P[i](1));
            int xf = _x(P[i](0)+V[i](0)/maxvnorm);
            int yf = _y(P[i](1)+V[i](1)/maxvnorm);
            _draw_arrow (xi, yi, xf, yf, Clr);
        }
    }
}

inline void DrawArea2D::_draw_arrow (int xi, int yi, int xf, int yf, Fl_Color Clr)
{
	fl_color      (Clr);
	fl_line_style (FL_SOLID, 1);
	fl_line       (xi, yi, xf, yf);
	fl_line_style (FL_SOLID, 2);
	fl_line       ((xi+xf)/2, (yi+yf)/2, xf, yf); 
}

inline void DrawArea2D::_draw_loaded_nodes ()
{
	// Check
	if (_g->nLd()<1) return;

	// Boundaries
	double maxfnx = _g->Ld(0)(0);
	double maxfny = _g->Ld(0)(1);
	for (size_t i=0; i<_g->nLd(); ++i)
	{
		if (_g->Ld(i)(0)>maxfnx) maxfnx = _g->Ld(i)(0);
		if (_g->Ld(i)(1)>maxfny) maxfny = _g->Ld(i)(1);
	}

	// Draw
	char buf[80];
	for (size_t i=0; i<_g->nLd(); ++i)
	{
		int lfx = 20*static_cast<int>(fabs(_g->Ld(i)(0)/maxfnx+1));
		int lfy = 20*static_cast<int>(fabs(_g->Ld(i)(1)/maxfny+1));
		if (lfx>0 || lfy>0)
		{
			int xi = _x(_g->LdN(i)(0));
			int yi = _y(_g->LdN(i)(1));
			int xf = (_g->Ld(i)(0)>0 ? xi+lfx : xi-lfx);
			int yf = (_g->Ld(i)(1)>0 ? yi-lfy : yi+lfy);
			_draw_arrow (xi, yi, xf, yf, FL_MAGENTA);
			// Write text
			sprintf(buf, "%g\n%g",_g->Ld(i)(0)*(*_M),_g->Ld(i)(1)*(*_M));
			fl_color (FL_BLACK);
			fl_font  (0, 10);
			fl_draw  (buf, xf, yf, 0, 0, FL_ALIGN_CENTER);
		}
	}
}

inline void DrawArea2D::_draw_loaded_points ()
{
	for (size_t i=0; i<_p->nPoints(); ++i)
	{
		if (_p->HasT(i))
		{
			// coordinates
			int xi = _x(_p->P(i)(0));
			int yi = _y(_p->P(i)(1));
			int xf = _x(_p->P(i)(0)+_p->t(i)(0));
			int yf = _y(_p->P(i)(1)+_p->t(i)(1));
			// icon
			int r = _cp.Psz/2;
			fl_color      (FL_BLACK);
			fl_line_style (FL_SOLID, 1);
			fl_pie        (xi-r, yi-r, 2*r+1, 2*r+1, 0,360);
			// arrow
			_draw_arrow (xi,yi, xf,yf, FL_RED);
			// text
			char buf[256];
			sprintf(buf, "(%.3g,%.3g)", _p->t(i)(0), _p->t(i)(1));
			fl_color (FL_BLACK);
			fl_font  (0, 9);
			fl_draw  (buf, xf, yf, 0, 0, FL_ALIGN_CENTER);
		}
	}
}

inline void DrawArea2D::_draw_disp_points ()
{
	for (size_t i=0; i<_p->nPoints(); ++i)
	{
		if (_p->HasU(i))
		{
			// coordinates
			int xi = _x(_p->P(i)(0));
			int yi = _y(_p->P(i)(1));
			// icon
			int r = _cp.Psz/2+2;
			fl_color      (FL_RED);
			fl_line_style (FL_SOLID, 1);
			fl_pie        (xi-r, yi-r, 2*r+1, 2*r+1, 0,360);
			// text
			char buf[256];
			sprintf(buf, "%g", _p->fu(i)(0));
			fl_color (FL_BLACK);
			fl_font  (0, 8);
			fl_draw  (buf, xi, yi, 0, 0, FL_ALIGN_CENTER);
		}
	}
}

inline void DrawArea2D::_draw_coords()
{
	// Coordinates as a string
	char s[80];
	sprintf(s, "x=%.4g y=%.4g", _X(Fl::event_x()), _Y(Fl::event_y()));
	fl_color (FL_WHITE);
	fl_rectf (x()+w()-151,y()+h()-19,150,18);
	fl_color (FL_BLACK);
	fl_font  (0, 12);
	fl_draw  (s, x()+w()-148, y()+h()-18, 0, 0, FL_ALIGN_LEFT_TOP);
	/*
	// Background
	fl_color (0xefebe700);
	fl_rectf (x()+1,y()+1,150,18);
	// Write text
	fl_color (FL_BLACK);
	fl_font  (0, 12);
	fl_draw  (s, x()+2, y()+1, 0, 0, FL_ALIGN_LEFT_TOP);
	*/
	// Frame
	/*
	fl_color      (FL_BLACK);
	fl_line_style (FL_SOLID, 1);
	fl_rect       (x(), y(), 150, 18);
	*/
}

inline void DrawArea2D::_draw_selected()
{
	// Mat points
	if (_p!=NULL)
	{
		if (_p->nPoints()>0 && _pselp>=0) // erase previous
		{
			int xx = _x(_p->P(_pselp)(0));
			int yy = _y(_p->P(_pselp)(1));
			int r  = _cp.Psz/2;
			// Clear background
			fl_color      (FL_WHITE);
			fl_line_style (FL_SOLID, 3);
			fl_circle     (xx, yy, r);
			// Draw previous
			if (_cp.Typ==CT_POINTS || _cp.Typ==CT_BOTH)
			{
				if (_cp.Pch==1) // Open circle
				{
					fl_color      (FL_BLACK);
					fl_line_style (FL_SOLID, 1);
					fl_circle     (xx, yy, r);
				}
				else if (_cp.Pch==16) // Filled circle
				{
					fl_color      (_cp.Clr);
					fl_line_style (FL_SOLID, 1);
					fl_pie        (xx-r, yy-r, 2*r+1, 2*r+1, 0,360);
				}
			}
		}
		if (_p->nPoints()>0 && _selp>=0) // draw selected
		{
			int r = _cp.Psz/2;
			fl_color      (FL_GREEN);
			fl_line_style (FL_SOLID, 3);
			fl_circle     (_x(_p->P(_selp)(0)), _y(_p->P(_selp)(1)), r);
		}
	}

	// Grid nodes
	if (_g!=NULL)
	{
		if (_g->nNodes()>0 && _pseln>=0) // erase previous
		{
			int xx = _x(_g->X(_pseln)(0));
			int yy = _y(_g->X(_pseln)(1));
			int r  = _cg.Psz/2;
			// Clear background
			fl_color      (FL_WHITE);
			fl_line_style (FL_SOLID, 3);
			fl_circle     (xx, yy, r);
			// Draw pieces of grid lines
			if (_cg.Typ==CT_LINES || _cg.Typ==CT_BOTH)
			{
				fl_color      (_cg.Clr);
				fl_line_style (_cg.Lty, _cg.Lwd);
				fl_line       (xx-r-2, yy, xx+r+2, yy);
				fl_line       (xx, yy-r-2, xx, yy+r+2);
			}
			// Draw previous
			if (_cg.Typ==CT_POINTS || _cg.Typ==CT_BOTH)
			{
				if (_cg.Pch==1) // Open circle
				{
					fl_color      (FL_BLACK);
					fl_line_style (FL_SOLID, 1);
					fl_circle     (xx, yy, r);
				}
				else if (_cg.Pch==16) // Filled circle
				{
					fl_color      (_cg.Clr);
					fl_line_style (FL_SOLID, 1);
					fl_pie        (xx-r, yy-r, 2*r+1, 2*r+1, 0,360);
				}
			}
		}
		if (_g->nNodes()>0 && _seln>=0) // draw selected
		{
			int r = _cg.Psz/2;
			fl_color      (FL_MAGENTA);
			fl_line_style (FL_SOLID, 3);
			fl_circle     (_x(_g->X(_seln)(0)), _y(_g->X(_seln)(1)), r);
		}
	}
}

inline void DrawArea2D::_draw_colorbar()
{
	size_t     nclr = 0;
	Fl_Color * map  = NULL;
	switch (_cmtype)
	{
		case CMT_BWP: { nclr=ClrMap::BWP_NCOLORS; map=ClrMap::BWP; break; }
		case CMT_HOT: { nclr=ClrMap::HOT_NCOLORS; map=ClrMap::HOT; break; }
		case CMT_JET: { nclr=ClrMap::JET_NCOLORS; map=ClrMap::JET; break; }
	}
	if (map==NULL) return;
	int cmbh = _cmbh/nclr;
	fl_line_style (FL_SOLID, 1);
	int xi,yi,xf;//,yf;
	for (size_t i=0; i<nclr; ++i)
	{
		xi = x()+w()-_cmbw-20;
		yi = y()+_ctrh+_cmbh-cmbh*(i+1);
		xf = xi+_cmbw;
		//yf = yi+cmbh;
		int idx=(_cmrinv?nclr-1-i:i);
		fl_color (map[idx]);
		fl_rectf (xi, yi, _cmbw, cmbh);
		fl_color (FL_BLACK);
		fl_line  (xi, yi, xf, yi);
	}
	//fl_line (xi, yf, xf, yf);
}

inline void DrawArea2D::_draw_colorbar_labels(double min, double max)
{
	size_t     nclr = 0;
	Fl_Color * map  = NULL;
	switch (_cmtype)
	{
		case CMT_BWP: { nclr=ClrMap::BWP_NCOLORS; map=ClrMap::BWP; break; }
		case CMT_HOT: { nclr=ClrMap::HOT_NCOLORS; map=ClrMap::HOT; break; }
		case CMT_JET: { nclr=ClrMap::JET_NCOLORS; map=ClrMap::JET; break; }
	}
	if (map==NULL) return;
	int cmbh = _cmbh/nclr;
	fl_line_style (FL_SOLID, 1);
	int xi,yi;//,xf,yf;
	char buf[256];
	for (size_t i=0; i<nclr+1; ++i)
	{
		xi = x()+w()-_cmbw-22;
		yi = y()+_ctrh+_cmbh-cmbh*i;
		//xf = xi+_cmbw;
		//yf = yi+cmbh;
		fl_color (FL_BLACK);
		snprintf (buf, 256, "%g", min+i*(max-min)/static_cast<double>(nclr));
		fl_draw  (buf, xi, yi, 0, 0, FL_ALIGN_RIGHT);
	}
	//fl_line (xi, yf, xf, yf);
}

inline void DrawArea2D::_set()
{
	// Input
	_cp.Psz = atoi(_w_psz->value());
	_cp.Pch = _w_pty->value();
	_cmtype = static_cast<ClrMapType>(_w_cmt->value());
	_field  = static_cast<FieldType> (_w_fld->value());
	_cmrinv = (_w_inv->value()==1 ? true : false);
	_loads  = (_w_ldg->value()==1 ? true : false);
	_veloc  = (_w_vel->value()==1 ? true : false);
	_selec  = (_w_sel->value()==1 ? true : false);
	_cbar   = (_w_cbr->value()==1 ? true : false);

	// Redraw
	redraw   ();
	Fl::wait (0);
}

}; // namespace MPM

#endif // MPM_DRAWAREA2D_H
