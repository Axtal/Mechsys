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

#ifndef MPM_INFOBOX_H
#define MPM_INFOBOX_H

#include <FL/Fl.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Double_Window.H>

namespace MPM {

void HideInfoBox (Fl_Widget * IB, void * Win) { ((Fl_Window*)Win)->hide(); }

void InfoBox (char const * Message)
{
	Fl_Window * w = new Fl_Window (420,103,"Error");
	w->begin();
		Fl_Box    * i = new Fl_Box   (10,10, 50,50,"!");
		Fl_Box    * o = new Fl_Box   (70,25,350,20,Message);
		Fl_Button * b = new Fl_Button(310,70,90,23,"Ok");
	w->end();
	i->box        (FL_THIN_UP_BOX);
	i->labelfont  (FL_TIMES_BOLD);
	i->labelsize  (34);
	i->color      (FL_WHITE);
	i->labelcolor (FL_BLUE);
	o->align      (FL_ALIGN_LEFT|FL_ALIGN_INSIDE|FL_ALIGN_WRAP);
	b->callback   (&HideInfoBox, w);
	b->color      (0xefebe700);
	w->color      (0xefebe700);
	w->set_modal  ();
	w->show       ();
}

}; // namespace MPM

#endif // MPM_INFOBOX_H
