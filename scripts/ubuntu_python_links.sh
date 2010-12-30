#!/bin/bash

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

echo Considering that MechSys is in $MECHSYS_ROOT

PYVER=2.6

MODULES="drawmesh invariants linfit plotter fig matvec fcrits"
for m in $MODULES; do
  sudo ln -s $MECHSYS_ROOT/mechsys/lib/python/msys_$m.py /usr/lib/python$PYVER/dist-packages
done
