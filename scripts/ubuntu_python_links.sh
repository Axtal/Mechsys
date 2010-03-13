#!/bin/bash

echo 'Considering that MechSys is in ~/mechsys'

PYVER=2.6

MODULES="drawmesh invariants linfit plotter readdata"
for m in $MODULES; do
	ln -s $HOME/mechsys/lib/python/msys_$m.py /usr/lib/python$PYVER/dist-packages
done
