#!/bin/bash

MODULES="3dlink cad dict fem gui main mesh shandler"
for m in $MODULES; do
	ln -s $HOME/mechsys/lib/blender/msys_$m.py $HOME/.blender/scripts
done
