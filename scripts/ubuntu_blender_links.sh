#!/bin/bash

echo 'Considering that MechSys is in ~/mechsys'

MODULES="3dlink cad dict fem dem gui main mesh shandler mexpt mex3dlink"
for m in $MODULES; do
	ln -s $HOME/mechsys/lib/blender/msys_$m.py $HOME/.blender/scripts
done
