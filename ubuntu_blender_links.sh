#!/bin/bash

MODULES="3dlink cad dict draw fem"
for m in $MODULES; do
	sudo ln -s $HOME/mechsys/lib/blender/msys_$m.py $HOME/.blender/scripts
done
