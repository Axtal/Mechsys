#!/bin/bash

MODULES="3dlink cad dict draw fem shandler"
for m in $MODULES; do
	ln -s $HOME/mechsys/lib/blender/msys_$m.py /Applications/blender.app/Contents/MacOS/.blender/scripts
done
