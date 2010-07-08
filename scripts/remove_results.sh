#!/bin/bash

FILES="*.vtk *.vtu *.cal *.res *.mpy *.stress *.disp"

for f in $FILES; do
	find . -name "$f" -exec rm {} \; > /dev/null 2>&1
done
