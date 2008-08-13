#!/bin/bash

PROGS="tpstrain01 tpstrain02 theat01 tquad4 tquad8 tsolid01 ttri6 ttruss01 tmesh01 tmesh02 ttree01"
for p in $PROGS; do
	./$p | grep "Errors"
done

PYS="tpstrain02.py ttri6.py ttruss01.py tmesh01.py"
for p in $PYS; do
	python $p | grep "Errors"
done
