#!/bin/bash

PROGS="tpstrain01 tquad4 tquad8 tsolid01 ttri6 ttruss01"

for p in $PROGS; do
	./$p | grep "Errors"
done
