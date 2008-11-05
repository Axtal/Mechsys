#!/bin/bash

FILES="*.vtu *.cal"

for f in $FILES; do
	find . -name "$f" -exec rm {} \; > /dev/null 2>&1
done
