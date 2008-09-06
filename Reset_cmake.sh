#!/bin/bash

make clean > /dev/null 2>&1

FILES="CMakeCache.txt cmake_install.cmake Makefile CTestTestfile.cmake"
SFILES="install_manifest.txt"
DIRS="CMakeFiles build"

for f in $FILES; do
	find . -name "$f" -exec rm {} \; > /dev/null 2>&1
done

for f in $SFILES; do
	find . -name "$f" -exec sudo rm {} \; > /dev/null 2>&1
done

for d in $DIRS; do
	find . -name "$d" -exec rm -rf {} \; > /dev/null 2>&1
done
