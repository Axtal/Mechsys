#!/bin/bash

make clean > /dev/null 2>&1

FILES="CMakeCache.txt cmake_install.cmake Makefile CTestTestfile.cmake"
SFILES="install_manifest.txt"
DIRS="CMakeFiles build Testing"

echo "Removing CMake files"
for f in $FILES; do
	find . -name "$f" -exec rm {} \; > /dev/null 2>&1
done

echo "Removing CMake directories"
for d in $DIRS; do
	find . -name "$d" -exec rm -rf {} \; > /dev/null 2>&1
done

echo "Removing installation log"
for f in $SFILES; do
	find . -name "$f" -exec sudo rm {} \; > /dev/null 2>&1
done

echo "Removing all output files (vtk, vtu, pvd, cal, mpy, ply, and draw)"
find . -iname "*.pvd" -exec rm {} \;
find . -iname "*.vtk" -exec rm {} \;
find . -iname "*.vtu" -exec rm {} \;
find . -iname "*.cal" -exec rm {} \;
find . -iname "*.mpy" -exec rm {} \;
find . -iname "*.ply" -exec rm {} \;
find . -iname "*.draw" -exec rm {} \;

echo "Removing all temporary files (pyc)"
find . -iname "*.pyc" -exec rm {} \;

exit 0
