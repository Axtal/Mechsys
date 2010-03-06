#!/bin/bash

set -e

if [ -d "$HOME/mechsys" ]; then
    echo
    echo "Found: $HOME/mechsys ==> OK"
    echo
else
    echo
    echo "Directory named 'mechsys' does not exist"
    echo "Please, download 'mechsys' first:"
    echo
    echo "   hg clone http://hg.savannah.nongnu.org/hgweb/mechsys/"
    echo
    exit 1
fi

TRIANGLE=triangle1.6
TETGEN=tetgen1.4.3
VORO=voro++0.3.1
HDF5=hdf5-1.8.4-patch1

test -d $HOME/pkg || mkdir $HOME/pkg

echo
echo "[1;31mDownloading MTL4   #########################################################################################################[0m"
echo
cd ~/pkg
svn co https://svn.osl.iu.edu/tlc/trunk/mtl4/trunk mtl4

echo
echo "[1;31mDownloading and Compiling Triangle   #######################################################################################[0m"
echo
cd ~/pkg
if [ -d "$HOME/pkg/$TRIANGLE" ]; then
    rm -rf $HOME/pkg/$TRIANGLE
fi
if [ -e $TRIANGLE.tar.gz ]; then
    rm $TRIANGLE.tar.gz;
fi
wget http://mechsys.nongnu.org/software/$TRIANGLE.tar.gz
tar xzvf $TRIANGLE.tar.gz
cd $TRIANGLE
sh ~/mechsys/patches/triangle/do_patch.sh
make

echo
echo "[1;31mDownloading and Compiling Tetgen   #########################################################################################[0m"
echo
cd ~/pkg
if [ -d "$HOME/pkg/$TETGEN" ]; then
    rm -rf $HOME/pkg/$TETGEN
fi
if [ -e $TETGEN.tar.gz ]; then
    rm $TETGEN.tar.gz;
fi
wget http://mechsys.nongnu.org/software/$TETGEN.tar.gz
tar xzvf $TETGEN.tar.gz
cd $TETGEN
sh ~/mechsys/patches/tetgen/do_patch.sh
make

echo
echo "[1;31mDownloading and Patching Voro++    #########################################################################################[0m"
echo
cd ~/pkg
if [ -d "$HOME/pkg/$VORO" ]; then
    rm -rf $HOME/pkg/$VORO
fi
if [ -e $VORO.tar.gz ]; then
    rm $VORO.tar.gz;
fi
wget http://mechsys.nongnu.org/software/$VORO.tar.gz
tar xzvf $VORO.tar.gz
cd $VORO
sh ~/mechsys/patches/voro/do_patch.sh

echo
echo "[1;31mDownloading HDF5    ########################################################################################################[0m"
echo
cd ~/pkg
if [ -d "$HOME/pkg/$HDF5" ]; then
    rm -rf $HOME/pkg/$HDF5
fi
if [ -e $HDF5.tar.gz ]; then
    rm $HDF5.tar.gz;
fi
wget http://www.hdfgroup.org/ftp/HDF5/current/src/$HDF5.tar.gz
tar xzvf $HDF5.tar.gz
cd $HDF5
./configure
make

echo
echo "[1;32mFinished   #################################################################################################################[0m"
echo
