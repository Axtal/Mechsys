#!/bin/bash

set -e

echo
echo "****************************************************************************"
echo "* You can call this script with an option to force recompiling everything  *"
echo "* and/or an option to also download packages                               *"
echo "*                                                                          *"
echo "* Example:                                                                 *"
echo "*   sh ~/mechsys/scripts/do_compile_deps.sh {0,1} {0,1}                    *"
echo "*                                                                          *"
echo "* By default, the code will not be compiled if this was done before.       *"
echo "****************************************************************************"

if [ -d "$HOME/mechsys" ]; then
    echo
    echo "Found: $HOME/mechsys ==> OK"
else
    echo
    echo "Directory named 'mechsys' does not exist"
    echo "Please, download 'mechsys' first:"
    echo
    echo "   hg clone http://hg.savannah.nongnu.org/hgweb/mechsys/"
    echo
    exit 1
fi

RECOMPILE=0
if [ "$#" -gt 0 ]; then
    RECOMPILE=$1
    if [ "$RECOMPILE" -lt 0 -o "$RECOMPILE" -gt 1 ]; then
        echo
        echo "The option for re-compilation must be either 1 or 0. ($1 is invalid)"
        echo
        exit 1
    fi
fi

FULL=0
if [ "$#" -gt 1 ]; then
    FULL=$2
    if [ "$FULL" -lt 0 -o "$FULL" -gt 1 ]; then
        echo
        echo "The option for downloading and compilation of additional packages must be either 1 or 0. ($2 is invalid)"
        echo
        exit 1
    fi
fi

TRIANGLE=triangle1.6
TETGEN=tetgen1.4.3
VORO=voro++0.3.1
HDF5=hdf5-1.8.4-patch1
MTL4=mtl4

test -d $HOME/pkg || mkdir $HOME/pkg

download_and_compile() {
    IS_SVN=0
    DO_PATCH=1
    DO_MAKE=1
    case "$1" in
        triangle)
            PKG=$TRIANGLE
            LOCATION=http://mechsys.nongnu.org/software/$PKG.tar.gz
            ;;
        tetgen)
            PKG=$TETGEN
            LOCATION=http://mechsys.nongnu.org/software/$PKG.tar.gz
            ;;
        voro)
            PKG=$VORO
            LOCATION=http://mechsys.nongnu.org/software/$PKG.tar.gz
            DO_MAKE=0
            ;;
        hdf5)
            PKG=$HDF5
            LOCATION=http://www.hdfgroup.org/ftp/HDF5/current/src/$HDF5.tar.gz
            DO_PATCH=0
            ;;
        mtl4)
            PKG=$MTL4
            LOCATION=https://svn.osl.iu.edu/tlc/trunk/mtl4/trunk
            IS_SVN=1
            DO_PATCH=0
            DO_MAKE=0
            ;;
        *)
            echo
            echo "download_and_compile_tar_gz: __Internal_error__"
            exit 1
            ;;
    esac
    echo
    echo "********************************** ${1} ********************************"
    cd ~/pkg
    if [ -d "$HOME/pkg/$PKG" ]; then
        if [ "$RECOMPILE" -eq 1 ]; then
            echo "    Recompiling $PKG"
            if [ "$IS_SVN" -eq 1 ]; then
                echo "    Updating $PKG SVN repository"
                svn co $LOCATION $PKG
            else
                rm -rf $HOME/pkg/$PKG
            fi
        else
            echo "    Using existing $PKG"
            return
        fi
    else
        echo "    New compilation of $PKG"
        if [ "$IS_SVN" -eq 1 ]; then
            echo "    Downloading $PKG SVN repository"
            svn co $LOCATION $PKG
        fi
    fi
    if [ "$IS_SVN" -eq 0 ]; then
        if [ -e $PKG.tar.gz ]; then
            echo "    Using local <$PKG.tar.gz>"
        else
            echo "    Downloading $PKG.tar.gz"
            wget $LOCATION
        fi
        echo "        . . . uncompressing . . ."
        tar xzf $PKG.tar.gz
    fi
    cd $PKG
    if [ "$DO_PATCH" -eq 1 ]; then
        echo "        . . . patching . . ."
        sh ~/mechsys/patches/${1}/do_patch.sh
    fi
    if [ "$DO_MAKE" -eq 1 ]; then
        echo "        . . . compiling . . ."
        make > /dev/null 2> /dev/null
    fi
}

download_and_compile triangle
download_and_compile tetgen
download_and_compile voro
download_and_compile hdf5
download_and_compile mtl4

echo
echo "Finished ###################################################################"
echo
