#!/bin/bash

set -e

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

echo
echo "****************************************************************************"
echo "* You can call this script with an option to force recompiling everything  *"
echo "* and/or an option to also download packages                               *"
echo "*                                                                          *"
echo "* Example:                                                                 *"
echo "*   sh $MECHSYS_ROOT/mechsys/scripts/do_compile_deps.sh {0,1} {0,1}        *"
echo "*                                                                          *"
echo "* By default, the code will not be compiled if this was done before.       *"
echo "****************************************************************************"

if [ -d "$MECHSYS_ROOT/mechsys" ]; then
    echo
    echo "Found: $MECHSYS_ROOT/mechsys ==> OK"
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

test -d $MECHSYS_ROOT/pkg || mkdir $MECHSYS_ROOT/pkg

TRIANGLE=triangle1.6
TETGEN=tetgen1.4.3
VORO=voro++0.3.1
MTL4=mtl4
SCALAPACK=scalapack_installer
SCALAPACK_DIR=scalapack_installer_0.96
MUMPS=MUMPS_4.9.2
PROC_VER=3.2.8
PROC=procps-$PROC_VER
SPHASH=sparsehash-1.8.1
IGRAPH=igraph-0.5.4
SOPLEX=soplex-1.5.0

compile_scalapack() {
    LDIR=$MECHSYS_ROOT/pkg/$SCALAPACK_DIR/lib
    python setup.py --notesting --mpiincdir=/usr/lib/openmpi/include/ --lapacklib=/usr/lib/liblapack.so --blaslib=/usr/lib/libblas.so
    #python setup.py --notesting --mpiincdir=/home/dorival/pkg/openmpi-1.4.2/ompi/include/ --lapacklib=/usr/lib/liblapack.so --blaslib=/usr/lib/libblas.so
    ln -s $LDIR/blacs.a    $LDIR/libblacs.a
    ln -s $LDIR/blacsC.a   $LDIR/libblacsC.a
    ln -s $LDIR/blacsF77.a $LDIR/libblacsF77.a
}

compile_mumps() {
    cp $MECHSYS_ROOT/mechsys/patches/mumps/Makefile.inc .
    make
}

proc_links() {
    LDIR=$MECHSYS_ROOT/pkg/$PROC/proc
    ln -s $LDIR/libproc-$PROC_VER.so $LDIR/libproc.so
}

download_and_compile() {
    IS_SVN=0
    DO_PATCH=1
    DO_CONF=0
    DO_MAKE=1
    CONF_PRMS=""
    CMD=""
    case "$1" in
        triangle)
            PKG=$TRIANGLE
            PKG_DIR=$PKG
            LOCATION=http://mechsys.nongnu.org/software/$PKG.tar.gz
            ;;
        tetgen)
            PKG=$TETGEN
            PKG_DIR=$PKG
            LOCATION=http://mechsys.nongnu.org/software/$PKG.tar.gz
            ;;
        voro)
            PKG=$VORO
            PKG_DIR=$PKG
            LOCATION=http://mechsys.nongnu.org/software/$PKG.tar.gz
            DO_MAKE=0
            ;;
        mtl4)
            PKG=$MTL4
            PKG_DIR=$PKG
            LOCATION=https://svn.osl.iu.edu/tlc/trunk/mtl4/trunk
            IS_SVN=1
            DO_PATCH=0
            DO_MAKE=0
            ;;
        scalapack)
            PKG=$SCALAPACK
            PKG_DIR=$SCALAPACK_DIR
            LOCATION=http://www.netlib.org/scalapack/$PKG.tgz
            DO_PATCH=0
            DO_MAKE=0
            CMD=compile_scalapack
            ;;
        mumps)
            PKG=$MUMPS
            PKG_DIR=$PKG
            LOCATION=""
            DO_PATCH=0
            DO_MAKE=0
            CMD=compile_mumps
            ;;
        proc)
            PKG=$PROC
            PKG_DIR=$PROC
            LOCATION=http://procps.sourceforge.net/$PKG.tar.gz
            DO_PATCH=0
            CMD=proc_links
            ;;
        sphash)
            PKG=$SPHASH
            PKG_DIR=$PKG
            LOCATION=http://google-sparsehash.googlecode.com/files/$PKG.tar.gz
            DO_PATCH=0
            DO_CONF=1
            DO_MAKE=1
            ;;
        igraph)
            PKG=$IGRAPH
            PKG_DIR=$PKG
            LOCATION=http://sourceforge.net/projects/igraph/files/C%20library/0.5.4/igraph-0.5.4.tar.gz
            DO_PATCH=0
            DO_MAKE=1
            DO_CONF=1
            ;;
        soplex)
            PKG=$SOPLEX
            PKG_DIR=$PKG
            LOCATION=http://soplex.zib.de/download/soplex-1.5.0.tgz
            DO_PATCH=0
            DO_MAKE=1
            ;;
        *)
            echo
            echo "download_and_compile_tar_gz: __Internal_error__"
            exit 1
            ;;
    esac
    echo
    echo "********************************** ${1} ********************************"

    # erase existent directory ?
    cd $MECHSYS_ROOT/pkg
    if [ -d "$MECHSYS_ROOT/pkg/$PKG_DIR" ]; then
        if [ "$RECOMPILE" -eq 1 ]; then
            echo "    Recompiling $PKG"
            if [ "$IS_SVN" -eq 1 ]; then
                echo "    Updating $PKG SVN repository"
                svn co $LOCATION $PKG
            else
                rm -rf $MECHSYS_ROOT/pkg/$PKG_DIR
            fi
        else
            echo "    Using existing $PKG_DIR"
            return
        fi

    # download SVN package
    else
        echo "    New compilation of $PKG"
        if [ "$IS_SVN" -eq 1 ]; then
            echo "    Downloading $PKG SVN repository"
            svn co $LOCATION $PKG
        fi
    fi

    # download and uncompressing package
    if [ "$IS_SVN" -eq 0 ]; then
        if [ -e $PKG.tar.gz ]; then
            echo "    Using local <$PKG.tar.gz>"
        else
            if [ -e $PKG.tgz ]; then
                echo "    Using local <$PKG.tgz>"
            else
                echo "    Downloading package"
                if [ -z "$LOCATION" ]; then
                    echo
                    echo
                    echo "    [1;31mPlease, download <$PKG.tar.gz> first[0m"
                    echo
                    echo
                    return
                else
                    wget $LOCATION
                fi
            fi
        fi
        echo "        . . . uncompressing . . ."
        if [ -e $PKG.tar.gz ]; then
            tar xzf $PKG.tar.gz
        else
            tar xzf $PKG.tgz
        fi
    fi

    # patch
    cd $PKG_DIR
    if [ "$DO_PATCH" -eq 1 ]; then
        echo "        . . . patching . . ."
        sh $MECHSYS_ROOT/mechsys/patches/${1}/do_patch.sh
    fi

    # configure
    if [ "$DO_CONF" -eq 1 ]; then
        echo "        . . . configuring . . ."
        ./configure $CONF_PRMS 2> /dev/null
    fi

    # compilation
    if [ "$DO_MAKE" -eq 1 ]; then
        echo "        . . . compiling . . ."
        make > /dev/null 2> /dev/null
    fi

    # given command
    if [ ! -z "$CMD" ]; then
        echo "        . . . command . . . . . ."
        $CMD
    fi
}

download_and_compile triangle
download_and_compile tetgen
download_and_compile voro
download_and_compile mtl4
download_and_compile scalapack
download_and_compile mumps
download_and_compile proc
download_and_compile sphash
download_and_compile igraph
download_and_compile soplex

echo
echo "Finished ###################################################################"
echo
