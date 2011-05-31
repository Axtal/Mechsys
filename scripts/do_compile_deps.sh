#!/bin/bash

set -e

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

echo
echo "****************************************************************************"
echo "* You can call this script with an option to force recompiling everything  *"
echo "* and/or an option to also download packages                               *"
echo "*                                                  recompile     download  *"
echo "*                                                          |     |         *"
echo "* Example:                                                 V     V         *"
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

FORCEDOWNLOAD=0
if [ "$#" -gt 1 ]; then
    FORCEDOWNLOAD=$2
    if [ "$FORCEDOWNLOAD" -lt 0 -o "$FORCEDOWNLOAD" -gt 1 ]; then
        echo
        echo "The option for downloading and compilation of additional packages must be either 1 or 0. ($2 is invalid)"
        echo
        exit 1
    fi
fi

test -d $MECHSYS_ROOT/pkg || mkdir $MECHSYS_ROOT/pkg

VER_TRIANGLE=1.6
VER_TETGEN=1.4.3
VER_VORO=0.3.1
VER_SCALAPACK=0.96
VER_MUMPS=4.10.0
VER_PROC=3.2.8
VER_SPHASH=1.8.1
VER_IGRAPH=0.5.4
VER_SOPLEX=1.5.0

INC_OPENMPI=/usr/lib/openmpi/include/  # $MECHSYS_ROOT/pkg/openmpi-1.4.2
LIB_LAPACK=/usr/lib/liblapack.so
LIB_BLAS=/usr/lib/libblas.so

compile_scalapack() {
    LDIR=$MECHSYS_ROOT/pkg/scalapack_installer_$VER_SCALAPACK/lib
    python setup.py --notesting --mpiincdir=$INC_OPENMPI --lapacklib=$LIB_LAPACK --blaslib=$LIB_BLAS
    ln -s $LDIR/blacs.a    $LDIR/libblacs.a
    ln -s $LDIR/blacsC.a   $LDIR/libblacsC.a
    ln -s $LDIR/blacsF77.a $LDIR/libblacsF77.a
}

compile_mumps() {
    cp $MECHSYS_ROOT/mechsys/patches/mumps/Makefile.inc .
    make clean
    make
}

proc_links() {
    LDIR=$MECHSYS_ROOT/pkg/procps-$VER_PROC/proc
    ln -s $LDIR/libproc-$VER_PROC.so $LDIR/libproc.so
}

error_message() {
    echo
    echo
    echo "    [1;31m Error: $1 [0m"
    echo
    echo
}

download_and_compile() {
    PKG=""
    PKG_DIR=""
    EXT=tar.gz
    LOCATION=""
    EXTRA_CMD=""
    CONF_PRMS=""
    IS_SVN=0
    DO_PATCH=0
    DO_CONF=0
    DO_MAKE=1
    case "$1" in
        triangle)
            PKG=triangle$VER_TRIANGLE
            LOCATION=http://mechsys.nongnu.org/software/$PKG.$EXT
            DO_PATCH=1
            ;;
        tetgen)
            PKG=tetgen$VER_TETGEN
            LOCATION=http://mechsys.nongnu.org/software/$PKG.$EXT
            DO_PATCH=1
            ;;
        voro)
            PKG=voro++$VER_VORO
            LOCATION=http://mechsys.nongnu.org/software/$PKG.$EXT
            DO_PATCH=1
            DO_MAKE=0
            ;;
        mtl4)
            PKG=mtl4
            LOCATION=https://svn.osl.iu.edu/tlc/trunk/mtl4/trunk
            IS_SVN=1
            DO_MAKE=0
            ;;
        scalapack)
            PKG=scalapack_installer
            PKG_DIR=scalapack_installer_$VER_SCALAPACK
            EXT=tgz
            LOCATION=http://www.netlib.org/scalapack/$PKG.$EXT
            DO_MAKE=0
            EXTRA_CMD=compile_scalapack
            ;;
        mumps)
            PKG=MUMPS_$VER_MUMPS
            #LOCATION=""
            LOCATION=http://mumps.enseeiht.fr/$PKG.$EXT
            DO_MAKE=0
            EXTRA_CMD=compile_mumps
            ;;
        proc)
            PKG=procps-$VER_PROC
            LOCATION=http://procps.sourceforge.net/$PKG.$EXT
            EXTRA_CMD=proc_links
            ;;
        sphash)
            PKG=sparsehash-$VER_SPHASH
            LOCATION=http://google-sparsehash.googlecode.com/files/$PKG.$EXT
            DO_CONF=1
            ;;
        igraph)
            PKG=igraph-$VER_IGRAPH
            LOCATION=http://sourceforge.net/projects/igraph/files/C%20library/$VER_IGRAPH/$PKG.$EXT
            DO_CONF=1
            ;;
        soplex)
            PKG=soplex-$VER_SOPLEX
            EXT=tgz
            LOCATION=http://soplex.zib.de/download/$PKG.$EXT
            ;;
        *)
            error_message "download_and_compile: __Internal_error__"
            exit 1
            ;;
    esac
    echo
    echo "********************************** ${1} ********************************"

    # change into the packages directory
    cd $MECHSYS_ROOT/pkg

    # package filename and directory
    PKG_FILENAME=$PKG.$EXT
    if [ -z "$PKG_DIR" ]; then PKG_DIR=$PKG; fi

    # check for package that must be existing (cannot be downloaded)
    if [ -z "$LOCATION" ]; then
        if [ ! -e "$PKG_FILENAME" ]; then
            error_message "Please download <$PKG_FILENAME> first"
            return
        fi
    fi

    # (re)compile or return (erasing existing package) ?
    if [ "$IS_SVN" -eq 0 ]; then
        if [ -d "$MECHSYS_ROOT/pkg/$PKG_DIR" ]; then
            if [ "$RECOMPILE" -eq 1   -o   "$FORCEDOWNLOAD" -eq 1 ]; then
                echo "    Erasing existing $PKG_DIR"
                rm -rf $MECHSYS_ROOT/pkg/$PKG_DIR
            else
                echo "    Using existing $PKG_DIR"
                return
            fi
        fi
    else
        if [ -d "$MECHSYS_ROOT/pkg/$PKG_DIR" ]; then
            if [ "$RECOMPILE" -eq 1   -o   "$FORCEDOWNLOAD" -eq 1 ]; then
                echo "    Updating existing $PKG SVN repository"
                cd $PKG_DIR
                svn up
                cd $MECHSYS_ROOT/pkg
            else
                echo "    Using existing $PKG SVN repository in $PKG_DIR"
                return
            fi
        fi
    fi

    # download package
    if [ "$IS_SVN" -eq 0 ]; then
        if [ "$FORCEDOWNLOAD" -eq 1   -o   ! -e "$PKG_FILENAME" ]; then
            if [ -e "$PKG_FILENAME" ]; then
                echo "    Removing existing <$PKG_FILENAME>"
                rm $PKG_FILENAME
            fi
            echo "    Downloading <$PKG_FILENAME>"
            wget $LOCATION
        fi
    else
        if [ ! -d "$MECHSYS_ROOT/pkg/$PKG_DIR" ]; then
            echo "    Downloading new $PKG SVN repository"
            svn co $LOCATION $PKG
        fi
    fi

    # uncompress package
    if [ "$IS_SVN" -eq 0 ]; then
        echo "        . . . uncompressing . . ."
        tar xzf $PKG_FILENAME
    fi

    # change into the package directory
    cd $PKG_DIR

    # patch
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

    # execute specific command
    if [ ! -z "$EXTRA_CMD" ]; then
        echo "        . . . command . . . . . ."
        $EXTRA_CMD
    fi
    echo "        . . . finished . . . . . "
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
