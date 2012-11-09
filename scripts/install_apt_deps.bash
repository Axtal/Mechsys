#!/bin/bash

set -e

sudo apt-get install \
    wget patch \
    g++ gfortran make cmake-curses-gui \
    libgsl0-dev libsuitesparse-dev \
    libboost-python-dev \
    python-tk python-numpy python-scipy python-matplotlib \
    libxml2-dev \
    libmumps-dev libparmetis-dev libvtk5-dev \
    libgtk2.0-dev libfltk1.3-dev libhdf5-serial-dev libxml2-dev \
    libcgal-dev

# note: libxml2-dev is for igraph
#       libmumps-dev will install libopenmpi-dev

# subversion 
# libblitz0-dev
# libgtk2.0-dev libfltk1.3-dev libhdf5-serial-dev libxml2-dev libvtk5-dev \
# mencoder pkg-config
