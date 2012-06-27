#!/bin/bash

set -e

sudo apt-get install \
    wget subversion patch \
    g++ gfortran make cmake-curses-gui \
    libgsl0-dev libblitz0-dev libsuitesparse-dev \
    libboost-python-dev python-tk python-numpy python-scipy python-matplotlib \
    libgtk2.0-dev libfltk1.3-dev libhdf5-serial-dev libxml2-dev libvtk5-dev \
    mencoder pkg-config
