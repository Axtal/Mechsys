#!/bin/bash

set -e

sudo apt-get install \
    wget subversion patch \
    g++ gfortran make cmake-curses-gui \
    libgsl0-dev libblitz0-dev libsuitesparse-dev \
    libboost-python-dev python-tk python-numpy python-scipy python-matplotlib \
    libgtk2.0-dev
