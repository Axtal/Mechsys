#!/bin/bash

patch makefile    $HOME/mechsys/patches/tetgen/makefile.diff
patch tetgen.cxx  $HOME/mechsys/patches/tetgen/tetgen1.4.3.cxx.diff
patch tetgen.h    $HOME/mechsys/patches/tetgen/tetgen1.4.3.h.diff
