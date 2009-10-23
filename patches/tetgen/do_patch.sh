#!/bin/bash

patch makefile         $HOME/mechsys/patches/tetgen/makefile.diff
patch tetgen1.4.3.cxx  $HOME/mechsys/patches/tetgen/tetgen1.4.3.cxx.diff
patch tetgen1.4.3.h    $HOME/mechsys/patches/tetgen/tetgen1.4.3.h.diff
