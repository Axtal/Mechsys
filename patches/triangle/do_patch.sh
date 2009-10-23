#!/bin/bash

patch makefile   $HOME/mechsys/patches/triangle/makefile.diff
patch triangle.c $HOME/mechsys/patches/triangle/triangle.c.diff
patch triangle.h $HOME/mechsys/patches/triangle/triangle.h.diff
