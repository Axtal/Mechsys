#!/usr/bin/env python
 
from distutils.core import setup
from distutils.extension import Extension

themodule = Extension("mechsys",
                      ["mechsys.cpp"],
                      include_dirs = ['../lib'],
                      libraries    = ['lapack','boost_python'])
 
setup(name        = "MechSys",
      description = "Open library for mechanics",
      version     =  "1.0",
      ext_modules = [themodule])
