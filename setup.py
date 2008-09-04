#!/usr/bin/env python
 
from distutils.core import setup
from distutils.extension import Extension

mechsys = Extension("mechsys",
                    ["src/mechsys.cpp"],
                    include_dirs = ['lib', '/usr/local/include'],
                    library_dirs = ['.'],
                    libraries    = ['lapack', 'boost_python', 'igraph'])

setup(name        = "MechSys",
      description = "Open library for mechanics",
      version     =  "1.0",
      ext_modules = [mechsys])
