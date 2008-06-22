#!/usr/bin/env python
 
from distutils.core import setup
from distutils.extension import Extension

msysfem_mod = Extension("msysfem",
                        ["msysfem.cpp"],
                        include_dirs = ['../lib'],
                        library_dirs = ['.'],
                        libraries    = ['lapack', 'boost_python'])
 
msysmesh_mod = Extension("msysmesh",
                         ["msysmesh.cpp"],
                         include_dirs = ['../lib'],
                         library_dirs = ['.'],
                         libraries    = ['lapack', 'boost_python'])

setup(name        = "MechSys",
      description = "Open library for mechanics",
      version     =  "1.0",
      ext_modules = [msysfem_mod,msysmesh_mod])
