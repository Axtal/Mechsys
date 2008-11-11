#!/usr/bin/env python
 
from distutils.core import setup
from distutils.extension import Extension

mechsys = Extension("mechsys",
                    ["lib/mesh/jrs_triangle.c", "src/mechsys.cpp"],
                    define_macros      = [('O3','1'), ('TRILIBRARY','1'), ('NDEBUG','1'), ('HAVE_UMFPACK','1')],
                    include_dirs       = ['lib', '/usr/include/suitesparse', '/usr/local/include', '/usr/local/include/suitesparse'],
                    library_dirs       = ['.'],
                    libraries          = ['lapack', 'amd', 'umfpack', 'boost_python', 'igraph', 'CGAL'],
                    extra_compile_args = ['-w'])

setup(name        = "MechSys",
      description = "Open library for mechanics",
      version     =  "1.0",
      ext_modules = [mechsys])
