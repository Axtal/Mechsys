#!/usr/bin/env python
 
from distutils.core import setup
from distutils.extension import Extension
 
setup(name="mechsys",
    ext_modules=[
        Extension("mechsys", ["mechsys.cpp"],
        libraries = ["boost_python"])
    ])
