#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension

setup(name="novatel",
  ext_modules=[
    Extension("novatel", sources = ["src/novatel_python.cpp",
                                    "src/novatel.cpp"],
    libraries = ["boost_python-py27", "boost_thread", "serial"],
    include_dirs = ["include/"],
    extra_compile_args=['-std=c++11'])
  ])

