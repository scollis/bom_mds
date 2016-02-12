#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


pyradar_module = Extension('_pyradar',
                           sources=['pyradar_wrap.c', 'pyradar.c'], include_dirs=['./include/'], library_dirs=['./lib/'], libraries=['rsl']
                           )

setup (name = 'pyradar',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       ext_modules = [pyradar_module],
       py_modules = ["pyradar"],
       )
