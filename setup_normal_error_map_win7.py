#from distutils.core import setup
#from distutils.extension import Extension
from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy as np


ext_modules = [Extension("cython_all_fit_error_map", ["cython_all_fit_error_map.pyx"],language='c++')]
               

setup(
  cmdclass = {'build_ext': build_ext},
  include_dirs = [np.get_include()],
  ext_modules = ext_modules,
)
