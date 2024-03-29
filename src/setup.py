# Setup for building C extensions.

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("dycom", ["dycom.pyx"])]

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
