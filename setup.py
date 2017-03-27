from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

cmorsesmale_module = Extension('_cpy_morse_smale',
                               include_dirs=[numpy_include, 'src'],
                               libraries=['blas', 'lapack'],
                               sources=['src/ANN.cpp',
                                        'src/bd_fix_rad_search.cpp',
                                        'src/bd_pr_search.cpp',
                                        'src/bd_search.cpp',
                                        'src/bd_tree.cpp',
                                        'src/brute.cpp',
                                        'src/kd_fix_rad_search.cpp',
                                        'src/kd_pr_search.cpp',
                                        'src/kd_search.cpp',
                                        'src/kd_split.cpp',
                                        'src/kd_tree.cpp',
                                        'src/kd_util.cpp',
                                        'src/cpy_morse_smale.cpp',
                                        'src/_cpy_morse_smale.pyx'],
                               language="c++")

setup(name='morsesmale',
      version='0.0.8',
      description='Morse Smale Complex approximation and visualization library providing a topological exploration of '
                  'manifolds',
      setup_requires=['numpy', 'cython'],
      py_modules=['morse_smale'],
      author="Emmett McQuinn",
      ext_modules=cythonize(cmorsesmale_module))
