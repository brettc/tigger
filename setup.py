from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import sys

# HACK for now
sys.argv = ['setup.py', 'build_ext', '--inplace']

extensions = [
    Extension(
        "tigger/_bounce", 
        ["tigger/_bounce.pyx"],
        extra_compile_args = [
            '-Wno-unused-function', 
        ],
        language = 'c++',
    )
]
setup(
    name='partfinder',
    ext_modules=cythonize(
        extensions,
        aliases={ 'NUMPY_PATH': numpy.get_include() },
        nthreads=4,
    ),
)
