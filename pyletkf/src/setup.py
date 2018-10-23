from Cython.Distutils import build_ext
from setuptools import setup, Extension
import numpy as np

ext_modules = [
    Extension('letkf',
              sources=['./letkf/letkf.pyx'],
              include_dirs=[np.get_include()],
              extra_compile_args=['-O3','-fopenmp'],
              extra_link_args=['-fopenmp'])
]

setup(
    name='letkf',
    packages=['letkf'],
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext}
)
