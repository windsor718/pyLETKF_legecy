from setuptools import setup
import argparse

mklDir = /opt/intel/mkl

# 20181015 need fortran compilation discription

setup(
    name='pyletkf',
    version='2.0.3',
    description='Python-Core-API for the Local Transformed Ensemble Kalman Filter',
    long_description="README.md",
    author='Yuta Ishitsuka',
    author_email='winzer718@gmail.com',
    install_requires=["numpy"],
    url='https://github.com/windsor718/pyLETKF'
)
