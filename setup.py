# -*- coding: utf-8 -*-
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy as np

setup(
    name='bio',
    author='Michael Sanders',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy>=1.9.1',
        'cython>=0.21.1',
    ],
    ext_modules=cythonize("csequencing.pyx"),
    include_dirs=[np.get_include()]
)
