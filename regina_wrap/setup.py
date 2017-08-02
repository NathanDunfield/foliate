from setuptools import setup, Extension
from Cython.Build import cythonize
import os, shutil

regina_root = '/usr/'
include_dirs = [regina_root + '/include/regina-normal', '/usr/include/libxml2']
extra_link_args = ['-lxml2', '-lgmp', '-lregina-engine']

pyregina_ext = Extension('pyregina',
                         sources=['pyregina.pyx'],
                         language='c++',
                         extra_compile_args=['-std=c++11'],
                         include_dirs=include_dirs,
                         extra_link_args=extra_link_args)

setup(
    name='pyregina',
    version='0.1',
    ext_modules=cythonize(pyregina_ext)
)

stupid_info = 'pyregina.egg-info'
if os.path.exists(stupid_info):
    shutil.rmtree(stupid_info)
