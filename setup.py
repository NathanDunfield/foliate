from setuptools import setup

setup(
    name = 'foliar',
    version = '0.2',
    install_requires = ['snappy'],
    packages = ['foliar'],
    package_dir = {'foliar':'src'},
)
