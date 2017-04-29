from distutils.core import setup
from os import path
from glob import glob

setup(
    name='SubtractSim',
    version='0.9.0',
    author='David Guevel',
    author_email='guevel.david@gmail.com',
    scripts=glob(path.join('bin/*')),
    license='LICENSE.txt',
    description='Interface to test image subtraction',
    requires=['numpy','astropy','scipy', 'PyZOGY'],
    packages=['SubtractSim'],
)