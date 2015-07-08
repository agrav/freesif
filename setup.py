# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
from codecs import open  # To use a consistent encoding
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.txt'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='freesif',
    version='0.1',
    description='Get data from Sesam Interface Files',
    long_description=long_description,
#    url='https://github.com/agrav/freesif',
    author='Audun Gravdal Johansen',
    author_email='audun.gravdal.johansen@gmail.com',
    license='MIT',
    classifiers=[
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
#        'Programming Language :: Python :: 3.4',
    ],
    keywords='sesam structural hydrodynamic',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=['tables', 'numpy'],
)
