# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
import os


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    # package data
    name='freesif',
    version='0.1.3',
    description='Get data from Sesam Interface Files',
    use_scm_version=True,
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    package_data=dict(),
    python_requires='~=3.7',
    setup_requires=['setuptools_scm'],
    install_requires=[
        'tables>=3.6,<4',
        'numpy>=1.17,<2'
    ],
    zip_safe=True,

    # meta data
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    keywords='sesam structural hydrodynamic',
    url='https://github.com/agrav/freesif',
    author='Audun Gravdal Johansen',
    author_email='audun.gravdal.johansen@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.4',
    ],
)
