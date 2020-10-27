from setuptools import setup, Extension
from setuptools_scm import get_version
import os

setup(name='variableapec',
      version=get_version(root=".", relative_to_"variableapec.py"),
      description='variableapec tool for AtomDB python library.',
      url='http://www.atomdb.org',
      author='Keri Heuer',
      author_email='kh3286@drexel.edu',
      license='Smithsonian',
      classifiers=['Development Status :: 4 - Beta',\
                   'Environment :: Console',\
                   'Intended Audience :: Developers',\
                   'Intended Audience :: Education',\
                   'Intended Audience :: End Users/Desktop',\
                   'Intended Audience :: Science/Research',\
                   'Topic :: Scientific/Engineering :: Astronomy',\
                   'Topic :: Scientific/Engineering :: Physics',\
                   'Programming Language :: Python :: 3',\
                   'Operating System :: POSIX'],
      zip_safe=False,
      long_description = README_TEXT,\
      install_requires=[
      "requests",\
      "wget",\
      "numpy>=1.9.0",\
      "scipy>=1.4.0",\
      "joblib",\
      "mock",\
      "astropy",\
      "pycurl"],
      ext_modules = extmos)
