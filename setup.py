#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 21:21:17 2020

@author: tamertevetoglu
"""
from setuptools import setup

with open("README.md", "r") as fh:
          long_description = fh.read()

setup(
      name='predictr',
      version='0.1.8',
      description='Weibull Analysis Utilities',
      author='Tamer Tevetoglu',
      author_email="tamer.tevetoglu@ima.uni-stuttgart.de",
      url="https://github.com/tvtoglu/predictr",
      py_modules=["predictr"],
      package_dir={'': 'src'},
      long_description=long_description,
      long_description_content_type='text/markdown',
      classifiers=[
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
          "Operating System :: OS Independent",],
      install_requires= [
          "numpy >= 1.16.0",
          "scipy >= 1.3.0",
          "pandas >= 1.0.0",
          ]
      )
