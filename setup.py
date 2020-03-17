#!/usr/bin/env python3

from setuptools import setup


"""Setup CLIPcontext"""

setup(
    name='clipcontext',
    version='0.1',
    description='Extract CLIP-seq binding regions with both genomic and transcript context',
    long_description=open('README.md').read(),
    author='Michael Uhl',
    author_email='uhlm@informatik.uni-freiburg.de',
    url='https://github.com/BackofenLab/CLIPcontext',
    license='MIT',
    scripts=['bin/clipcontext'],
    packages=['clipcontext'],
    zip_safe=False,
)

