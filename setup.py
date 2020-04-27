#!/usr/bin/env python3

from setuptools import setup


"""Setup CLIPcontext"""

setup(
    name='clipcontext',
    version='0.4',
    description='Extract CLIP-seq binding regions with both genomic and transcript context',
    long_description=open('README.md').read(),
    author='Michael Uhl',
    author_email='uhlm@informatik.uni-freiburg.de',
    url='https://github.com/BackofenLab/CLIPcontext',
    license='MIT',
    install_requires=["seaborn>=0.10.0", "matplotlib>=3.1.3", "markdown>=3.2.1", "pandas>=1.0.3"],
    scripts=['bin/clipcontext'],
    packages=['clipcontext'],
    package_data={'clipcontext': ['content/*']},
    zip_safe=False,
)

