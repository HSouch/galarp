
#!/usr/bin/env python
# Licensed under an MIT license - see LICENSE

import os
import sys

from galarp import __version__


from setuptools import setup, find_packages


from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.rst").read_text()


setup(
    name='galarp',
    version=__version__,
    author='Harrison Souchereau',
    author_email='harrison.souchereau@yale.edu',
    description='A ram pressure add-on for Gala numerical integration of gravitational potentials',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/HSouch/galarp',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    keywords='galaxies, gravity, gravitational potentials',
    install_requires=[
        'gala',
        'astropy',
        'numpy',
        'matplotlib'
    ],
)