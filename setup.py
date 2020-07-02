"""
"""

from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
with open(path.join(here, 'requirements.txt')) as f:
    requires = f.read().strip().split('\n')

setup(
    name='htp-md',
    version='0.0.1',
    description='A library to analyze trajectory data from Molecular Dynamics Simulations',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='',
    package_dir={'': 'htp_md'},
    packages=find_packages(where='htp_md'),
    python_requires='>=3.7',
    install_requires=requires,
    entry_points={
        'console_scripts': [
            'htp-md=htp_md.main:main',
        ],
    }
)