import os
from setuptools import setup, find_packages
from pip._internal.req import parse_requirements
from pip._internal.network.session import PipSession

this_dir = os.path.dirname(os.path.abspath(__file__))
pip_requirements = parse_requirements(
    os.path.join(this_dir, "requirements.txt"), PipSession())
pip_requirements_test = parse_requirements(
    os.path.join(this_dir, "requirements-test.txt"), PipSession())

reqs = [pii.requirement for pii in pip_requirements]
reqs_test = [pii.requirement for pii in pip_requirements_test]

readme_path = os.path.join(this_dir, "README.md")

with open(readme_path, "r") as f:
    long_description = f.read()

setup(
    name='htpmd',
    version='0.1.4',
    description='A library to analyze trajectory data from Molecular Dynamics Simulations',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='',
    packages=find_packages(),
    python_requires='>=3.7',
    install_requires=reqs,
    extras_require={"tests": reqs_test},
    entry_points={
        'console_scripts': [
            'htpmd=htpmd.main:main',
        ],
    }
)
