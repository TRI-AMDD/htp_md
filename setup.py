import os
from os.path import join
from setuptools import setup, find_packages
from pip._internal.req import parse_requirements
from pip._internal.network.session import PipSession

this_dir = os.path.dirname(os.path.abspath('setup.py'))
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
    version='1.0.0',
    description='A library to analyze trajectory data from Molecular Dynamics Simulations',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/TRI-AMDD/htp_md',
    packages=find_packages(),
    python_requires='>=3.9',
    install_requires=reqs,
    extras_require={
        "tests": reqs_test
    },
    entry_points={
        'console_scripts': [
            'htpmd=htpmd.main:main',
        ],
    },
    data_files=[(join('htpmd', 'ml_models', 'pre_trained_gnns'),
                 [join('htpmd', 'ml_models', 'pre_trained_gnns', 'conductivity.pth'),
                  join('htpmd', 'ml_models', 'pre_trained_gnns', 'li_diffusivity.pth'),
                  join('htpmd', 'ml_models', 'pre_trained_gnns', 'poly_diffusivity.pth'),
                  join('htpmd', 'ml_models', 'pre_trained_gnns', 'tfsi_diffusivity.pth'),
                  join('htpmd', 'ml_models', 'pre_trained_gnns', 'transference_number.pth'),
                  ]),
                (join('htpmd', 'ml_models', 'pre_trained_rfs'),
                 [join('htpmd', 'ml_models', 'pre_trained_rfs', 'rf_conductivity.sav'),
                  join('htpmd', 'ml_models', 'pre_trained_rfs', 'rf_li_diffusivity.sav'),
                  join('htpmd', 'ml_models', 'pre_trained_rfs', 'rf_poly_diffusivity.sav'),
                  join('htpmd', 'ml_models', 'pre_trained_rfs', 'rf_tfsi_diffusivity.sav'),
                  join('htpmd', 'ml_models', 'pre_trained_rfs', 'rf_transference_number.sav'),
                  ])
                ],
)
