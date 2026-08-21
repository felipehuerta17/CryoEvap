from setuptools import setup, find_packages

setup(
    name='cryoevap',
    version='0.1',
    packages=find_packages(),
    package_data={'cryoevap': ['cryogens/Coeffs/*.csv']},
    include_package_data=True,
    description='Simulation suite for the evaporation of cryogenic liquids in storage tanks',
    author='Felipe Huerta',
    author_email=' fnhuerta@uc.cl',
)
