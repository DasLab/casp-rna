from setuptools import setup, find_packages

setup(
    name='casp_rna',
    version='0.1',
    packages=find_packages(include=['casp_rna', 'casp_rna.*'])
)