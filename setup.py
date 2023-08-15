from setuptools import find_packages
from setuptools import setup

setup(
    name='FingerprintFinder',
    version='2.4.3',
    description='Takes one pdb, a .smi database, finds and stores a set of similar fingerprints in an excel files',
    author='Ludovico Pipitò',
    author_email='pipitol@uni.coventry.ac.uk',
    python_requires=">=3.6.6",
    packages=find_packages(),
    install_requires=[
        'rdkit',
        'pandas',
        'setuptools'
    ],
)
