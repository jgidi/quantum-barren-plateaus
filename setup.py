#!/usr/bin/env python3

from setuptools import setup, find_packages
setup(
    name = 'qubap',
    packages = find_packages(),
    install_requires = [
        "numpy",
        "qiskit",
        "qiskit_aer",
    ],
)
