#!/usr/bin/env python3
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.opflow.primitive_ops import PauliSumOp

def parse_hamiltonian(paulis, coefs):
    if isinstance(paulis, str):
        paulis = [paulis]
        coefs  = [coefs]

    return PauliSumOp( SparsePauliOp(paulis, coefs) )

def paulistrings2hamiltonian(pauli_strings, coeffs):
    return PauliSumOp( SparsePauliOp( [Pauli(string) for string in pauli_strings], coeffs ) )
