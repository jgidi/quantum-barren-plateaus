#!/usr/bin/env python3

import numpy as np

import qiskit.opflow as of
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.opflow.primitive_ops import PauliSumOp

from .tools import parse_hamiltonian

def ladder_hamiltonian(num_qubits, transverse_field_intensity=0):
    def interactions(i, j):
        paulis = ['I'] * num_qubits
        for l in [i, j]:
            if l < num_qubits:
                paulis[l] = 'Z'
        return ''.join(paulis)

    # Add a term for each interacting pair
    operators = []
    for n in range(num_qubits):
        # Even n only, with upper bound
        if (not n%2) and (n < num_qubits-1):
            operators.append( interactions(n, n+1) )
        # Even and odd, with upper bound
        if n < num_qubits-2:
            operators.append( interactions(n, n+2) )
    operator_coeffs = [1]*len(operators)

    # If there is magnetic field
    if transverse_field_intensity:
        for i in range(num_qubits):
            paulis = ['I'] * num_qubits
            paulis[i] = 'X'
            operators.append(''.join(paulis))
            operator_coeffs.append(transverse_field_intensity)

    return parse_hamiltonian(operators, operator_coeffs)

def test_hamiltonian( num_qubits, coeff ):
    ops = []

    z = np.zeros(num_qubits)
    x = np.ones(num_qubits)
    ops.append( Pauli((z,x)) )

    z = np.ones(num_qubits)
    ops.append( Pauli((z,x)) )

    x = np.zeros(num_qubits)
    ops.append( Pauli((z,x)) )

    hamiltonian = PauliSumOp( SparsePauliOp( ops, coeff ) )

    return hamiltonian.reduce()

def test_hamiltonian_2( num_qubits ):

    Z = of.Z #Pauli Z
    I = of.I #Identidad
    Zero = 0.5*( I + Z )

    hamiltonian = eval( '('+(num_qubits-1)*'I^'+'I)-('+(num_qubits-1)*'Zero^'+'Zero)' )

    return hamiltonian.reduce()
