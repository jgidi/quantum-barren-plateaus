#!/usr/bin/env python3
import numpy as np

import qiskit.opflow as of
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.opflow.primitive_ops import PauliSumOp, PauliOp

from .tools import parse_hamiltonian

def global2local( hamiltoniano, reduce=True ):
    """
    Descriptiom

    Input:
        hamiltonian (PauliSumOp):
        reduce(, optional):

    Output:
        ():
    """
    num_qubits = hamiltoniano.num_qubits

    ops_local   = []
    coeff_local = []

    paulis_global = hamiltoniano.to_pauli_op()
    if isinstance( paulis_global, PauliOp ):
        paulis_global = [paulis_global]

    for pauli_global, coeff in zip( paulis_global, hamiltoniano.coeffs ):
        pauli_label = pauli_global.primitive.to_label()

        for qb in range(num_qubits):
            pauli_local = pauli_label[qb]
            x = np.zeros(num_qubits)
            z = np.zeros(num_qubits)

            if pauli_local == 'X':
                x[qb] = 1
            elif pauli_local == 'Z':
                z[qb] = 1
            elif pauli_local == 'Y':
                x[qb] = 1
                z[qb] = 1

            ops_local.append( Pauli((z,x)) )
            coeff_local.append( coeff/num_qubits )

    hamiltoniano_local = PauliSumOp( SparsePauliOp( ops_local, coeff_local ) )

    return hamiltoniano_local.reduce() if reduce else hamiltoniano_local


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
    """
    Descriptiom

    Input:
        num_qubits (int): numbers of qubits
        reduce(, optional):

    Output:
        ():
    """
    Z = of.Z #Pauli Z
    I = of.I #Identidad
    Zero = 0.5*( I + Z )
    hamiltonian = eval( '('+(num_qubits-1)*'I^'+'I)-('+(num_qubits-1)*'Zero^'+'Zero)' )

    return hamiltonian.reduce()
