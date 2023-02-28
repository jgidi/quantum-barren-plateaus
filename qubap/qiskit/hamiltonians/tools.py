#!/usr/bin/env python3
from numbers import Number
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.opflow.primitive_ops import PauliOp, PauliSumOp

def parse_hamiltonian(paulis, coefs):
    if isinstance(paulis, str):
        paulis = [paulis]
        coefs  = [coefs]

    return PauliSumOp( SparsePauliOp(paulis, coefs) )

def global2local( hamiltoniano, reduce=True ):

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
            # elif pauli_local == 'I':
            #     continue

            ops_local.append( Pauli((z,x)) )
            coeff_local.append( coeff/num_qubits )

    hamiltoniano_local = PauliSumOp( SparsePauliOp( ops_local, coeff_local ) )

    return hamiltoniano_local.reduce() if reduce else hamiltoniano_local

def paulistrings2hamiltonian(pauli_strings, coeffs):
    return PauliSumOp( SparsePauliOp( [Pauli(string) for string in pauli_strings], coeffs ) )
