import numpy as np
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.opflow.primitive_ops import PauliSumOp, PauliOp
from qiskit.algorithms.optimizers import SPSA
from qubap.qiskit.luciano.variational_algorithms import VQE, energy_evaluation
from qubap.qiskit.jorge.tools import make_data_and_callback
import qiskit.opflow as of


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
