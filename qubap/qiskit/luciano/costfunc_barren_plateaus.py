import qiskit.opflow as of
import numpy as np
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.opflow.primitive_ops import  PauliSumOp

def global2local( hamiltoniano, reduce=True ):

    num_qubits = hamiltoniano.num_qubits

    ops_local   = []
    coeff_local = []

    for pauli_global, coeff in zip( hamiltoniano.to_pauli_op(), hamiltoniano.coeffs ):
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
            coeff_local.append( coeff )

    hamiltoniano_local = PauliSumOp( SparsePauliOp( ops_local, coeff_local ) )

    return hamiltoniano_local.reduce() if reduce else hamiltoniano_local
