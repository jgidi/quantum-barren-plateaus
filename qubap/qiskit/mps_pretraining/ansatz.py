#!/usr/bin/env python3

import numpy as np

from qiskit.circuit import ParameterVector
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

from .gates import U, W
from qubap.qiskit.luciano.variational_algorithms import VQE

def Ansatz(num_qubits,
           diag_gate=U,
           offdiag_gate=W,
           diag_num_params=15,
           offdiag_num_params=3,
           # reps = 1
           diagonal=True):

    qc = QuantumCircuit(num_qubits)

    # Starting parameter indices per gate
    diag_start = 0
    offdiag_start = (num_qubits - 1)*diag_num_params # Total params in the diag

    num_params = offdiag_start + offdiag_num_params * num_qubits * num_qubits//2
    params = ParameterVector("θ", num_params)
    reps = 1                                  # TODO
    for _ in range(reps):
        for i in range(num_qubits-1):
            # Below diagonal
            if not diagonal:
                for j in range(i%2, i, 2):
                    offdiag_end = offdiag_start + offdiag_num_params
                    offdiag_gate(qc, params[offdiag_start:offdiag_end], qc.qubits[j:j+2])
                    offdiag_start = offdiag_end

            # Diagonal
            diag_end = diag_start + diag_num_params
            diag_gate(qc, params[diag_start:diag_end], qc.qubits[i:i+2])
            diag_start = diag_end

            # Above diagonal
            if not diagonal:
                for j in range(i+2, num_qubits-1, 2):
                    offdiag_end = offdiag_start + offdiag_num_params
                    offdiag_gate(qc, params[offdiag_start:offdiag_end], qc.qubits[j:j+2])
                    offdiag_start = offdiag_end

        if not diagonal:
            for j in range(1, num_qubits-1, 2):
                offdiag_end = offdiag_start + offdiag_num_params
                offdiag_gate(qc, params[offdiag_start:offdiag_end], qc.qubits[j:j+2])
                offdiag_start = offdiag_end

    return qc

def VQE_pretrained(hamiltonian, quantum_instance, iters_vqe, iters_train, returns=['x', 'fx']):

    returns_mps = np.unique(np.append(['x'], returns))

    qc_mps = Ansatz(hamiltonian.num_qubits, diagonal=True)
    backend_mps = AerSimulator(method='matrix_product_state',
                               matrix_product_state_max_bond_dimension=2,
                               shots=2**13)
    guess_mps = np.random.rand(qc_mps.num_parameters) * np.pi
    results_mps = VQE(hamiltonian, qc_mps, guess_mps, iters_train, backend_mps, returns=returns_mps)
    params_mps = np.mean(results_mps['x'][-10:], axis=0)

    qc_full = Ansatz(hamiltonian.num_qubits, diagonal=False)
    num_params_full = qc_full.num_parameters

    guess_full = np.append(params_mps, np.zeros(num_params_full - len(params_mps)))
    results_full = VQE(hamiltonian, qc_full, guess_full, iters_vqe, quantum_instance,
                       iter_start=iters_train, returns=returns)

    return results_full, results_mps