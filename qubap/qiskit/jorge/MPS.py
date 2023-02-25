#!/usr/bin/env python3

from qiskit.circuit import ParameterVector
from qiskit import QuantumCircuit
from .gates import U, W

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
