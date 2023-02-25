#!/usr/bin/env python3

from qiskit import QuantumCircuit

def _GateW(params):
    qc = QuantumCircuit(2, name="W")
    qc.ry( params[0], qc.qubits[0])
    qc.ry( params[1], qc.qubits[1])
    qc.crx(params[2], qc.qubits[1], qc.qubits[0])

    return qc.to_gate()

def _GateU(params):

    qc = QuantumCircuit(2, name="U")
    qc.u(params[0], params[1], params[2], qc.qubits[0])
    qc.u(params[3], params[4], params[5], qc.qubits[1])

    qc.rxx(params[6], qc.qubits[0], qc.qubits[1])
    qc.ryy(params[7], qc.qubits[0], qc.qubits[1])
    qc.rzz(params[8], qc.qubits[0], qc.qubits[1])

    qc.u(params[ 9], params[10], params[11], qc.qubits[0])
    qc.u(params[12], params[13], params[14], qc.qubits[1])

    return qc.to_gate()

def W(self, params, qubits):
    self.append(_GateW(params), qubits)

def U(self, params, qubits):
    self.append(_GateU(params), qubits)

# TODO keep? To make U/W available as QC methods
QuantumCircuit.W = W
QuantumCircuit.U = U
