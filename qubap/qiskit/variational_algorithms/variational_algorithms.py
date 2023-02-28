#!/usr/bin/env python3

import qiskit.opflow as of
from qiskit.algorithms import NumPyMinimumEigensolver

from qubap.qiskit.jorge.tools import SPSA_calibrated
from qubap.qiskit.jorge.tools import make_data_and_callback

def energy_evaluation(hamiltonian, ansatz, parameters , quantum_instance, callback=None):
    """
    Evaluate the energy given an ansatz and a Hamiltonian

    Input:
    hamiltonian (PauliSumOp): Hamiltonian of the system
    ansatz (QuantumCircuit): 
    initial_guess (ndarray):
    quantum_instance (QuantumInstance):

    Output:
        (dict): 
    """
    ansatz_state = of.StateFn( ansatz.bind_parameters(parameters) )  
    measurement = of.StateFn( hamiltonian ).adjoint() @ ansatz_state
    pauli_circs = of.PauliExpectation().convert(measurement)
    sampler = of.CircuitSampler(quantum_instance).convert(pauli_circs)
    evaluation = sampler.eval().real
    if callback is not None:
        callback( parameters, evaluation )
    return evaluation
