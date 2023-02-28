#!/usr/bin/env python3

from .tools import energy_evaluation, make_adiabatic_cost_and_callback
from qubap.qiskit.tools import make_data_and_callback, SPSA_calibrated
from qubap.qiskit.hamiltonians.tools import global2local
from qiskit.algorithms.optimizers import SPSA

def VQE(hamiltonian, ansatz, initial_guess, num_iters, quantum_instance,
        returns=['x', 'fx'], iter_start=0):
    """
    Standard VQE

    Input:
        hamiltonian (PauliSumOp):
        ansatz (QuantumCircuit):
        initial_guess (ndarray):
        num_iters (int): number of iteration of the VQE algorithm
        quantum_instance (QuantumInstance):
    Output:
        (dict):
    """
    results, callback = make_data_and_callback(save=returns)

    energy_hamiltonian = lambda params : energy_evaluation(hamiltonian, ansatz, params, quantum_instance)

    optimizer = SPSA_calibrated(energy_hamiltonian, initial_guess,
                                iter_start=iter_start, maxiter=num_iters,
                                callback=callback)

    optimizer.minimize(energy_hamiltonian , initial_guess)

    return results

def VQE_adiabatic(hamiltonian, ansatz, initial_guess, num_iters,
                  quantum_instance, transition_lims=(0.0, 1.0),
                  returns=['x', 'fx']):
    """
    VQE in the adiabatic regime

    Input:
        hamiltonian (PauliSumOp):
        ansatz (QuantumCircuit):
        initial_guess (ndarray):
        num_iters (int): number of iteration of the VQE algorithm
        quantum_instance (QuantumInstance):

    Output:
        ():
    """
    hamiltonian_local = global2local( hamiltonian )
    acc_adiabatic, cb = make_data_and_callback(save=returns)
    cost, cb = make_adiabatic_cost_and_callback(Hglobal = hamiltonian,
                                                Hlocal  = hamiltonian_local,
                                                circ    = ansatz,
                                                backend = quantum_instance,
                                                niters  = num_iters,
                                                transition_lims = transition_lims,
                                                callback = cb)
    optimizer = SPSA(maxiter=num_iters, callback=cb)
    optimizer.minimize(cost, initial_guess)

    return acc_adiabatic

def VQE_shift( hamiltonian, ansatz, initial_guess, max_iter, shift_iter, quantum_instance, iter_start=0, returns=['x', 'fx']):
    """
    Description

    Input:
        hamiltonian (PauliSumOp):
        ansatz (QuantumCircuit):
        initial_guess (ndarray):
        num_iters (int): number of iteration of the VQE algorithm
        quantum_instance (QuantumInstance):

    Output:
        ():
    """
    hamiltonian_local = global2local( hamiltonian )

    results_local  = VQE(hamiltonian_local, ansatz, initial_guess, shift_iter,
                         quantum_instance, returns, iter_start=iter_start)
    results_global = VQE(hamiltonian, ansatz, results_local['x'][-1], max_iter-shift_iter,
                         quantum_instance, iter_start=shift_iter+iter_start)

    results  = {'local'  : results_local,
                'global' : results_global}

    return results
