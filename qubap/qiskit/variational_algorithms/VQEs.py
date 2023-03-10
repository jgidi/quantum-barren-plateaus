#!/usr/bin/env python3

from .tools import energy_evaluation, make_adiabatic_cost_and_callback
from qiskit.algorithms.optimizers import SPSA
from .tools import make_data_and_callback, SPSA_calibrated
# from qubap.qiskit.tools import make_data_and_callback, SPSA_calibrated

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

def VQE_adiabatic(hamiltonian_in, hamiltonian_out, ansatz, initial_guess, num_iters,
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

    acc_adiabatic, cb = make_data_and_callback(save=returns)
    cost, cb = make_adiabatic_cost_and_callback(Hglobal = hamiltonian_out,
                                                Hlocal  = hamiltonian_in,
                                                circ    = ansatz,
                                                backend = quantum_instance,
                                                niters  = num_iters,
                                                transition_lims = transition_lims,
                                                callback = cb)
    optimizer = SPSA(maxiter=num_iters, callback=cb)
    optimizer.minimize(cost, initial_guess)

    return acc_adiabatic


def VQE_shift( hamiltonian_in, hamiltonian_out, ansatz, initial_guess, max_iter, shift_iter, quantum_instance, iter_start=0, returns=['x', 'fx']):
    """
    Description

    Input:
        hamiltonian_in (PauliSumOp):
        hamiltonian_out (PauliSumOp):
        ansatz (QuantumCircuit):
        initial_guess (ndarray):
        num_iters (int): number of iteration of the VQE algorithm
        quantum_instance (QuantumInstance):

    Output:
        ():
    """

    results_in  = VQE(hamiltonian_in, ansatz, initial_guess, shift_iter,
                         quantum_instance, returns, iter_start=iter_start)
    results_out = VQE(hamiltonian_out, ansatz, results_in['x'][-1], max_iter-shift_iter,
                         quantum_instance, iter_start=shift_iter+iter_start)

    results  = {'in'  : results_in,
                'out' : results_out}

    return results
