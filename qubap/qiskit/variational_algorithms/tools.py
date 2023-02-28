#!/usr/bin/env python3

import numpy as np
import qiskit.opflow as of
from qiskit.algorithms import NumPyMinimumEigensolver
from qiskit.algorithms.optimizers import SPSA

def classical_solver(hamiltonian):
    eig = NumPyMinimumEigensolver()
    results = eig.compute_minimum_eigenvalue( hamiltonian )
    return results


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

def make_adiabatic_cost_and_callback(Hlocal, Hglobal, circ, backend, niters,
                                     transition_lims=(0.0, 1.0), callback=None):
    s = [0]
    a1, a2 = np.min(transition_lims), np.max(transition_lims)
    def get_a(x):
        if a1 < x < a2:
            return float((x - a1) / (a2 - a1))
        else:
            return int( x > a1 )
    def cost(x):
        # 'a' is linearly increasing.
        # Starts at 0 for a1% of the iterations, and reaches 1.0 at a2%
        a = get_a( s[0] / niters )
        if a <= 0:
            H = Hlocal
        elif a >= 1:
            H = Hglobal
        else:
            H = (1 - a)*Hlocal + a*Hglobal
        return energy_evaluation(H, circ, x, backend)
    def update(i=None):
        if i is None:
            s[0] += 1
        else:
            s[0] = i
    def cb_wrapper(nfev, x, fx, dx, is_accepted=True):
        update()
        if callback is not None:
            callback(nfev, x, fx, dx, is_accepted)
    return cost, cb_wrapper

def make_data_and_callback(save=['x', 'fx']):
    if isinstance(save, str): save = [save]
    data = {key : [] for key in save}
    def cb(nfev, x, fx, dx, is_accepted=True):
        values = locals()
        for key in save:
            data[key].append(values[key])
    return data, cb

def SPSA_calibrated(fun, x0, iter_start=1, maxiter=100, **spsa_args):

    lr, pert = SPSA(**spsa_args).calibrate(fun, np.asarray(x0))
    ak, bk = lr(), pert()

    for _ in range(iter_start-1):
        next(ak)
        next(bk)

    ak = [next(ak) for _ in range(maxiter)]
    bk = [next(bk) for _ in range(maxiter)]

    return SPSA(learning_rate=ak, perturbation=bk, maxiter=maxiter, **spsa_args)