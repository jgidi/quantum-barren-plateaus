import qiskit.opflow as of
from qiskit.algorithms import NumPyMinimumEigensolver

from qubap.qiskit.jorge.tools import SPSA_calibrated
from qubap.qiskit.jorge.tools import make_data_and_callback

def energy_evaluation(hamiltonian, ansatz, parameters , quantum_instance, callback=None):
    
    ansatz_state = of.StateFn( ansatz.bind_parameters(parameters) )
    measurement = of.StateFn( hamiltonian ).adjoint() @ ansatz_state
    
    pauli_circs = of.PauliExpectation().convert(measurement)
    sampler = of.CircuitSampler(quantum_instance).convert(pauli_circs)
    
    evaluation = sampler.eval().real

    if callback is not None:
        callback( parameters, evaluation )

    return evaluation
    
def VQE(hamiltonian, ansatz, initial_guess, num_iters, quantum_instance,
        returns='x', iter_start=1):

    results, callback = make_data_and_callback(save=returns)
    
    energy_hamiltonian = lambda params : energy_evaluation(hamiltonian, ansatz, params, quantum_instance)

    optimizer = SPSA_calibrated(energy_hamiltonian, initial_guess,
                                iter_start=iter_start, maxiter=num_iters,
                                callback=callback)

    optimizer.minimize(energy_hamiltonian , initial_guess)

    return results 


def classical_solver(hamiltonian):
    eig = NumPyMinimumEigensolver()
    results = eig.compute_minimum_eigenvalue( hamiltonian )
    return results
