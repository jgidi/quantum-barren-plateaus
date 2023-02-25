import qiskit.opflow as of

def energy_evaluation( hamiltonian , ansatz, parameters , quantum_instance, callback=None  ):
    
    ansatz_state = of.StateFn( ansatz.bind_parameters( parameters ) )
    measurement = of.StateFn( hamiltonian ).adjoint() @ ansatz_state
    
    pauli_circs = of.PauliExpectation().convert( measurement )
    sampler = of.CircuitSampler(quantum_instance).convert( pauli_circs )
    
    evaluation = sampler.eval().real

    if callback is not None:
        callback( parameters, evaluation )

    return evaluation
    
def VQE( hamiltonian , ansatz , initial_guess, optimizator , quantum_instance, callback=None  ):
    
    energy_hamiltonian = lambda params : energy_evaluation( hamiltonian, ansatz, params , quantum_instance, callback  )
    result = optimizator.minimize( energy_hamiltonian , initial_guess )
    
    return result

from qiskit.algorithms import NumPyMinimumEigensolver

def classical_solver(hamiltonian):
    eig = NumPyMinimumEigensolver()
    results = eig.compute_minimum_eigenvalue( hamiltonian )
    return results