import qiskit.opflow as of
from qubap.qiskit.jorge.tools import SPSA_calibrated
from qubap.qiskit.jorge.tools import make_array_and_callback

def energy_evaluation( hamiltonian , ansatz, parameters , quantum_instance, callback=None  ):
    
    ansatz_state = of.StateFn( ansatz.bind_parameters( parameters ) )
    measurement = of.StateFn( hamiltonian ).adjoint() @ ansatz_state
    
    pauli_circs = of.PauliExpectation().convert( measurement )
    sampler = of.CircuitSampler(quantum_instance).convert( pauli_circs )
    
    evaluation = sampler.eval().real

    if callback is not None:
        callback( parameters, evaluation )

    return evaluation
    
# def VQE( hamiltonian , ansatz , initial_guess, optimizator , quantum_instance, callback=None  ):
    
#     energy_hamiltonian = lambda params : energy_evaluation( hamiltonian, ansatz, params , quantum_instance, callback  )
#     result = optimizator.minimize( energy_hamiltonian , initial_guess )
    
#     return result

def VQE( hamiltonian, ansatz, initial_guess, num_iters, quantum_instance, iter_start=1 ):

    num_params = len( initial_guess )
    results, callback = make_array_and_callback( num_iters, num_params )
    
    energy_hamiltonian = lambda params : energy_evaluation( hamiltonian, ansatz, params , quantum_instance  )

    optimizer = SPSA_calibrated( energy_hamiltonian, initial_guess, iter_start=iter_start, 
                                maxiter=num_iters, callback=callback )

    optimizer.minimize( energy_hamiltonian , initial_guess )

    return results 


from qiskit.algorithms import NumPyMinimumEigensolver

def classical_solver(hamiltonian):
    eig = NumPyMinimumEigensolver()
    results = eig.compute_minimum_eigenvalue( hamiltonian )
    return results