import numpy as np
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.opflow.primitive_ops import PauliSumOp, PauliOp
from qiskit.algorithms.optimizers import SPSA
from qubap.qiskit.luciano.variational_algorithms import VQE, energy_evaluation
from qubap.qiskit.jorge.tools import make_data_and_callback
import qiskit.opflow as of


def global2local( hamiltoniano, reduce=True ):

    num_qubits = hamiltoniano.num_qubits

    ops_local   = []
    coeff_local = []

    paulis_global = hamiltoniano.to_pauli_op()
    if isinstance( paulis_global, PauliOp ):
        paulis_global = [paulis_global]

    for pauli_global, coeff in zip( paulis_global, hamiltoniano.coeffs ):
        pauli_label = pauli_global.primitive.to_label()

        for qb in range(num_qubits):
            pauli_local = pauli_label[qb]
            x = np.zeros(num_qubits)
            z = np.zeros(num_qubits)

            if pauli_local == 'X':
                x[qb] = 1
            elif pauli_local == 'Z':
                z[qb] = 1
            elif pauli_local == 'Y':
                x[qb] = 1
                z[qb] = 1
            # elif pauli_local == 'I':
            #     continue
            
            ops_local.append( Pauli((z,x)) )
            coeff_local.append( coeff/num_qubits )

    hamiltoniano_local = PauliSumOp( SparsePauliOp( ops_local, coeff_local ) )

    return hamiltoniano_local.reduce() if reduce else hamiltoniano_local


def test_hamiltonian( num_qubits, coeff ):

    ops = []

    z = np.zeros(num_qubits)
    x = np.ones(num_qubits)
    ops.append( Pauli((z,x)) )

    z = np.ones(num_qubits)
    ops.append( Pauli((z,x)) )

    x = np.zeros(num_qubits)
    ops.append( Pauli((z,x)) )

    hamiltonian = PauliSumOp( SparsePauliOp( ops, coeff ) )

    return hamiltonian.reduce()

def test_hamiltonian_2( num_qubits ):

    Z = of.Z #Pauli Z
    I = of.I #Identidad
    Zero = 0.5*( I + Z )

    hamiltonian = eval( '('+(num_qubits-1)*'I^'+'I)-('+(num_qubits-1)*'Zero^'+'Zero)' )

    return hamiltonian.reduce()  

def make_adiabatic_cost_and_callback(Hlocal, Hglobal, circ, backend, niters, callback=None):
    s = [0]
    def cost(x):
        # Linearly increasing. Reaches 1.0 at 80% of the iterations.
        a = (s[0] / niters) * 5/4
        a = min(a, 1.0)
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


def VQE_adiabatic( hamiltonian, ansatz, initial_guess, num_iters, quantum_instance, returns='x'):

    hamiltonian_local = global2local( hamiltonian )
    acc_adiabatic, cb = make_data_and_callback(save=returns)
    cost, cb = make_adiabatic_cost_and_callback(Hglobal = hamiltonian, 
                                                Hlocal  = hamiltonian_local, 
                                                circ    = ansatz, 
                                                backend = quantum_instance, 
                                                niters  = num_iters, 
                                                callback = cb)
    optimizer = SPSA(maxiter=num_iters, callback=cb)
    optimizer.minimize(cost, initial_guess)

    return acc_adiabatic


def VQE_shift( hamiltonian, ansatz, initial_guess, max_iter, shift_iter, quantum_instance, returns=['x', 'fx']):

    hamiltonian_local = global2local( hamiltonian )

    results_local  = VQE(hamiltonian_local, ansatz, initial_guess, shift_iter,
                         quantum_instance, returns)
    results_global = VQE(hamiltonian, ansatz, results_local['x'][-1], max_iter-shift_iter, 
                         quantum_instance, iter_start=shift_iter) 

    return np.append( results_local['x'], results_global['x'], axis=0 )
