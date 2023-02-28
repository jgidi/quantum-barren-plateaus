import numpy as np
from qiskit.algorithms.optimizers import SPSA
from qiskit.quantum_info import SparsePauliOp
from qiskit.opflow.primitive_ops import PauliSumOp

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


def parse_hamiltonian(paulis, coefs):
    if isinstance(paulis, str): paulis = [paulis]
    if isinstance(coefs, str): coefs = [coefs]

    return PauliSumOp( SparsePauliOp(paulis, coefs) )

def ladder_hamiltonian(num_qubits, transverse_field_intensity=0):
    def interactions(i, j):
        paulis = ['I'] * num_qubits
        for l in [i, j]:
            if l < num_qubits:
                paulis[l] = 'Z'
        return ''.join(paulis)
    
    # Add a term for each interacting pair
    operators = []
    for n in range(num_qubits):
        # Even n only, with upper bound
        if (not n%2) and (n < num_qubits-1):
            operators.append( interactions(n, n+1) )
        # Even and odd, with upper bound
        if n < num_qubits-2:
            operators.append( interactions(n, n+2) )
    operator_coeffs = [1]*len(operators)
    
    # If there is magnetic field
    if transverse_field_intensity:
        for i in range(num_qubits):
            paulis = ['I'] * num_qubits
            paulis[i] = 'X'
            operators.append(''.join(paulis))
            operator_coeffs.append(transverse_field_intensity)
        

    
    return parse_hamiltonian(operators, operator_coeffs)