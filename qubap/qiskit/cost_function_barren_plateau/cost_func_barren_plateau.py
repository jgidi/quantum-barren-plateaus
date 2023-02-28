
import qiskit.opflow as of
from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector
import numpy as np
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.opflow.primitive_ops import PauliOp, PauliSumOp

""""
In this program we are trying to replicate the numerical 
simulations of a quantum autoencoder from the paper

Cerezo, M., Sone, A., Volkoff, T. et al. Cost function dependent
barren plateaus in shallow parametrized quantum circuits. 
Nat Commun 12, 1791 (2021). https://doi.org/10.1038/s41467-021-21728-w
"""

def global2local( hamiltoniano, reduce=True ):
    """
    Take a global Hamiltonian and reduce it to a local Hamiltonian

    Input:
        hamiltonian (PauliSumOp): Given Hamiltonian
        reduce (Boolean): Use or not reduce() in Local Hamiltonian

    Output:
        (PauliSumOp): Local Hamiltonian
    """
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




def global_observable(n_qbitsB, n_qbitsA=1):
    """
    Global observable (or Hamiltonian) of the form
    O_G = I_{AB} - I_{A} \otimes \ket{0}\bra{0}

    Input:
        n_qbitsB (int): number of qubits of Bob
        n_qbitsA (int, optional): number of qubits of Alice

    Output:
        (PauliSumOp): global observable
    """
    Z = of.Z
    I = of.I
    Zero = 0.5*( I + Z )
    I_AB = I
    IA_ZeroB = I
    for _ in range(n_qbitsA+n_qbitsB):
        I_AB = I_AB^I
    for _ in range(n_qbitsA):
        IA_ZeroB = IA_ZeroB^I
    for _ in range(n_qbitsB):
        IA_ZeroB = IA_ZeroB^Zero
    OG = I_AB - IA_ZeroB
    # OG = OG.to_pauli_op()
    return OG

def initial_state_ex(n_qbitsB, n_qbitsA=1):
    """
    Initial state given the states in eq. (25) and eq. (26) for
    the example in numerical simmulations

    Input:
        n_qbitsB (int): number of qubits of Bob
        n_qbitsA (int, optional): number of qubits of Alice

    Output:
        (QuantumCircuit): initial state
    """
    qbt_ancilla = 1
    num_total = n_qbitsA+n_qbitsB+qbt_ancilla
    circuit = QuantumCircuit(num_total)
    theta = 2*np.arccos(np.sqrt(2/3))
    circuit.ry(theta,1)
    circuit.cnot(1,0)
    circuit.cnot(1,2)
    circuit.cnot(2,3)
    return circuit

def variational_circuit(n_qbitsB, n_qbitsA=1, layers=1):
    """
    Variational circuit V(theta) following Fig. 4. Note that for the example
    the number of Alice's qubits are fixed n_A=1

    Input:
        n_qbitsB (int): number of qubits of Bob
        n_qbitsA (int, optional): number of qubits of Alice
        layers (int, optional): number of layers. Correspond to the deep of the circuit

    Output:
        (QuantumCircuit): variational quantum circuit
    """
    qbt_ancilla = 1
    n_total = n_qbitsA+n_qbitsB+qbt_ancilla
    circuit = QuantumCircuit(n_total)
    n_params = 2*n_qbitsB*layers +n_qbitsA+n_qbitsB
    params = ParameterVector(r"$\theta$", n_params)
    n = n_qbitsA+n_qbitsB-1
    p_ry2 = 2*n_qbitsB-1 # variable to change when the second row of ry start
    p_lay = 2*n_qbitsB # when we add layers, this variable we will help us to
                    # give continuity to the parameters
    for i in range(1,n_total):
        circuit.ry(params[i-1], i)
    circuit.barrier()
    for i in range(layers):
        for k in range(1,n_total-2):
            circuit.cz( k, k+1 )
        for k in range(1, n_total-1):
            circuit.ry(params[n+k+i*p_lay],k)
        for k in range(2,n_total-1):
            circuit.cz( k, k+1)
        for k in range(2,n_total):
            circuit.ry(params[p_ry2+k+i*p_lay],k)
        circuit.barrier()
    return circuit

def ansatz_numerical(n_qbitsB, n_qbitsA=1, layers=1):
    """"
    Circuit made by composing the initial state and the variational quantum circuit

    Input:
        n_qbitsB (int): number of qubits of Bob
        n_qbitsA (int, optional): number of qubits of Alice
        layers (int, optional): number of layers. Correspond to the deep of the circuit

    Output:
        (QuantumCircuit): Composed circuit
    """
    circuit= initial_state_ex(n_qbitsB, n_qbitsA)
    circuit.barrier()
    circuit.compose(variational_circuit(n_qbitsB, n_qbitsA, layers), inplace=True)
    return circuit

def local_observable(n_qbitsB, n_qbitsA=1):
    """
    Local Observable

    Input:
        n_qbitsB (int): number of qubits of Bob
        n_qbitsA (int, optional): number of qubits of Alice

    Output:
        (PauliSumOp): local observable
    """""
    local_observable = global2local(global_observable(n_qbitsB, n_qbitsA))
    return local_observable
