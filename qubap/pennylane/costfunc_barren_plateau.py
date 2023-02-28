import pennylane as qml
import pennylane.numpy as np

def global2local( Hamiltonian ):

    num_wires = len(Hamiltonian.wires)

    coeff = []
    ops   = []

    coeff_local = Hamiltonian.coeffs

    for i, hamiltonian in enumerate(Hamiltonian.ops):
        
        if int(hamiltonian.num_wires)== -1 or int(hamiltonian.num_wires)==1:
            ops.append( qml.Identity(range(num_wires)) )
            coeff.append( coeff_local[i] )

        else:   
            pauli_local = hamiltonian.obs

            for pauli in pauli_local:
                ops.append( pauli )
                coeff.append( coeff_local[i]/num_wires  )

            ops.append( qml.Identity(range(num_wires)) )
            coeff.append( coeff_local[i]*(num_wires - hamiltonian.num_wires)/num_wires )

    return qml.Hamiltonian( np.array(coeff), ops).simplify()

def test_hamiltonian( num_wires ):

    Paulis = qml.Hamiltonian( [0.5, 0.5], [qml.Identity(0), qml.PauliZ(0) ] )
    for i in range(1,num_wires):
        zero = qml.Hamiltonian( [0.5, 0.5], [qml.Identity(i), qml.PauliZ(i) ] )
        Paulis = Paulis@zero

    Hamiltonian = qml.Hamiltonian( list(-Paulis.coeffs)+[1], Paulis.ops+[ qml.Identity(0) ] )

    return Hamiltonian.simplify()


def test_hamiltonian_2( num_wires, coeff ):

    obs = [ qml.pauli.string_to_pauli_word(num_wires*'X'),
            qml.pauli.string_to_pauli_word(num_wires*'Y') ,
            qml.pauli.string_to_pauli_word(num_wires*'Z') ] 

    return qml.Hamiltonian(coeff, obs)
