import pennylane as qml
from pennylane import numpy as np

def cnot_layer( n_qbits, n_cnot='Full_connect' ):
    """
    Create an entangling layer circuit.

    Input:
    n_qbits (int): Number of qubits of the circuit
    n_cnot (int, list, str):  Specify the type of entangling layer:
        str = 'Full_connect' ; that means the circuit is full fill with CNOT's gates. 
        int = Number of CNOT's to implement from the first qubit to the last of the first half of qubits in descendet order.
        list= List of pair values e.g. [[1,2], [7,6]], first element of each list is control qubit e.g 1, 7 and the second the target qubit e.g 2, 6.

    Output:
    (QuantumCircuit): Entangling layer 
    """ 
    sysA = int( len(n_qbits)/2 )
        
    if isinstance( n_cnot, int ):
        for i in range( n_cnot ):
            qml.CNOT( wires=[ n_qbits[i] , n_qbits[i + sysA] ] )

    elif n_cnot == 'Full_connect':        
        for i in range( sysA ):
            qml.CNOT( wires=[ n_qbits[i] , n_qbits[i + sysA] ] )

    else:
        n_cnot = np.shape( n_cnot )[0]
        for i in range( n_cnot ):
            qml.CNOT( wires=[ n_cnot[i][0] , n_cnot[i][1] ] )


def SCL( params, n_qbits, n_qbits_crz=2, deep=1 ):
    """
    Schmidt Coefficient Layer (SCL) that performs a Schmidth decomposition or a basis change.

    Input:
    params (list):      Parameters to initialize the circuit.
    n_qbits (iterable): Qubits of the layer.
    n_qbits_crz (list): Qubits in the control Z gate. Default (2) is for 1 control qubit and 1 target qubit.    
    deep (int):         Times that the circuit is repeated.
    name (str):         Name of the circuit.

    Output:
    (QuantumCircuit): A Parametric Quantum Circuit.
    """

    qubits= len( n_qbits )
    par = -(qubits % 2) # parity test: -1 for odd, 0 for even
        
    for i in range( qubits ):
        qml.RY( params[i], wires=[ n_qbits[i] ] )

    cc = -1
    for _ in range( 1, deep+1 ):

        for i in range( 0 , qubits + par, n_qbits_crz ):
            for l in range( 1, n_qbits_crz ):
                cc += 1
                qml.CZ( wires=[ n_qbits[i], n_qbits[i+l] ] )
        
        for i in range( qubits ):
            cc += 1
            qml.RY( params[ qubits + cc ], wires=[ n_qbits[i] ] )

        for i in range( 1, qubits-1, n_qbits_crz ):
            for l in range( 1, n_qbits_crz ):
                cc += 1
                qml.CZ( wires=[ n_qbits[i], n_qbits[i+l] ] )
                qml.RY( params[ qubits + cc ], wires=[ n_qbits[i] ] )
                qml.RY( params[ qubits + cc ], wires=[ n_qbits[i + l] ] )


def ansatz_constructor( n_qbits, unitaries=[SCL, SCL, SCL], n_qb_crz=[2,2,2], deep= [1,1,1], n_cnot='Full_connect' ):
    """
    The State Efficient Ansatz parametric quantum circuit (PQC) that perform a schmith decomposition, entanglement and a basis change.
    The ansatz circuit is made of three PQC's and one entangling layer. 
    The input parameters will define this four gates of the ansatz circuit U_1, V, U_2 and U_3.

    Input:
    n_qbits (int):            Number of qubits.    
    unitaries (list):         Containing the three PQC's that conform the SEA.
    n_qbits_crz (list):       Qubits in the control Z gate for each PQC. Default ([2,2,2]) is for 1 control qubit and 1 target qubit for each PQC.    
    deep (list):              Times that each circuit is repeated.
    n_cnot (int, list, str):  Specify the type of entangling layer:
        str = 'Full_connect' ; that means the circuit is full fill with CNOT's gates. 
        int = Number of CNOT's to implement from the first qubit to the last of the first half of qubits in descendet order.
        list= List of pair values e.g. [[1,2], [7,6]], first element of each list is control qubit e.g 1, 7 and the second the target qubit e.g 2, 6.

    Output:
    (QuantumCircuit): The State Efficient Ansatz parametric quantum circuit (PQC).
    """
   
    def circuit(params):
        unitaries[0]( params[0], range( int(n_qbits/2) ), n_qbits_crz=int(n_qb_crz[0]), deep=deep[0] )
        cnot_layer( range(n_qbits), n_cnot=n_cnot )
        unitaries[1]( params[1], range( int(n_qbits/2) ), n_qbits_crz=int(n_qb_crz[1]), deep=deep[1] )
        unitaries[2]( params[2], range( int(n_qbits/2), n_qbits ), n_qbits_crz=int(n_qb_crz[2]), deep=deep[2] )

    return circuit