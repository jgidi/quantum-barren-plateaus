import numpy as np

from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector

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
    sysA = int( n_qbits/2 )
    circ = QuantumCircuit( n_qbits )
        
    if isinstance( n_cnot, int ):
        for indx in range( n_cnot ):
            circ.cnot( indx , indx + sysA )

    elif n_cnot == 'Full_connect':        
        for indx in range( sysA ):
            circ.cnot( indx , indx + sysA )

    else:
        n_cnot = np.shape( n_cnot )[0]
        for indx in range( n_cnot ):
            circ.cnot( n_cnot[indx][0] , n_cnot[indx][1] )
            
    return circ


def SCL( params, n_qbits, n_qbits_crz=2, deep=1, name=None ):
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
    qc = QuantumCircuit( qubits, name=name )
        
    for i in range( qubits ):
        qc.ry( params[i], i )

    cc = -1
    for _ in range( 1, deep+1 ):

        for i in range( 0 , qubits + par, n_qbits_crz ):
            for l in range( 1, n_qbits_crz ):
                cc += 1
                qc.cz( i, i+l )
        
        for i in range( qubits ):
            cc += 1
            qc.ry( params[ qubits + cc ], i )

        for i in range( 1, qubits-1, n_qbits_crz ):
            for l in range( 1, n_qbits_crz ):
                cc += 1
                qc.cz( i, i+l )
                qc.ry( params[ qubits + cc ], i )
                qc.ry( params[ qubits + cc ], i + l )


    return qc.to_gate( )


def ansatz_constructor( n_qbits, unitaries=[SCL, SCL, SCL], n_qb_crz=[2,2,2], deep= [1,1,1], n_cnot='Full_connect', set_barrier=False ):
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

    num_params_SCL = ( deep[0] +1 )*n_qbits + deep[0]*(n_qbits - 2)
      
    params_1 = ParameterVector(r"$\theta$", num_params_SCL)
    params_2 = ParameterVector(r"$\phi$", num_params_SCL)
    params_3 = ParameterVector(r"$\omega$", num_params_SCL)
    qc     = QuantumCircuit( n_qbits )
    

    U_1   = unitaries[0]( params_1, range( int(n_qbits/2) ), n_qbits_crz=int(n_qb_crz[0]), deep=deep[0], name="U1" )
    U_2   = unitaries[1]( params_2, range( int(n_qbits/2) ), n_qbits_crz=int(n_qb_crz[1]), deep=deep[1], name="U2" )
    U_3   = unitaries[2]( params_3, range( int(n_qbits/2) ), n_qbits_crz=int(n_qb_crz[2]), deep=deep[2], name="U3" )
    ent_l = cnot_layer( n_qbits, n_cnot=n_cnot )


    qc.compose( U_1 , qubits=range( int(n_qbits/2) ), clbits=None, inplace=True )         # Compose entangled layer cnots
    if set_barrier== True:
        qc.barrier()

    qc.compose( ent_l , qubits=range( n_qbits ), clbits=None, inplace=True )              # Compose entangled layer cnots
    if set_barrier== True:
        qc.barrier()

    qc.compose( U_2 , qubits=range( int(n_qbits/2) ), clbits=None, inplace=True )         # Compose entangled layer cnots
    qc.compose( U_3 , qubits=range( int(n_qbits/2), n_qbits ), clbits=None, inplace=True ) # Compose entangled layer cnots
    
    return qc
