from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector

import numpy as np

def cnot_layer( n_qbits, n_cnot='Full_connect' ):
    """
    Create an entangling layer circuit.

    Input:
        -n_qbits: Number of qubits of the circuit
        
        -n_cnot:  Specify the type of entangling layer:
                    -Could be 'Full_connect' that means the circuit is full fill with CNOT's gates. 
                    -A intenger number that indicates the number of CNOT's to implement from the first qubit to the last of the first half of qubits in descendet order.
                    -A list of pair list e.g. [[1,2], [7,6]] where the first element of the array is the control qubit e.g 1, 7 and the second qubit of the list the target qubit e.g 2, 6.

    Output:
        - The entangling layer as a Quantum Circuit.
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


def parity_test( n ):
    """
    Check if the number is even or odd.

    Input:
        -n: A number.

    Output:
        -Return 0 if the number is even or -1 if the number is odd.
    """
    if n%2 == 0:
        p = 0
    else:
        p = -1
    
    return p


def SCL( params, n_qbits, n_qbits_crz=2, deep=1, name=None ):
    """
    Schmith Coefficient Layer.

    Input:
        -params:      The parameters to initialize the circuit.

        -n_qbits:     The number of qubits of the circuit.
        
        -n_qbits_crz: How many qubits are used in the control Z gate. Default is two; one control and one target.
        
        -deep:        The deep of the circuit or how many times the circuit is repeated. Default is one.
        
        -name:        The name of the circuit. Default None.

    Output:
        -A Parametric Quantum Circuit (PQC) that in our case could perform a schmith decomposition or a basis change.
    """

    qubits= len( n_qbits )
    par = parity_test( qubits )
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


    return qc.to_gate( )


def ansatz_constructor( n_qbits, unitaries=[SCL, SCL, SCL], n_qb_crz=[2,2,2], deep= [1,1,1], n_cnot='Full_connect' ):
    """
    The ansatz circuit that are made of three Parametrized Quantum Circuits (PQC's) and one entangling layer. The input parameters will define this four parts of the ansatz circuit.

    Input:
        -n_qbits:     The number of qubits of the circuit.
        
        -unitaries:   A list of three PQC's
        
        -n_qbits_crz: A list containing how many qubits are used in the control Z gate for each PQC's. Default are [2,2,2].
        
        -deep:        A list containing the deep of each PQC. Default is [1,1,1].

        -n_cnot:  Specify the type of entangling layer:
                    -Could be 'Full_connect' that means the circuit is full fill with CNOT's gates. 
                    -A intenger number that indicates the number of CNOT's to implement from the first qubit to the last of the first half of qubits in descendet order.
                    -A list of pair list e.g. [[1,2], [7,6]] where the first element of the array is the control qubit e.g 1, 7 and the second qubit of the list the target qubit e.g 2, 6.


    Output:
        -A Parametric Quantum Circuit (PQC) that in our case could perform a schmith decomposition or a basis change.
    """

    num_params_SCL = ( deep[0] +1 )*n_qbits
      
    params_1 = ParameterVector(r"$\theta$", num_params_SCL)
    params_2 = ParameterVector(r"$\phi$", num_params_SCL)
    params_3 = ParameterVector(r"$\omega$", num_params_SCL)
    qc     = QuantumCircuit( n_qbits )
    

    U_1   = unitaries[0]( params_1, range( int(n_qbits/2) ), n_qbits_crz=int(n_qb_crz[0]), deep=deep[0], name="U1" )
    U_2   = unitaries[1]( params_2, range( int(n_qbits/2) ), n_qbits_crz=int(n_qb_crz[1]), deep=deep[1], name="U2" )
    U_3   = unitaries[2]( params_3, range( int(n_qbits/2) ), n_qbits_crz=int(n_qb_crz[2]), deep=deep[2], name="U3" )
    ent_l = cnot_layer( n_qbits, n_cnot=n_cnot )


    qc.compose( U_1 , qubits=range( int(n_qbits/2) ), clbits=None, inplace=True )         # Compose entangled layer cnots
    qc.barrier()
    
    qc.compose( ent_l , qubits=range( n_qbits ), clbits=None, inplace=True )              # Compose entangled layer cnots
    qc.barrier()
    
    qc.compose( U_2 , qubits=range( int(n_qbits/2) ), clbits=None, inplace=True )         # Compose entangled layer cnots
    qc.compose( U_3 , qubits=range( int(n_qbits/2), n_qbits ), clbits=None, inplace=True ) # Compose entangled layer cnots
    
    return qc