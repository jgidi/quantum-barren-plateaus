from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector

import numpy as np

def cnot_layer( n_qbits, n_cnot='Full_connect' ):   
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
    if n%2 == 0:
        p = 0
    else:
        p = -1
    
    return p


def _GateSmith( params, qubits, n_qbits_crz=2, deep=1, name=None ):
      
    qubits= len( qubits )
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



def SCL( params, qubits, n_qbits_crz=2, deep=1, name=None ):
    "Schmidt Coefficient Layer"
    #assert ( (qubits[-1]+1)%2 == 0 ) and ( type(qubits) == list ), "The qubits must be given in a list and the total number of qubits must be an even integer."

    #self.append( _GateSmith( params, qubits, n_qbits_crz=n_qbits_crz, deep=deep, name=name ), qubits )
    qc = _GateSmith( params, qubits, n_qbits_crz=n_qbits_crz, deep=deep, name=name )

    return qc 



def ansatz_constructor( qubits, unitaries=[SCL, SCL, SCL], n_qb_crz=[2,2,2], deep= [1,1,1], n_cnot='Full_connect' ):
    
    #assert isinstance( n_cnot, int ) or isinstance( n_cnot, list ) or ( n_cnot == 'Full_connect' ), "The variable n_cnot must be of type integer or a numpy.ndarray"
    #assert isinstance( unitaries, list ), "The unitaries must be given in a numpy.ndarray variable"
    
    num_params_SCL = ( deep[0] +1 )*qubits
    
    
    params_1 = ParameterVector(r"$\theta$", num_params_SCL)
    params_2 = ParameterVector(r"$\phi$", num_params_SCL)
    params_3 = ParameterVector(r"$\omega$", num_params_SCL)
    qc     = QuantumCircuit( qubits )
    

    U_1   = unitaries[0]( params_1, range( int(qubits/2) ), n_qbits_crz=int(n_qb_crz[0]), deep=deep[0], name="U1" )
    U_2   = unitaries[1]( params_2, range( int(qubits/2) ), n_qbits_crz=int(n_qb_crz[1]), deep=deep[1], name="U2" )
    U_3   = unitaries[2]( params_3, range( int(qubits/2) ), n_qbits_crz=int(n_qb_crz[2]), deep=deep[2], name="U3" )
    ent_l = cnot_layer( qubits, n_cnot=n_cnot )

    qc.compose( U_1 , qubits=range( int(qubits/2) ), clbits=None, inplace=True )         # Compose entangled layer cnots
    qc.barrier()
    
    qc.compose( ent_l , qubits=range( qubits ), clbits=None, inplace=True )              # Compose entangled layer cnots
    qc.barrier()
    
    qc.compose( U_2 , qubits=range( int(qubits/2) ), clbits=None, inplace=True )         # Compose entangled layer cnots
    qc.compose( U_3 , qubits=range( int(qubits/2), qubits ), clbits=None, inplace=True ) # Compose entangled layer cnots
    
    return qc