
import numpy as np
from qiskit_nature.circuit.library import HartreeFock
from qiskit_nature.transformers.second_quantization.electronic import FreezeCoreTransformer
from qiskit_nature.problems.second_quantization.electronic import ElectronicStructureProblem
from qiskit_nature.mappers.second_quantization import ParityMapper, JordanWignerMapper, BravyiKitaevMapper
from qiskit_nature.converters.second_quantization.qubit_converter import QubitConverter
from qiskit_nature.drivers.second_quantization import PySCFDriver

def HeisenbergHamiltonian(J=1, H=1, num_qubits=2, neighbours=None):
    """
    Qiskit operator of the 3-D Heisenberg Hamiltonian of a lattice of spins.
    H = - J Σ_j ( X_j X_{j+1} + Y_j Y_{j+1} + Z_j Z_{j+1} ) - H Σ_j Z_j
    Parameters
    ----------
    J: float
        Coupling constant.
    H: float
        External magnetic field.
    num_qubits: int.
        Number of qubits.
    neighbours: list(tuples).
        Coupling between the spins.
    Return
    ------
    Hamiltonian: SummedOp
        Heisenberg Hamiltonian of the system.
    """

    if neighbours is None:
        neighbours = [(0, 1)]

    num_op = num_qubits + 3 * len(neighbours)
    Hamiltonian_op_x = []
    Hamiltonian_op_z = []
    Hamiltonian_coeff = num_qubits * [-H] + num_op * [-J]

    for idx in range(num_qubits):
        op_x = np.zeros(num_qubits)
        op_z = np.zeros(num_qubits)
        op_z[idx] = 1
        Hamiltonian_op_x.append(op_x.copy())
        Hamiltonian_op_z.append(op_z.copy())

    for idx in neighbours:
        op_x = np.zeros(num_qubits)
        op_z = np.zeros(num_qubits)
        op_x[idx[0]] = 1
        op_x[idx[1]] = 1
        Hamiltonian_op_x.append(op_x.copy())
        Hamiltonian_op_z.append(op_z.copy())
        op_z[idx[0]] = 1
        op_z[idx[1]] = 1
        Hamiltonian_op_x.append(op_x.copy())
        Hamiltonian_op_z.append(op_z.copy())
        op_x[idx[0]] = 0
        op_x[idx[1]] = 0
        Hamiltonian_op_x.append(op_x.copy())
        Hamiltonian_op_z.append(op_z.copy())

    Hamiltonian = SummedOp(
        [PauliOp(Pauli((Hamiltonian_op_z[j], Hamiltonian_op_x[j])), Hamiltonian_coeff[j]) for j in range(num_op)])

    return Hamiltonian


def H2(distance=None, freeze_core=True, remove_orbitals=False, initial_state=False, operator=True,
       mapper_type='ParityMapper'):
    """
    Qiskit operator of the LiH
    Parameters
    ----------
    distance: float (optional)
        Distance between atoms of Li and H
    freeze_core: Bool (optional)
        If freeze some cores that do highly impact in the energy
    remove_orbitals: Bool (optional)
        Remove some orbitals that do no impact in the energy
    initial_state: Bool (optional)
        Return the initial Hartree Fock state
    operator: Bool (optional)
    mapper_type: str (optional)
        Type of mapping between orbitals and qubits. Available options:
            'ParityMapper'
            'JordanWignerMapper'
            'BravyiKitaevMapper'
    Returns
    -------
    qubit_op: SummedOp
        Pauli strings and coefficients for the Hamiltonian
    init_state: QuantumCircuit (if initial_state=True)
        Quantum Circuit with the initial state given by Hartree Fock
    """

    if distance is None:
        distance = .761

    molecule = 'H .0 .0 .0; H .0 .0 ' + str(distance)

    try:
        driver = PySCFDriver(molecule)
    except Exception:
        from qiskit_nature.drivers.second_quantization.pyquanted import PyQuanteDriver
        driver = PyQuanteDriver(molecule)

    # qmolecule = driver.run()

    if remove_orbitals is False:
        Transformer = FreezeCoreTransformer(freeze_core=freeze_core)
    else:
        Transformer = FreezeCoreTransformer(freeze_core=freeze_core, remove_orbitals=remove_orbitals)

    problem = ElectronicStructureProblem(driver, transformers=[Transformer])

    # Generate the second-quantized operators
    second_q_ops = problem.second_q_ops()

    # Hamiltonian
    main_op = second_q_ops[0]

    # Setup the mapper and qubit converter
    if mapper_type == 'ParityMapper':
        mapper = ParityMapper()
    elif mapper_type == 'JordanWignerMapper':
        mapper = JordanWignerMapper()
    elif mapper_type == 'BravyiKitaevMapper':
        mapper = BravyiKitaevMapper()
    else:
        # TODO: Raise an error
        return None

    # The fermionic operators are mapped
    converter = QubitConverter(mapper=mapper, two_qubit_reduction=True)

    if operator is False:
        return converter, problem
    else:
        particle_number = problem.grouped_property_transformed.get_property("ParticleNumber")
        num_particles = (particle_number.num_alpha, particle_number.num_beta)
        num_spin_orbitals = particle_number.num_spin_orbitals
        qubit_op = converter.convert(main_op, num_particles=num_particles)
        if initial_state is False:
            return qubit_op
        else:
            init_state = HartreeFock(num_spin_orbitals, num_particles, converter)
            return qubit_op, init_state


def LiH(distance=None, freeze_core=True, remove_orbitals=None, initial_state=False, operator=True,
        mapper_type='ParityMapper'):
    """
    Qiskit operator of the LiH
    Parameters
    ----------
    distance: float (optional)
        Distance between atoms of Li and H
    freeze_core: Bool (optional)
        If freeze some cores that do highly impact in the energy
    remove_orbitals: Bool (optional)
        Remove some orbitals that do no impact in the energy
    initial_state: Bool (optional)
        Return the initial Hartree Fock state
    operator: Bool (optional)
    mapper_type: str (optional)
        Type of mapping between orbitals and qubits. Available options:
            'ParityMapper'
            'JordanWignerMapper'
            'BravyiKitaevMapper'
    Returns
    -------
    qubit_op: SummedOp
        Pauli strings and coefficients for the Hamiltonian
    init_state: QuantumCircuit (if initial_state=True)
        Quantum Circuit with the initial state given by Hartree Fock
    """

    if distance is None:
        distance = 1.5474

    if remove_orbitals is None:
        remove_orbitals = [3, 4]

    molecule = 'Li 0.0 0.0 0.0; H 0.0 0.0 ' + str(distance)

    try:
        driver = PySCFDriver(molecule)
    except Exception:
        from qiskit_nature.drivers.second_quantization.pyquanted import PyQuanteDriver
        driver = PyQuanteDriver(molecule)

    # qmolecule = driver.run()

    if remove_orbitals is False:
        Transformer = FreezeCoreTransformer(freeze_core=freeze_core)
    else:
        Transformer = FreezeCoreTransformer(freeze_core=freeze_core, remove_orbitals=remove_orbitals)

    problem = ElectronicStructureProblem(driver, transformers=[Transformer])

    # Generate the second-quantized operators
    second_q_ops = problem.second_q_ops()

    # Hamiltonian
    main_op = second_q_ops[0]

    # Setup the mapper and qubit converter
    if mapper_type == 'ParityMapper':
        mapper = ParityMapper()
    elif mapper_type == 'JordanWignerMapper':
        mapper = JordanWignerMapper()
    elif mapper_type == 'BravyiKitaevMapper':
        mapper = BravyiKitaevMapper()
    else:
        return None

    # The fermionic operators are mapped
    converter = QubitConverter(mapper=mapper, two_qubit_reduction=True)

    if operator is False:
        return converter, problem
    else:
        # # The fermionic operators are mapped to qubit operators
        # num_particles = (problem.grouped_property_transformed.get_property("ParticleNumber").num_alpha,
        #                  problem.grouped_property_transformed.get_property("ParticleNumber").num_beta)
        # num_spin_orbitals = 2 * problem.molecule_data_transformed.num_molecular_orbitals

        particle_number = problem.grouped_property_transformed.get_property("ParticleNumber")
        num_particles = (particle_number.num_alpha, particle_number.num_beta)
        num_spin_orbitals = particle_number.num_spin_orbitals
        qubit_op = converter.convert(main_op, num_particles=num_particles)
        if initial_state is False:
            return qubit_op
        else:
            init_state = HartreeFock(num_spin_orbitals, num_particles, converter)
            return qubit_op, init_state


def BeH2(distance=None, freeze_core=True, remove_orbitals=None, initial_state=False, operator=True,
         mapper_type='ParityMapper'):
    """
    Qiskit operator of the BeH2
    Parameters
    ----------
    distance: float (optional)
        Distance between atoms of Be and H
    freeze_core: Bool (optional)
        If freeze some cores that do highly impact in the energy
    remove_orbitals: Bool (optional)
        Remove some orbitals that do no impact in the energy
    initial_state: Bool (optional)
        Return the initial Hartree Fock state
    operator: Bool (optional)
    mapper_type: str (optional)
        Type of mapping between orbitals and qubits. Available options:
            'ParityMapper'
            'JordanWignerMapper'
            'BravyiKitaevMapper'
    Returns
    -------
    qubit_op: SummedOp
        Pauli strings and coefficients for the Hamiltonian
    init_state: QuantumCircuit (if initial_state=True)
        Quantum Circuit with the initial state given by Hartree Fock
    """

    if distance is None:
        distance = 1.339

    if remove_orbitals is None:
        remove_orbitals = [3, 6]

    molecule = 'H 0.0 0.0 -' + str(distance) + '; Be 0.0 0.0 0.0; H 0.0 0.0 ' + str(distance)

    try:
        driver = PySCFDriver(molecule)
    except Exception:
        from qiskit_nature.drivers.second_quantization.pyquanted import PyQuanteDriver
        driver = PyQuanteDriver(molecule)

    # qmolecule = driver.run()
    if remove_orbitals is False:
        Transformer = FreezeCoreTransformer(freeze_core=freeze_core)
    else:
        Transformer = FreezeCoreTransformer(freeze_core=freeze_core, remove_orbitals=remove_orbitals)

    problem = ElectronicStructureProblem(driver, transformers=[Transformer])

    # Generate the second-quantized operators
    second_q_ops = problem.second_q_ops()

    # Hamiltonian
    main_op = second_q_ops[0]

    # Setup the mapper and qubit converter
    if mapper_type == 'ParityMapper':
        mapper = ParityMapper()
    elif mapper_type == 'JordanWignerMapper':
        mapper = JordanWignerMapper()
    elif mapper_type == 'BravyiKitaevMapper':
        mapper = BravyiKitaevMapper()
    else:
        return None

    # The fermionic operators are mapped
    converter = QubitConverter(mapper=mapper, two_qubit_reduction=True)

    if operator is False:
        return converter, problem
    else:
        # num_particles = (problem.grouped_property_transformed.get_property("ParticleNumber").num_alpha,
        #                  problem.grouped_property_transformed.get_property("ParticleNumber").num_beta)

        particle_number = problem.grouped_property_transformed.get_property("ParticleNumber")
        num_particles = (particle_number.num_alpha, particle_number.num_beta)
        num_spin_orbitals = particle_number.num_spin_orbitals
        qubit_op = converter.convert(main_op, num_particles=num_particles)
        if initial_state is False:
            return qubit_op
        else:
            init_state = HartreeFock(num_spin_orbitals, num_particles, converter)
            return qubit_op, init_state


def H2O(distance=None, freeze_core=True, remove_orbitals=None, initial_state=False, operator=True,
        mapper_type='ParityMapper'):
    """
    Qiskit operator of the BeH2
    Parameters
    ----------
    distance: float (optional)
        Distance between atoms of Be and H
    freeze_core: Bool (optional)
        If freeze some cores that do highly impact in the energy
    remove_orbitals: Bool (optional)
        Remove some orbitals that do no impact in the energy
    initial_state: Bool (optional)
        Return the initial Hartree Fock state
    operator: Bool (optional)
    mapper_type: str (optional)
        Type of mapping between orbitals and qubits. Available options:
            'ParityMapper'
            'JordanWignerMapper'
            'BravyiKitaevMapper'
    Returns
    -------
    qubit_op: SummedOp
        Pauli strings and coefficients for the Hamiltonian
    init_state: QuantumCircuit (if initial_state=True)
        Quantum Circuit with the initial state given by Hartree Fock
    """

    if distance is None:
        distance = 0.9584

    if remove_orbitals is None:
        remove_orbitals = [4]

    x = distance * np.sin(np.deg2rad(104.45 / 2))
    y = distance * np.cos(np.deg2rad(104.45 / 2))

    molecule = 'O 0.0 0.0 0.0; H ' + str(x) + ' ' + str(y) + ' 0.0; H -' + str(x) + ' ' + str(y) + ' 0.0'

    try:
        driver = PySCFDriver(molecule)
    except Exception:
        from qiskit_nature.drivers.second_quantization.pyquanted import PyQuanteDriver
        driver = PyQuanteDriver(molecule)

    # qmolecule = driver.run()
    if remove_orbitals is False:
        Transformer = FreezeCoreTransformer(freeze_core=freeze_core)
    else:
        Transformer = FreezeCoreTransformer(freeze_core=freeze_core, remove_orbitals=remove_orbitals)

    problem = ElectronicStructureProblem(driver, transformers=[Transformer])

    # Generate the second-quantized operators
    second_q_ops = problem.second_q_ops()

    # Hamiltonian
    main_op = second_q_ops[0]

    # Setup the mapper and qubit converter
    if mapper_type == 'ParityMapper':
        mapper = ParityMapper()
    elif mapper_type == 'JordanWignerMapper':
        mapper = JordanWignerMapper()
    elif mapper_type == 'BravyiKitaevMapper':
        mapper = BravyiKitaevMapper()
    else:
        return None

    # The fermionic operators are mapped
    converter = QubitConverter(mapper=mapper, two_qubit_reduction=True)

    if not operator:
        return converter, problem
    else:
        # num_particles = (problem.grouped_property_transformed.get_property("ParticleNumber").num_alpha,
        #                 problem.grouped_property_transformed.get_property("ParticleNumber").num_beta)

        particle_number = problem.grouped_property_transformed.get_property("ParticleNumber")
        num_particles = (particle_number.num_alpha, particle_number.num_beta)
        num_spin_orbitals = particle_number.num_spin_orbitals
        qubit_op = converter.convert(main_op, num_particles=num_particles)
        if initial_state is False:
            return qubit_op
        else:
            init_state = HartreeFock(num_spin_orbitals, num_particles, converter)
            return qubit_op, init_state


def CH4(distance=None, freeze_core=True, remove_orbitals=None, initial_state=False, operator=True,
        mapper_type='ParityMapper'):
    """
    Qiskit operator of the CH4
    Parameters
    ----------
    distance: float (optional)
        Distance between atoms of Be and H
    freeze_core: Bool (optional)
        If freeze some cores that do highly impact in the energy
    remove_orbitals: Bool (optional)
        Remove some orbitals that do no impact in the energy
    initial_state: Bool (optional)
        Return the initial Hartree Fock state
    operator: Bool (optional)
    mapper_type: str (optional)
        Type of mapping between orbitals and qubits. Available options:
            'ParityMapper'
            'JordanWignerMapper'
            'BravyiKitaevMapper'
    Returns
    -------
    qubit_op: SummedOp
        Pauli strings and coefficients for the Hamiltonian
    init_state: QuantumCircuit (if initial_state=True)
        Quantum Circuit with the initial state given by Hartree Fock
    """

    if distance is None:
        distance = 0.9573

    if remove_orbitals is None:
        remove_orbitals = [7, 8]

    #          H(1)
    #          O
    #   H(2)      H(3)   H(4)

    theta = 109.5
    r_inf = distance * np.cos(np.deg2rad(theta - 90))
    height_low = distance * np.sin(np.deg2rad(theta - 90))

    H1 = np.array([0, 0, distance])
    H2 = np.array([r_inf, 0, -height_low])
    H3 = np.array([-r_inf * np.cos(np.pi / 3), r_inf * np.sin(np.pi / 3), -height_low])
    H4 = np.array([-r_inf * np.cos(np.pi / 3), -r_inf * np.sin(np.pi / 3), -height_low])

    molecule = 'O 0 0 0; H {}; H {}; H {}; H {}'.format(str(H1)[1:-1], str(H2)[1:-1], str(H3)[1:-1], str(H4)[1:-1])

    try:
        driver = PySCFDriver(molecule)
    except Exception:
        from qiskit_nature.drivers.second_quantization.pyquanted import PyQuanteDriver
        driver = PyQuanteDriver(molecule)

    if remove_orbitals is False:
        Transformer = FreezeCoreTransformer(freeze_core=freeze_core)
    else:
        Transformer = FreezeCoreTransformer(freeze_core=freeze_core, remove_orbitals=remove_orbitals)

    problem = ElectronicStructureProblem(driver, transformers=[Transformer])

    # Generate the second-quantized operators
    second_q_ops = problem.second_q_ops()

    # Hamiltonian
    main_op = second_q_ops[0]

    # Setup the mapper and qubit converter
    if mapper_type == 'ParityMapper':
        mapper = ParityMapper()
    elif mapper_type == 'JordanWignerMapper':
        mapper = JordanWignerMapper()
    elif mapper_type == 'BravyiKitaevMapper':
        mapper = BravyiKitaevMapper()
    else:
        return None

    # The fermionic operators are mapped
    converter = QubitConverter(mapper=mapper, two_qubit_reduction=True)

    if not operator:
        return converter, problem
    else:
        particle_number = problem.grouped_property_transformed.get_property("ParticleNumber")
        num_particles = (particle_number.num_alpha, particle_number.num_beta)
        num_spin_orbitals = particle_number.num_spin_orbitals
        qubit_op = converter.convert(main_op, num_particles=num_particles)
        if initial_state is False:
            return qubit_op
        else:
            init_state = HartreeFock(num_spin_orbitals, num_particles, converter)
            return qubit_op, init_state


def C2H2(distance=None, freeze_core=True, remove_orbitals=None, initial_state=False, operator=True,
         mapper_type='ParityMapper'):
    """
    Qiskit operator of the C2H2
    Parameters
    ----------
    distance: float (optional)
        Distance between atoms of Be and H
    freeze_core: Bool (optional)
        If freeze some cores that do highly impact in the energy
    remove_orbitals: Bool (optional)
        Remove some orbitals that do no impact in the energy
    initial_state: Bool (optional)
        Return the initial Hartree Fock state
    operator: Bool (optional)
    mapper_type: str (optional)
        Type of mapping between orbitals and qubits. Available options:
            'ParityMapper'
            'JordanWignerMapper'
            'BravyiKitaevMapper'
    Returns
    -------
    qubit_op: SummedOp
        Pauli strings and coefficients for the Hamiltonian
    init_state: QuantumCircuit (if initial_state=True)
        Quantum Circuit with the initial state given by Hartree Fock
    """

    if distance is None:
        distance = [1.2, 1.06]

    if remove_orbitals is None:
        remove_orbitals = [11]

    #   H(1)  C(1)  C(2)  H(2)

    H1 = str(np.array([0, 0, 0]))[1:-1]
    C1 = str(np.array([0, 0, distance[1]]))[1:-1]
    C2 = str(np.array([0, 0, distance[1] + distance[0]]))[1:-1]
    H2 = str(np.array([0, 0, 2 * distance[1] + distance[0]]))[1:-1]

    molecule = 'H {}; C {}; C {}; H {}'.format(H1, C1, C2, H2)

    try:
        driver = PySCFDriver(molecule)
    except Exception:
        from qiskit_nature.drivers.second_quantization.pyquanted import PyQuanteDriver
        driver = PyQuanteDriver(molecule)

    if remove_orbitals is False:
        Transformer = FreezeCoreTransformer(freeze_core=freeze_core)
    else:
        Transformer = FreezeCoreTransformer(freeze_core=freeze_core, remove_orbitals=remove_orbitals)

    problem = ElectronicStructureProblem(driver, transformers=[Transformer])

    # Generate the second-quantized operators
    second_q_ops = problem.second_q_ops()

    # Hamiltonian
    main_op = second_q_ops[0]

    # Setup the mapper and qubit converter
    if mapper_type == 'ParityMapper':
        mapper = ParityMapper()
    elif mapper_type == 'JordanWignerMapper':
        mapper = JordanWignerMapper()
    elif mapper_type == 'BravyiKitaevMapper':
        mapper = BravyiKitaevMapper()
    else:
        return None

    # The fermionic operators are mapped
    converter = QubitConverter(mapper=mapper, two_qubit_reduction=True)

    if not operator:
        return converter, problem
    else:
        particle_number = problem.grouped_property_transformed.get_property("ParticleNumber")
        num_particles = (particle_number.num_alpha, particle_number.num_beta)
        num_spin_orbitals = particle_number.num_spin_orbitals
        qubit_op = converter.convert(main_op, num_particles=num_particles)
        if initial_state is False:
            return qubit_op
        else:
            init_state = HartreeFock(num_spin_orbitals, num_particles, converter)
            return qubit_op, init_state


def molecules(molecule_name, distance=None, freeze_core=True, remove_orbitals=None, operator=True, initial_state=False,
              mapper_type='ParityMapper', load=False):
    if load:
        try:
            qubit_op = np.load('../data/molecules_qubitop.npy', allow_pickle=True).item()[molecule_name]
            print('Molecule loaded')
            return qubit_op
        except KeyError:
            print('Computing molecule')

    molecule_name = molecule_name.lower()
    if molecule_name == 'h2':
        return H2(distance=distance, freeze_core=freeze_core, remove_orbitals=remove_orbitals,
                  initial_state=initial_state, operator=operator, mapper_type=mapper_type)
    elif molecule_name == 'lih':
        return LiH(distance=distance, freeze_core=freeze_core, remove_orbitals=remove_orbitals,
                   initial_state=initial_state, operator=operator, mapper_type=mapper_type)
    elif molecule_name == 'beh2':
        return BeH2(distance=distance, freeze_core=freeze_core, remove_orbitals=remove_orbitals,
                    initial_state=initial_state, operator=operator, mapper_type=mapper_type)
    elif molecule_name == 'h2o':
        return H2O(distance=distance, freeze_core=freeze_core, remove_orbitals=remove_orbitals,
                   initial_state=initial_state, operator=operator, mapper_type=mapper_type)
    elif molecule_name == 'ch4':
        return CH4(distance=distance, freeze_core=freeze_core, remove_orbitals=remove_orbitals,
                   initial_state=initial_state, operator=operator, mapper_type=mapper_type)
    elif molecule_name == 'c2h2':
        return C2H2(distance=distance, freeze_core=freeze_core, remove_orbitals=remove_orbitals,
                    initial_state=initial_state, operator=operator, mapper_type=mapper_type)
    else:
        raise Exception('The molecule {} is not implemented.'.format(molecule_name))
