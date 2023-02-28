import pennylane as qml
import pennylane.numpy as np


def VQE( Hamiltonian, 
            ansatz, 
            params, 
            optimizer, 
            max_iterations = 100,
            conv_tol       = None, 
            device         = "default.qubit",
            shots          = None,
            callback       = None
             ):

    num_wires = len(Hamiltonian.wires)
    dev = qml.device(device, wires=num_wires, shots=shots)
    @qml.qnode(dev)
    def cost_func(params):
        ansatz(params)
        return qml.expval(Hamiltonian)

    # store the values of the cost function
    energy_k = [cost_func(params)]

    # store the values of the circuit parameter
    params_k = [params]

    for n in range(max_iterations):
        params, prev_energy = optimizer.step_and_cost(cost_func, params)

        energy_k.append(cost_func(params))
        params_k.append(params)

        if callback is not None:
            callback(params)

        # if n % 2 == 0:
            # print(f"Step = {n},  Energy = {energy_k[-1]:.8f} Ha")

        if conv_tol is not None:
            conv = np.abs(energy_k[-1] - prev_energy)
            if conv <= conv_tol:
                break

    return params_k, energy_k