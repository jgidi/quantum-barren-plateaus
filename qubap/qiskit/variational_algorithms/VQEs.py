#!/usr/bin/env python3

def VQE(hamiltonian, ansatz, initial_guess, num_iters, quantum_instance,
        returns=['x', 'fx'], iter_start=0):
    """
    Standard VQE

    Input:
        hamiltonian (PauliSumOp):
        ansatz (QuantumCircuit):
        initial_guess (ndarray):
        num_iters (int): number of iteration of the VQE algorithm
        quantum_instance ()

    Output:
        (dict):
    """
    results, callback = make_data_and_callback(save=returns)

    energy_hamiltonian = lambda params : energy_evaluation(hamiltonian, ansatz, params, quantum_instance)

    optimizer = SPSA_calibrated(energy_hamiltonian, initial_guess,
                                iter_start=iter_start, maxiter=num_iters,
                                callback=callback)

    optimizer.minimize(energy_hamiltonian , initial_guess)

    return results


def classical_solver(hamiltonian):
    eig = NumPyMinimumEigensolver()
    results = eig.compute_minimum_eigenvalue( hamiltonian )
    return results

def make_adiabatic_cost_and_callback(Hlocal, Hglobal, circ, backend, niters,
                                     transition_lims=(0.0, 1.0), callback=None):
    s = [0]
    a1, a2 = np.min(transition_lims), np.max(transition_lims)
    def get_a(x):
        if a1 < x < a2:
            return float((x - a1) / (a2 - a1))
        else:
            return int( x > a1 )
    def cost(x):
        # 'a' is linearly increasing.
        # Starts at 0 for a1% of the iterations, and reaches 1.0 at a2%
        a = get_a( s[0] / niters )
        if a <= 0:
            H = Hlocal
        elif a >= 1:
            H = Hglobal
        else:
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

def VQE_adiabatic(hamiltonian, ansatz, initial_guess, num_iters,
                  quantum_instance, transition_lims=(0.0, 1.0),
                  returns=['x', 'fx']):

    hamiltonian_local = global2local( hamiltonian )
    acc_adiabatic, cb = make_data_and_callback(save=returns)
    cost, cb = make_adiabatic_cost_and_callback(Hglobal = hamiltonian,
                                                Hlocal  = hamiltonian_local,
                                                circ    = ansatz,
                                                backend = quantum_instance,
                                                niters  = num_iters,
                                                transition_lims = transition_lims,
                                                callback = cb)
    optimizer = SPSA(maxiter=num_iters, callback=cb)
    optimizer.minimize(cost, initial_guess)

    return acc_adiabatic

def VQE_shift( hamiltonian, ansatz, initial_guess, max_iter, shift_iter, quantum_instance, iter_start=0, returns=['x', 'fx']):

    hamiltonian_local = global2local( hamiltonian )

    results_local  = VQE(hamiltonian_local, ansatz, initial_guess, shift_iter,
                         quantum_instance, returns, iter_start=iter_start)
    results_global = VQE(hamiltonian, ansatz, results_local['x'][-1], max_iter-shift_iter,
                         quantum_instance, iter_start=shift_iter+iter_start)

    results  = {'local'  : results_local,
                'global' : results_global}

    return results
