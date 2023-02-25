import numpy as np
from qiskit.algorithms.optimizers import SPSA
from qubap.qiskit.luciano.variational_algorithms import energy_evaluation

def make_array_and_callback(niters, nx):
    xacc = np.empty((niters, nx))
    i = np.zeros(1, dtype=int)
    def cb(nfev, x, fx, dx, is_accepted=True):
        xacc[i[0], :] = x
        i[0] += 1

    return xacc, cb

def SPSA_calibrated(fun, x0, iter_start=1, maxiter=100, **spsa_args):

    lr, pert = SPSA(**spsa_args).calibrate(fun, np.asarray(x0))
    ak, bk = lr(), pert()

    for _ in range(iter_start-1):
        next(ak)
        next(bk)

    ak = np.array([next(ak) for _ in range(maxiter)])
    bk = np.array([next(bk) for _ in range(maxiter)])

    return SPSA(learning_rate=ak, perturbation=bk, maxiter=maxiter, **spsa_args)

def make_adiabatic_cost_and_update(Hlocal, Hglobal, circ, backend, niters):
    s = np.zeros(1)
    def update(i=None):
        if i is None:
            s[0] += 1
        else:
            s[0] = i
    def cost(x):
        # Linearly increasing. Reaches 1.0 at 80% of the iterations.
        a = float(s[0] / niters) * 5/4
        a = max(a, 1.0)

        # if a > 0:
        #     print(f'a = {a}, type of a: {type(a)}')
        H = (1 - a)*Hlocal + a*Hglobal
        return energy_evaluation(H, circ, x, backend)
    return cost, update

def make_adiabatic_cost_and_callback(Hlocal, Hglobal, circ, backend, niters, callback=None):
    cost, update = make_adiabatic_cost_and_update(Hlocal, Hglobal, circ, backend, niters)
    def cb_wrapper(nfev, x, fx, dx, is_accepted=True):
        update()
        if callback is not None:
            callback(nfev, x, fx, dx, is_accepted)
    return cost, cb_wrapper