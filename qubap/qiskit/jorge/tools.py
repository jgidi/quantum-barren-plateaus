import numpy as np
from qiskit.algorithms.optimizers import SPSA

def make_list_and_callback(save='x'): # or save='fx'
    acc = []
    if save == 'x':
        def cb(nfev, x, fx, dx, is_accepted=True):
            acc.append(x)
    else:
        def cb(nfev, x, fx, dx, is_accepted=True):
            acc.append(fx)

    return acc, cb

def SPSA_calibrated(fun, x0, iter_start=1, maxiter=100, **spsa_args):

    lr, pert = SPSA(**spsa_args).calibrate(fun, np.asarray(x0))
    ak, bk = lr(), pert()

    for _ in range(iter_start-1):
        next(ak)
        next(bk)

    ak = np.array([next(ak) for _ in range(maxiter)])
    bk = np.array([next(bk) for _ in range(maxiter)])

    return SPSA(learning_rate=ak, perturbation=bk, maxiter=maxiter, **spsa_args)

