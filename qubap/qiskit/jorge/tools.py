import numpy as np
from qiskit.algorithms.optimizers import SPSA

def make_data_and_callback(save=['x', 'fx']):
    if isinstance(save, str): save = [save]
    data = {key : [] for key in save}
    def cb(nfev, x, fx, dx, is_accepted=True):
        values = locals()
        for key in save:
            data[key].append(values[key])
    return data, cb

def SPSA_calibrated(fun, x0, iter_start=1, maxiter=100, **spsa_args):

    lr, pert = SPSA(**spsa_args).calibrate(fun, np.asarray(x0))
    ak, bk = lr(), pert()

    for _ in range(iter_start-1):
        next(ak)
        next(bk)

    ak = [next(ak) for _ in range(maxiter)]
    bk = [next(bk) for _ in range(maxiter)]

    return SPSA(learning_rate=ak, perturbation=bk, maxiter=maxiter, **spsa_args)

