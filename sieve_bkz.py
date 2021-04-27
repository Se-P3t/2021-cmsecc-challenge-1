"""
BKZ

References:
    g6k/bkz.py
"""

from __future__ import absolute_import
from __future__ import print_function
import os
import re
import sys
import time
from collections import OrderedDict

from fpylll import BKZ as BKZ_FPYLLL, GSO, IntegerMatrix
from fpylll.tools.quality import basis_quality

from g6k.algorithms.bkz import pump_n_jump_bkz_tour
from g6k.siever import Siever
from g6k import SieverParams
from g6k.utils.cli import run_all, pop_prefixed_params
from g6k.utils.stats import SieveTreeTracer
from g6k.utils.util import sanitize_params_names, db_stats
import numpy as np

import six
from six.moves import range

from util import str_mat, print_stats, load_matrix_file


__all__ = ['solve_bkz']


def bkz_kernel(arg0, params=None, seed=None):
    """
    Run the BKZ algorithm with different parameters.
    :param d: the dimension of the lattices to BKZ reduce
    :param params: parameters for BKZ:
        - bkz/blocksizes: given as low:high:inc perform BKZ reduction
          with blocksizes in range(low, high, inc) (after some light)
          prereduction
        - bkz/pre_blocksize: prereduce lattice with fpylll BKZ up
          to this blocksize
        - bkz/tours: the number of tours to do for each blocksize
        - bkz/extra_dim4free: lift to indices extra_dim4free earlier in
          the lattice than the currently sieved block
        - bkz/jump: the number of blocks to jump in a BKZ tour after
          each pump
        - bkz/dim4free_fun: in blocksize x, try f(x) dimensions for free,
          give as 'lambda x: f(x)', e.g. 'lambda x: 11.5 + 0.075*x'
        - pump/down_sieve: sieve after each insert in the pump-down
          phase of the pump
    """
    # Pool.map only supports a single parameter
    if params is None and seed is None:
        d, params, seed = arg0
    else:
        d = arg0

    # params for underlying BKZ/workout/pump
    dim4free_fun = params.pop("bkz/dim4free_fun")
    extra_dim4free = params.pop("bkz/extra_dim4free")
    jump = params.pop("bkz/jump")
    pump_params = pop_prefixed_params("pump", params)

    # flow of the bkz experiment
    algbkz = params.pop("bkz/alg")
    blocksizes = params.pop("bkz/blocksizes")
    blocksizes = eval("range(%s)" % re.sub(":", ",", blocksizes))
    pre_blocksize = params.pop("bkz/pre_blocksize")
    tours = params.pop("bkz/tours")

    load_matrix = params.pop("load_matrix")

    # misc
    verbose = params.pop("verbose")

    if blocksizes[-1] > d:
        print('set a smaller maximum blocksize with --blocksizes')
        return

    A, bkz = load_matrix_file(load_matrix, randomize=(seed != -1), seed=seed, int_type="mpz", float_type="mpfr")
    if verbose:
        print(("Loaded file '%s'" % load_matrix))
    for b in range(min(pre_blocksize, 10), pre_blocksize + 1):
        print("\r created, LLLed, BKZed %d" % b, end=" ")
        sys.stdout.flush()
        par = BKZ_FPYLLL.Param(
            b, strategies=BKZ_FPYLLL.DEFAULT_STRATEGY, max_loops=1, flags=BKZ_FPYLLL.MAX_LOOPS
        )
        bkz(par)

    MM = GSO.Mat(A, float_type="double",
                 U=IntegerMatrix.identity(A.nrows, int_type=A.int_type),
                 UinvT=IntegerMatrix.identity(A.nrows, int_type=A.int_type))

    g6k = Siever(MM, params, seed=seed)
    tracer = SieveTreeTracer(g6k, root_label=("bkz", d), start_clocks=True)

    M = g6k.M

    T0 = time.time()
    for blocksize in blocksizes:

        for t in range(tours):
            with tracer.context("tour", t, dump_gso=True):
                pump_n_jump_bkz_tour(g6k, tracer, blocksize, jump=jump,
                                     dim4free_fun=dim4free_fun,
                                     extra_dim4free=extra_dim4free,
                                     pump_params=pump_params)

            if verbose:
                slope = basis_quality(M)["/"]
                fmt = "{'alg': '%10s', 'jump':%2d, 'pds':%d, 'extra_d4f': %2d, 'beta': %2d, 'slope': %.5f, 'total walltime': %.3f}" # noqa
                print(fmt % (algbkz + "+" + g6k.params.default_sieve,
                             jump, pump_params["down_sieve"], extra_dim4free,
                             blocksize, slope, time.time() - T0))

    tracer.exit()
    slope = basis_quality(M)["/"]
    stat = tracer.trace
    stat.data['res'] = A
    stat.data["slope"] = np.array(slope)
    return stat


def bkz_tour(n, params, threads, **kwds):
    """
    Run bkz tours.
    ..  note :: that by default no information is printed.
        To enable set ``--dummy-tracer False`` and ``--verbose``.
    """

    lower_bound = kwds.get("lower_bound", n) # lowest lattice dimension to consider (inclusive)
    upper_bound = kwds.get("upper_bound", 0) # upper bound on lattice dimension to consider (exclusive)
    step_size = kwds.get("step_size", 2) # increment lattice dimension in these steps
    trials = kwds.get("trials", 1) # number of experiments to run per dimension
    workers = kwds.get("workers", 1) # number of parallel experiments to run
    seed = kwds.get("seed", 0) # randomness seed

    all_params = OrderedDict({f"'threads': {threads}, ": params})

    stats = run_all(bkz_kernel, list(all_params.values()),
                    lower_bound=lower_bound,
                    upper_bound=upper_bound,
                    step_size=step_size,
                    trials=trials,
                    workers=workers,
                    seed=seed)

    inverse_all_params = OrderedDict([(v, k) for (k, v) in six.iteritems(all_params)])
    stats = sanitize_params_names(stats, inverse_all_params)

    fmt = "{name:20s} :: n: {n:2d}, cputime {cputime:7.4f}s, walltime: {walltime:7.4f}s, "\
          "slope: {slope:1.5f}, |db|: 2^{avg_max:.2f}"
    print_stats(fmt, stats, ("cputime", "walltime", "slope", "avg_max"),
                extractf={"avg_max": lambda n, params, stat: db_stats(stat)[0]})

    res = list(stats.values())[0][0].data['res']

    return res

def solve_bkz(A, **kwds):
    """
    """
    if isinstance(A, list):
        return_list = True
        n = len(A)
    elif isinstance(A, IntegerMatrix):
        return_list = False
        n = A.nrows
    else:
        raise TypeError("Matrix `A` type '%s' unknown." % type(A))

    keep_tmpfile = kwds.pop('keep_tmpfile', False)
    threads = kwds.pop('threads', 4)
    load_matrix = kwds.pop('load_matrix', f'bkz-{n}.txt')
    verbose = kwds.pop('verbose', True)

    params = SieverParams(threads=threads,
                          load_matrix=load_matrix,
                          verbose=verbose)
    params["bkz/alg"] = "pump_and_jump"
    params["bkz/blocksizes"] = kwds.pop("blocksizes", "40:51:2")
    params["bkz/pre_blocksize"] = int(params["bkz/blocksizes"].split(':')[0])-1
    params["bkz/tours"] = kwds.pop("tours", 1)
    params["bkz/extra_dim4free"] = kwds.pop("extra_dim4free", 0)
    params["bkz/jump"] = kwds.pop("jump", 1)
    params["bkz/dim4free_fun"] = "default_dim4free_fun"
    params["pump/down_sieve"] = True

    with open(load_matrix, 'w') as f:
        f.write(str_mat(A))

    res = bkz_tour(n, params, threads, **kwds)

    if keep_tmpfile:
        with open(load_matrix, 'w') as f:
            f.write(str_mat(res))
    else:
        os.system(f'rm -f {load_matrix}')

    if return_list:
        B = [[res[i,j] for j in range(res.ncols)] for i in range(res.nrows)]
        return B
    else:
        return res
