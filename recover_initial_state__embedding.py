# -*- coding: utf-8 -*-
"""
recover initial state
"""
import os
import sys
import random
import argparse

import IPython
from fpylll import FPLLL

from util import read_data, save_solution, matrix_overview
from mrg import MRG, MRGSolver


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# chall params
parser.add_argument('m', type=int, help="modulus")
parser.add_argument('coeffs', type=str, help="coefficients `c_0,c_1,...,c_{n-1}`"
                    "of connection polynomial (sep by comma)")
parser.add_argument('d', type=int, help="number of outputs")
parser.add_argument('zbits', type=int, help="number of bits of unknowns `z_i`")
# basic args
parser.add_argument('--category', type=int, dest="category", help="challenge category")
parser.add_argument('--level', type=int, dest="level", help="challenge level")
parser.add_argument('--experiment', action='store_true', dest="experiment",
                    help="auto generate data and test the code")
parser.add_argument('-v', '--verbose', type=int, dest="verbose", default=0,
                    help="verbose level")
parser.add_argument('--threads', type=int, dest="threads", default=4, help="threads")
parser.add_argument('--sieve', action='store_true', dest="sieve",
                    help="sieve for asvp when solution is not found after BKZ reduction")
parser.add_argument('--seed',  type=int, dest="seed", default=0,
                    help="randomness seed (0 for auto choose)")
parser.add_argument('--debug', action='store_true', dest="debug",
                    help="call IPython.embed when running")
# args for BKZ
parser.add_argument('--block-size', type=int, dest="block_size", default=20,
                    help="BKZ block size")
# args for sieving
parser.add_argument('--workout/dim4free-dec', type=int, dest="workout__dim4free_dec", default=3,
                    help="By how much do we decreaseee dim4free at each iteration")

args, _ = parser.parse_known_args()


DEBUG = args.debug

m = args.m
d = args.d
zbits = args.zbits
mbits = m.bit_length()
coeffs = list(map(int, args.coeffs.replace(' ', '').split(',')))
n = len(coeffs)

""" verbose level
0: None
1: basic param & solution
2: verbose when sieving
3: matrix_overview of reduced basis
4: matrix_overview of input basis
5: verbose when BKZ
"""
VERBOSE = args.verbose
THREADS = args.threads
SIEVE = args.sieve
SEED = args.seed or int.from_bytes(os.urandom(8), 'big')
if VERBOSE >= 1:
    print(f"SEED: {SEED}\n")

random.seed(SEED)
FPLLL.set_threads(THREADS)

sieve_param = {
    'threads': THREADS,
    'verbose': (DEBUG or VERBOSE >= 2),
    'seed': SEED,
    #'dry_run': True,
    'workout/dim4free_dec': args.workout__dim4free_dec,
    'workout/save_prefix': f'logs/sieve-{args.category}-{args.level}',
    'workout/start_n': args.block_size+1,
    'pump/down_sieve': True,
}


if (args.category is None) or (args.level is None):
    if not args.experiment:
        raise ValueError("Undefined behavior")

    mrg = MRG.random(m, 16)
    if VERBOSE >= 1:
        print(f"initial state: {mrg.initial_state}")
        print()

    y_, z_ = mrg.output(d, zbits, return_z=True)

else:
    y_ = read_data(args.category, args.level)
    if len(y_) < d:
        raise ValueError(f"outputs is not enough (got {len(y_)}, expect >={d})")

    ybits = max(y_).bit_length()
    if ybits != mbits - zbits:
        raise ValueError(f"bit length of `y_i` is wrong (got {ybits}, expect {mbits-zbits})")

    mrg = MRG(m, coeffs)

solver = MRGSolver(mrg, zbits, y_)
solver.gen_lattice(d)
if DEBUG or VERBOSE >= 4:
    matrix_overview(solver.L)

solver.L.randomize_block()


solver.L.run_bkz(args.block_size, verbose=VERBOSE >= 5)
if DEBUG or VERBOSE >= 3:
    matrix_overview(solver.L)

if DEBUG and input('embed? '):
    IPython.embed()

for idx, row in enumerate(solver.L):
    init = solver.recover_init(list(row))
    if init is not None:
        if VERBOSE >= 1:
            print(f"solution: {init}")
        break
if init is not None:
    if mrg.initial_state is None:
        solution = {
            'modulus': m,
            'zbits': zbits,
            'coefficients': coeffs,
            'initial_state': init,
        }
        save_solution(args.category, args.level, solution)
    sys.exit(0)
elif not SIEVE:
    sys.exit(-1) # cannot find solution after BKZ


if VERBOSE >= 2:
    vol = solver.volf(m, n, zbits, d)
    expected_target = solver.mvf(m, zbits, d)
    expected_b0 = solver.ghf(vol, d + 1)
    print(f"E|v|         = {float(expected_target):.1f}")
    print(f"E|b[0]|      = {float(expected_b0):.1f}")
    print(f"E|v|/E|b[0]| = {float(expected_target/expected_b0):.3f}")
    print()

# NOTE: r0 in `workout` is squared
sieve_param['workout/goal_r0'] = (d + 1) * m ** 2
solver.L.solve_asvp(**sieve_param)
if DEBUG or VERBOSE >= 3:
    matrix_overview(solver.L)

if DEBUG and input('embed? '):
    IPython.embed()

for idx, row in enumerate(solver.L):
    init = solver.recover_init(list(row))
    if init is not None:
        if VERBOSE >= 1:
            print(f"solution: {init}")
        break
if init is not None:
    if mrg.initial_state is None:
        solution = {
            'modulus': m,
            'zbits': zbits,
            'coefficients': coeffs,
            'initial_state': init,
        }
        save_solution(args.category, args.level, solution)
else:
    sys.exit(-2) # cannot find solution after SIEVE
