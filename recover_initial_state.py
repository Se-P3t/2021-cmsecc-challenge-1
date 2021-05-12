# -*- coding: utf-8 -*-
"""
recover initial state
"""
import os
import argparse
import IPython
import random
from copy import copy

from sympy.matrices import Matrix, zeros, eye
from fpylll import FPLLL, IntegerMatrix, BKZ

from util import read_data, save_solution, matrix_overview
from sieve_asvp import solve_asvp


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
parser.add_argument('-s', '--step-size', type=int, dest="step_size", default=2,
                    help="increment lattice dimension in these steps")
parser.add_argument('--workout/dim4free-dec', type=int, dest="workout__dim4free_dec", default=3,
                    help="By how much do we decreaseee dim4free at each iteration")
parser.add_argument('--goal-r0/gh', type=float, dest="goal_r0__gh", default=1.05,
                    help="Quit when this is reached")

args, _ = parser.parse_known_args()


DEBUG = args.debug

m = args.m
d = args.d
zbits = args.zbits
mbits = m.bit_length()
coeffs = list(map(int, args.coeffs.replace(' ', '').split(',')))
n = len(coeffs)

""" verbose level
0: print warnning when no solution was found
1: basic param
2: verbose when sieving
3: matrix_overview of reduced basis
4: matrix_overview of input basis
5: verbose when BKZ
10: print all state
100: 
"""
VERBOSE = args.verbose
THREADS = args.threads
SIEVE = args.sieve
SEED = args.seed or int.from_bytes(os.urandom(8), 'big')
if VERBOSE >= 1:
    print(f"SEED: {SEED}\n")

FPLLL.set_threads(THREADS)


if (args.category is None) or (args.level is None):
    if not args.experiment:
        raise ValueError("Undefined behavior")
    # TODO: test the code
    raise NotImplementedError
else:
    y_ = read_data(args.category, args.level)
    if len(y_) < d:
        raise ValueError(f"outputs is not enough (got {len(y_)}, expect >={d})")
    #assert max(y_).bit_length() <= mbits - zbits
    SOL = None

Q = zeros(n)
for i in range(n):
    Q[i, n-1] = coeffs[i]
    if i == 0:
        continue
    Q[i, i-1] = 1

def mat_pow_mod(A, e, m):
    B = eye(A.shape[0])
    while e:
        if e & 1:
            B *= A
            B %= m
        A *= A
        A %= m
        e >>= 1
    return B

Qn = mat_pow_mod(Q, n, m)

Q_ = [None] * n
Qj = Qn
for j in range(n, d):
    Q_.extend(Qj[:, 0].transpose().tolist())
    Qj *= Q
    Qj %= m


def check(initial_state):
    assert len(initial_state) == n
    for a_i, y_i in zip(initial_state, y_):
        if a_i >> zbits != y_i:
            return False
    state = copy(initial_state)
    for y_j in y_[n:]:
        a_j = sum(a_i*c_i for a_i, c_i in zip(state[-n:], coeffs)) % m
        if a_j >> zbits != y_j:
            return False
        state.append(a_j)

    if DEBUG or VERBOSE >= 10:
        print(f"\nstate: {state}")
    return True


def find_solution(B):
    nrows = B.nrows
    if SOL is not None:
        print(f"SOL: {SOL}")
    for idx, b in enumerate(B):
        b = list(b)
        if abs(b[0]) != bias*scale:
            if DEBUG:
                print(f"row {idx}:: b[0] != bias")
            continue
        if b[0] == bias*scale:
            b = [-bi for bi in b]

        z_ = list(int(bi/scale) + 2**(zbits-1) for bi in b[1:n+1])
        if max(z_) >= 2**zbits:
            if DEBUG:
                print(f"row {idx}:: b[i] too big")
            continue

        a_ = list(2**zbits*yi + zi for yi, zi in zip(y_, z_))

        if SOL is None:
            res = check(a_)
        else:
            res = tuple(z_) == SOL
        if VERBOSE >= 1:
            print(f"row {idx+1}/{nrows}:: {res}")
        if res:
            if VERBOSE >= 1:
                print(f"\nsolution: {a_}")
            return a_

    print("not found\n")
    return None


########## build basis matrix ##########
scale = m >> (zbits-1)
bias = 1 << (zbits-1)

M = [[0]*(d+1) for _ in range(d+1)]
for j in range(d+1):
    if j < 1:
        M[j][j] = bias * scale
    elif j < n+1:
        M[0][j] = bias * scale
        M[j][j] = 1 * scale
    elif j < d+1:
        c_j = 2**zbits * (y_[j-1] - sum(int(Q_[j-1][l])*y_[l] for l in range(n)))
        c_j %= m
        M[j][j] = m * scale
        M[0][j] = (c_j + 2**(zbits-1)) * scale
        for i in range(n):
            M[i+1][j] = int(Q_[j-1][i]) * scale

if DEBUG or VERBOSE >= 4:
    matrix_overview(M)

random.seed(SEED)
random.shuffle(M)

########## BKZ ##########
B = IntegerMatrix.from_matrix(M)
flags = BKZ.DEFAULT | BKZ.AUTO_ABORT
if DEBUG or VERBOSE >= 5:
    flags |= BKZ.VERBOSE
BKZ.reduction(B, BKZ.EasyParam(block_size=args.block_size, flags=flags))
#BKZ2(B)(BKZ.EasyParam(block_size=args.block_size, flags=flags))

if DEBUG or VERBOSE >= 3:
    matrix_overview(B)

if DEBUG and input('embed? '):
    IPython.embed()

a_ = find_solution(B)
if SOL is None and a_:
    solution = {
        'modulus': m,
        'zbits': zbits,
        'coefficients': coeffs,
        'initial_state': a_,
    }
    save_solution(args.category, args.level, solution)
    exit(0)

########## sieve ##########
if not SIEVE:
    exit(-1)


sieve_param = {
    'workout/dim4free_dec': args.workout__dim4free_dec,
    'workout/save_prefix': f'logs/sieve-{args.category}-{args.level}',
    'workout/start_n': args.block_size+1,
    'pump/down_sieve': True,
}
B = solve_asvp(
    B,
    threads=THREADS,
    verbose=(DEBUG or VERBOSE >= 2),
    seed=SEED,
    step_size=args.step_size,
    goal_r0__gh=args.goal_r0__gh,
    #dry_run=True,
    **sieve_param
)

if DEBUG or VERBOSE >= 3:
    matrix_overview(B)

if DEBUG and input('embed? '):
    IPython.embed()

a_ = find_solution(B)
if SOL is None and a_:
    solution = {
        'modulus': m,
        'zbits': zbits,
        'coefficients': coeffs,
        'initial_state': a_,
    }
    save_solution(args.category, args.level, solution)
