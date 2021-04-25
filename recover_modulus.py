# -*- coding: utf-8 -*-
"""
recover modulus
"""

import os
import math
import random
import argparse
import IPython
from copy import copy

from Crypto.Util.number import getPrime
from sympy.matrices import Matrix, zeros, eye
from fpylll import FPLLL, IntegerMatrix, BKZ, LLL

from util import read_data, save_solution, matrix_overview, str_mat
from sieve_asvp import solve_asvp


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# chall params
parser.add_argument('mbits', type=int, help="number of bits of the modulus")
parser.add_argument('n', type=int, help="degree of feedback polynomial")
parser.add_argument('r', type=int)
parser.add_argument('t', type=int)
parser.add_argument('zbits', type=int, help="number of bits of unknowns `z_i`")
parser.add_argument('--N', type=int, dest='N', default=0, help="number of outputs")
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

mbits = args.mbits
n = args.n
r = args.r
t = args.t
zbits = args.zbits
N = args.N


if not (r > t > n):
    raise ValueError("r > t > n")


""" verbose level
0: print warnning when no solution was found
1: basic param
2: verbose when sieving
3: matrix_overview of reduced basis
4: matrix_overview of input basis
5: verbose when BKZ
10: print all state
100: print shffule infos
"""
VERBOSE = args.verbose
THREADS = args.threads
SIEVE = args.sieve
SEED = args.seed or int.from_bytes(os.urandom(8), 'big')
if VERBOSE >= 1:
    print(f"SEED: {SEED}\n")
    print()
random.seed(SEED)

idx = 0

# def is_linear_independent(vectors, vec):
#     _, indexes = Matrix(list(vectors)+[vec]).T.rref()
#     return len(indexes) == len(vectors)+1
# 
# def rank(vectors):
#     _, indexes = Matrix(vectors).T.rref()
#     return len(indexes)

def det(vectors):
    A = Matrix(vectors)
    return abs(A.det())

if (args.category is None) or (args.level is None):
    if not args.experiment:
        raise ValueError("Undefined behavior")
    if N == 0:
        raise ValueError("number of outputs is not specified in experiment")
    if N < r+t:
        raise ValueError("too little output")
    m = getPrime(mbits)
    COEFFS = [random.randint(0, m-1) for _ in range(n)]
    #COEFFS = [509369830, 464884724]
    INIT_STATE = [random.randint(0, m-1) for _ in range(n)]
    if VERBOSE >= 1:
        print(f"modulus: {m}")
        print(f"coefficients: {COEFFS}")
        print(f"initial state: {INIT_STATE}")
        print()
    STATE = copy(INIT_STATE)
    for j in range(n, N):
        a_j = sum(c_i*a_i for c_i, a_i in zip(COEFFS, STATE[-n:])) % m
        STATE.append(a_j)
    y_ = [a_i >> zbits for a_i in STATE]
    SOL = tuple(COEFFS)

    Q = zeros(n)
    for i in range(n):
        Q[i, n-1] = COEFFS[i]
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

    q_ = [None] * n
    Qj = Qn
    for j in range(n, r):
        q_.extend(Qj[:, 0].transpose().tolist())
        Qj *= Q
        Qj %= m

else:
    y_ = read_data(args.category, args.level)
    N = len(y_)
    if N < r+t-1:
        raise ValueError(f"outputs is not enough (got {N}, expect >={r+t-1})")
    SOL = None
    

    #raise NotImplementedError




MM = 1 << (mbits - zbits)
BB = math.ceil((2*MM*r)**(1.0*t/(r-t)))
KK = math.ceil(math.sqrt(r)*2**((r-1.0)/2) * BB)

groups = N-(r+t-2)
full_rank = r-n+1
vectors_per_group = 3

if VERBOSE >= 1:
    print(f"groups: {groups}")
    print(f"full_rank: {full_rank}")
    print()

ETA = []

while True:
    offset = random.randint(0, groups-1)

    # build basis matrix
    M = [[0]*(t+r) for _ in range(r)]
    for i in range(r):
        for j in range(t):
            M[i][j] = y_[offset+i+j] *KK
        M[i][t+i] = 1

    random.seed(SEED)
    random.shuffle(M)

    # BKZ
    B = IntegerMatrix.from_matrix(M)
    flags = BKZ.DEFAULT | BKZ.AUTO_ABORT
    if VERBOSE >= 5:
        flags |= BKZ.VERBOSE
    BKZ.reduction(B, BKZ.EasyParam(block_size=min(B.nrows, args.block_size), flags=flags))

    #if DEBUG or VERBOSE >= 3:
    #    matrix_overview(B)

    count = 0
    for i in range(B.nrows):
        b = list(B[i])
        if any(b_i != 0 for b_i in b[:t]):
            break
        eta = list(b[t:])

        if VERBOSE >= 100 or DEBUG:
            if SOL is not None:
                print(f"{i+1}/{B.nrows}", sum(e*e for e in eta), sum(e*a for e, a in zip(eta, STATE[offset:])), (eta[idx]+sum(eta[j]*q_[j][idx] for j in range(n, r)))%m)
            else:
                print(f"{i+1}/{B.nrows}", sum(e*e for e in eta))

        if SOL is not None:
            if (eta[idx]+sum(eta[j]*q_[j][idx] for j in range(n, r)))%m == 0:
                if 1:# offset == 0 or is_linear_independent(ETA, eta):
                    ETA.append(eta)
        else:
            if count < vectors_per_group:# and is_linear_independent(ETA, eta):
                ETA.append(eta)
                count += 1

    if ETA:
        ETA_m = IntegerMatrix.from_matrix(ETA, int_type='mpz')
        LLL.reduction(ETA_m)
        ETA = [list(b) for b in ETA_m if any(list(b))]

    if VERBOSE < 5:
        print(f"\roffset = {offset} :: rank(ETA) = {len(ETA)}/{full_rank}", end='')

    if len(ETA) >= full_rank:
        break

if VERBOSE < 5:
    print()



KK = math.ceil(100 * m * 2**((r-n)/2))

M = []
for eta in ETA:
    row = [eta[idx], eta[n]] + [e *KK for e in eta[n+1:r]]
    M.append(row)


B = IntegerMatrix.from_matrix(M)
flags = BKZ.DEFAULT | BKZ.AUTO_ABORT
if DEBUG or VERBOSE >= 5:
    flags |= BKZ.VERBOSE
BKZ.reduction(B, BKZ.EasyParam(block_size=args.block_size, flags=flags))

M = []
for b in B:
    b = list(b)
    if not any(b):
        continue
    row = b[:2] + [bi//KK for bi in b[2:]]
    M.append(row)

if DEBUG or VERBOSE >= 3:
    matrix_overview(M)

if len(M) == full_rank:
    detM = det(M)
    print(f"modulus = {detM}")
    if not args.experiment:
        with open(f"ETA-{args.category}-{args.level}.mat", 'w') as f:
            f.write(str_mat(ETA))

if DEBUG and input('embed? '):
    IPython.embed()
