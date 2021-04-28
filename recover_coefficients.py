# -*- coding: utf-8 -*-
"""
recover initial state
"""
import os
import math
import random
import argparse
import IPython
from copy import copy

from sympy.matrices import Matrix, zeros, eye
from fpylll import FPLLL, IntegerMatrix, BKZ, LLL

from util import read_data, save_solution, matrix_overview
from sieve_asvp import solve_asvp


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# chall params
parser.add_argument('m', type=int, help="modulus")
parser.add_argument('n', type=int, help="degree of connection polynomial")
parser.add_argument('r', type=int)
parser.add_argument('t', type=int)
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
n = args.n
r = args.r
t = args.t
zbits = args.zbits
mbits = m.bit_length()

if not (r > t > n):
    raise ValueError("r > t > n")


""" verbose level
0: no output
1: basic param
2: print norm of first two vector
3: matrix_overview of reduced basis
4: 
5: verbose when BKZ
10: 
100: 
"""
VERBOSE = args.verbose
THREADS = args.threads
SIEVE = args.sieve
SEED = args.seed or int.from_bytes(os.urandom(8), 'big')
if VERBOSE >= 1:
    print(f"SEED: {SEED}\n")
random.seed(SEED)


def is_linear_independent(vectors, vec):
    _, indexes = Matrix(list(vectors)+[vec]).T.rref()
    return len(indexes) == len(vectors)+1



if (args.category is None) or (args.level is None):
    if not args.experiment:
        raise ValueError("Undefined behavior")
    COEFFS = [random.randint(0, m-1) for _ in range(n)]
    INIT_STATE = [random.randint(0, m-1) for _ in range(n)]
    if VERBOSE >= 1:
        print(f"coefficients: {COEFFS}")
        print(f"initial state: {INIT_STATE}")
    STATE = copy(INIT_STATE)
    for j in range(n, r+t-1):
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
    if len(y_) < r+t-1:
        raise ValueError(f"outputs is not enough (got {len(y_)}, expect >={r+t-1})")
    SOL = None


MM = 1 << (mbits - zbits)
BB = math.ceil((2*MM*r)**(1.0*t/(r-t)))
KK = math.ceil(math.sqrt(r)*2**((r-1.0)/2) * BB)

ETA = []

M = [[0]*(t+r) for _ in range(r)]
for i in range(r):
    for j in range(t):
        #M[i][j] = y_[i+j] *KK
        M[i][j] = (2*y_[i+j]+1) *KK
    M[i][t+i] = 1

random.shuffle(M)

B = IntegerMatrix.from_matrix(M)
flags = BKZ.DEFAULT | BKZ.AUTO_ABORT
if DEBUG or VERBOSE >= 5:
    flags |= BKZ.VERBOSE
BKZ.reduction(B, BKZ.EasyParam(block_size=min(B.nrows, args.block_size), flags=flags))

if DEBUG or VERBOSE >= 3:
    matrix_overview(B)

ETA = []
for i in range(2):
    b = list(B[i])
    if any(b[:t]):
        raise ValueError("we need `\sum_i \eta_i Y_i = 0`")
    eta = list(b[t:])

    if VERBOSE >= 2:
        if SOL is not None:
            print(i, sum(e*e for e in eta), sum(e*a for e, a in zip(eta, STATE)))
        else:
            print(i, sum(e*e for e in eta))

    if SOL is not None:
        if sum(e*a for e, a in zip(eta, STATE)) != 0:
            raise ValueError("we need `\sum_i \eta_i A_i = 0`")

    ETA.append(eta)

if VERBOSE >= 2:
    print()

ETA_m = IntegerMatrix.from_matrix(ETA, int_type='mpz')
LLL.reduction(ETA_m)
ETA = [list(b) for b in ETA_m if any(list(b))]



from sage.all import PolynomialRing, ZZ, Matrix
#from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence

from resultant import solve_with_resultant

PR = PolynomialRing(ZZ, names=[f"C{i}" for i in range(n)])
varC_ = PR.gens()

varQ = Matrix(PR, n)
for i in range(n):
    varQ[i, n-1] = varC_[i]
    if i == 0:
        continue
    varQ[i, i-1] = 1

varq_ = [None] * n

varQj = varQ**n
for j in range(n, r):
    varq_.append(varQj[:, 0].T[0].list())
    varQj *= varQ

polys = []
for eta in ETA:
    for i in range(n):
        gi = PR(eta[i] + sum(eta[j]*varq_[j][i] for j in range(n, r)))
        polys.append(gi)


for root in solve_with_resultant(polys[:n+1], m, verbose=VERBOSE):
    if all(poly(*root) % m == 0 for poly in polys):
        if VERBOSE >= 1:
            print(f"\ncoefficients: {root}")
        # only one solution
        break

if not args.experiment:
    solution = {
        'modulus': m,
        'coefficients': root,
    }
    save_solution(args.category, args.level, solution)



if DEBUG and input("embed? "):
    IPython.embed()




"""level 1
 % sage -python recover_coefficients.py 2147483647 2 30 8 17 --category 2 --level 1 --verbose 1 --block-size 2
SEED: 14823415482135119576

known roots: []
known roots: [1596998372]

coefficients: (1596998372, 913674193)
"""

"""level 2
 % sage -python recover_coefficients.py 2147483647 2 60 15 23 --category 2 --level 2 --verbose 1 --block-size 2
SEED: 16807067413575740166

known roots: []
known roots: [423368878]

coefficients: (423368878, 1375517413)
"""








