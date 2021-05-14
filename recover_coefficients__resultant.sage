# -*- coding: utf-8 -*-
"""
recover coefficients
    solve `g_i` with resultant
"""
import os
import argparse
import IPython

from fpylll import FPLLL, IntegerMatrix, BKZ, LLL

from util import read_data, save_solution, matrix_overview
from resultant import solve_with_resultant


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
parser.add_argument('--seed',  type=int, dest="seed", default=0,
                    help="randomness seed (0 for auto choose)")
parser.add_argument('--debug', action='store_true', dest="debug",
                    help="call IPython.embed when running")
parser.add_argument('--check', action='store_true', dest="check",
                    help="check the solution")
# args for BKZ
parser.add_argument('--block-size', type=int, dest="block_size", default=20,
                    help="BKZ block size")

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
0: none
1: basic param
2: print norm of first two vector
3: matrix_overview of reduced basis
4: 
5: verbose when BKZ
"""
VERBOSE = args.verbose
SEED = args.seed or int.from_bytes(os.urandom(8), 'big')
if VERBOSE >= 1:
    print(f"SEED: {SEED}\n")
set_random_seed(SEED)


if (args.category is None) or (args.level is None):
    if not args.experiment:
        raise ValueError("Undefined behavior")
    COEFFS = [randint(0, m-1) for _ in range(n)]
    INIT_STATE = [randint(0, m-1) for _ in range(n)]
    if VERBOSE >= 1:
        print(f"coefficients: {COEFFS}")
        print(f"initial state: {INIT_STATE}")
        print()
    STATE = copy(INIT_STATE)
    N = r + t - 1
    for j in range(n, N):
        a_j = sum(c_i*a_i for c_i, a_i in zip(COEFFS, STATE[-n:])) % m
        STATE.append(a_j)
    y_ = [a_i >> zbits for a_i in STATE]
    z_ = [a_i % 2**zbits for a_i in STATE]
    SOL = tuple(COEFFS)

    Q = Matrix(Zmod(m), n)
    for i in range(n):
        Q[i, n-1] = COEFFS[i]
        if i == 0:
            continue
        Q[i, i-1] = 1

    Q_ = Q.powers(r)
    q_ = [Qj.change_ring(ZZ)[:, 0].T[0].list() for Qj in Q_]

else:
    y_ = read_data(args.category, args.level)
    N = len(y_)
    if N < r+t-1:
        raise ValueError(f"outputs is not enough (got {N}, expect >={r+t-1})")
    SOL = None


MM = 1 << (mbits - zbits)
BB = ceil((2*MM*r)**(t/(r-t)))
KK = ceil(sqrt(r)*2**((r-1)/2) * BB)

bkz_flags = BKZ.DEFAULT | BKZ.AUTO_ABORT
if VERBOSE >= 5:
    bkz_flags |= BKZ.VERBOSE

M = [[0]*(t+r) for _ in range(r)]
for i in range(r):
    for j in range(t):
        M[i][j] = y_[i+j] *KK
    M[i][t+i] = 1

if DEBUG and VERBOSE >= 3:
    matrix_overview(M)


B = IntegerMatrix.from_matrix(M)
BKZ.reduction(B, BKZ.EasyParam(block_size=min(B.nrows, args.block_size), flags=bkz_flags))
if VERBOSE >= 3:
    matrix_overview(B)


expect_vectors = 2
ETA = []

for i in range(expect_vectors):
    b = list(B[i])
    if any(b[:t]):
        raise ValueError(r"we need `\sum_i \eta_i Y_i = 0`")
    eta = list(b[t:])

    if VERBOSE >= 2:
        if SOL is not None:
            print(i, sum(e*e for e in eta)//r, sum(e*a for e, a in zip(eta, STATE)),
                  sum(e*y for e, y in zip(eta, y_)), sum(e*z for e, z in zip(eta, z_)))
        else:
            print(i, sum(e*e for e in eta)//r)

    if SOL is not None:
        if sum(e*a for e, a in zip(eta, STATE)) != 0:
            raise ValueError(r"we need `\sum_i \eta_i A_i = 0`")

    if SOL is not None:
        ETA.append(eta)

if VERBOSE >= 2:
    print()



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


if DEBUG and input("embed? "):
    IPython.embed()


"""level 1
 % sage recover_coefficients__resultant.sage 2147483647 2 30 8 17 --experiment --verbose 2 --block-size 10
SEED: 9714827701253331061

coefficients: [2109157988, 1364229973]
initial state: [234746902, 1509028418]

0 38 0 0 0
1 40 0 0 0

known roots: []
compute resultant of polys[2] (28) and polys[0] (28) with respect to `C1`
compute resultant of polys[2] (28) and polys[1] (28) with respect to `C1`
compute GDD
known roots: [2109157988]
compute GDD

coefficients: (2109157988, 1364229973)
"""

"""level 2
 % sage recover_coefficients__resultant.sage 2147483647 2 60 15 23 --experiment --verbose 2 --block-size 10
SEED: 7783158384627216085

coefficients: [1393533710, 1194216880]
initial state: [1864857680, 298943375]

0 3 0 0 0
1 3 0 0 0

known roots: []
compute resultant of polys[2] (56) and polys[0] (56) with respect to `C1`
compute resultant of polys[2] (56) and polys[1] (56) with respect to `C1`
compute GDD
known roots: [1393533710]
compute GDD

coefficients: (1393533710, 1194216880)
"""
