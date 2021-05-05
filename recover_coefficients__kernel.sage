# -*- coding: utf-8 -*-
"""
recover coefficients
    find `z_i` in kernel(ETA)
"""
import os
import json
import argparse
import IPython

from fpylll import FPLLL, IntegerMatrix, BKZ, LLL

from util import read_data, save_solution, matrix_overview


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
0: no output
1: basic param
2: print norm
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

ETA = []

USE_MODULUS = False
if USE_MODULUS:
    M = [[0]*(t+r) for _ in range(r)]
else:
    M = [[0]*(t+r) for _ in range(r)]
for i in range(r):
    for j in range(t):
        M[i][j] = y_[i+j] *KK
    M[i][t+i] = 1
if USE_MODULUS:
    for i in range(t):
        M[r+i][i] = m

if DEBUG and VERBOSE >= 3:
    matrix_overview(M)


B = IntegerMatrix.from_matrix(M)
BKZ.reduction(B, BKZ.EasyParam(block_size=min(B.nrows, args.block_size), flags=bkz_flags))
if VERBOSE >= 3:
    matrix_overview(B)

expect_vectors = ceil((r+t-1)/t)#+1
if VERBOSE >= 1:
    print(f"expect_vectors: {expect_vectors}")
    print()

for i in range(B.nrows):
    b = list(B[i])
    if any(b[:t]):
        if i >= expect_vectors:
            break
        raise ValueError("we need `\sum_i \eta_i Y_i = 0`")
    eta = list(b[t:])

    if VERBOSE >= 2:
        if SOL is not None:
            print(i, sum(e*e for e in eta)//r, sum(e*a for e, a in zip(eta, STATE)),
                  sum(e*y for e, y in zip(eta, y_)), sum(e*z for e, z in zip(eta, z_)))
        else:
            print(i, sum(e*e for e in eta)//r)

    if SOL is not None:
        if sum(e*a for e, a in zip(eta, STATE)) != 0:
            if i >= expect_vectors:
                break
            raise ValueError("we need `\sum_i \eta_i A_i = 0`")

    if SOL is not None or i < expect_vectors:
        ETA.append(eta)

if VERBOSE >= 2:
    print()

ETA_m = IntegerMatrix.from_matrix(ETA, int_type='mpz')
LLL.reduction(ETA_m)
ETA = [list(b) for b in ETA_m if any(list(b))]

if SOL is not None:
    assert len(ETA) >= expect_vectors, f"rank ETA: {len(ETA)}"

M = Matrix(ZZ, r+t-1, t*expect_vectors)
for j in range(t):
    for i in range(r):
        for a in range(expect_vectors):
            M[j+i, expect_vectors*j+a] = ETA[a][i]
B = M.left_kernel(basis="LLL").basis_matrix()

nrows = B.nrows()
if VERBOSE >= 1:
    print(f"kernel rank: {nrows}")
    print()
assert 0 < nrows <= 2

if DEBUG and input('embed? '):
    IPython.embed()
    if input("exit? "):
        exit(0)

b = B[0]
if b[0] < 0:
    b *= -1
assert min(b) >= 0 and max(b) < 2**zbits, f"range b: {min(b)}, {max(b)}"

a_ = [2**zbits*y + int(z) for y, z in zip(y_, b)]

if SOL is not None:
    res = (tuple(STATE) == tuple(a_))
    print("find z_i:", res)
    if not res:
        exit(-1)

N = N - 1 # leave one pair to check the solution

M = [[0]*(N+1) for _ in range(N+1)]
for j in range(N-n):
    for i in range(n+1):
        M[i][j] = a_[j+i] *KK
    M[n+1+j][j] = m *KK
for i in range(n+1):
    M[i][N-n+i] = 1
M[n][N] = m >> 1

if DEBUG and VERBOSE >= 3:
    matrix_overview(M)

B = IntegerMatrix.from_matrix(M)
#BKZ.reduction(B, BKZ.EasyParam(block_size=min(20, args.block_size), flags=bkz_flags))
LLL.reduction(B)

if DEBUG or VERBOSE >= 3:
    matrix_overview(B)

for b in B:
    b = list(b)
    if any(b[:N-n]):
        continue
    if abs(b[N]) != m >> 1:
        continue
    if b[N] == m >> 1:
        c_ = [int(-c_i)%m for c_i in b[N-n:N]]
    else:
        c_ = [int(c_i)%m for c_i in b[N-n:N]]

    if SOL is not None:
        if tuple(SOL) == tuple(COEFFS):
            print("found the coefficients")
            break
    else:
        if sum(c*a for c, a in zip(c_, a_[-n-1:-1]))%m == a_[-1]:
            print(c_)
            break
else:
    raise ValueError("not found")

if SOL is None:
    solution = {
        'modulus': m,
        'coefficients': c_,
        'initial_state': tuple(a_[:n]),
    }
    save_solution(args.category, args.level, solution)


if DEBUG and input('embed? '):
    IPython.embed()
