# -*- coding: utf-8 -*-
"""
recover coefficients
    find `z_i` in kernel(ETA)
"""
import os
import sys
import random
import argparse
import itertools
import subprocess

import IPython

from util import read_data, save_solution, matrix_overview
from lattice import Lattice


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

USE_SUBS = False
if 2*zbits > mbits + 2: # roughly check
    USE_SUBS = True

if not r > t > n:
    raise ValueError("r > t > n")


""" verbose level
0: solution bound check
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
random.seed(SEED)


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

Ls = []
L = Lattice(r, t+r)
for i in range(r):
    for j in range(t):
        L[i, j] = y_[i+j] *KK
    L[i, t+i] = 1
Ls.append(L)
if USE_SUBS:
    L1 = Lattice(r, t+r)
    for i in range(r):
        for j in range(t):
            L1[i, j] = (y_[i+j]*2+1) *KK
        L1[i, t+i] = 1
    Ls.append(L1)

if DEBUG and VERBOSE >= 3:
    matrix_overview(L)


for L in Ls:
    L.randomize_block()
    L.run_bkz(block_size=args.block_size, verbose=VERBOSE >= 5)
if VERBOSE >= 3:
    matrix_overview(L)

expect_vectors = ceil((r+t-1)/t)#+1
if VERBOSE >= 1:
    print(f"expect_vectors: {expect_vectors}")
    print()

ETAs = []
if SOL is None:
    for L in Ls:
        ETA = []

        for i in range(expect_vectors):
            b = list(L[i])
            if any(b[:t]):
                raise ValueError(r"we need `\sum_i \eta_i Y_i = 0`")

            eta = list(b[t:])

            if VERBOSE >= 2:
                print(i, sqrt(sum(e*e for e in eta)/r).n())

            ETA.append(eta)

        ETAs.append(ETA)

else: # experiment
    for L in Ls:
        ETA = []

        for i, b in enumerate(L):
            b = list(b)

            if any(b[:t]):
                if i >= expect_vectors:
                    break
                raise ValueError(r"we need `\sum_i \eta_i Y_i = 0`")

            eta = list(b[t:])
            sum_e_a = sum(e*a for e, a in zip(eta, STATE))

            if VERBOSE >= 2:
                print(
                    i,
                    sqrt(sum(e*e for e in eta)/r).n(),
                    sum_e_a,
                    sum(e*y for e, y in zip(eta, y_)),
                    sum(e*z for e, z in zip(eta, z_))
                )

            if sum_e_a != 0:
                if i >= expect_vectors:
                    break
                raise ValueError(r"we need `\sum_i \eta_i A_i = 0`")
            ETA.append(eta)

        ETAs.append(ETA)

if VERBOSE >= 2:
    print()

if SOL is not None:
    for ETA in ETAs:
        if len(ETA) < expect_vectors:
            raise RuntimeError(f"rank ETA: {len(ETA)}")


M = Matrix(ZZ, r+t-1, t*expect_vectors)
if USE_SUBS:
    M1 = Matrix(ZZ, r+t-1, t*expect_vectors)
for j, i, a in itertools.product(
        range(t), range(r), range(expect_vectors)):
    M[j+i, expect_vectors*j+a] = ETAs[0][a][i]
    if USE_SUBS:
        M1[j+i, expect_vectors*j+a] = ETAs[1][a][i]
B = M.left_kernel(basis="LLL").basis_matrix()
if USE_SUBS:
    B1 = M1.left_kernel(basis="LLL").basis_matrix()

nrows = B.nrows()
if USE_SUBS:
    nrows1 = B1.nrows()
if VERBOSE >= 1:
    print(f"kernel rank: {nrows}")
    if USE_SUBS:
        print(f"kernel' rank: {nrows1}")
    print()

if DEBUG and input('embed? '):
    IPython.embed()
    if input("exit? "):
        sys.exit(int(0))

assert nrows == 2, r"some `\eta` are wrong"
if USE_SUBS:
    assert nrows1 == 2, r"some `\eta` are wrong"

if USE_SUBS:
    for code in range(1 << 3):
        M = Matrix(ZZ, [
            list(B[0] if B[0,0] > 0 else -B[0]), # y
            list(-B[1] if (code >> 1)&1 else B[1]), # \pm (z + k_1 y)
            list(e-2**(zbits-1) if code & 1 else e+2**(zbits-1)
                for e in B1[1]), # \pm (z + k_2(2y+1))
            [1]*B.ncols() # 1
        ])
        v = M.left_kernel().basis_matrix()
        assert v.nrows() == 1
        v = copy(v[0])
        assert abs(v[1]) == abs(v[2]) == 1
        if v[1] == v[2]:
            continue
        v *= v[1]
        k2 = v[3]
        k1 = -2*k2 + v[0] if (code >> 2)&1 else 2*k2 - v[0]
        b = B[1] - k1*B[0]
        if b[0] < 0:
            b *= -1
        if min(b) < 0 or max(b) >= 2**zbits:
            continue
        try:
            _ = B1.solve_left(vector(ZZ, [zi-2**(zbits-1) for zi in b]))
            print("find k:", k1, k2)
            break
        except ValueError:
            continue
    else:
        print("cannot find `z_i`")
        if DEBUG and input('embed? '):
            IPython.embed()
        sys.exit(int(-1))
else:
    b = B[0]
    if b[0] < 0:
        b *= -1
print("checking z_i:", end=' ')
if SOL is not None:
    res = (tuple(z_) == tuple(b))
    print(res)
    if not res:
        if DEBUG and input('embed? '):
            IPython.embed()
        sys.exit(int(-1))
else:
    if min(b) >= 0 and max(b) < 2**zbits:
        print("maybe")
    else:
        print(f"False :: range {min(b)}, {max(b)}")
        sys.exit(int(-1))

a_ = [int(2**zbits*y + z) for y, z in zip(y_, b)]

A0 = Matrix(Zmod(m), n)
A1 = Matrix(Zmod(m), n)
for i in range(n):
    for j in range(n):
        A0[i, j] = a_[i+j]
        A1[i, j] = a_[i+j+1]
Q1 = A0.solve_right(A1)
print("checking c_i:", end=' ')
if SOL is not None:
    res = (Q1 == Q)
    print(res)
    if DEBUG and input('embed? '):
        IPython.embed()
    if not res:
        sys.exit(int(-1))
else:
    Q = Q1.change_ring(ZZ)
    c_ = []
    for i in range(n):
        for j in range(n-1):
            entry = 1 if i == j+1 else 0
            if Q[i, j] != entry:
                print("False")
                sys.exit(int(-1))
        c_.append(int(Q[i, n-1]))
    print("True")

    if SOL is None:
        solution = {
            'modulus': int(m),
            'zbits': zbits,
            'coefficients': c_,
            'initial_state': tuple(a_[:n]),
        }
        save_solution(args.category, args.level, solution)

    if args.check:
        _ = subprocess.check_call(f"python check.py {args.category} {args.level}", shell=True)
