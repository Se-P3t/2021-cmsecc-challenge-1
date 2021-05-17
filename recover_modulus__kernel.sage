# -*- coding: utf-8 -*-
"""
recover modulus
    find `z_i` in kernel(ETA)
"""
import os
import argparse
import subprocess
import IPython

from fpylll import FPLLL, IntegerMatrix, BKZ

from util import read_data, save_solution, matrix_overview, str_mat
from hlll import hlll_wrapper
from sieve_asvp import solve_asvp


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# chall params
parser.add_argument('mbits', type=int, help="number of bits of the modulus")
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
parser.add_argument('--check', action='store_true', dest="check",
                    help="check the solution")
# args for hpLLL
parser.add_argument('--threshold', type=int, dest="threshold", default=40,
                    help="HLLL: block size threshold for parallelism")
# args for BKZ
parser.add_argument('--block-size', type=int, dest="block_size", default=20,
                    help="BKZ: block size")
# args for sieving
parser.add_argument('--max-dim', type=int, dest="max_dim", default=128,
                    help="sieve max dim")
parser.add_argument('--timeout', type=int, dest="timeout", default=300,
                    help="sieve timeout")
parser.add_argument('--workout/dim4free-dec', type=int, dest="workout__dim4free_dec", default=3,
                    help="By how much do we decreaseee dim4free at each iteration")

args, _ = parser.parse_known_args()


DEBUG = args.debug

mbits = args.mbits
n = args.n
r = args.r
t = args.t
zbits = args.zbits

USE_SUBS = False
if 2*zbits > mbits + 2: # roughly check
    USE_SUBS = True

if not (r > t > n):
    raise ValueError("r > t > n")


""" verbose level
0: solution bound check
1: basic param
2: verbose for sieve
3: print norm
4: 
5: verbose when BKZ
"""
VERBOSE = args.verbose
THREADS = args.threads
SEED = args.seed or int.from_bytes(os.urandom(8), 'big')
if VERBOSE >= 1:
    print(f"SEED: {SEED}")
    print(f"USE_SUBS: {USE_SUBS}")
    print()
set_random_seed(SEED)
FPLLL.set_threads(THREADS)


if (args.category is None) or (args.level is None):
    if not args.experiment:
        raise ValueError("Undefined behavior")
    m = random_prime(2**mbits, lbound=2**(mbits-1))
    COEFFS = [randint(0, m-1) for _ in range(n)]
    INIT_STATE = [randint(0, m-1) for _ in range(n)]
    if VERBOSE >= 1:
        print(f"modulus: {m}")
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


if DEBUG and args.experiment and input('embed? '):
    tmpfile = f"exp-{SEED}-{mbits}-{n}-{r}-{t}-{zbits}.mat"
    if os.path.exists(tmpfile):
        A = IntegerMatrix.from_file(tmpfile)
        print(f"Loaded tmpfile {tmpfile} :: {A.nrows} x {A.ncols}")
    IPython.embed()
    if input("exit? "):
        exit(int(0))



bkz_flags = BKZ.DEFAULT | BKZ.AUTO_ABORT
if VERBOSE >= 5:
    bkz_flags |= BKZ.VERBOSE

sieve_param = {
    'threads': THREADS,
    'verbose': (DEBUG or VERBOSE >= 2),
    'seed': SEED,
    #'dry_run': True,
    'workout/dim4free_dec': args.workout__dim4free_dec,
    'workout/save_prefix': f'logs/sieve-{args.category}-{args.level}',
    'workout/start_n': args.block_size+1,
    'workout/goal_r0': 0,
    'workout/dim4free_min': r-t-args.max_dim,
    'workout/timeout': args.timeout,
    'pump/down_sieve': True,
}

def left_kernel_lll(M:list, expect_rank) -> list:
    r = len(M)
    t = len(M[0])
    MM = max(M[0])
    BB = ceil((2*MM*r)**(t/(r-t)))
    KK = ceil(sqrt(r)*2**((r-1)/2) * BB)

    M = block_matrix(ZZ,
        [[Matrix(ZZ, M)*KK, identity_matrix(r)]]
    )
    M = [[M[i,j] for j in range(M.ncols())] for i in range(M.nrows())]

    shuffle(M)
    B = hlll_wrapper(M, threads=THREADS, threshold=args.threshold)
    M = []
    for i, row in enumerate(B):
        if any(row[:t]):
            if i < expect_rank:
                # LLL not enough
                #matrix_overview(B)
                if VERBOSE >= 1:
                    print(f"got non-zero at row {i}")
                M = B
                break
            if VERBOSE >= 1:
                print(f"find the kernel (rank {expect_rank})\n")
            if expect_rank != 2:
                tmpfile = f"eta-{args.category}-{args.level}.mat"
                if args.experiment:
                    tmpfile = f"exp-{SEED}-{mbits}-{n}-{r}-{t}-{zbits}.mat"
                open(tmpfile, "w").write(
                    f"[{str_mat(M)}]"
                )
            return M
        M.append(list(row[t:]))

    # more precisely, there is at least one \eta not work
    raise RuntimeError("cannot find the kernel, expect rank too large")



M = [[0]*t for _ in range(r)]
if USE_SUBS:
    M1 = [[0]*t for _ in range(r)]
for i in range(r):
    for j in range(t):
        M[i][j] = y_[i+j]
        if USE_SUBS:
            M1[i][j] = (y_[i+j]*2+1)

B = left_kernel_lll(M, r-t)
if USE_SUBS:
    B1 = left_kernel_lll(M1, r-t)

try:
    B = IntegerMatrix.from_matrix(B)
    BKZ.reduction(B,BKZ.EasyParam(block_size=args.block_size, flags=bkz_flags))
    if args.sieve:
        B = solve_asvp(B, **sieve_param)
    if USE_SUBS:
        B1 = IntegerMatrix.from_matrix(M1)
        BKZ.reduction(B1, BKZ.EasyParam(block_size=args.block_size, flags=bkz_flags))
        if args.sieve:
            B1 = solve_asvp(B1, **sieve_param)
except KeyboardInterrupt:
    IPython.embed()
    exit(int(-1))

expect_vectors = ceil((r+t-1)/t)#+1
if VERBOSE >= 1:
    print(f"expect_vectors: {expect_vectors}")
    print()

ETA = []
if USE_SUBS:
    ETA1 = []
for i in range(B.nrows):
    if SOL is None and i >= expect_vectors:
        break
    eta = list(B[i])
    if USE_SUBS:
        eta1 = list(B1[i])

    if VERBOSE >= 3:
        if SOL is not None:
            print(i, sum(e*e for e in eta)//r, sum(e*a for e, a in zip(eta, STATE)),
                  sum(e*y for e, y in zip(eta, y_)), sum(e*z for e, z in zip(eta, z_)))
        else:
            print(i, sum(e*e for e in eta)//r)

    if SOL is not None:
        if (sum(e*a for e, a in zip(eta, STATE)) != 0 or (
            USE_SUBS and sum(e*a for e, a in zip(eta1, STATE)) != 0
        )):
            if i > expect_vectors:
                break
            raise ValueError(r"we need `\sum_i \eta_i A_i = 0`")

    ETA.append(eta)
    if USE_SUBS:
        ETA1.append(eta1)

if VERBOSE >= 3:
    print()

if SOL is not None:
    assert len(ETA) >= expect_vectors, f"rank ETA: {len(ETA)}"
    if USE_SUBS:
        assert len(ETA1) >= expect_vectors, f"rank ETA': {len(ETA1)}"


M = [[0]*t*expect_vectors for _ in range(r+t-1)]
if USE_SUBS:
    M1 = [[0]*t*expect_vectors for _ in range(r+t-1)]
for j in range(t):
    for i in range(r):
        for a in range(expect_vectors):
            M[j+i][expect_vectors*j+a] = ETA[a][i]
            if USE_SUBS:
                M1[j+i][expect_vectors*j+a] = ETA1[a][i]
B = left_kernel_lll(M, 2)
B = Matrix(ZZ, B).LLL()
if USE_SUBS:
    B1 = left_kernel_lll(M1, 2)
    B1 = Matrix(ZZ, B1).LLL()

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
        exit(int(0))

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
        exit(int(-1))
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
        exit(int(-1))
else:
    if min(b) >= 0 and max(b) < 2**zbits:
        print("maybe")
    else:
        print(f"False :: range {min(b)}, {max(b)}")
        exit(int(-1))
         

a_ = [int(2**zbits*y + z) for y, z in zip(y_, b)]

print("finding modulus")
km = 0
for offset in range(N-2*n+2):
    A0 = Matrix(ZZ, n)
    A1 = Matrix(ZZ, n)
    A2 = Matrix(ZZ, n)
    for i in range(n):
        for j in range(n):
            A0[i, j] = a_[offset+i+j]
            A1[i, j] = a_[offset+i+j+1]
            A2[i, j] = a_[offset+i+j+2]
    tmp = (A1 * A0.inverse() * A1)[n-1, n-1]
    km = gcd(km, tmp.numerator() - A2[n-1, n-1]*tmp.denominator())
    bits = int(km).bit_length()
    print(f"\r{bits:3d} bits", end=' ')
    if bits == mbits:
        if SOL is None:
            print("maybe")
        else:
            res = (km == m)
            print(res)
            if not res:
                if DEBUG and input('embed? '):
                    IPython.embed()
                exit(int(-1))
        break
else:
    print("\rcannot find modulus")
    exit(int(-1))

m = int(km)

print("checking c_i:", end=' ')
A0 = Matrix(Zmod(m), n)
A1 = Matrix(Zmod(m), n)
for i in range(n):
    for j in range(n):
        A0[i, j] = a_[i+j]
        A1[i, j] = a_[i+j+1]
Q1 = A0.solve_right(A1)
if SOL is not None:
    res = (Q1 == Q)
    print(res)
    if DEBUG and input('embed? '):
        IPython.embed()
    if not res:
        exit(int(-1))

else:
    Q = Q1.change_ring(ZZ)
    c_ = []
    for i in range(n):
        for j in range(n-1):
            entry = 1 if i == j+1 else 0
            if Q[i, j] != entry:
                print("False")
                exit(int(-1))
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
