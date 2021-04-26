# -*- coding: utf-8 -*-
"""
recover initial state
"""
import os
import argparse
import IPython
import random
import subprocess
from copy import copy

from sympy.matrices import Matrix, zeros, eye
from fpylll import FPLLL, IntegerMatrix, BKZ

from util import read_data, save_solution, matrix_overview
from sieve_asvp import solve_asvp

DEBUG = 0


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# chall params
parser.add_argument('m', type=int, help="modulus")
parser.add_argument('coeffs', type=str, help="coefficients `c_0,c_1,...,c_{n-1}`"
                    "of feedback polynomial (sep by comma)")
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

# FPLLL.set_threads(args.threads)

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
100: print shffule infos
"""
VERBOSE = args.verbose
THREADS = args.threads
SIEVE = args.sieve
SEED = args.seed or int.from_bytes(os.urandom(8), 'big')
if VERBOSE >= 1:
    print(f"SEED: {SEED}\n")

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
if VERBOSE >= 100: # rm
    random.seed(SEED)
    index = list(range(d+1))
    random.shuffle(index)
    print(f"shuffle: {index}")
    matrix_overview(M)


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
    data = {
        'modulus': m,
        'coefficients': coeffs,
        'initial_state': a_,
    }
    save_solution(args.category, args.level, data)
    exit(0)

########## sieve ##########
if not SIEVE:
    exit(-1)



B = solve_asvp(
    B,
    threads=THREADS,
    verbose=(DEBUG or VERBOSE >= 2),
    seed=SEED,
    step_size=args.step_size,
    trials=args.trials,
    workers=args.workers,
    workout__dim4free_dec=args.workout__dim4free_dec,
    goal_r0__gh=args.goal_r0__gh
)

if DEBUG or VERBOSE >= 3:
    matrix_overview(B)

if DEBUG and input('embed? '):
    IPython.embed()

a_ = find_solution(B)
if SOL is None and a_:
    data = {
        'modulus': m,
        'coefficients': coeffs,
        'initial_state': a_,
    }
    save_solution(args.category, args.level, data)
    exit(0)




"""
$ time python recover_initial_state.py 2147483647 "257, 0, 0, 0, 1048576, 0, 0, 0, 0, 0, 2097152, 0, 0, 131072, 0, 32768" 22 2 --category 1 --level 1 --verbose 1
SEED: 6082509956042779895

row 2/23:: True

solution: [257122259, 561754033, 340567466, 370127561, 1603890151, 789710361, 438276282, 1205614745, 929435387, 273101854, 1330188497, 1927005651, 70974738, 512638222, 655376420, 1799812840]

real    0m1.625s
user    0m2.578s
sys     0m1.797s


$ time python recover_initial_state.py 2147483647 "257, 0, 0, 0, 1048576, 0, 0, 0, 0, 0, 2097152, 0, 0, 131072, 0, 32768" 23 3 --category 1 --level 2 --verbose 1
SEED: 13355306573255156583

row 1/24:: True

solution: [720759561, 571519362, 1471239584, 2146374199, 1943406434, 1021387208, 911847040, 972284090, 1269846527, 2125700056, 582856260, 1326807684, 820744516, 83875554, 1033926054, 974403091]

real    0m1.653s
user    0m2.359s
sys     0m1.859s


$ time python recover_initial_state.py 2147483647 "257, 0, 0, 0, 1048576, 0, 0, 0, 0, 0, 2097152, 0, 0, 131072, 0, 32768" 30 14 --category 1 --level 3 --verbose 1
SEED: 3103757694693629835

row 1/31:: True

solution: [2056319062, 1175325906, 34344908, 300877541, 871395921, 953051611, 1276817066, 429446330, 716236050, 1498148665, 419137803, 983513800, 877501144, 162784581, 615817441, 1532450981]

real    0m1.665s
user    0m2.422s
sys     0m1.938s


$ time python recover_initial_state.py 2147483647 "257, 0, 0, 0, 1048576, 0, 0, 0, 0, 0, 2097152, 0, 0, 131072, 0, 32768" 51 21 --category 1 --level 4 --verbose 1 --block-size 29
SEED: 14561653990545178586

row 2/52:: True

solution: [146123803, 1660954690, 553686861, 1592631770, 2039784960, 874444650, 1462760700, 1629573947, 927239148, 2020986341, 2134682761, 1440980008, 214415113, 823589071, 1178840115, 237668181]

real    0m2.012s
user    0m2.828s
sys     0m1.906s


$ time python recover_initial_state.py 2147483647 "257, 0, 0, 0, 1048576, 0, 0, 0, 0, 0, 2097152, 0, 0, 131072, 0, 32768" 75 24 --category 1 --level 5 --verbose 1 --block-size 38
SEED: 17451625552935355143

row 1/76:: True

solution: [1683538635, 247852287, 1110454234, 2134379965, 524624671, 1135659173, 611030817, 917303344, 67431860, 844532212, 2121384016, 1630029172, 1311537205, 1289944654, 358213366, 2024311471]

real    0m5.382s
user    0m6.016s
sys     0m1.719s


$ time python recover_initial_state.py 2147483647 "257, 0, 0, 0, 1048576, 0, 0, 0, 0, 0, 2097152, 0, 0, 131072, 0, 32768" 128 26 --category 1 --level 6 --verbose 2 --threads 6 --sieve
SEED: 3039381073685026565

not found

Loaded file 'svpchallenge-129.txt'
gh = 657627253642393944064.000000, goal_r0/gh = 1.102500, r0/gh = 4.671924
  42: ↑ 42 ↓  8  T:   0.39783s, TT:   0.39788s, q:  14.51358 r0/gh:   3.39586
  45: ↑ 45 ↓ 12  T:   0.48302s, TT:   0.88092s, q:  13.60473 r0/gh:   3.39586
  48: ↑ 48 ↓ 24  T:   0.62368s, TT:   1.50462s, q:  11.94789 r0/gh:   3.32193
  51: ↑ 51 ↓ 20  T:   0.80179s, TT:   2.30643s, q:  11.22410 r0/gh:   3.32193
  54: ↑ 54 ↓ 21  T:   0.92295s, TT:   3.22939s, q:  10.11974 r0/gh:   3.32193
  57: ↑ 57 ↓ 22  T:   1.13283s, TT:   4.36224s, q:   8.95521 r0/gh:   3.32193
  60: ↑ 60 ↓ 34  T:   1.34190s, TT:   5.70416s, q:   7.57498 r0/gh:   3.32193
  63: ↑ 63 ↓ 35  T:   1.95294s, TT:   7.65712s, q:   6.86204 r0/gh:   2.90969
  66: ↑ 66 ↓ 39  T:   2.91808s, TT:  10.57523s, q:   6.02354 r0/gh:   2.90969
  69: ↑ 69 ↓ 41  T:   4.60971s, TT:  15.18496s, q:   5.13003 r0/gh:   2.59308
  72: ↑ 72 ↓ 46  T:   9.33832s, TT:  24.52330s, q:   4.57562 r0/gh:   2.59308
  75: ↑ 75 ↓ 44  T:  16.49760s, TT:  41.02092s, q:   4.10923 r0/gh:   2.07538
  78: ↑ 74       T:  11.29234s, TT:  52.31329s, q:   3.72726 r0/gh:   0.31697
svp: norm 14437807996.7 ,hf 0.56300
'threads': 6,        :: n: 129, cputime 258.4375s, walltime: 52.3138s, flast: 51.00, |db|: 2^17.24
row 1/129:: True

solution: [512741665, 1096178369, 808807049, 608491186, 420891056, 1682771835, 358966452, 55989687, 890238631, 2137448551, 2058494244, 38743896, 96170410, 49379854, 1551639435, 1878181314]

real    0m58.058s
user    4m23.297s
sys     0m3.313s

"""
