import re
import IPython
from random import SystemRandom

from fpylll import IntegerMatrix

from pylattice.algorithms.sieve_asvp import solve_asvp



RNG = SystemRandom()
set_verbose(2)
DEBUG = 0

with open('环上序列的截位还原问题——数据.txt', mode='r', encoding='GB2312') as f:
    DATA = f.read()

def read_data(category, level):
    assert 1 <= category <= 3
    assert 1 <= level <= 9
    pattern = r'挑战{}\.{}数据[^\n]+\n([^\n]+)\n'.format(category, level)
    match_text = re.search(pattern, DATA).groups()[0]
    data = list(map(int, match_text.split()))
    return data


def save_solution(category, level, data):
    with open(f"solutions/sol-{category}-{level}.txt", 'w') as f:
        for d in data:
            f.write(f"{d}\n")

def matrix_overview(BB, level=2, bound = None):
    """
    Print a matrix overview if the current verbosity is at least level.
    """
    if get_verbose() >= level:
        mid = ''
        if BB.nrows() < 60:
            mid = ' '
        print()
        for ii in range(BB.nrows()):
            a = f'{ii:03d} '
            for jj in range(BB.ncols()):
                B_ij = abs(BB[ii, jj])
                if B_ij == 0:
                    a += ' '
                elif B_ij == 1:
                    a += '1'
                else:
                    a += 'X'
                a += mid
            if bound is not None and BB[ii,ii] > bound:
                a += '~'
            print(a)
        print()


def shuffle_matrix(M):
    M = list(M)
    shuffle(M, random=RNG.random)
    M = Matrix(ZZ, M)
    return M


def sieve(B, threads=4, goal_r0__gh=1.05):
    A = IntegerMatrix.from_matrix(B)
    AA = solve_asvp(A, threads=threads, goal_r0__gh=goal_r0__gh)
    B = Matrix(ZZ, AA)
    return B


############## 第一类 ##############

category = 1

m = 2^31-1
mbits = m.nbits()
PR.<x> = PolynomialRing(Zmod(m))
f = x^16 - 2^15*x^15 - 2^17*x^13 - 2^21*x^10 - 2^20*x^4 - 2^8 - 1
n = f.degree()
coeffs = (x^n - f).coefficients(sparse=False)

Q = Matrix(Zmod(m), n)
for i in range(n):
    Q[i, n-1] = coeffs[i]
    if i == 0:
        continue
    Q[i, i-1] = 1

def q_jl(Q, j):
    return (Q^j).change_ring(ZZ).column(0)


def check(Q, initial_state, outputs, zbits):
    n = Q.nrows()
    m = Q.base_ring().characteristic()
    assert len(initial_state) >= n
    a = vector(Zmod(m), initial_state[:n])
    for yj in outputs[n:]:
        a *= Q
        yj_ = int(a[-1]) >> zbits
        if yj != yj_:
            return False
    return True

method = 0

scale = 1
bais = 1

def recover_initial_state(Q, y_, zbits, d = None, block_size=20):
    global scale, bais
    m = Q.base_ring().characteristic()
    n = Q.nrows()
    if d is None:
        d = len(y_)
    bound = floor(m^((d-n)/(d+1)))
    verbose(f"d = {d} :: bound = {bound}", level=1)
    #if 2^zbits-1 >= bound:
    #    print("unsolvable")
    #    return None

    if d > len(y_):
        raise ValueError("not enough")

    if method == 0:
        scale = m >> (zbits-1)
        bais = 2^(zbits-1)
    elif method == 1:
        scale = m // (2^zbits-1)
        bais = 2^(zbits-1)

    M = Matrix(ZZ, d+1)
    for j in range(d+1):
        if j < 1:
            if method == 0:
                M[j, j] = bais * scale
            elif method == 1:
                for i in range(n+1):
                    if i == j:
                        M[i, j] = -bais*(n+1) * scale
                        continue
                    M[i, j] = -1 * scale
        elif j < n+1:
            if method == 0:
                M[0, j] = bais * scale
                M[j, j] = 1 * scale
            elif method == 1:
                for i in range(n+1):
                    if i == j:
                        M[i, j] = (n+1) * scale
                        continue
                    if i == 0:
                        M[i, j] = bais * scale
                    else:
                        M[i, j] = -1 * scale
        elif j < d+1:
            M[j, j] = m * scale
            q_j = q_jl(Q, j-1)
            c_j = 2^zbits*(y_[j-1] - sum(q_j[l]*y_[l] for l in range(n)))
            c_j %= m
            if method == 0:
                M[0, j] = (c_j + 2^(zbits-1)) * scale
            elif method == 1:
                M[0, j] = c_j * scale
            for i in range(n):
                M[i+1, j] = q_j[i] * scale

    matrix_overview(M, level=2)

    if not DEBUG:
        M = shuffle_matrix(M)
        matrix_overview(M, level=3)

    #B = M.LLL()
    B = M.BKZ(block_size=block_size)
    matrix_overview(B, level=2)

    if DEBUG:
        IPython.embed()

    return B


def find_solution(B, Q, y_, zbits):
    nrows = B.nrows()
    d = B.ncols() - 2
    m = Q.base_ring().characteristic()
    n = Q.nrows()
    for idx, b in enumerate(B):
        if method == 1:
            b = copy(b)
            s = sum(b[:n+1])
            for i in range(n+1):
                b[i] = (b[i]+s) / (n+2)
            b[0] *= -1
        if abs(b[0]) != bais*scale:
            #print(idx, "b0 =", b[0])
            continue
        if b[0] == bais*scale:
            b *= -1
        if method == 0:
            z_ = list(int(bi/scale)+2^(zbits-1) for bi in b[1:n+1])
        else:
            z_ = list(map(int, b[1:n+1]/scale))
        if max(z_) >= 2^zbits:
            #print(idx, "big")
            continue
        a_ = tuple(2^zbits*yi + zi for yi, zi in zip(y_, z_))
        verbose(f"{idx+1}/{nrows}, {check(Q, a_, y_, zbits)}", level=1, caller_name="checking")
        if check(Q, a_, y_, zbits):
            verbose(str(z_), level=1, caller_name="solutions")
            return a_

    print("not found")
    return None





############## 1.1 ##############
"""
level = 1
ybits = 29
zbits = mbits - ybits # 2
y_ = read_data(category, level)
B = recover_initial_state(Q, y_, zbits)
a_ = find_solution(B, Q, y_, zbits)
if a_:
    save_solution(category, level, [a_])
"""

############## 1.2 ##############
"""
level = 2
ybits = 28
zbits = mbits - ybits # 3
y_ = read_data(category, level)
B = recover_initial_state(Q, y_, zbits)
a_ = find_solution(B, Q, y_, zbits)
if a_:
    save_solution(category, level, [a_])
"""

############## 1.3 ##############
"""
level = 3
ybits = 17
zbits = mbits - ybits # 14
y_ = read_data(category, level)
B = recover_initial_state(Q, y_, zbits)
a_ = find_solution(B, Q, y_, zbits)
if a_:
    save_solution(category, level, [a_])
"""

############## 1.4 ##############
"""
level = 4
ybits = 10
zbits = mbits - ybits # 21
y_ = read_data(category, level)
B = recover_initial_state(Q, y_, zbits)
a_ = find_solution(B, Q, y_, zbits)
if a_:
    save_solution(category, level, [a_])
"""

############## 1.5 ##############
"""
level = 5
ybits = 7
zbits = mbits - ybits # 24
y_ = read_data(category, level)
B = recover_initial_state(Q, y_, zbits)
a_ = find_solution(B, Q, y_, zbits)
if a_:
    save_solution(category, level, [a_])
else:
    B = sieve(B, 4, 0.95)
    a_ = find_solution(B, Q, y_, zbits)
    if a_:
        save_solution(category, level, [a_])
"""

############## 1.6 ##############
"""
level = 6
ybits = 5
zbits = mbits - ybits # 26
y_ = read_data(category, level)
B = recover_initial_state(Q, y_, zbits, d=128)
a_ = find_solution(B, Q, y_, zbits)
if a_:
    save_solution(category, level, [a_])
else:
    B = sieve(B, 4)
    a_ = find_solution(B, Q, y_, zbits)
    if a_:
        save_solution(category, level, [a_])
"""

############## 1.7 ##############
"""
level = 7
ybits = 4
zbits = mbits - ybits # 27
y_ = read_data(category, level)
#B = recover_initial_state(Q, y_, zbits, d=150)
#a_ = find_solution(B, Q, y_, zbits)
#if a_:
#    save_solution(category, level, [a_])
#else:
#    B = sieve(B, 6)
#    a_ = find_solution(B, Q, y_, zbits)
#    if a_:
#        save_solution(category, level, [a_])
"""
















"""
d = len(y_)

PR = PolynomialRing(Zmod(m), d, 'z')
PRZ = PR.change_ring(ZZ)
z_ = PR.gens()

monos = list(z_) + [PR(1)]
polys = []
for j in range(n, d):
    q_j = q_jl(Q, j)
    c_j = 2^zbits*(y_[j] - sum(q_j[l]*y_[l] for l in range(n)))
    f = sum(q_j[l]*z_[l] for l in range(n)) - z_[j] - c_j
    polys.append(f)


nrows = len(polys)
ncols = len(monos)
verbose(f"dim = ({nrows}, {ncols})", level=0, caller_name='param')

M = Matrix(ZZ, nrows, ncols)
for i in range(M.nrows()):
    for j in range(M.ncols()):
        M[i, j] = polys[i](*(z*(2^zbits-1) for z in z_)).monomial_coefficient(monos[j])


matrix_overview(M, level=0)

B = M.LLL()

matrix_overview(B, level=0)


v = vector([(mono)(*(z/(2^zbits-1) for z in z_)) for mono in monos])
hs = [PRZ(row*v) for row in B]

for i, h in enumerate(hs):
    coeffs = h(*(z*(2^zbits-1) for z in z_)).coefficients()
    w = len(coeffs)
    h_norm = sum(c*c for c in coeffs)
    print(i, w * h_norm < m**2)
"""
