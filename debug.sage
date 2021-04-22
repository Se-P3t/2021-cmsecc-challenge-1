

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


m = 2^31-1
mbits = m.nbits()
PR.<x> = PolynomialRing(Zmod(m))
f = x^16 - 2^15*x^15 - 2^17*x^13 - 2^21*x^10 - 2^20*x^4 - 2^8 - 1
n = f.degree()
coeffs = (x^n - f).coefficients(sparse=False)
coeffs = vector(Zmod(m), coeffs)

Q = Matrix(Zmod(m), n)
for i in range(n):
    Q[i, n-1] = coeffs[i]
    if i == 0:
        continue
    Q[i, i-1] = 1

def q_jl(Q, j):
    return (Q^j).change_ring(ZZ).column(0)


d = 75
ybits = 7
zbits = mbits - ybits # 24

state = [randint(0, m-1) for _ in range(n)]
for _ in range(d-n):
    out = ZZ(coeffs * vector(Zmod(m), state[-n:]))
    state.append(out)

y_ = [a>>zbits for a in state]
z_ = [a % (2^zbits) for a in state]



M = Matrix(ZZ, d-n+2)
M[0, 0] = 1 #* m // 2
M[1, 1] = 1 #* (m // (2^zbits-1) // n)
for j in range(n, d):
    col = j - (n-2)
    q_j = q_jl(Q, j)
    c_j = 2^zbits*(y_[j] - sum(q_j[l]*y_[l] for l in range(n)))
    M[0, col] = (c_j % m) #* (m // (2^zbits-1))
    M[1, col] = (sum(q_j) % m) #* (m // (2^zbits-1))
    M[col, col] = m #* (m // (2^zbits-1))

matrix_overview(M, level=0)

B = M.BKZ(block_size=20)

matrix_overview(B, level=0)


s = sum(z_[:n])
for i in range(B.nrows()):
    if abs(B[i, 0]) == 1: #m // 2:
        print(i, abs(B[i, 1]) == s) #* (m // (2^zbits-1) // n))
