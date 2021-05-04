N = 40

m = 2**31 - 1
n = 16
PR = PolynomialRing(ZZ, names=[f'c{i}' for i in range(n)]+[f"a{i}" for i in range(n)])
var = PR.gens()
c_, a_ = var[:n], var[n:]

"""
Q = Matrix(PR, n)
for i in range(n):
    Q[i, n-1] = c_[i]
    if i == 0:
        continue
    Q[i, i-1] = 1

q_ = [None] * n
Qj = Q**n
for j in range(n, N):
    print(f"\rQ^{j}", end='')
    q_.append(Qj[:, 0].T[0].list())
    Qj *= Q
print()
"""

state = list(a_)
for j in range(n, N):
    print(f"\ra_{j}", end='')
    a_j = sum(c*a for c, a in zip(c_, state[-n:]))
    state.append(a_j)
print()
