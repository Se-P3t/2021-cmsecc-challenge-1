# -*- coding: utf-8 -*-
"""
"""
from copy import copy
from sage.all import Zmod, PolynomialRing

def solve_with_resultant(polys, m, verbose=0, root = None):
    """
    solve multivariate polynomial system with resultant

    :param polys: polynomials that have roots over `Zmod(m)`,
        while this variable should be defined over integer ring
    :param m: modulus
    :param verbose: verbose level (default: 0)
    :param root: do NOT change this parameter
    """
    poly = polys[0]
    n = poly.nvariables()
    x_ = poly.parent().gens()

    if root is None:
        root = []
        if len(polys) != n + 1:
            # we need one more polynomial to caculate the GCD
            raise ValueError("polys must have the length of `pol.nvariables() + 1`")

    if verbose >= 1:
        print(f"known roots: {root}")

    len_known = len(root)
    polys = copy(polys)
    polys = [h(*(root + list(x_[len_known:]))) for h in polys]

    PRm = PolynomialRing(Zmod(m), 'x')

    if len_known == n - 1:

        if verbose >= 2:
            print(f"compute GDD")

        f1 = polys[0].univariate_polynomial(PRm)
        f2 = polys[1].univariate_polynomial(PRm)
        g = f1.gcd(f2)

        if g.degree() == 1:
            x0 = int(-g.constant_coefficient())
            yield tuple([*root, x0])

        else:
            if verbose >= 3:
                print(f"g = {g}")
            raise RuntimeError("no solution")

    else:

        polys0 = copy(polys)

        for j in range(n - len_known, 1, -1):
            new_polys = []
            hj = polys[j]

            for i in range(j):

                if verbose >= 2:
                    print(f"compute resultant of polys[{j}] and polys[{i}] with respect to `{x_[j-1]}`")

                new_hi = polys[i].resultant(hj, x_[j-1])
                if new_hi.is_constant():
                    raise RuntimeError("polys must be algebraically independent")
                new_polys.append(new_hi)

            polys = new_polys

        if verbose >= 2:
            print(f"compute GDD")

        f1 = polys[0].univariate_polynomial(PRm)
        f2 = polys[1].univariate_polynomial(PRm)
        g = f1.gcd(f2)

        if g.degree() == 1:
            x0 = int(-g.constant_coefficient())
            yield from solve_with_resultant(
                polys=polys0[:-1],
                m=m,
                verbose=verbose,
                root=(root + [x0]),
            )

        else:
            if verbose >= 3:
                print(f"g = {g}")
            raise RuntimeError("no solution")

