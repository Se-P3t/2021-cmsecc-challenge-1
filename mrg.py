# pylint: disable=C0103
"""
recover initial state with truncated outputs from MRG using lattice attack
"""
import random
from functools import lru_cache

import mpmath
from sympy.matrices import zeros, eye
from fpylll import IntegerMatrix, LLL, BKZ

from hlll import hlll_wrapper


class MRG:
    """
    multiple recursive generator
    """

    def __init__(self, modulus, coeffs, initial_state=None):
        """
        Args:
            modulus (int): modulus
            coeffs (list-like): coefficients `c_0,c_1,...,c_{n-1}` of connection polynomial
            initial_state (list-like, optional): initial state of MRG
        """
        self.modulus = modulus
        self.mbits = int(modulus).bit_length()
        self.coeffs = tuple(c % modulus for c in coeffs)
        self.degree = len(coeffs)
        if initial_state is None:
            self.initial_state = None
        else:
            if len(initial_state) != self.degree:
                raise ValueError(
                    "length of initial_state must be "
                    f"{self.degree}, got {len(initial_state)}"
                )
            self.initial_state = tuple(initial_state)

    @classmethod
    def random(cls, modulus=0x7fffffff, degree=16):
        """
        Construct new random MRG

        Args:
            modulus (int, optional): Defaults to 0x7fffffff.
            degree (int, optional): Defaults to 16.

        Returns:
            `MRG`: a MRG object
        """
        if degree == 16:
            coeffs = [257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768]
        else:
            # TODO: generate irreducible polynomial
            pass

        initial_state = [random.randint(0, modulus-1) for _ in range(degree)]

        return cls(modulus, coeffs, initial_state)

    @property
    @lru_cache()
    def characteristic(self):
        """
        characteristic polynomial
        """
        poly = ''

        for i, coeff in enumerate(self.coeffs):
            if i == 0 and coeff:
                poly = f' - {coeff}'
                continue
            if coeff:
                poly = f' - {coeff}*x^{i}' + poly

        poly = f'f(x) = x^{self.degree}' + poly

        return poly

    def __repr__(self):
        name = "MRG"
        if self.initial_state is None:
            name += " with unknown state"
        return f"{name} defined by {self.characteristic} over Zmod({self.modulus})"

    @staticmethod
    def mat_pow_mod(A, e, m):
        """
        compute `A^e % m` for a given IntegerMatrix

        Args:
            A (IntegerMatrix): an IntegerMatrix
            e (int): power
            m (int): modulus

        Returns:
            IntegerMatrix:
        """
        #assert A.nrows == A.ncols
        B = eye(A.shape[0])
        while e:
            if e & 1:
                B *= A
                B %= m
            A *= A
            A %= m
            e >>= 1
        return B

    @lru_cache()
    def gen_qji(self, d):
        """
        generate `q_{j,i}`

        Args:
            d (int): number of MRG outputs
            coeffs (list-like): coefficients
            m (int): modulus

        Returns:
            tuple: a tuple of `q_j`
        """
        n = self.degree
        m = self.modulus

        Q = zeros(n)
        for i in range(n):
            Q[i, n-1] = self.coeffs[i]
            if i == 0:
                continue
            Q[i, i-1] = 1

        Qn = self.mat_pow_mod(Q, n, m)

        Q_ = [None] * n
        Qj = Qn
        for _ in range(n, d):
            Q_.extend(Qj[:, 0].transpose().tolist())
            Qj *= Q
            Qj %= m

        return tuple(Q_)

    def output(self, length, zbits=0, return_z=False):
        """
        Args:
            length (int): a positive integer
            zbits (int, optional): truncated bits. Defaults to 0.
            return_z (bool, optional): Defaults to False.

        Returns:
            tuple: MRG output `y_0, y_1, ..., y_{length-1}`
        """
        assert self.initial_state is not None
        state = list(self.initial_state)
        n = self.degree

        for _ in range(n, length):
            a_j = sum(c_i * a_i for c_i, a_i in zip(self.coeffs, state[-n:]))
            state.append(a_j)

        y_list = tuple(a_i >> zbits for a_i in state)
        z_list = tuple(a_i % 2**zbits for a_i in state)

        if return_z:
            return y_list, z_list
        return y_list


class MRGSolver:
    """
    Solve MRG with truncated outputs
    """

    def __init__(self, mrg, zbits, outputs):
        """
        Args:
            mrg (`MRG`): a MRG object
            zbits (int): truncated bits
            outputs (list-like): truncated outputs
        """
        self.mrg = mrg
        self.zbits = zbits
        self.outputs = tuple(outputs)
        self.length = len(outputs)
        self.L = None

    def gen_lattice(self, d):
        """
        Args:
            d (int): number of samples

        Raises:
            ValueError: outputs not enough

        Returns:
            IntegerMatrix:
        """
        if d > self.length:
            raise ValueError(f"outputs is not enough (got {self.length}, expect >={d})")

        y_ = self.outputs
        size = 1 << (self.zbits - 1)
        scale = self.mrg.modulus >> (self.zbits - 1)
        Q_ = self.mrg.gen_qji(self.length)

        indices = random.sample(range(self.mrg.degree, self.length), d - self.mrg.degree)

        A = IntegerMatrix(d + 1, d + 1, int_type="mpz")

        A[0, 0] = size * scale
        for j in range(self.mrg.degree):
            A[0, j+1] = size * scale
            A[j+1, j+1] = 1 * scale

        for j in range(d - self.mrg.degree):
            idx = indices[j]
            col = j + self.mrg.degree + 1
            c_j = 2**self.zbits * (y_[idx] - sum(
                int(Q_[idx][l])*y_[l] for l in range(self.mrg.degree)
            )) % self.mrg.modulus
            A[col, col] = self.mrg.modulus * scale
            A[0, col] = (c_j + 2**(self.zbits-1)) * scale
            for i in range(self.mrg.degree):
                A[i+1, col] = int(Q_[idx][i]) * scale

        return A

    def have_construct_lattice(self):
        if self.L is None:
            raise ValueError("lattice not construct, call `gen_lattice` first")

    def run_lll(self):
        """
        run LLL over lattice
        """
        self.have_construct_lattice()

        LLL.reduction(self.L)

    def run_hlll(self, threads=4, threshold=40, **kwds):
        """
        run hLLL over lattice

        Args:
            threads (int, optional): Defaults to 4.
            threshold (int, optional): Defaults to 40.
            precision (int, optional):
            keeptmp (bool, optional): Defaults to False.
        """
        self.have_construct_lattice()

        self.L = hlll_wrapper(self.L, threads=threads, threshold=threshold, **kwds)

    def run_bkz(self, block_size=20, verbose=False):
        """
        run BKZ over lattice

        Args:
            block_size (int, optional): Defaults to 20.
            verbose (bool, optional): Defaults to False.
        """
        self.have_construct_lattice()

        bkz_flags = BKZ.DEFAULT | BKZ.AUTO_ABORT
        if verbose:
            bkz_flags |= BKZ.VERBOSE
        par = BKZ.EasyParam(block_size=block_size, flags=bkz_flags)

        BKZ.reduction(self.L, par)

    def check_init(self, init):
        """
        check the initial state

        Args:
            init (list-like): the state to be checked

        Returns:
            bool:
        """
        n = self.mrg.degree
        assert len(init) == n

        if self.mrg.initial_state is not None:
            return self.mrg.initial_state == tuple(init)

        zbits = self.zbits
        y_ = self.outputs

        for a_i, y_i in zip(init, y_):
            if a_i >> zbits != y_i:
                return False

        state = list(init)
        for y_j in y_[n:]:
            a_j = (
                sum(a_i*c_i for a_i, c_i in zip(state[-n:], self.mrg.coeffs))
                % self.mrg.modulus
            )
            if a_j >> zbits != y_j:
                return False
            state.append(a_j)

        return True

    def recover_init(self, sol):
        """
        recover initial state from the solution vector

        Args:
            sol (list-like): solution vector

        Returns:
            tuple: the initial state, None for wrong vector
        """
        zbits = self.zbits
        size = 1 << (zbits - 1)
        scale = self.mrg.modulus >> (zbits - 1)

        if abs(sol[0]) != size * scale:
            return None

        if sol[0] == size * scale:
            sol = [-sol_i for sol_i in sol]

        z_ = tuple(int(sol_i/scale) + size for sol_i in sol[1:self.mrg.degree+1])
        if max(z_) >= 2**zbits:
            return None

        a_ = tuple(2**zbits*yi + zi for yi, zi in zip(self.outputs, z_))

        if self.check_init(a_):
            return a_
        return None

    @classmethod
    def volf(cls, m, n, zbits, d):
        """
        lattice volume

        Args:
            m (int): modulus
            n (int): degree
            zbits (int): truncated bits
            d (int): number of samples

        Returns:
            int:
        """
        size = 1 << (zbits - 1)
        scale = m >> (zbits - 1)

        vol = size * scale
        vol *= scale ** n
        vol *= (m * scale) ** (d - n)

        return vol

    @classmethod
    def ghf(cls, vol, dim, prec=53):
        """
        Estimate norm of shortest vector according to Gaussian Heuristic

        Args:
            vol (int): lattice volume
            dim (int): dimension
            prec (int, optional): precision to use. Defaults to 53.

        Returns:
            mpf:
        """
        # NOTE: The Gaussian Heuristic does not hold in small dimensions
        mpmath.mp.prec = prec

        dim = mpmath.mpf(dim)
        vol = mpmath.mpf(vol)

        gh = (mpmath.gamma(1 + dim / 2) * vol) ** (1 / dim) / mpmath.sqrt(mpmath.pi())

        return gh

    @classmethod
    def evf(cls, m, zbits, d, prec=53):
        """
        Estimate norm of target vector

        Args:
            m (int): modulus
            zbits (int): truncated bits
            d (int): number of samples
            prec (int, optional): precision to use. Defaults to 53.

        Returns:
            mpf:
        """
        mpmath.mp.prec = prec

        size = 1 << (zbits - 1)
        scale = m >> (zbits - 1)

        return mpmath.sqrt((size * scale)**2 + d * m ** 2 / 12)

    @classmethod
    def mvf(cls, m, zbits, d, prec=53):
        """
        Maximal norm of target vector

        Args:
            m (int): modulus
            zbits (int): truncated bits
            d (int): number of samples
            prec (int, optional): precision to use. Defaults to 53.

        Returns:
            mpf:
        """
        mpmath.mp.prec = prec

        size = 1 << (zbits - 1)
        scale = m >> (zbits - 1)

        return mpmath.sqrt((size * scale)**2 + d * m ** 2 / 12)
