# pylint: disable=C0103
"""
a Lattice object with reduction algorithms
"""
import random

from fpylll import IntegerMatrix, LLL, BKZ


class Lattice(IntegerMatrix):
    """
    IntegerMatrix with some reduction algorithms
    """

    def randomize_block(self, min_row=0, max_row=None, density=0):
        """
        Randomize basis between from ``min_row`` and ``max_row`` (exclusive)

        Args:
            min_row (int, optional): start in this row. Defaults to 0.
            max_row (int, optional): stop at this row (exclusive). Defaults to None.
            density (int, optional): number of non-zero coefficients in lower triangular
                transformation matrix. Defaults to 0.
        """
        if max_row is None:
            max_row = self.nrows

        # 1. permute rows
        niter = 4 * (max_row-min_row)  # some guestimate
        for _ in range(niter):
            b = a = random.randint(min_row, max_row-1)
            while b == a:
                b = random.randint(min_row, max_row-1)
            self.swap_rows(a, b)

        # 2. triangular transformation matrix with coefficients in -1,0,1
        for a in range(min_row, max_row-2):
            for _ in range(density):
                b = random.randint(a+1, max_row-1)
                s = random.randint(0, 1)
                self[a].addmul(self[b], x=2*s-1)

    def run_lll(self, delta=0.999, eta=0.501, verbose=False):
        """
        Args:
            delta (float, optional): Defaults to 0.999.
            eta (float, optional): Defaults to 0.501.
            verbose (bool, optional): Defaults to False.
        """
        flags = LLL.DEFAULT
        if verbose:
            flags |= LLL.VERBOSE

        LLL.reduction(self, delta=delta, eta=eta, flags=flags)

    def run_bkz(self, block_size=20, verbose=False, float_type=None, prec=0):
        """
        Args:
            block_size (int, optional): Defaults to 20.
            verbose (bool, optional): Defaults to False.
            float_type (str, optional): Defaults to None.
            prec (int, optional): Defaults to 0.
        """
        bkz_flags = BKZ.DEFAULT | BKZ.AUTO_ABORT
        if verbose:
            bkz_flags |= BKZ.VERBOSE
        par = BKZ.EasyParam(block_size=block_size, flags=bkz_flags)

        BKZ.reduction(self, par, float_type=float_type, precision=prec)

    def solve_asvp(self, **lwds):
        # TODO
        raise NotImplementedError
