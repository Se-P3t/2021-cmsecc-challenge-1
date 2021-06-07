# pylint: disable=C0103
"""
a Lattice object with reduction algorithms
"""
import sys
import random

from fpylll import IntegerMatrix, LLL, BKZ
from fpylll.util import gaussian_heuristic

from g6k.algorithms.workout import workout
from g6k.siever import Siever
from g6k import SieverParams
from g6k.utils.cli import pop_prefixed_params
from g6k.utils.stats import SieveTreeTracer
from g6k.utils.util import db_stats

import six
from six.moves import range


class Lattice(IntegerMatrix):
    """
    IntegerMatrix with some reduction algorithms
    """

    def randomize_block(self, min_row=0, max_row=None, density=None):
        """
        Randomize basis between from ``min_row`` and ``max_row`` (exclusive)

        Args:
            min_row (int, optional): start in this row. Defaults to 0.
            max_row (int, optional): stop at this row (exclusive). Defaults to None.
            density (int, optional): number of non-zero coefficients in lower triangular
                transformation matrix. Defaults to None.
        """
        if max_row is None:
            max_row = self.nrows
        if density is None:
            density = self.nrows // 4

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

    def __asvp_kernel(self, params):
        """
        TODO
        """
        seed = params.pop("seed", None)

        goal_r0__gh = params.pop('goal_r0__gh')
        pump_params = pop_prefixed_params("pump", params)
        workout_params = pop_prefixed_params("workout", params)
        verbose = params.pop("verbose")
        if verbose:
            workout_params["verbose"] = True

        g6k = Siever(self, params, seed=seed)
        tracer = SieveTreeTracer(g6k, root_label=("asvp", self.nrows), start_clocks=True)

        gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(self.nrows)])
        goal_r0 = workout_params.pop('goal_r0', (goal_r0__gh**2) * gh)
        r0 = sum([x * x for x in self[0]])
        if verbose:
            print(
                (
                    "gh = %f, goal_r0/gh = %f, r0/gh = %f"
                    % (gh, goal_r0 / gh, r0 / gh)
                )
            )
        if goal_r0__gh * gh > r0:
            tracer.exit()
            tracer.trace.data["flast"] = -1
            return tracer.trace

        flast = workout(
            g6k, tracer, 0, self.nrows, goal_r0=goal_r0, pump_params=pump_params, **workout_params
        )

        tracer.exit()
        stat = tracer.trace
        stat.data["flast"] = flast

        if verbose:
            norm = self[0].norm()
            print("svp: norm %.1f ,hf %.5f" % (norm**.5, (norm/gh)**.5))

        return tracer.trace

    def solve_asvp(self, threads=1, seed=None, **kwds):
        """
        A G6K Approx-SVP Solver

        Args:
            threads (int, optional): Defaults to 1.
            seed (int, optional): Defaults to None.
            verbose (bool, optional): Defaults to True.
            goal_r0__gh (int, optional): Defaults to 1.05.
        """
        dry_run = kwds.pop('dry_run', False)
        verbose = kwds.pop('verbose', True)
        goal_r0__gh = kwds.pop('goal_r0__gh', 1.05)

        params = SieverParams(
            threads=threads,
            seed=seed,
            verbose=verbose,
            goal_r0__gh=goal_r0__gh,
        )

        for k, v in six.iteritems(kwds):
            params[k] = v
        if dry_run:
            print(params)
            sys.exit(0)

        stats = self.__asvp_kernel(params)

        if verbose:
            data = stats.data
            avg_max, _ = db_stats(stats)
            print(
                f"{stats.label} :: "
                f"cputime {data['cputime'].sum:7.4f}s, "
                f"walltime: {data['walltime'].sum:7.4f}s, "
                f"flast: {data['flast']:3.2f}, "
                f"|db|: 2^{avg_max:.2f}"
            )
