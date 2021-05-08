# -*- coding: utf-8 -*-
"""
hlll wrapper
"""
import os
import subprocess

from fpylll import IntegerMatrix

from util import str_mat


def hlll_wrapper(A, threads=4, threshold=40, **kwds):
    r"""
    Run hlll for a given matrix

    :param A: an `IntegerMatrix` or `list`
    :param threads: ... (default: 4)
    :param threshold: block size threshold for parallelism (default: 40)
    :param delta: LLL parameter delta, 0 < \delta < 1 (default: 0.999)
    :param verbose: print verbose info (default: False)
    """
    if not os.path.exists("hlll"):
        raise RuntimeError("cannot find `hlll`")

    if isinstance(A, list):
        return_list = True
        m, n = len(A), len(A[0])
    elif isinstance(A, IntegerMatrix):
        return_list = False
        m, n = A.nrows, A.ncols
    else:
        raise TypeError("Matrix `A` type '%s' unknown." % type(A))
    tmpfile = f"hlll_tmp_matrix-{m}-{n}.mat"

    command = f"./hlll {tmpfile} "
    command += f"--threads {threads} "
    command += f"--threshold {threshold} "
    command += f"--delta {kwds.pop('delta', 0.999)} "
    if kwds.pop('verbose', False):
        command += "--verbose "

    with open(tmpfile, "w") as f:
        f.write(f"[{str_mat(A)}]")

    try:
        subprocess.check_call(command, shell=True)
        A = IntegerMatrix.from_file(tmpfile)
    finally:
        os.remove(tmpfile)

    if return_list:
        A = [[A[i, j] for j in range(A.ncols)] for i in range(A.nrows)]
    return A
