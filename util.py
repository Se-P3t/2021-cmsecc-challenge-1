# -*- coding: utf-8 -*-
import sys
import re
import json
import subprocess
from random import SystemRandom
from collections import OrderedDict

from fpylll import IntegerMatrix


RNG = SystemRandom()


def asvp_warpper(A, *cmd_args, verbose=False):
    """
    wraper of script `svp_challenge.py`

    EXAMPLE::

        >>> from fpylll import IntegerMatrix
        >>> A = IntegerMatrix.random(80, 'uniform', bits=30)
        >>> sol = asvp_warpper(A, '--seed 123', '--threads 1', '--dry-run')
        >>> sol
        ...
    """
    try:
        dim = A.nrows
    except:
        dim = len(A)
    open('/tmp/asvp-{}.mat'.format(dim), 'w').write(str_mat(A))
    cmd = 'svp_challenge.py 40 --load-matrix /tmp/asvp-{}.mat '.format(dim) + \
        ' '.join(cmd_args)
    if verbose:
        print(cmd)

    with subprocess.Popen(cmd, shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as proc, \
            open('/tmp/asvp-{}.log'.format(dim), 'bw') as logfile:
        while True:
            byte = proc.stdout.read(1)
            if byte:
                if verbose:
                    sys.stdout.buffer.write(byte)
                    sys.stdout.flush()
                logfile.write(byte)
                logfile.flush()
            else:
                break

    exit_status = proc.returncode
    if exit_status != 0:
        raise RuntimeError("exit({})".format(exit_status))

    with open('/tmp/asvp-{}.log'.format(dim), 'r') as logfile:
        logs = logfile.read()
        pattern = r'sol\s\d+,\s\(([-,\s\d]+)\)'
        sol = re.search(pattern, logs).groups()
        if not sol:
            return None
        sol = list(map(int, sol.split(', ')))
        return sol




def read_data(category, level):
    assert 1 <= category <= 3
    assert 1 <= level <= 9

    with open('环上序列的截位还原问题——数据.txt', mode='r', encoding='GB2312') as f:
        DATA = f.read()

    pattern = r'挑战{}\.{}数据[^\n]+\n([^\n]+)\n'.format(category, level)
    match_text = re.search(pattern, DATA).groups()[0]
    data = list(map(int, match_text.split()))
    return data


def save_solution(category, level, data):
    with open(f"solutions/sol-{category}-{level}.json", 'w') as f:
        f.write(json.dumps(data))


def matrix_overview(BB):
    try:
        for ii in range(BB.nrows):
            a = ('%03d ' % ii)
            for jj in range(BB.ncols):
                if BB[ii, jj] == 0:
                    a += ' '
                elif abs(BB[ii, jj]) == 1:
                    a += '1'
                #elif abs(BB[ii, jj]) < 1:
                #    a += '.'
                else:
                    a += 'X'
                # if len(BB) < 60:
                #     a += ' '
            print(a, flush=True)
        print('', flush=True)
    except:
        for ii in range(len(BB)):
            a = ('%03d ' % ii)
            for jj in range(len(BB[0])):
                if BB[ii][jj] == 0:
                    a += ' '
                elif abs(BB[ii][jj]) == 1:
                    a += '1'
                #elif abs(BB[ii][jj]) < 1:
                #    a += '.'
                else:
                    a += 'X'
                # if len(BB) < 60:
                #     a += ' '
            print(a, flush=True)
        print('', flush=True)


def str_mat(m):
    if isinstance(m, IntegerMatrix):
        return str(m)
    elif isinstance(m, list):
        s = ''
        for v in m:
            s += str(v).replace(',', '')
            s += '\n'
        return s
    raise TypeError(f"unknown type ({type(m)}) of input")


def print_stats(fmt, stats, keys, extractf=None):
    if extractf is None:
        extractf = {}
    for (n, params), stat in stats.items():
        kv = OrderedDict()
        for key in keys:
            if key in extractf:
                value = extractf[key](n, params, stat)
            else:
                value = sum([float(node[key]) for node in stat]) / len(stat)
            kv[key] = value

        print(fmt.format(name=params, n=n, **kv))
