# -*- coding: utf-8 -*-
import re
import json

from fpylll import IntegerMatrix


def read_data(category, level):
    assert 1 <= category <= 3
    assert 1 <= level <= 9

    with open('data.txt', mode='r', encoding='GB2312') as f:
        DATA = f.read()

    pattern = r'挑战{}\.{}数据[^\n]+\n([^\n]+)\n'.format(category, level)
    match_text = re.search(pattern, DATA).groups()[0]
    data = list(map(int, match_text.split()))
    return data


def save_solution(category, level, data):
    with open(f"solutions/sol-{category}-{level}.json", 'w') as f:
        f.write(json.dumps(data))

def load_solution(category, level):
    with open(f"solutions/sol-{category}-{level}.json", 'r') as f:
        data = json.loads(f.read())
    return data


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
