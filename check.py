# -*- coding: utf-8 -*-
"""
check the solution
"""
import argparse
from copy import copy

from util import load_solution, read_data


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('category', type=int, help="challenge category")
parser.add_argument('level', type=int, help="challenge level")

args, _ = parser.parse_known_args()

try:
    data = load_solution(args.category, args.level)
except FileNotFoundError:
    print(f"sol-{args.category}-{args.level}.json not found, "
        "solve it first")
    exit(-1)
y_ = read_data(args.category, args.level)

globals().update(data)

state = copy(initial_state)
while initial_state:
    if initial_state.pop(0) >> zbits != y_.pop(0):
        exit(-1)
while y_:
    a_j = sum(c*a for c, a in zip(coefficients, state)) % modulus
    _ = state.pop(0)
    state.append(a_j)
    if a_j >> zbits != y_.pop(0):
        exit(-1)
