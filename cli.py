from __future__ import absolute_import
from __future__ import print_function
import argparse
import copy
import datetime
import logging
import multiprocessing_logging
import os
import re
import socket
import subprocess
import sys
from collections import OrderedDict
from multiprocessing import Pool



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
    with open(f"solutions/sol-{category}-{level}.txt", 'w') as f:
        for d in data:
            f.write(f"{d}\n")


def log_filenamef(category, level):
    date = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")
    hostname = socket.gethostname()
    log_filename = "chall-{category}-{level}-{date}-{hostname}.log".format(category=category,
                                                                           level=level,
                                                                           date=date,
                                                                           hostname=hostname)
    log_filename = os.path.join("logs", log_filename)
    return log_filename


