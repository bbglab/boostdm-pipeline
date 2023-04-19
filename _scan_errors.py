import glob
import pandas as pd
import sys

input = sys.stdin.read().rstrip().split('\n')
for item in input:
    for fn in glob.glob(f'./work/{item}*/.command.log'):
        with open(fn, 'rt') as f:
            lines = list(map(str.rstrip, f.readlines()))
            if len(lines) > 0:
                print(fn, '>>>', lines[-1])
