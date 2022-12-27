#!/usr/bin/env python3

import hashlib
from itertools import permutations
from hashlib import sha256
from multiprocessing import Pool
import time

print(time.time())

def hash(s):
    m = hashlib.sha256()
    m.update(s.encode('ascii'))
    return m.digest()

def find_least_prefix(args):
    (letters, prefix) = args
    remaining_letters  = [l for l in letters if l not in prefix]
    least_hash = b"\x7f"
    least_message = ""
    for perm in permutations(remaining_letters):
        m = ''.join(prefix + perm)
        h = hash(m)
        if h < least_hash:
            least_hash = h
            least_message = m
    return (least_hash, least_message)
        
def find_least_linear(letters):
    r = []
    for prefix in permutations(letters, 2):
        r.append(find_least_prefix((letters, prefix)))
    return min(r)

def find_least_parallel(letters):
    prefixes = [(letters, prefix) for prefix in permutations(letters, 2)] 
    with Pool(processes=30) as pool:
        return min(pool.imap_unordered(find_least_prefix, prefixes))


fourteen = [str(digit) for digit in range(0, 10)] + ['a', 'b', 'c', 'd']
ten = [str(digit) for digit in range(0, 10)]

print("Checking on 14 permuted thingies!")

print(find_least_parallel(fourteen))
print(time.time())