#!/usr/bin/env python3
# Parse tab-delimited penalty matrices
# author:   Kevin Chau


blosum62_penalties_raw = []
blosum62 = {}
with open('blosum62.txt', 'r') as blos:
    b_lines = blos.readlines()
    letters = b_lines[0].strip().split()
    for pen_set in b_lines[1:]:
        blosum62_penalties_raw.append(pen_set.strip().split()[1:])
for v, let_v in enumerate(letters):
    for h, let_h in enumerate(letters):
        blosum62[(let_v, let_h)] = int(blosum62_penalties_raw[v][h])

pam250_penalties_raw = []
pam250 = {}
with open('pam250.txt', 'r') as blos:
    b_lines = blos.readlines()
    letters = b_lines[0].strip().split()
    for pen_set in b_lines[1:]:
        pam250_penalties_raw.append(pen_set.strip().split()[1:])
for v, let_v in enumerate(letters):
    for h, let_h in enumerate(letters):
        pam250[(let_v, let_h)] = int(pam250_penalties_raw[v][h])

basic_penalties_raw = []
basic = {}
with open('BASIC.txt', 'r') as blos:
    b_lines = blos.readlines()
    letters = b_lines[0].strip().split()
    for pen_set in b_lines[1:]:
        basic_penalties_raw.append(pen_set.strip().split()[1:])
for v, let_v in enumerate(letters):
    for h, let_h in enumerate(letters):
        basic[(let_v, let_h)] = int(basic_penalties_raw[v][h])

