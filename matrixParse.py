#!/usr/bin/env python3
# Parse tab-delimited penalty matrices and output to JSON
# author:   Kevin Chau

import json


def parse_matrix(pen_mat):
    """Parse penalty matrix from file

    Keyword arguments:
        pen_mat --  Penalty matrix file

    Returns:
        Dictionary of amino acid single-letter code pairs to penalty values
    """
    penalties_raw = []
    penalties = {}
    matrix_lines = [_.strip().split() for _ in pen_mat.readlines()]
    letters = matrix_lines[0]
    for pen_set in matrix_lines[1:]:
        penalties_raw.append(pen_set[1:])
    for v, let_v in enumerate(letters):
        for h, let_h in enumerate(letters):
            penalties["{}.{}".format(let_v, let_h)] = int(penalties_raw[v][h])
    return penalties

json.dump(parse_matrix(open('blosum62.txt', 'r')), open('blosum62.json', 'w'))
json.dump(parse_matrix(open('PAM250.txt', 'r')) , open('pam250.json', 'w'))
json.dump(parse_matrix(open('BASIC.txt', 'r')), open('basic.json', 'w'))
