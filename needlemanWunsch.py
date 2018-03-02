#!/usr/bin/env python3
# Implementation of the Needleman-Wunsch Global Alignment Algorithm
# author:   Kevin Chau

import sys

def global_alignment(string1, string2, mu_mat, sigma):
    """ Global alignment of two amino acid strings with penalties

    Keyword arguments:
    string1 --  First string
    string2 --  Second string
    mu_mat  --  Mismatch penalty matrix
    sigma   --  Indel pentalty
    """
    sigma = int(sigma)
    alignment_score = 0
    top_string = ''
    bot_string = ''

    # backtrack 2d array
    backtrack = [[None for _ in range(len(string2) + 1)]
                 for _ in range(len(string1) + 1)
                 ]

    # initialize value of paths to node (i, j)
    s = [[0 for _ in range(len(string2) + 1)] for _ in
         range(len(string1) + 1)]

    for _ in range(1, len(string2) + 1):
        backtrack[0][_] = 'RIGHT'
        s[0][_] = s[0][_ - 1] - 5
    for _ in range(1, len(string1) + 1):
        backtrack[_][0] = 'DOWN'
        s[_][0] = s[_ - 1][0] - 5
    s[0][0] = 0

    for j in range(1, len(string2) + 1):
        for i in range(1, len(string1) + 1):

            # which node led to current node?
            match = s[i - 1][j - 1] + mu_mat[
                (string1[i - 1], string2[j - 1])]

            s[i][j] = max([s[i - 1][j] - sigma, s[i][j - 1] - sigma, match])

            # backtrack assignments based on possible pre-nodes
            if s[i][j] == s[i - 1][j] - sigma:
                backtrack[i][j] = "DOWN"
            elif s[i][j] == s[i][j - 1] - sigma:
                backtrack[i][j] = "RIGHT"
            elif s[i][j] == match:
                backtrack[i][j] = "DOWN-RIGHT"

    i = len(string1)
    j = len(string2)

    while i > 0 or j > 0:
        if backtrack[i][j] == "DOWN":
            top_string += string1[i - 1]
            bot_string += '-'
            i -= 1
            alignment_score -= sigma
        elif backtrack[i][j] == "RIGHT":
            top_string += '-'
            bot_string += string2[j - 1]
            j -= 1
            alignment_score -= sigma
        elif backtrack[i][j] == "DOWN-RIGHT":
            top_string += string1[i - 1]
            bot_string += string2[j - 1]
            alignment_score += mu_mat[(string1[i - 1], string2[j - 1])]
            i -= 1
            j -= 1

    return alignment_score, top_string[::-1], bot_string[::-1]


if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit("Insufficient arguments!\nRequire <stringsfile> ('blosum62' OR 'pam250') <sigma>")

    filename, penalty_matrix, sigma_val = sys.argv[1:]

    with open(filename, 'r') as input_file:
        lines = input_file.readlines()

    score, s1, s2 = global_alignment(lines[0].strip(), 
                                     lines[1].strip(), 
                                     penalty_matrix, 
                                     int(sig_val))

    with open('result.txt', 'w') as result:
        result.write('\n'.join([str(score), str(s1), str(s2)]))
