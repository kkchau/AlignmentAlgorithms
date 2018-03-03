#!/usr/bin/env python3
# Implementation of the Needleman-Wunsch Global Alignment Algorithm
# author:   Kevin Chau

import sys, json

penalty_matrices = {
        'blosum62': 'blosum62.json',
        'pam250': 'pam250.json',
        'basic': 'basic.json'
        }

def global_alignment(string1, string2, pen_mat, sigma):
    """ Global alignment of two amino acid strings with penalties

    Keyword arguments:
    string1 --  First string
    string2 --  Second string
    mu_mat  --  Mismatch penalty matrix
    sigma   --  Indel penalty
    """

    sigma = int(sigma)
    alignment_score = 0
    top_string = ''
    bot_string = ''

    # case switching for penalty matrix
    mu_mat = json.load(open(penalty_matrices[pen_mat], 'r')) 

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
            matchup = "{}.{}".format(string1[i - 1], string2[j - 1])
            match = s[i - 1][j - 1] + mu_mat[matchup]

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
            matchup = "{}.{}".format(string1[i - 1], string2[j - 1])
            alignment_score += mu_mat[matchup]
            i -= 1
            j -= 1

    return alignment_score, top_string[::-1], bot_string[::-1]


def affine_alignment(string1, string2, mu_mat, sigma, epsilon):
    """ Global alignment of two amino acid strings with affine gap penalties

    Keyword arguments:
    string1 --  First string
    string2 --  Second string
    mu_mat  --  Mismatch penalty matrix
    sigma   --  Indel opening penalty
    epsilon --  Gap extension penalty
    """

    sigma = int(sigma)
    epsilon = int(epsilon)
    top_string = ''
    bot_string = ''

    # backtrack matrix
    backtrack = [[None for _ in range(len(string2) + 1)]
                 for _ in range(len(string1) + 1)]
    for i in range(1, len(string1) + 1):
        backtrack[i][0] = 0
    for j in range(1, len(string2) + 1):
        backtrack[0][j] = 2

    # dags for lower, middle, and upper scores
    l = [[0 for _ in range(len(string2) + 1)] for _ in range(len(string1) + 1)]
    m = [[0 for _ in range(len(string2) + 1)] for _ in range(len(string1) + 1)]
    u = [[0 for _ in range(len(string2) + 1)] for _ in range(len(string1) + 1)]

    for i in range(1, len(string1) + 1):
        l[i][0] = -sigma + i*(-epsilon)
        m[i][0] = -float('inf')
        u[i][0] = -float('inf')

    for j in range(1, len(string2) + 1):
        l[0][j] = -float('inf')
        m[0][j] = -float('inf')
        u[0][j] = -sigma + j*(-epsilon)

    for i in range(1, len(string1) + 1):
        for j in range(1, len(string2) + 1):

            matchup = "{}.{}".format(string1[i - 1], string2[j - 1])

            l[i][j] = max([l[i-1][j] - epsilon, m[i-1][j] - sigma])

            u[i][j] = max([u[i][j-1] - epsilon, m[i][j-1] - sigma])
            m[i][j] = max(mu_mat[matchup] + m[i-1][j-1], l[i][j], [i][j])

            pre = max(l[i][j], m[i][j], u[i][j])

            if pre == u[i][j]:
                backtrack[i][j] = 2
            elif pre == m[i][j]:
                backtrack[i][j] = 1
            elif pre == l[i][j]:
                backtrack[i][j] = 0

    print('\n'.join(str(x) for x in backtrack) + '\n')

    curr_node = (len(string1), len(string2))
    while curr_node[0] > 0 or curr_node[1] > 0:
        if backtrack[curr_node[0]][curr_node[1]] == 0:
            top_string += string1[curr_node[0] - 1]
            bot_string += '-'
            curr_node = (curr_node[0] - 1, curr_node[1])
        elif backtrack[curr_node[0]][curr_node[1]] == 1:
            top_string += string1[curr_node[0] - 1]
            bot_string += string2[curr_node[1] - 1]
            curr_node = (curr_node[0] - 1, curr_node[1] - 1)
        elif backtrack[curr_node[0]][curr_node[1]] == 2:
            top_string += '-'
            bot_string += string2[curr_node[1] - 1]
            curr_node = (curr_node[0], curr_node[1] - 1)
    print(u[3][4])

    return (str(max(l[len(string1)][len(string2)],
                    m[len(string1)][len(string2)],
                    u[len(string1)][len(string2)])),
            top_string[::-1],
            bot_string[::-1])


if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit("Insufficient arguments!\nRequire <stringsfile> ('blosum62' OR 'pam250') <sigma>")

    filename, penalty_matrix, sigma_val = sys.argv[1:]

    with open(filename, 'r') as input_file:
        lines = input_file.readlines()

    score, s1, s2 = global_alignment(lines[0].strip(), 
                                     lines[1].strip(), 
                                     penalty_matrix, 
                                     int(sigma_val))

    with open('result.txt', 'w') as result:
        result.write('\n'.join([str(score), str(s1), str(s2)]))
