#!/usr/bin/env python3
# Implementation of the Smith-Waterman Local Alignment Algorithm
# author:   Kevin Chau

import sys, json

penalty_matrices = {
        'blosum62': 'blosum62.json',
        'pam250': 'pam250.json',
        'basic': 'basic.json'
        }

def local_ali_dag(string1, string2, pen_mat, sigma):
    """Local alignment of two amino acid strings with penalties

    Keyword arguments:
    string1 --  First sting
    string2 --  Second string
    mu_mat  --  Mismatch penalty matrix
    sigma   --  Indel penalty
    """

    # case switching for penalty matrix
    mu_mat = json.load(open(penalty_matrices[pen_mat], 'r')) 

    sigma = int(sigma)
    node_properties = [[{} for _ in range(len(string2) + 1)]
                       for _ in range(len(string1) + 1)]
    top_string = ''
    bot_string = ''
    ali_score = 0
    total_max_path_weight = 0
    max_i = 0
    max_j = 0

    print("Initializing nodes...")

    # Initialize nodes; row-by-row traversal is topologically ordered
    for i in range(len(string1) + 1):
        for j in range(len(string2) + 1):
            node_properties[i][j] = {'PREVIOUS': None, 'P_WEIGHT': None}
            if j == 0 and i >= 1:
                node_properties[i][j]['PREVIOUS'] = (i - 1, j)
                node_properties[i][j]['P_WEIGHT'] = i * (-sigma)
            if i == 0 and j >= 1:
                node_properties[i][j]['PREVIOUS'] = (i, j - 1)
                node_properties[i][j]['P_WEIGHT'] = j * (-sigma)
            if i == 0 and j == 0:
                node_properties[i][j]['P_WEIGHT'] = 0

    print("Computing paths for all nodes...")

    # traverse nodes and update properties based on longest paths
    for i in range(1, len(string1) + 1):
        for j in range(1, len(string2) + 1):

            matchup = "{}.{}".format(string1[i - 1], string2[j - 1])
            ver = node_properties[i - 1][j]['P_WEIGHT'] - sigma
            hor = node_properties[i][j - 1]['P_WEIGHT'] - sigma
            dia = (node_properties[i - 1][j - 1]['P_WEIGHT']
                   + mu_mat[matchup])

            node_properties[i][j]['P_WEIGHT'] = max([0, ver, hor, dia])
            this_weight = node_properties[i][j]['P_WEIGHT']

            if node_properties[i][j]['P_WEIGHT'] >= total_max_path_weight:
                max_i = i
                max_j = j
                total_max_path_weight = this_weight

            if this_weight == 0:
                node_properties [i][j]['PREVIOUS'] = (0, 0)
            elif this_weight == ver:
                node_properties[i][j]['PREVIOUS'] = (i - 1, j)
            elif this_weight == hor:
                node_properties[i][j]['PREVIOUS'] = (i, j - 1)
            elif this_weight == dia:
                node_properties[i][j]['PREVIOUS'] = (i - 1, j - 1)

    if (node_properties[len(string1)][len(string2)]['P_WEIGHT']
            < total_max_path_weight):
        node_properties[len(string1)][len(string2)]['PREVIOUS'] = (max_i, max_j)

    print("Backtracking...")

    curr_node = (len(string1), len(string2))
    while node_properties[curr_node[0]][curr_node[1]]['P_WEIGHT'] != 0:
        pre = node_properties[curr_node[0]][curr_node[1]]['PREVIOUS']

        if pre[0] < curr_node[0] - 1 or pre[1] < curr_node[1] - 1:
            curr_node = (pre[0], pre[1])
            continue
        elif pre[0] == curr_node[0] - 1 and pre[1] == curr_node[1] - 1:
            top_string += string1[pre[0]]
            bot_string += string2[pre[1]]
        elif pre[0] == curr_node[0] - 1:
            top_string += string1[pre[0]]
            bot_string += '-'
        elif pre[1] == curr_node[1] - 1:
            top_string += '-'
            bot_string += string2[pre[1]]

        curr_node = (pre[0], pre[1])

    for p in range(len(top_string)):
        if top_string[p] == '-' or bot_string[p] == '-':
            ali_score -= sigma
        else:
            matchup = "{}.{}".format(top_string[p], bot_string[p])
            ali_score += mu_mat[matchup]

    return ali_score, top_string[::-1], bot_string[::-1]


if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit("Insufficient arguments!\nRequire <stringsfile> ('blosum62' OR 'pam250') <sigma>")

    filename, penalty_matrix, sigma_val = sys.argv[1:]

    with open(filename, 'r') as input_file:
        lines = input_file.readlines()

    score, s1, s2 = local_ali_dag(lines[0].strip(),
                                  lines[1].strip(), 
                                  penalty_matrix, 
                                  int(sigma_val))

    with open('result.txt', 'w') as result:
        result.write('\n'.join([str(score), str(s1), str(s2)]))
