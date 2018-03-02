"""
    Kevin Chau
    BIMM 181
    Chapter 4 Coding Challenges
"""

import Chapter1
import Chapter2
import Chapter3
import sys
from copy import deepcopy

"""
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
"""

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


# Input: target value money and denominations coins
# Output: fewest number of coins to reach target value
def dyn_prog_change(money, coins):
    min_count_list = [0]    # target val 0 takes 0 coins
    for m in range(1, money + 1):
        min_count_list.append(sys.maxsize)
        for i in range(0, len(coins)):
            if m >= coins[i]:
                if min_count_list[m - coins[i]] + 1 < min_count_list[m]:
                    min_count_list[m] = min_count_list[m - coins[i]] + 1
    return str(min_count_list[money])


# Input: ints n and m, down grid and right grid
# Output: length of longest path from (0, 0) to (n, m)
# Note: n is vert axis, m is horiz axis
def manhattan_tourist(n, m, down, right):
    n = int(n)
    m = int(m)
    l = []
    for _ in range(n+1):
        l.append([0] * (m+1))

    for i in range(1, n + 1):
        l[i][0] = l[i-1][0] + int(down[i-1][0])

    for j in range(1, m + 1):
        l[0][j] = l[0][j-1] + int(right[0][j-1])

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            print(l[i][j], i, j)
            l[i][j] = max([l[i-1][j] + int(down[i-1][j]),
                           l[i][j-1] + int(right[i][j-1])])

    return str(l[n][m])


# Input: two strings string1 and string2
# Output: assigns backtracking pointers to each node
def lcs_backtrack(string1, string2):

    # backtrack 2d array
    backtrack = [[None for _ in range(len(string2) + 1)]
                 for _ in range(len(string1) + 1)
                 ]

    # initialize value of paths to node (i, j)
    s = [[0 for _ in range(len(string2) + 1)] for _ in range(len(string1) + 1)]

    for j in range(1, len(string2) + 1):
        for i in range(1, len(string1) + 1):

            # which node led to current node?
            if string1[i-1] == string2[j-1]:
                match = s[i-1][j-1] + 1
            else:
                match = s[i-1][j-1]
            s[i][j] = max([s[i-1][j], s[i][j-1], match])

            # backtrack assignments based on possible pre-nodes
            if s[i][j] == s[i-1][j]:
                backtrack[i][j] = "DOWN"
            elif s[i][j] == s[i][j-1]:
                backtrack[i][j] = "RIGHT"
            elif s[i][j] == (s[i-1][j-1] + 1) and string1[i-1] == string2[j-1]:
                backtrack[i][j] = "DOWN-RIGHT"

    return backtrack


# Input: array of backtracking pointers, a string, ints i, j, and empty string
#        i, j == len string1 and len string2
# Output: updates empty string to longest common subsequence
# ISSUE: too many recursive calls
def ret_lcs(backtrack_list, string, i, j, file_obj):
    i = int(i)
    j = int(j)
    if i == 0 or j == 0:
        return
    if backtrack_list[i][j] == "DOWN":
        ret_lcs(backtrack_list, string, i-1, j, file_obj)
    elif backtrack_list[i][j] == "RIGHT":
        ret_lcs(backtrack_list, string, i, j-1, file_obj)
    elif backtrack_list[i][j] == "DOWN-RIGHT":
        ret_lcs(backtrack_list, string, i-1, j-1, file_obj)
        file_obj.write(string[i-1])
        return string[i-1]

    return


# Iterative implementation of ret_lcs
def iter_ret_lcs(backtrack_list, string, i , j):
    i = int(i)
    j = int(j)
    output_string = []
    while i > 0 and j > 0:
        if backtrack_list[i][j] == "DOWN":
            i -= 1
        elif backtrack_list[i][j] == "RIGHT":
            j -= 1
        elif backtrack_list[i][j] == "DOWN-RIGHT":
            output_string.append(string[i-1])
            i -= 1
            j -= 1
    return ''.join(reversed(output_string))


# Input: source node, sink node, set of nodes, set of unweighted edges
# Output: list of nodes in topological order
def top_order(incoming_edges_orig, source=None, sink= None):
    incoming_edges = deepcopy(incoming_edges_orig)
    topological = []
    
    if source:
        n = source
        
    while len(incoming_edges) > 1:
        if not source:
            for e in incoming_edges:
                if not incoming_edges[e]:
                    n = e

        topological.append(n)
        del incoming_edges[n]

        for e in incoming_edges:
            incoming_edges[e].pop(n, None)

        for e in incoming_edges:
            if not incoming_edges[e] and e != sink:
                n = e
                
    topological.append(n) if not sink else topological.append(sink)
    return topological


# Input: adjacency list of a dag with weighted edges
#        adj_list format: <source> -> <target>:<weight>
# Output: heaviest path in the dag
def dag_longest(adj_list, source, sink):
    dag_path = []
    from_source = []
    from_target = {}
    node_props = {}

    for edge in adj_list:
        edge = edge.strip().split('->')
        target, weight = edge[1].strip().split(':')

        if target not in from_target:
            from_target[target] = {}

        if edge[0] not in from_target:
            from_target[edge[0]] = {}

        from_target[target][edge[0]] = weight

    top_dag = top_order(from_target, source, sink)

    for n in top_dag:
        node_props[n] = {'BACK_NODE': None, 'CURR_LEN': 0}

    from_source.append(source)

    for n in top_dag[1:]:

        for pre_node in from_target[n]:

            # only choose nodes in path from source
            if pre_node not in from_source:
                continue

            from_pre_to_tar = (int(node_props[pre_node]['CURR_LEN'])
                               + int(from_target[n][pre_node]))

            if from_pre_to_tar > node_props[n]['CURR_LEN']:
                node_props[n]['CURR_LEN'] = from_pre_to_tar
                node_props[n]['BACK_NODE'] = pre_node
                from_source.append(n)

    curr_node = sink
    while node_props[curr_node]['BACK_NODE']:
        dag_path.append(curr_node)
        curr_node = node_props[curr_node]['BACK_NODE']
    dag_path.append(curr_node)

    print(top_dag)

    return node_props[sink]['CURR_LEN'], reversed(dag_path)


# Input: two strings
# Output: value of max alignment along with alignmnet
# Essentially a copy of ret_lcs + lcs_backtrack but with penalties
def global_alignment(string1, string2, mu_mat, sigma):
    sigma = int(sigma)
    alignment_score = 0
    top_string = ''
    bot_string = ''

    # backtrack 2d array
    backtrack = [[None for _ in range(len(string2) + 1)]
                 for _ in range(len(string1) + 1)
                 ]

    # initialize value of paths to node (i, j)
    s = [[0 for _ in range(len(string2) + 1)] for _ in range(len(string1) + 1)]

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
            match = s[i-1][j-1] + mu_mat[(string1[i-1], string2[j-1])]

            s[i][j] = max([s[i-1][j] - sigma, s[i][j-1] - sigma, match])

            # backtrack assignments based on possible pre-nodes
            if s[i][j] == s[i-1][j] - sigma:
                backtrack[i][j] = "DOWN"
            elif s[i][j] == s[i][j-1] - sigma:
                backtrack[i][j] = "RIGHT"
            elif s[i][j] == match:
                backtrack[i][j] = "DOWN-RIGHT"

    # print('\n'.join([str(line) for line in s]))

    i = len(string1)
    j = len(string2)

    while i > 0 or j > 0:
        if backtrack[i][j] == "DOWN":
            top_string += string1[i-1]
            bot_string += '-'
            i -= 1
            alignment_score -= sigma
        elif backtrack[i][j] == "RIGHT":
            top_string += '-'
            bot_string += string2[j-1]
            j -= 1
            alignment_score -= sigma
        elif backtrack[i][j] == "DOWN-RIGHT":
            top_string += string1[i-1]
            bot_string += string2[j-1]
            alignment_score += mu_mat[(string1[i-1], string2[j-1])]
            i -= 1
            j -= 1

    return alignment_score, top_string[::-1], bot_string[::-1]


# Input: two strings
# Output: value of max local alignment along with alignments
# very similar to global_alignment, with caveat that source targets all nodes
# and sink is targetted by all nodes
# BASED ON MANHATTAN-STRUCTURE
def local_alignment(string1, string2, mu_mat, sigma):
    sigma = int(sigma)
    alignment_score = 0
    top_string = ''
    bot_string = ''

    # backtrack 2d array
    backtrack = [[None for _ in range(len(string2) + 1)]
                 for _ in range(len(string1) + 1)
                 ]

    # initialize value of paths to node (i, j)
    s = [[0 for _ in range(len(string2) + 1)] for _ in range(len(string1) + 1)]

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
            match = s[i-1][j-1] + mu_mat[(string1[i-1], string2[j-1])]

            s[i][j] = max([0, s[i-1][j] - sigma, s[i][j-1] - sigma, match])

            # backtrack assignments based on possible pre-nodes
            if s[i][j] == 0:
                backtrack[i][j] = "SOURCE"
            elif s[i][j] == s[i-1][j] - sigma:
                backtrack[i][j] = "DOWN"
            elif s[i][j] == s[i][j-1] - sigma:
                backtrack[i][j] = "RIGHT"
            elif s[i][j] == match:
                backtrack[i][j] = "DOWN-RIGHT"

    max_v_index = 0
    max_h_index = 0
    max_s = 0
    for j in range(1, len(string2)):
        for i in range(1, len(string1)):
            if s[i][j] > max_s:
                max_v_index = i
                max_h_index = j
                max_s = s[i][j]

    # print('\n'.join([str(line) for line in s]))

    if backtrack[len(string1)][len(string2)] == "DOWN-RIGHT":
        i = len(string1)
        j = len(string2)
    else:
        i = max_v_index
        j = max_h_index
    s[i][j] = max_s

    while i > 0 or j > 0:
        if s[i][j] == 0:
            break
        elif backtrack[i][j] == "DOWN":
            top_string += string1[i-1]
            bot_string += '-'
            i -= 1
        elif backtrack[i][j] == "RIGHT":
            top_string += '-'
            bot_string += string2[j-1]
            j -= 1
        elif backtrack[i][j] == "DOWN-RIGHT":
            top_string += string1[i-1]
            bot_string += string2[j-1]
            i -= 1
            j -= 1

    for p in range(len(top_string)):
        if top_string[p] == '-' or bot_string[p] == '-':
            alignment_score -= sigma
        else:
            alignment_score += mu_mat[(top_string[p], bot_string[p])]

    return alignment_score, top_string[::-1], bot_string[::-1]


# DAG principles to solve local_alignment problem
def local_ali_dag(string1, string2, mu_mat, sigma):
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

            # print(i, j)

            ver = node_properties[i - 1][j]['P_WEIGHT'] - sigma
            hor = node_properties[i][j - 1]['P_WEIGHT'] - sigma
            dia = (node_properties[i - 1][j - 1]['P_WEIGHT']
                   + mu_mat[string1[i - 1], string2[j - 1]])

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
            ali_score += mu_mat[top_string[p], bot_string[p]]

    return ali_score, top_string[::-1], bot_string[::-1]


if __name__ == '__main__':
    filename = sys.argv[1]

    with open(filename, 'r') as input_file:
        lines = input_file.readlines()

    score, s1, s2 = local_ali_dag(lines[0].strip(), lines[1].strip(), pam250, 5)

    with open('result.txt', 'w') as result:
        result.write('\n'.join([str(score), str(s1), str(s2)]))