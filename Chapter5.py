"""
    Kevin Chau
    BIMM 181
    Chapter 5 Coding Challenges
"""

# import Chapter1
# import Chapter2
# import Chapter3
import sys
import math
from copy import deepcopy


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

    # print('\n'.join([str(line) for line in s]))

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


# Input: two strings
# Output: value of max alignment along with alignmnet
# Essentially a copy of ret_lcs + lcs_backtrack but with penalties
def multi_alignment(s1, s2, s3):

    def basic_score(a, b, c):
        return 1 if a == b and b == c else 0

    alignment_score = 0
    f_s = ''
    s_s = ''
    t_s = ''

    # backtrack 2d array
    backtrack = [[[None for _ in range(len(s3) + 1)]
                  for _ in range(len(s2) + 1)]
                 for _ in range(len(s1) + 1)]

    # initialize value of paths to node (i, j)
    s = [[[0 for _ in range(len(s3) + 1)]
          for _ in range(len(s2) + 1)]
         for _ in range(len(s1) + 1)]

    for i in range(len(s1) + 1):
        for j in range(len(s2) + 1):
            for k in range(len(s3) + 1):

                if 0 in [i, j, k]:
                    if i > 0 and j == 0 and k == 0:
                        s[i][j][k] = s[i-1][j][k]
                        backtrack[i][j][k] = 0
                    elif i == 0 and j > 0 and k == 0:
                        s[i][j][k] = s[i][j-1][k]
                        backtrack[i][j][k] = 1
                    elif i == 0 and j == 0 and k > 0:
                        s[i][j][k] = s[i][j][k-1]
                        backtrack[i][j][k] = 2
                    elif i > 0 and j > 0 and k == 0:
                        scores = [s[i-1][j][k],
                                  s[i][j-1][k],
                                  -sys.maxsize,
                                  s[i-1][j-1][k]]
                        s[i][j][k] = max(scores)
                        backtrack[i][j][k] = scores.index(s[i][j][k])
                    elif i > 0 and j == 0 and k > 0:
                        scores = [s[i-1][j][k],
                                  -sys.maxsize,
                                  s[i][j][k-1],
                                  -sys.maxsize,
                                  s[i-1][j][k-1]]
                        s[i][j][k] = max(scores)
                        backtrack[i][j][k] = scores.index(s[i][j][k])
                    elif i == 0 and j > 0 and k > 0:
                        scores = [-sys.maxsize,
                                  s[i][j-1][k],
                                  s[i][j][k-1],
                                  -sys.maxsize,
                                  -sys.maxsize,
                                  s[i][j-1][k-1]]
                        s[i][j][k] = max(scores)
                        backtrack[i][j][k] = scores.index(s[i][j][k])
                    continue

                full_align = basic_score(s1[i-1], s2[j-1], s3[k-1])
                scores = [s[i-1][j][k],
                          s[i][j-1][k],
                          s[i][j][k-1],
                          s[i-1][j-1][k],
                          s[i-1][j][k-1],
                          s[i][j-1][k-1],
                          s[i-1][j-1][k-1] + full_align
                          ]

                s[i][j][k] = max(scores)
                backtrack[i][j][k] = scores.index(s[i][j][k])

    i, j, k = [len(s1), len(s2), len(s3)]
    b_switch = {None: (0, 0, 0),
                0: (-1, 0, 0),
                1: (0, -1, 0),
                2: (0, 0, -1),
                3: (-1, -1, 0),
                4: (-1, 0, -1),
                5: (0, -1, -1),
                6: (-1, -1, -1)
                }
    while i > 0 or j > 0 or k > 0:
        b = b_switch[backtrack[i][j][k]]

        f_s = f_s + s1[i-1] if b[0] else f_s + '-'
        s_s = s_s + s2[j-1] if b[1] else s_s + '-'
        t_s = t_s + s3[k-1] if b[2] else t_s + '-'

        i += b[0]
        j += b[1]
        k += b[2]

    ali_score = sum(1 if triplet[0] == triplet[1] and triplet[1] == triplet[2]
                    else 0 for triplet in list(zip(f_s, s_s, t_s)))

    return str(ali_score), f_s[::-1], s_s[::-1], t_s[::-1]


# Input: two strings, scoring matrix, sigma and epsilon
# Output: global alignment implementing affine gap penalties
def affine_alignment(string1, string2, mu_mat, sigma, epsilon):
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

            l[i][j] = max([l[i-1][j] - epsilon, m[i-1][j] - sigma])

            u[i][j] = max([u[i][j-1] - epsilon, m[i][j-1] - sigma])
            m[i][j] = max(mu_mat[(string1[i-1], string2[j-1])] + m[i-1][j-1],
                                                                 l[i][j],
                                                                 u[i][j])

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


def to_middle(t1, t2, m_score, sig):
    # initialize first two columns
    s = [[i * -sig, 0] for i in range(len(t1) + 1)]
    b = [None for _ in range(len(t1) + 1)]

    # build columns by pairs
    for j in range(1, (len(t2) // 2) + 1):
        for i in range(len(t1) + 1):

            if i == 0:
                s[i][1] = j*-sig
                continue

            scores = [s[i - 1][1] - sig,
                      s[i][0] - sig,
                      s[i - 1][0] + m_score[(t1[i - 1], t2[j - 1])]]
            s[i][1] = max(scores)
            b[i] = scores.index(s[i][1])

        for i in range(len(t1) + 1):
            s[i][0] = s[i][1]

    return [[r[0]] for r in s], b


def middle_edge(t1, t2, m_score, sig):
    mid_edge = [None, None]

    # source_to_middle and middle_to_sink calculations
    s_to_m = to_middle(t1, t2, m_score, sig)[0]
    m_to_s, back = to_middle(t1[::-1],
                             t2[::-1] + ['', ' '][len(t2) % 2 == 1
                                                  and len(t2) > 1],
                             m_score,
                             sig)
    m_to_s = m_to_s[::-1]
    back = back[::-1]

    # calculate i_paths
    i_scores = []
    for i in range(len(s_to_m)):
        i_scores.append((s_to_m[i] + m_to_s[i]))

    # index of max i_path
    max_i = max(range(len(s_to_m)), key=lambda n: i_scores[n])

    # construct middle edge
    mid_edge[0] = (max_i, len(t2) // 2)
    mid_edge[1] = [(max_i + 1, len(t2) // 2),
                   (max_i, (len(t2) // 2) + 1),
                   (max_i + 1, (len(t2) // 2) + 1)][back[max_i]]

    return mid_edge


def lin_space_ali(t1, t2, score, sig):

    def lin_space_rec(top, bot, lef, rig):
        if lef == rig:
            return [t1[top:bot], '-'*(bot-top)]
        elif top == bot:
            return['-'*(rig-lef), t2[lef:rig]]
        elif 1 in [bot - top, rig - lef]:
            return global_alignment(t1[top:bot], t2[lef:rig], score, sig)[1:]

        mid_edge = middle_edge(t1[top:bot], t2[lef:rig], score, sig)
        mid_node, next_node = mid_edge

        mid_node = (mid_node[0] + top, mid_node[1] + lef)
        next_node = (next_node[0] + top, next_node[1] + lef)

        m_top = ['-', t1[mid_node[0] % len(t1)]][next_node[0] - mid_node[0]]
        m_bot = ['-', t2[mid_node[1] % len(t2)]][next_node[1] - mid_node[1]]

        first_top, first_bot = lin_space_rec(top, mid_node[0], lef, mid_node[1])
        last_top, last_bot = lin_space_rec(next_node[0], bot, next_node[1], rig)

        t_string = first_top + m_top + last_top
        b_string = first_bot + m_bot + last_bot

        return t_string, b_string

    top_string, bot_string = lin_space_rec(0, len(t1), 0, len(t2))
    ali_score = sum(-sig if '-' in alignment
                    else score[(alignment[0], alignment[1])]
                    for alignment in list(zip(top_string, bot_string)))

    return str(ali_score), top_string, bot_string


if __name__ == '__main__':
    filename = sys.argv[1]

    with open(filename, 'r') as input_file:
        lines = input_file.readlines()

    with open('result.txt', 'w') as result:
        result.write('\n'.join(lin_space_ali(lines[0].strip(),
                                             lines[1].strip(),
                                             blosum62,
                                             5)))
