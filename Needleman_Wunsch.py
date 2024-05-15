import numpy as np

X = np.full((3, 4, 2), (0, 0))
X[0][0][1] = 2
X[1][0][0] = -1
X[0][0] = (99, 99)
X[0][1] = (2, max(2, 5))



print(X)
print(X.shape[0])
print(X.shape[1])


def generateTable(S, T, file, match, mismatch, slippage, gap):
    # read FASTA into S and T sequences
    # in each table coordinate, index 0 refers to direction
    # index 1 refers to the score

    # direction:
    # 0 = diagonal
    # 1 = up
    # 2 = left
    n = len(S)
    m = len(T)

    X = np.full((n+1, m+1, 2), (0, 0))
    print('table before NW \n', X)
    # this generates an (n+1)x(m+1) sized matrix with each entry
    # holding a direction and a best score

    for i in range(n+1):
        for j in range(m+1):

            if i == 0:
                if j == 0:
                    X[i][j] = (2, 0)
                elif j == 1:
                    X[i][j] = (2, gap)
                else:
                    gap_cost = X[i][j-1][1] + gap
                    if T[j - 1] == T[j - 2]:
                        slip_cost = X[i][j - 1][1] + slippage
                        X[i][j] = (2, max(gap_cost, slip_cost))
                    else:
                        X[i][j] = (2, gap_cost)

            elif j == 0:
                if i == 0:
                    X[i][j] = (1, 0)
                elif i == 1:
                    X[i][j] = (1, gap)
                else:
                    gap_cost = X[i-1][j][1] + gap
                    if S[i - 1] == S[i - 2]:
                        slip_cost = X[i - 1][j][1] + slippage
                        X[i][j] = (1, max(gap_cost, slip_cost))
                    else:
                        X[i][j] = (1, gap_cost)

            else:
                dir_score = optDirAndScore(S, T, i, j, X, gap, slippage, match, mismatch)
                X[i][j] = dir_score
    return X


def optDirAndScore(S, T, i, j, X, gap, slippage, match, mismatch):
    # assume diagonal has best score
    diagonal = X[i-1][j-1][1] + scoringFunc(S[i-1], T[j-1], match, mismatch)
    bestDirection = 0
    bestScore = diagonal


    up = X[i-1][j][1]
    left = X[i][j-1][1]

    # if i = 1 and j = 1, then we're at the beginning of both sequences so there can't be a slippage
    if i == 1 and j == 1:
        up += gap
        if up > bestScore:
            bestScore = up
            bestDirection = 1

        left += gap
        if left > bestScore:
            bestScore = left
            bestDirection = 2

        return bestDirection, bestScore

    elif i == 1 and j != 1:
        # we're in the second row
        up += gap
        if up > bestScore:
            bestScore = up
            bestDirection = 1

        if T[j-1] == T[j-2]:
            left += slippage
        else:
            left += gap

        if left > bestScore:
            bestScore = left
            bestDirection = 2

        return bestDirection, bestScore


    elif i != 1 and j == 1:
        if S[i-1] == S[i-2]:
            up += slippage
        else:
            up += gap

        if up > bestScore:
            bestScore = up
            bestDirection = 1

        left += gap
        if left > bestScore:
            bestScore = left
            bestDirection = 2

        return bestDirection, bestScore

    else:
        # up puts a gap in T so we analyze S
        if S[i - 1] == S[i - 2]:
            up += slippage
            if up > bestScore:
                bestScore = up
                bestDirection = 1
        else:
            up += gap
            if up > bestScore:
                bestScore = up
                bestDirection = 1

        # left puts a gap  in S so we analyze T
        if T[j - 1] == T[j - 2]:
            left += slippage
            if left > bestScore:
                bestScore = left
                bestDirection = 2
        else:
            left += gap
            if left > bestScore:
                bestScore = left
                bestDirection = 2

    return (bestDirection, bestScore)


def scoringFunc(n1, n2, match, mismatch):
    if n1 == n2:
        return match
    else:
        return mismatch


def traceBack(X, S, T):
    aligned_S = ""
    aligned_T = ""
    i = X.shape[0] - 1
    j = X.shape[1] - 1

    while i != 0 or j != 0:
        cur_direction = X[i][j][0]
        if cur_direction == 0:
            aligned_S = S[i - 1] + aligned_S
            aligned_T = T[j - 1] + aligned_T
            i = i - 1
            j = j - 1
        elif cur_direction == 1:
            aligned_S = S[i - 1] + aligned_S
            aligned_T = "-" + aligned_T
            i = i - 1
        else:
            aligned_S = "-" + aligned_S
            aligned_T = T[j - 1] + aligned_T
            j = j - 1

    return aligned_S, aligned_T


short_s = "TCATAAATAATTTTGCTTGCTGAAGGAAGAAAAAGTGTTTTTCATAAACCCATTATCCAGGACTGTTTATAGCTGTTGGAAGGACTAGGTCTTCCCTAGCCCCCCCAGTGTGCAAGGGCAGTGAAGACTTGATTGTACAAAATACGTTTTGTAAATGTTGTGCTGTTAACACTGCAAATAAACTTGGTAGCAAACACTTC"
short_t = "TACTTCATAAATAATTTTGCTTGCTGAAGGAAGAGAACATATTTTTCATGAACTCATAATCTTGGACTGTTTATAGATATTGGGGGGAGTAGGTCTTCTCTGGCCCCCCCCAGTATGCAAGGGCAATGGAGACCTGATTATATAAAGTATGTTTATAAATGCTGTGCTGTTAATACTGCAAATAAACTAAATAGCAAACA"
long_s = "TTCAGGCTGTTGTTGGCTTAGGGCTGGAAGCACAGAGTGGCTTGGCCTCAAGAGAATAGCTGGTTTCCCTAAGTTTACTTCTCTAAAACCCTGTGTTCACAAAGGCAGAGAGTCAGACCCTTCAATGGAAGGAGAGTGCTTGGGATCGATTATGTGACTTAAAGTCAGAATAGTCCTTGGGCAGTTCTCAAATGTTGGAGTGGAACATTGGGGAGGAAATTCTGAGGCAGGTATTAGAAATGAAAAGGAAACTTGAAACCTGGGCATGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGTGGGCAGATCACTGGAGGTCAGGAGTTCGAAACCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAATACAGAAATTAGCCGGTCATGGTGGTGGACACCTGTAATCCCAGCTACTCAGGTGGCTAAGGCAGGAGAATCACTTCAGCCCGGGAGGTGGAGGTTGCAGTGAGCCAAGATCATACCACGGCACTCCAGCCTGGGTGACAGTGAGACTGTGGCTCAAAAAAAAAAAAAAAAAAAGGAAAATGAAACTAGAAGAGATTTCTAAAAGTCTGAGATATATTTGCTAGATTTCTAAAGAATGTGTTCTAAAACAGCAGAAGATTTTCAAGAACCGGTTTCCAAAGACAGTCTTCTAATTCCTCATTAGTAATAAGTAAAATGTTTATTGTTGTAGCTCTGGTATATAATCCATTCCTCTTAAAATATAAGACCTCTGGCATGAATATTTCATATCTATAAAATGACAGATCCCACCAGGAAGGAAGCTGTTGCTTTCTTTGAGGTGATTTTTTTCCTTTGCTCCCTGTTGCTGAAACCATACAGCTTCATAAATAATTTTGCTTGCTGAAGGAAGAAAAAGTGTTTTTCATAAACCCATTATCCAGGACTGTTTATAGCTGTTGGAAGGACTAGGTCTTCCCTAGCCCCCCCAGTGTGCAAGGGCAGTGAAGACTTGATTGTACAAAATACGTTTTGTAAATGTTGTGCTGTTAACACTGCAAATAAACTTGGTAGCAAACACTTC"
long_t = "TAAGAAATGGTCTTCCATGCCCGGGGAGCACCTCTCAGTCTTTGGTGCTTCCACACCCTACTTACTAAGTATTTTTCATGCATCCGCCTGGAAAAACCAGCTTCTACCTATGTGAGGGTGGGTGCCTTAAAGGTTTTCTAACCAAGTCCCCCTTTGAAGTCTGTCTTGAACACAAAAAGTTTTCACCTGAGAAAGCTTTATCCCACCAATTGAGCAAGATGCTAATCTGTTGCATATCGTTCTTGGAGAGACAGTAGAGTGGCTTGGCTGATCTGAGCTGGGCATGGTGGTACTTGCCTATAATCACAGCATTTGGTGGAAAGAGATGGGCATGTCGTGAGTTGTATATGACTGAAGAACAAATAGCTAGTGGTCCTGATTTGTTCCTTATAAAATTGCATGTTCATAAAGGCCAAGAGTCAAGTCCTTCAAGGAGATGGCTTGAGATCTGTGGTTTGACTTAAAATCAGTGCTCCTCAAGGCTCCAGGTTTAGTAGCCAGCACCAAACAACACCCAGAATAGGGTTTGTGCAGGACACTGAGGAAATTCAGAGCCAGATTATGTATTTGATGTGAATGTGAGGTTTTTCAAGTGCCAGTGTCAAGGAGAGCCTTCTAATTCCCCATAGCTTATGAGCAGTGTATCTGTTGTAGCAGCTCTAAGTAACCCATTCTTCTTTAAATACACGGACTCTGATATGAATATTTCATGTCTGTACAATGATGGATTCCACCAGGAAAGGAGTTGTTGCTTTCTTTGAAGTAATTTCTTCCCTTTTGCTCCCTGTTGCTGAGGCATACTACTTCATAAATAATTTTGCTTGCTGAAGGAAGAGAACATATTTTTCATGAACTCATAATCTTGGACTGTTTATAGATATTGGGGGGAGTAGGTCTTCTCTGGCCCCCCCCAGTATGCAAGGGCAATGGAGACCTGATTATATAAAGTATGTTTATAAATGCTGTGCTGTTAATACTGCAAATAAACTAAATAGCAAACA"


mismatch = -1
match = 1
slippage = -1
gap = -2

S = long_s
T = long_t

testTable = generateTable(S, T, None, match, mismatch, slippage, gap)
print(testTable)
newS, newT = traceBack(testTable, S, T)
print(newS)
print(newT)
final_score = testTable[len(S)][len(T)][1]
print("The Optimal Score of this Alignment is:", final_score)
