from collections import Counter


def ham_dis(x, y):
    """ Calculate the edit distance between s1 to s2 """
    dis = 0
    for i in range(0, len(x)):
        if x[i] != y[i]:
            dis += 1
    return dis


def edit_dis(s1, s2):
    """ Calculate the edit distance between s1 to s2 """
    m = len(s1) + 1
    n = len(s2) + 1

    tbl = {}
    for i in range(m):
        tbl[i, 0] = i
    for j in range(n):
        tbl[0, j] = j
    for i in range(1, m):
        for j in range(1, n):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            tbl[i, j] = min(tbl[i, j - 1] + 1, tbl[i - 1, j] + 1, tbl[i - 1, j - 1] + cost)

    return tbl[i, j]


def jaccard(s1, s2, n):
    """ Calculate the Jaccard index for string similarity between s1 to s2 \
    where the tokens are all the substrings in"""
    n_gram_dict = {}
    counter = 0
    for i in range(len(s1)-(n-1)):
        n_gram_dict[s1[i:i + n]] = 1
    for i in range(len(s2)-(n-1)):
        if s2[i:i + n] in n_gram_dict:
            n_gram_dict[s2[i:i + n]] = -1
        else:
            counter += 1
    # print(n_gram_dict)
    # print([key for key, item in n_gram_dict.items() if item < 1])
    return 100*len([key for key, item in n_gram_dict.items() if item < 1])/(len(n_gram_dict)+counter)


def GPM_quick_ratio(s1: str, s2: str) -> float:
    """Return an upper bound on Gestalt Pattern Matching relatively quickly."""
    length = len(s1) + len(s2)

    if not length:
        return 1.0

    intersect = Counter(s1) & Counter(s2)
    matches = sum(intersect.values())
    return (2.0 * matches / length) * 100