import argparse
import os
import shutil
import threading
import random
from datetime import datetime

import numpy as np

from simulator import *
import time
import math
# random.seed(datetime.now())

Number_of_steps = 220
Windows_size = 4
Similarity_threshold = 15


def file_to_clustering(file_path):
    """
    Return a tuple of all the attributes for ClusteringInfo object from a givan file in the evyat format.
    The result tuple contains:
    |   clustering: A dict in the form: {cluster_id (int): list of the cluster reads (list of strings)}
    |   original_strand_dict: {strand_id (int): the actual read (string)}
    |   reads_err_original_strand_dict: {read_id (int): the id of the origin strand of the read (int)}
    |   reads_err: A list of all the reads. the i'th element is the read with read_id = i.
    |   reads_err_dict: {read_id (int): the read itself (string)}
    |   C_dict: { cluster rep (full str) : List of all the reads that belong to that cluster }
    |   C_reps: [(Read, Cluster rep of the cluster to which the read belongs to)] must be sorted lexicographic.
    """
    reads_err = []  # will have all the reads from the sequencing stage
    reads_err_dict = {}
    strand_id = 0
    clustering = {}
    original_strand_dict = {}  # map from orig strand id to the actual strand
    reads_err_original_strand_dict = {}  # map from read_err to it's orig strand id
    C_dict = {}  # { cluster rep (full str) : List of all the reads that belong to that cluster }
    C_reps = []  # [(Read, Cluster rep of the cluster to which the read belongs to)] must be sorted lexicographic by read.
    with open(file_path, 'r') as evyat_f:
        line = "j"
        while line:
            line = evyat_f.readline()
            if not line:
                break
            original_strand_dict.update({strand_id: line.strip()})
            clustering[strand_id] = []
            line = evyat_f.readline()  # line == '*****************************\n':
            line = evyat_f.readline()
            cluster = []
            while line != '\n':
                striped_line = line.strip()
                reads_err.append(striped_line)
                clustering[strand_id].append(striped_line)
                reads_err_dict.update({len(reads_err)-1: striped_line})
                reads_err_original_strand_dict.update({len(reads_err)-1: strand_id})
                C_reps.append((striped_line, original_strand_dict[strand_id]))
                cluster.append(striped_line)
                line = evyat_f.readline()
            C_dict.update({original_strand_dict[strand_id]: cluster})
            strand_id = strand_id + 1
            line = evyat_f.readline()
    C_reps = sorted(C_reps, key=lambda x: x[0])
    return clustering, original_strand_dict, reads_err_original_strand_dict, reads_err, reads_err_dict, C_dict, C_reps


def hash_fun(x, a, w, l):
    """ An implementation of the hash function from the article section 3.2"""
    # assume that x has substring a
    # but what happens if it return -1???
    ind = x.find(a)
    return x[ind:min(len(x), ind + w + l)]


def dna_str_to_decimal(st):
    """
    return the value of the DNA string in base 10.
    We interpreted DNA string as a number in base 4
    where A = 0 | C = 1 | G = 2 | T = 3
    """
    # All 3-grams: AAA,AAG,AAC,AAT,AGA,AGC,...,TTA,TTG,TTC,TTT
    # A-0, G-1, C-2, T-3
    # index of CAT = Decimal representaion of (203)
    # it's actually 2*4^2 + 0*4^1 + 3*4^0 its quad not decimal
    N_q = {"A": 0, "C": 1, "G": 2, "T": 3}
    dec = 0
    for i in range(0, len(st)):
        dec += N_q[st[i]] * (4 ** (len(st) - i - 1))
    return dec


def decimal_to_dna_str(num):
    """
    return the value of a decimal number in base 4 in which:
            A = 0 | C = 1 | G = 2 | T = 3
    """
    res = ''
    map_dict = {0: "A", 1: "C", 2: "G", 3: "T"}
    i = num
    while i != 0:
        res = map_dict[i % 4] + res
        i = math.floor(i / 4)
    return res


def bin_sig(x, q):
    """ An implementation of the binary signature from the article section 3.3"""
    bs = [0] * (4 ** q)  # if bs[i] = 1 it means that there is at least 1 substring
                         # of x that in quad base means i
                         # for example: if x is CAT... then CAT is substring in x
                         # and it's value is 203(in base 4) = i(in base 10)
    for i in range(0, len(x) - q + 1):
        st = x[i:i + q]
        bs[dna_str_to_decimal(st)] = 1
    bs_str = ''.join(str(e) for e in bs)
    return bs_str


def ham_dis(x, y):
    """ Calculate the edit distance between s1 to s2 """
    dis = 0
    for i in range(0, len(x)):
        if x[i] != y[i]:
            dis += 1
    return dis


def rand_perm(w):
    """ return random DNA strand of length w"""
    return ''.join(random.choice('ACGT') for _ in range(w))


def gen_rand_input(strand_len, num_of_strands, file_path=None):
    """
    Generate num_of_strands random strands each in strand_len length.
    Returns (list of this strands, string of all the strands).
    optionally also write is to a given file.
    """
    strands = []
    res_str = ''
    for i in range(0, num_of_strands):
        strand = rand_perm(strand_len)
        res_str += strand + "\n"
        strands.append(strand)
    if file_path is not None:
        with open(file_path, 'w', newline='\n') as f:
            f.write(res_str)
    return strands, res_str


def rep_find(inp, reads_leaders):
    """
    | Args:
    |    inp is a read id
    |    parent[i] = the minimum id of all the read in the cluster that contains read i
    |    e.g. if read i is in the cluster {k, i, j} where k < i < j
    |    than parent[k], parent[i], parent[j] = k
    | return:
    |    the id of the cluster that contains the read with id=ind
    """
    temp = inp
    cnt = 0
    while reads_leaders[temp] != temp and cnt < 10:
        cnt += 1
        temp = reads_leaders[temp]
    return temp


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


def create_clustering(reads_leaders):
    """
    | Args:
    |     reads_leaders[i] = the minimum id of all the read in the cluster that contains read i
    |     e.g. if read i is in the cluster {k, i, j} where k < i < j
    |     than reads_leaders[k], reads_leaders[i], reads_leaders[j] = k

    | return:
    |    Clustering that is construct from the parent list.
    |    The output is a dict where the keys are the id of a cluster
    |    and the values are the cluster - list of reads id.
    """
    clstr = {}
    for i in range(0, len(reads_leaders)):
        clstr[i] = []
    for i in range(0, len(reads_leaders)):
        clstr[rep_find(i, reads_leaders)].append(i)
    return clstr


def min_max(val1, val2):
    """ return the tuple: (min{val1, val2}, max{val1, val2})"""
    min_val = min(val1, val2)
    max_val = max(val1, val2)
    return min_val, max_val


def condition0(bin_sig_arr, hash_C_til, id1, read1, read2, th_low, th_high, r, index_size=0):
    return (((ham_dis(bin_sig_arr[hash_C_til[id1][0]], bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_low) or
             ((ham_dis(bin_sig_arr[hash_C_til[id1][0]],
                       bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_high) and
              edit_dis(read1, read2) <= r)))


def condition1(bin_sig_arr, hash_C_til, id1, read1, read2, th_low, th_high, r, index_size):
    if index_size == 0:
        return condition0(bin_sig_arr, hash_C_til, id1, read1, read2, th_low, th_high, r)
    return ((index_size == 0 or (edit_dis(read1[:index_size], read2[:index_size]) <= 3)) and
            ((ham_dis(bin_sig_arr[hash_C_til[id1][0]], bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_low) or
             ((ham_dis(bin_sig_arr[hash_C_til[id1][0]],
                       bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_high) and
              edit_dis(read1, read2) <= r)))


def condition2(bin_sig_arr, hash_C_til, id1, read1, read2, th_low, th_high, r, index_size):
    if index_size == 0:
        return condition0(bin_sig_arr, hash_C_til, id1, read1, read2, th_low, th_high, r)
    index_ed = edit_dis(read1[:index_size], read2[:index_size])
    th_low_new = th_low + (index_size - index_ed)
    th_high_new = th_high + (index_size - index_ed)
    return ((ham_dis(bin_sig_arr[hash_C_til[id1][0]], bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_low_new) or
            ((ham_dis(bin_sig_arr[hash_C_til[id1][0]], bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_high_new) and
            edit_dis(read1, read2) <= r))


def condition3(bin_sig_arr, hash_C_til, id1, read1, read2, th_low, th_high, r, index_size):
    if index_size == 0:
        return condition0(bin_sig_arr, hash_C_til, id1, read1, read2, th_low, th_high, r)
    index_ed = edit_dis(read1[:index_size], read2[:index_size])
    th_low_new = th_low - index_ed
    th_high_new = th_high - index_ed
    return ((ham_dis(bin_sig_arr[hash_C_til[id1][0]], bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_low_new) or
            ((ham_dis(bin_sig_arr[hash_C_til[id1][0]], bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_high_new) and
            edit_dis(read1, read2) <= r))


def condition4(bin_sig_arr, hash_C_til, id1, read1, read2, th_low, th_high, r, index_size):
    if index_size == 0:
        return condition0(bin_sig_arr, hash_C_til, id1, read1, read2, th_low, th_high, r)
    index_ed = edit_dis(read1[:index_size], read2[:index_size])
    th_low_new = th_low + (math.floor(index_size/2) - index_ed)
    th_high_new = th_high - index_ed
    return ((ham_dis(bin_sig_arr[hash_C_til[id1][0]], bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_low_new) or
            ((ham_dis(bin_sig_arr[hash_C_til[id1][0]], bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_high_new) and
            edit_dis(read1, read2) <= r))


def condition5(bin_sig_arr, hash_C_til, id1, read1, read2, th_low, th_high, r, index_size):
    if index_size == 0:
        return condition0(bin_sig_arr, hash_C_til, id1, read1, read2, th_low, th_high, r)
    index_ed = edit_dis(read1[:index_size], read2[:index_size])
    th_low_new = th_low - index_ed
    th_high_new = th_high - index_ed
    if ham_dis(bin_sig_arr[hash_C_til[id1][0]], bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_low_new:
        return True
    if ham_dis(bin_sig_arr[hash_C_til[id1][0]], bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_high_new:
        data_ed = edit_dis(read1, read2)
        score = alpha * index_ed + beta * data_ed
        if score < threshold:
            return True
    return False


def hash_based_cluster(reads_err, number_of_steps=Number_of_steps, windows_size=Windows_size,
                       similarity_threshold=Similarity_threshold, index_size=0, cond_func=condition0):
    """
    Implementation of the microsoft hash based clustering algorithm.
    | Args:
    |   reads_err: reads_err: A list in which the i'th element is the string DNA of the read with id = i
    |   number_of_steps: The number of iterations of the algorithm. Default is 220.
    |   windows_size: A parameter for calculate the binary signature. Default is 4.
    |   similarity_threshold: A bound for the algorithm. Default is 15.
    |   index_size: the number of chars in each index. e.g. for AAA,AAC,AAG,AAT,.... the size is 3.
    Returns a clustering. A list of clusters, each cluster in a sorted list if all the cluster reads ids.
    """
    # reads_err = clustering_info.reads_err  # will have all the reads from the sequencing stage

    reads_err_ind = [0] * (len(reads_err))
    read_leaders = [0] * (len(reads_err))
    bin_sig_arr = []

    local_comm = number_of_steps

    q = windows_size
    # computing the binary signature for all the reads
    for i in range(0, len(reads_err)):
        reads_err_ind[i] = (i, reads_err[i][index_size:])
        read_leaders[i] = i
        bin_sig_arr.append(bin_sig(reads_err[i][index_size:], q))

    C_til = create_clustering(read_leaders)  # short for C^(~) i.e. C tilda
    dist_arr = []

    # th_low, th_high:
    # finding Hamming distance from the first read
    for i in range(1, len(reads_err)):
        dist_arr.append(ham_dis(bin_sig_arr[i], bin_sig_arr[0]))

    # sort the dist array
    dist_arr.sort()
    # the weird definition of theta_low and theta_high
    # I should probably ask why this specific theta's
    for i in range(0, 1000):
        if dist_arr[i + 1] - dist_arr[i] > 10:
            k = i
            break
    th_low = min(dist_arr[i] + 5, dist_arr[i + 1] - 8)
    th_high = min(dist_arr[i + 1] - 5, th_low + 10)

    r = similarity_threshold
    # set parameters from the article section 5.1:
    w = math.ceil(math.log(len(reads_err[0][index_size:]), 4))
    l = math.ceil(math.log(len(reads_err), 4))
    for lcl_step in range(0, local_comm):
        # report_func(total_bar_size, (2 * tmp_bar) + lcl_step)
        hash_select = [0] * (len(reads_err))
        # picking random representatives samples from each cluster:
        for i in range(0, len(C_til)):
            if len(C_til[i]) != 0:
                hash_select[random.choice(C_til[i])] = 1

        # computing the hash value for each representative:
        a = rand_perm(w)
        hash_C_til = [0] * (len(reads_err))
        for i in range(0, len(reads_err_ind)):
            if hash_select[i] == 0:
                hash_C_til[i] = (i, "")
            else:
                hash_C_til[i] = (i, hash_fun(reads_err_ind[i][1], a, w, l))
        # sort the hash values by the strings values
        hash_C_til.sort(key=lambda x: x[1])

        cnt = 0
        # implementation of lines 8-10 in the algorithm in the article:
        for i in range(0, len(hash_C_til) - 1):
            if hash_C_til[i][1] == "":
                continue
            else:
                if hash_C_til[i][1] == hash_C_til[i + 1][1]:
                    x = reads_err[hash_C_til[i][0]][index_size:]
                    y = reads_err[hash_C_til[i + 1][0]][index_size:]
                    # (edit_dis(x[:5], y[:5]) <= 3) and
                    # if ((index_size == 0 or (edit_dis(x[:index_size], y[:index_size]) <= 3)) and
                    #     ((ham_dis(bin_sig_arr[hash_C_til[i][0]], bin_sig_arr[hash_C_til[i + 1][0]]) <= th_low) or
                    #         ((ham_dis(bin_sig_arr[hash_C_til[i][0]],
                    #                        bin_sig_arr[hash_C_til[i + 1][0]]) <= th_high) and
                    #          edit_dis(x, y) <= r))):
                    if cond_func(bin_sig_arr, hash_C_til, i, x, y, th_low, th_high, r, index_size):
                        cnt += 1
                        min_temp, max_temp = min_max(rep_find(hash_C_til[i][0], read_leaders),
                                                          rep_find(hash_C_til[i + 1][0], read_leaders))
                        # merging x and y clusters
                        C_til[min_temp].extend(C_til[max_temp])
                        C_til[max_temp] = []
                        read_leaders[max_temp] = min_temp

    return [sorted(x) for x in list(C_til.values()) if x != []]


def str_clustering_to_ids(cluster_info):
    """
    Make a the evyat file like the output of the algorithm.
    | Args:
    |   cluster_info: An instance of class ClusteringInfo
    Return:
        A list of clusters. Each cluster is a list of reads_ids
    """
    clustering = []
    for i in range(len(cluster_info.original_strand_dict())):
        clustering.append([])
    for key, value in cluster_info.reads_err_original_strand_dict.items():
        clustering[value].append(key)
    return [sorted(cluster) for cluster in clustering]


def find_unwanted_rebellious_reads(stats, gamma=0.5):
    """
    | Args:
    |   stats:  {index of the algo cluster: value}
    |           Where value is a dict:
    |           {id_of_origin_strand: tuple
    |           #element1: The ratio between the next two
    |           #element2: Number of reads that originated from this origin strand
    |           #element3: Size of the origin strand true cluster
    |   gamma:  A param that tell what we consider as a rebellious_reads.
    |           Default is 0.5.
    Returns:
        A dict similar to the arg stats. But contains only algo clusters that are considered rebellious_reads
        and contains only the information about them.
    """
    unwanted_rebellious_reads = {}
    for index, algo_cluster_stat in stats.items():
        is_rebellious = False
        if len(algo_cluster_stat.keys()) > 1:
            temp = {}
            for orig_cluster_id, stats_ in algo_cluster_stat.items():
                # temp = {}
                if stats_[0] < gamma and stats_[1] < 5:
                    temp[orig_cluster_id] = stats_
                    is_rebellious = True
            if is_rebellious:
                unwanted_rebellious_reads[index] = temp
    return unwanted_rebellious_reads


def find_unwanted_singletons(stats, gamma=0.5):
    """
    | Args:
    |   stats:  {index of the algo cluster: value}
    |           Where value is a dict:
    |           {id_of_origin_strand: tuple
    |           #element1: The ratio between the next two
    |           #element2: Number of reads that originated from this origin strand
    |           #element3: Size of the origin strand true cluster
    |   gamma:  A param that tell what we consider as a singleton.
    |           Default is 0.5.
    Returns:
        A dict similar to the arg stats. But contains only algo clusters that are considered singletons.
    """
    unwanted_singletons = {}
    for index, algo_cluster_stat in stats.items():
        if len(algo_cluster_stat.keys()) == 1:
            for orig_cluster_id, stats_ in algo_cluster_stat.items():
                if stats_[0] <= gamma and stats_[1] == 1:
                    unwanted_singletons[index] = algo_cluster_stat
    return unwanted_singletons


def find_unwanted_unions(stats, gamma=0.5):
    """
    | Args:
    |   stats:  {index of the algo cluster: value}
    |           Where value is a dict:
    |           {id_of_origin_strand: tuple
    |           #element1: The ratio between the next two
    |           #element2: Number of reads that originated from this origin strand
    |           #element3: Size of the origin strand true cluster
    |   gamma:  A param that tell us how much we consider as a union between two clusters.
    |           Default is 0.5.
    Returns:
        A dict similar to the arg stats. But contains only algo clusters that are unwanted unions.
    """
    unwanted_unions = {}
    for index, algo_cluster_stat in stats.items():
        if len(algo_cluster_stat.keys()) > 1:
            count = 0
            for orig_cluster_id, stats_ in algo_cluster_stat.items():
                if stats_[0] >= gamma:
                    count += 1
            if count >= 2:
                unwanted_unions[index] = algo_cluster_stat
    return unwanted_unions


def find_clusters_stats(algo_clustering, orig_cluster_info):
    """
    Returns tuple.
    The first element is a dict:
    {index of the algo cluster: value}
    Where value is a dict:
    {id_of_origin_strand: tuple
    #element1: The ratio between the next two
    #element2: Number of reads that originated from this origin strand
    #element3: Size of the origin strand true cluster

    The second element is the same as the first except for the keys.
    Now the first keys are id_of_origin_strand and the second is the index of the algo cluster.
    | Args:
    |   algo_clustering: The output of the hash based clustering algorithm - list of clusters.
    |                    Each cluster is a list of reads ids.
    |   orig_cluster_info: an object of class ClusteringInfo
    """
    clusters_sizes = {}
    res_stat = {}
    res_stat2 = {}
    cluster_index_dict = find_clusters_origins(algo_clustering, orig_cluster_info)
    for key, value in orig_cluster_info.clustering.items():
        clusters_sizes[key] = len(value)
    for index, cluster_stat in cluster_index_dict.items():
        res_stat[index] = {}
        for orig_id, algo_size in cluster_stat.items():
            percentage = algo_size/clusters_sizes[orig_id]
            res_stat[index][orig_id] = (percentage, algo_size, clusters_sizes[orig_id])
            if orig_id in res_stat2:
                res_stat2[orig_id][index] = (percentage, algo_size, clusters_sizes[orig_id])
            else:
                res_stat2[orig_id] = {index: (percentage, algo_size, clusters_sizes[orig_id])}
    return res_stat, res_stat2


def find_clusters_origins(algo_clustering, orig_cluster_info):
    """
    returns a dict:
    {index of the algo cluster: value}
    Where value is a dict:
    {id_of_origin_strand: the number of reads that originated from
    that original strand and in that specific algo cluster}

    | Args:
    |   algo_clustering: The output of the hash based clustering algorithm - list of clusters.
    |                    Each cluster is a list of reads ids.
    |   orig_cluster_info: an object of class ClusteringInfo
    """
    cluster_index_dict = {}
    # "clean" the output of the algo_clustering
    for index in range(len(algo_clustering)):
        if len(algo_clustering[index]) == 0:
            algo_clustering.pop(index)
    # for each cluster in the algo output check how many reads
    # are origin from the same strand.
    for index in range(len(algo_clustering)):
        cluster_index_dict[index] = {}
        for read_id in algo_clustering[index]:
            orig_id = orig_cluster_info.reads_err_original_strand_dict[read_id]
            if orig_id in cluster_index_dict[index]:
                cluster_index_dict[index][orig_id] += 1
            else:
                cluster_index_dict[index][orig_id] = 1
    return cluster_index_dict


def handle_singletons_with_index(algo_clustering, orig_cluster_info, index_size, threshold=6):
    reads_err = orig_cluster_info.reads_err
    # find singletons
    singletons = []
    for cluster_id, cluster in enumerate(algo_clustering):
        if len(cluster) == 1:
            singletons.append((cluster[0], cluster_id))
    # for each singleton search for his cluster's mates
    for singleton, singleton_id in singletons:
        largest_cluster = singleton_id
        for cluster_id, cluster in enumerate(algo_clustering):
            if cluster[0] != singleton:
                representatives = random.sample(cluster, min([threshold, len(cluster)]))
                identical_index_count = 0
                for read in representatives:
                    if reads_err[read][: index_size] == reads_err[singleton][: index_size]:
                        identical_index_count += 1
                if identical_index_count >= (1/2)*max([math.floor(len(representatives)), 1]):
                    if len(algo_clustering[largest_cluster]) > len(algo_clustering[cluster_id]):
                        largest_cluster = cluster_id
        if largest_cluster != singleton_id:
            algo_clustering[largest_cluster].extend(algo_clustering[singleton_id])
            algo_clustering[singleton_id] = []
    return [sorted(x) for x in algo_clustering if x != []]


def handle_singletons_with_index_ver2(algo_clustering, orig_cluster_info, index_size, threshold=6):
    reads_err = orig_cluster_info.reads_err
    # find singletons
    singletons = []
    for cluster_id, cluster in enumerate(algo_clustering):
        if len(cluster) == 1:
            singletons.append((cluster[0], cluster_id))
    # for each singleton search for his cluster's mates
    for singleton, singleton_id in singletons:
        largest_cluster = singleton_id
        for cluster_id, cluster in enumerate(algo_clustering):
            if len(cluster) >= len(algo_clustering[largest_cluster]):
                true_index = []
                if cluster[0] != singleton:
                    representatives = random.sample(cluster, min([threshold, len(cluster)]))
                    index_mat = np.array([list(reads_err[read][:index_size]) for read in representatives]).T
                    for row in index_mat:
                        row_list = list(row)
                        true_index.append(max(row_list, key=row_list.count))
                    str_index = ''.join(true_index)
                    # orig_id = orig_cluster_info.reads_err_original_strand_dict[representatives[0]]
                    # orig_index = orig_cluster_info.original_strand_dict[orig_id][: index_size]
                    # print(f'{str_index=} ; {orig_index=}')
                    ed = edit_dis(str_index, reads_err[singleton][: index_size])
                    if ed <= 1:
                        largest_cluster = cluster_id
        if largest_cluster != singleton_id:
            algo_clustering[largest_cluster].extend(algo_clustering[singleton_id])
            algo_clustering[singleton_id] = []
    return [sorted(x) for x in algo_clustering if x != []]


def arrange_clustering(algo_clustering, orig_cluster_info):
    """
    Arrange the output of the clustering algorithm so that the order of the clusters
    and the original strands will be similar to the original (true) clustering.
    | Args:
    |   algo_clustering: The output from the clustering algorithm.
    |   orig_cluster_info: all the known information on the true clustering.
    |                      It should be ClusteringInfo object.
    Returns the evyat format string of the arranged clustering.
    """
    res = ''
    for cluster in algo_clustering:
        orig_strand_candidates = []
        for cluster_element in cluster:
            orig_strand_candidates.append(orig_cluster_info.reads_err_original_strand_dict.get(cluster_element))
        orig_strand_id = max(orig_strand_candidates, key=orig_strand_candidates.count)
        res += str(orig_cluster_info.original_strand_dict.get(orig_strand_id)) + '\n'
        res += '*****************************\n'
        for cluster_element in cluster:
            res += orig_cluster_info.read_err_dict.get(cluster_element) + '\n'
        res += '\n\n'
    return res


def file_to_cluster(file_path):
    """ return a clustering from a file as a dict
    |   where the keys are cluster id and the value in a list of all the reads (strings) that are in the cluster.
    """
    strand_id = 0
    cluster = {}
    with open(file_path, 'r') as evyat_f:
        line = "j"
        while line:
            line = evyat_f.readline()
            if not line:
                break
            # cluster[strand_id] = [line.strip()]
            cluster[strand_id] = []
            line = evyat_f.readline()  # line == '*****************************\n':
            line = evyat_f.readline()
            while line != '\n':
                cluster[strand_id].append(line.strip())
                line = evyat_f.readline()
            strand_id = strand_id + 1
            line = evyat_f.readline()
    return cluster


def comp_clusters(cluster1, cluster2, gamma):
    """
    compute the indicator condition that is a part of the
    accuracy definition in the article section 2 definition 2.1
    | Args:
    |    cluster1: cluster from the algorithm clustering
    |    cluster2: cluster from the true clustering
    |   gamma: A parameter for tuning the 'strength' of the subset. Should be between 0.5 to 1

    Returns: 1 if cluster1 is subset of cluster2
             and the size of the Intersection between them is larger than gamma*|cluster2|
             otherwise return 0.
    """
    # if gamma < 0.5 or gamma > 1:
    #     print(f'gamma can only be value between 0.5 to 1. you gave gamma = {gamma}')
    #     return 0
    # if len(cluster1) < gamma * len(cluster2):
    #     return 0
    # for read1 in cluster1:
    #     exist_in_cluster2 = False
    #     for read2 in cluster2:
    #         if read1 == read2:
    #             exist_in_cluster2 = True
    #             break
    #     if not exist_in_cluster2:
    #         # means that cluster1 is not a subset of cluster2
    #         return 0
    # return 1

    # alternative implementation
    if gamma < 0.5 or gamma > 1:
        print(f'gamma can only be value between 0.5 to 1. you gave gamma = {gamma}')
        return 0
    if len(cluster1) < gamma * len(cluster2):
        return 0
    for read in cluster1:
        if read not in cluster2:
            return 0
    return 1


def calc_accuracy(algo_clustering, true_clustering, gamma):
    """
    calculate the accuracy of the algorithm output in the same way as in section 2.1 in the article.
    | Args:
    |   algo_clustering: A dict where the keys are cluster id and the value in a list of all the reads (strings)
    |                    that are in the cluster.
    |   true_clustering: Same as the above. This parameter is like reference point.
    |   gamma: A parameter for tuning the 'strength' of the subset. Should be between 0.5 to 1
    """
    if gamma < 0.5 or gamma > 1:
        print(f'gamma can only be value between 0.5 to 1. you gave gamma = {gamma}')
        return -1
    accuracy = 0

    for true_cluster in true_clustering.values():
        for algo_cluster in algo_clustering.values():
            if len(algo_cluster) >= 1:
                res = comp_clusters(algo_cluster, true_cluster, gamma)
                accuracy += res
                if res == 1:
                    break
    return accuracy/len(true_clustering)


#######################################
""" functions from the Clustering-algorithm-single-core-(Erlich-sampled dataset).ipynb file """
#######################################


def rep_in_C(read, C_reps):
    """
    return the representative of the given read from the given C_reps
    | Args:
    |   read - A strand of DNA (string)
    |   C_reps - A list of tuples in which the first element is a read (string) and the
    |            second element is the read representative (string) of the cluster that the first element belongs to.
    """
    lower = 0
    upper = len(C_reps) - 1
    while lower <= upper:
        mid = lower + int((upper - lower) / 2)
        #         print(upper,mid)
        res = -1
        if read == (C_reps[mid][0]):
            return C_reps[mid][1]
        if read > (C_reps[mid][0]):
            lower = mid + 1
        else:
            upper = mid - 1
    return -1


def comp_clstrs(alg_clstr, org_clstr, reads_err, gamma):
    """
    check if alg_cluster is at least gamma subset of the org_cluster.
    i.e there are at least gamma*|org_clstr| reads in org_clstr that are also in alg_clstr.
    | Args:
    |   alg_clstr: A list of reads ids
    |   org_clstr: A list of actual reads (i.e. the strings of DNA)
    |   reads_err: A list in which the i'th element is the string DNA of the read with id = i
    |   gamma: A parameter for tuning the 'strength' of the subset. Should be between 0.5 to 1
    """
    num_exist = 0
    if len(alg_clstr) > len(org_clstr):
        #         print(alg_clstr)
        return 0
    else:
        for i in range(0, len(alg_clstr)):
            flg_exist = 0
            for j in range(0, len(org_clstr)):
                if reads_err[alg_clstr[i]] == org_clstr[j]:
                # TODO: replace the ifs if you want to read from file: option2
                # if alg_clstr[i] == org_clstr[j]:
                    flg_exist = 1
                    num_exist += 1
                    break
            if flg_exist == 0:
                return 0
        if num_exist < gamma * len(org_clstr):
            return 0

        return 1


def calc_acrcy(clustering, reads_err, C_dict, C_reps, gamma):
    #     clustering = display_parent(parent)
    """
    calculate the accuracy of the algorithm output in a similar way to the article section 2.1
    | Args:
    |   clustering: A list of clusters. Each cluster is a sorted list of reads ids.
    |   read_err: a list of all the reads. reads_err[i] = the read with id i.
    |   C_dict: { cluster rep (full str) : List of all the reads that belong to that cluster }
    |   C_reps: [(Read, Cluster rep of the cluster to which the read belongs to)] must be sorted lexicographic by read.
    |   gamma: A parameter for tuning the 'strength' of the subset. Should be between 0.5 to 1
    """
    acrcy = 0
    for i in range(0, len(clustering)):
        if len(clustering[i]) >= 1:
            acrcy += comp_clstrs(clustering[i], C_dict[rep_in_C(reads_err[clustering[i][0]], C_reps)],
                                 reads_err, gamma)
            # TODO: replace the above line with the 2 lines below if you want to read from file: option2
            # c_til = clustering[i]
            # acrcy += comp_clstrs(c_til, C_dict[rep_in_C(c_til[0], C_reps)], reads_err, gamma)
    return acrcy/len(C_dict.keys())


class ClusteringInfo:
    """
    | Attributes:
    |   clustering: A dict in the form: {cluster_id (int): list of the cluster reads (list of strings)}
    |   original_strand_dict: {strand_id (int): the actual read (string)}
    |   reads_err_original_strand_dict: {read_id (int): the id of the origin strand of the read (int)}
    |   reads_err: A list of all the reads. the i'th element is the read with read_id = i.
    |   reads_err_dict: {read_id (int): the read itself (string)}
    |   C_dict: { cluster rep (full str) : List of all the reads that belong to that cluster }
    |   C_reps: [(Read, Cluster rep of the cluster to which the read belongs to)] must be sorted lexicographic.
    | methods:
    |   __init__: Initialize a class object via a given object of the same class
    |             or via a file in the format of evyat file.
    |   __str__: Return an evyat format string contains the output of the clustering algorithm.
    |   str_to_file: Writes the output of __str__ to the given file.
    """
    def __init__(self, clustering_info=None, file_path=None):
        clusters_info = None
        if clustering_info is not None:
            clusters_info = clustering_info
        elif file_path is not None:
            clusters_info = file_to_clustering(file_path)
        if clusters_info is not None:
            self.clustering = clusters_info[0]
            self.original_strand_dict = clusters_info[1]
            self.reads_err_original_strand_dict = clusters_info[2]
            self.reads_err = clusters_info[3]
            self.reads_err_dict = clusters_info[4]
            self.C_dict = clusters_info[5]
            self.C_reps = clusters_info[6]
        else:
            self.clustering = None
            self.original_strand_dict = None
            self.reads_err_original_strand_dict = None
            self.reads_err = None
            self.reads_err_dict = None
            self.C_dict = None
            self.C_reps = None

    def __str__(self, number_of_steps=Number_of_steps, windows_size=Windows_size,
                similarity_threshold=Similarity_threshold):
        """
        return a string that represent the output of the clustering algorithm.
        | Args for the clustering algorithm:
        |   number_of_steps: The number of iterations of the algorithm. Default is 220.
        |   windows_size: A parameter for calculate the binary signature. Default is 4.
        |   similarity_threshold: A bound for the algorithm. Default is 15.

        """
        algo_clustering = hash_based_cluster(self.reads_err, number_of_steps, windows_size, similarity_threshold)
        return arrange_clustering(algo_clustering, self)

    def str_to_file(self, file_path, number_of_steps=Number_of_steps, windows_size=Windows_size,
                    similarity_threshold=Similarity_threshold):
        """
        Write the output of the __str__ method to a given file.
        Meaning, write the result of the clustering algorithm to the given file.
        | Args for the clustering algorithm:
        |   number_of_steps: The number of iterations of the algorithm. Default is 220.
        |   windows_size: A parameter for calculate the binary signature. Default is 4.
        |   similarity_threshold: A bound for the algorithm. Default is 15.
        """
        with open(file_path, 'w', newline='\n') as f:
            f.write(self.__str__(number_of_steps, windows_size, similarity_threshold))
