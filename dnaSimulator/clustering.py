import argparse
import os
import shutil
import threading
import random
from datetime import datetime
from collections import Counter

import numpy as np

from simulator import *
from metrics import *
from evaluation import *
from gen_input import *
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
    #######################################
    # Add shuffle on reads_err_ind
    reads_err_ind_shuffle = random.sample(reads_err_ind, k=len(reads_err_ind))
    reads_err_ind_shuffle_ids = [0] * len(reads_err_ind_shuffle)
    for index, (read_id, read) in enumerate(reads_err_ind_shuffle):
        reads_err_ind_shuffle_ids[read_id] = index
    # print(reads_err_ind[0:10])
    # print(reads_err_ind_shuffle[0:10])
    #######################################
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
        for i in range(0, len(reads_err_ind_shuffle)):
            j = reads_err_ind_shuffle_ids[i]
            if hash_select[i] == 0:
                hash_C_til[i] = (reads_err_ind_shuffle[j][0], "")
            else:
                hash_C_til[i] = (reads_err_ind_shuffle[j][0], hash_fun(reads_err_ind_shuffle[j][1], a, w, l))
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

    return [sorted(x) for x in list(C_til.values()) if x != []], bin_sig_arr


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


def handle_singletons_with_index(algo_clustering, orig_cluster_info, index_size, threshold=100):
    reads_err = orig_cluster_info.reads_err
    # find singletons
    singletons = []
    for cluster_id, cluster in enumerate(algo_clustering):
        if len(cluster) == 1:
            singletons.append((cluster[0], cluster_id))
    # for each singleton search for his cluster's mates
    for singleton, singleton_cluster_id in singletons:
        largest_cluster = singleton_cluster_id
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
        if largest_cluster != singleton_cluster_id:
            algo_clustering[largest_cluster].extend(algo_clustering[singleton_cluster_id])
            algo_clustering[singleton_cluster_id] = []
    return [sorted(x) for x in algo_clustering if x != []]


def handle_singletons_with_index_ver2(algo_clustering, orig_cluster_info, bin_sign_arr, index_size, threshold=100):
    algo_clustering_copy = copy.deepcopy(algo_clustering)
    reads_err = orig_cluster_info.reads_err
    stat1, stat2 = find_clusters_stats(algo_clustering, orig_cluster_info)
    count_how_many_candidates_we_missed = 0
    count_true_unwanted_singletons = 0
    ham_dist_candidates_and_origins = {}
    ham_dist_from_true_algo_cluster = {}
    # find singletons
    unwanted_singletons = find_unwanted_singletons(stat1)
    singletons = []
    for cluster_id, cluster in enumerate(algo_clustering):
        if len(cluster) == 1:
            if cluster_id in unwanted_singletons.keys():
                count_true_unwanted_singletons += 1
            singletons.append((cluster[0], cluster_id))
    print(f'{count_true_unwanted_singletons=}')
    print(f'{len(singletons)=}')
    # for each singleton search for his cluster's mates
    for singleton, singleton_cluster_id in singletons:
        largest_cluster = singleton_cluster_id
        singleton_bin_sign = bin_sign_arr[singleton]  # bin_sig(reads_err[singleton][index_size:], 4)
        candidates = []
        for cluster_id, cluster in enumerate(algo_clustering):
            if len(cluster) >= len(algo_clustering[largest_cluster]) and cluster_id != singleton_cluster_id:
                true_index = []
                if cluster[0] != singleton: #and len(cluster) >= 35:
                    representatives = random.sample(cluster, min([threshold, len(cluster)]))
                    index_mat = np.array([list(reads_err[read][:index_size]) for read in representatives]).T
                    for row in index_mat:
                        row_list = list(row)
                        true_index.append(max(row_list, key=row_list.count))
                    str_index = ''.join(true_index)
                    # orig_id = orig_cluster_info.reads_err_original_strand_dict[representatives[0]]
                    # orig_index = orig_cluster_info.original_strand_dict[orig_id][: index_size]
                    # print(f'{str_index=} ; {orig_index=}')
                    # ed = edit_dis(str_index, reads_err[singleton][: index_size])
                    # if ed <= 1:
                    #     largest_cluster = cluster_id
                    ham_reps = []
                    for rep in representatives:
                        ham_reps.append(ham_dis(bin_sign_arr[rep], singleton_bin_sign))
                    if str_index == reads_err[singleton][: index_size] or \
                            edit_dis(reads_err[singleton][: index_size], str_index) <= 2 or \
                            np.mean(ham_reps) <= 95:
                        candidates.append(cluster_id)
        compare = []
        for cluster_id in candidates:
            representatives = random.sample(algo_clustering[cluster_id], min([len(algo_clustering[cluster_id]), 10]))
            ham = 0
            for i, rep in enumerate(representatives):
                ham += ham_dis(bin_sign_arr[rep], singleton_bin_sign)
            compare.append(ham/len(representatives))
        if len(compare) != 0:
            largest_cluster = candidates[np.argmin(compare)]
            origin_id = orig_cluster_info.reads_err_original_strand_dict[singleton]
            for key, value in stat2[origin_id].items():
                if value[0] >= 0.5:
                    reps_from_origin = random.sample(algo_clustering_copy[key],
                                                     min([len(algo_clustering_copy[key]), 10]))
                    ham = 0
                    for i, rep in enumerate(reps_from_origin):
                        ham += ham_dis(bin_sign_arr[rep], singleton_bin_sign)
                    ham_dist_candidates_and_origins[singleton] = (ham / len(reps_from_origin), min(compare),
                                                                  largest_cluster == key)
                    if key not in candidates:
                        ham = 0
                        count_how_many_candidates_we_missed += 1
                        true_algo_reps = random.sample(algo_clustering_copy[key],
                                                       min([len(algo_clustering_copy[key]), 10]))
                        index_mat = np.array([list(reads_err[read][:index_size]) for read in true_algo_reps]).T
                        true_index = []
                        for row in index_mat:
                            row_list = list(row)
                            true_index.append(max(row_list, key=row_list.count))
                        str_index = ''.join(true_index)
                        for true_rep in true_algo_reps:
                            ham += ham_dis(bin_sign_arr[true_rep], singleton_bin_sign)
                        ham_dist_from_true_algo_cluster[singleton] = (ham / len(true_algo_reps),
                                                                      ham_dist_candidates_and_origins[singleton][0],
                                                                      str_index, reads_err[singleton][: index_size])

            if largest_cluster != singleton_cluster_id:
                if min(compare) <= 90:
                    algo_clustering[largest_cluster].extend(algo_clustering[singleton_cluster_id])
                    algo_clustering[singleton_cluster_id] = []
        # print(f'singleton origin cluster in algo clustering = {stat2[origin_id]} ; {singleton_id =}; {candidates =}')
    print(f'{count_how_many_candidates_we_missed=}')
    print(f'{ham_dist_candidates_and_origins=}')
    print(f'{len(ham_dist_candidates_and_origins)=}')
    num_of_corrects = sum([int(x[2]) for x in ham_dist_candidates_and_origins.values()])
    sum_of_ham_dists_of_correct_reps = sum([int(x[1]) for x in ham_dist_candidates_and_origins.values() if x[2]])
    max_ham_dists_of_correct_reps = max(ham_dist_candidates_and_origins.values(), key=lambda x: x[1] if x[2] else 0)[1]
    num_of_incorrects = sum([1 for x in ham_dist_candidates_and_origins.values() if not x[2]])
    sum_of_ham_dists_of_incorrect_reps = sum([int(x[1]) for x in ham_dist_candidates_and_origins.values() if not x[2]])
    min_ham_dists_of_incorrect_reps = min(ham_dist_candidates_and_origins.values(), key=lambda x: x[1] if not x[2] else 1000)[1]
    sum_of_origin_ham_dists_of_incorrect_reps = sum([int(x[0]) for x in ham_dist_candidates_and_origins.values() if not x[2]])
    print(f'number of corrects: {num_of_corrects}')
    print(f'mean of corrects ham dist: {sum_of_ham_dists_of_correct_reps/num_of_corrects}')
    print(f'max ham dist of corrects: {max_ham_dists_of_correct_reps}')
    print(f'number of incorrects: {num_of_incorrects}')
    print(f'mean of incorrects ham dist: {sum_of_ham_dists_of_incorrect_reps / num_of_incorrects}')
    print(f'min ham dist of incorrects: {min_ham_dists_of_incorrect_reps}')
    print(f'mean of origin incorrects ham dist: {sum_of_origin_ham_dists_of_incorrect_reps / num_of_incorrects}')

    for singleton, stat in ham_dist_from_true_algo_cluster.items():
        print(f"singleton id: {singleton} ; true_algo_ham_mean: {stat[0]} ; origin_ham_mean: {stat[1]} ; "
              f"guessed index: {stat[2]} ; singleton index: {stat[3]}\n")

    return [sorted(x) for x in algo_clustering if x != []]


def handle_singletons_with_index_ver2_5(algo_clustering, orig_cluster_info, bin_sign_arr, index_size, threshold=100, log=True):
    # if log_path is None:
    #     log_path = "files/minion_idt/3000 strands in size 150 with x2 errors and cluster avg of 40/stats files/stats00_index.txt"
    # with open(log_path, 'w', newline='\n'):
    log_str = f'#################################################\n' \
              f'               handle_singletons_log             \n' \
              f'#################################################\n'
    algo_clustering_copy = copy.deepcopy(algo_clustering)
    reads_err = orig_cluster_info.reads_err
    stat1, stat2 = find_clusters_stats(algo_clustering, orig_cluster_info)
    count_how_many_candidates_we_missed = 0
    count_true_unwanted_singletons = 0
    ham_dist_candidates_and_origins = {}
    ham_dist_from_true_algo_cluster = {}
    wrong_right_dict = {}
    true_candidates_len_dict = {}
    count_missed_in_true_candidates = 0
    # find singletons
    unwanted_singletons = find_unwanted_singletons(stat1)
    singletons = []
    for cluster_id, cluster in enumerate(algo_clustering):
        if len(cluster) == 1:
            if cluster_id in unwanted_singletons.keys():
                count_true_unwanted_singletons += 1
            singletons.append((cluster[0], cluster_id))
    log_str += f'{count_true_unwanted_singletons=}\n'
    log_str += f'{len(singletons)=}\n'

    # find clusters reps and their index:
    clusters_reps = {}
    for cluster_id, cluster in enumerate(algo_clustering_copy):
        true_index = []
        representatives = random.sample(cluster, min([threshold, len(cluster)]))
        index_mat = np.array([list(reads_err[read][:index_size]) for read in representatives]).T
        for row in index_mat:
            row_list = list(row)
            true_index.append(max(row_list, key=row_list.count))
        str_index = ''.join(true_index)
        clusters_reps[cluster_id] = (representatives, str_index)

    # for each singleton search for his cluster's mates
    for singleton, singleton_cluster_id in singletons:
        largest_cluster = singleton_cluster_id
        singleton_index = reads_err[singleton][:index_size]
        singleton_bin_sign = bin_sign_arr[singleton]
        candidates = []
        for cluster_id, cluster in enumerate(algo_clustering_copy):
            if cluster_id != singleton_cluster_id:
                if cluster[0] != singleton:
                    representatives = clusters_reps[cluster_id][0]
                    str_index = clusters_reps[cluster_id][1]
                    ham_reps = []
                    for rep in representatives:
                        ham_reps.append(ham_dis(bin_sign_arr[rep], singleton_bin_sign))
                    if str_index == reads_err[singleton][: index_size] or \
                            edit_dis(reads_err[singleton][: index_size], str_index) <= 2 or \
                            np.mean(ham_reps) <= 95:
                        candidates.append(cluster_id)
        compare = []
        for cluster_id in candidates:
            representatives = clusters_reps[cluster_id][0]
            ham = 0
            for i, rep in enumerate(representatives):
                ham += ham_dis(bin_sign_arr[rep], singleton_bin_sign)
            compare.append(ham/len(representatives))
        if len(compare) != 0:
            # true candidates:
            min_ham = min(compare)
            min_ham_cluster_id = candidates[np.argmin(compare)]
            largest_cluster = min_ham_cluster_id
            true_candidates = [cluster_id for index, cluster_id in enumerate(candidates) if compare[index] <= min_ham + 6]
            true_candidates_compare = [compare[index] for index, cluster_id in enumerate(candidates) if compare[index] <= min_ham + 6]
            true_candidates_len_dict[singleton] = len(true_candidates)
            # TODO: continue with this
            best_cluster_index = clusters_reps[min_ham_cluster_id][1]
            diffs_dict = {}
            inside_ham_dist = {}
            if len(true_candidates) > 2:
                abs_diffs = {}
                for index, true_can in enumerate(true_candidates):
                    guessed_index = clusters_reps[true_can][1]
                    diff = edit_dis(best_cluster_index, singleton_index) - edit_dis(guessed_index, singleton_index)
                    diffs_dict[true_can] = diff
                    # TODO: check ham_dist inside each candidate to see if it closer to the singleton ham dist
                    ham = 0
                    for i, rep1 in enumerate(clusters_reps[true_can][0]):
                        for j, rep2 in enumerate(clusters_reps[true_can][0]):
                            if i < j:
                                ham += ham_dis(bin_sign_arr[rep1], bin_sign_arr[rep2])
                    inside_ham_dist[true_can] = 0 if ham == 0 else ham/((len(clusters_reps[true_can][0])**2 - len(clusters_reps[true_can][0]))/2)
                    abs_diffs[true_can] = abs(inside_ham_dist[true_can] - true_candidates_compare[index])
                largest_cluster = min(abs_diffs.items(), key=lambda pair: pair[1])[0]


            # best_cluster_id, max_diff = max(diffs_dict.items(), key=lambda pair: pair[1])
            # x = [(cluster_id, true_candidates_compare[i]) for i, cluster_id in enumerate(true_candidates) if diffs_dict[cluster_id] == max_diff]
            # x_cluster_id = min(x, key=lambda pair: pair[1])[0]
            # # if x_cluster_id != best_cluster_id:
            # #     print("what!!!!!!")
            # #     exit()
            # if max_diff > 0:
            #     largest_cluster = x_cluster_id # best_cluster_id
            # else:
            #     largest_cluster = min_ham_cluster_id
            # largest_cluster = candidates[np.argmin(compare)]
            origin_id = orig_cluster_info.reads_err_original_strand_dict[singleton]
            for key, value in stat2[origin_id].items():
                if value[0] >= 0.5 and key != singleton_cluster_id:
                    reps_from_origin = clusters_reps[key][0]
                    ham = 0
                    for i, rep in enumerate(reps_from_origin):
                        ham += ham_dis(bin_sign_arr[rep], singleton_bin_sign)
                    ham_dist_candidates_and_origins[singleton] = (ham / len(reps_from_origin), min(compare),
                                                                  largest_cluster == key)
                    if key not in candidates:
                        log_str += f"{singleton_cluster_id=} ; {value=}\n"
                        log_str += f"true algo cluster = {key} ; len = {len(algo_clustering_copy[key])} ; candidates = {candidates}\n"
                        ham1 = 0
                        ham2 = 0
                        count_how_many_candidates_we_missed += 1
                        true_algo_reps = clusters_reps[key]
                        for i, true_rep in enumerate(true_algo_reps[0]):
                            ham1 += ham_dis(bin_sign_arr[true_rep], singleton_bin_sign)
                            if i != 0:
                                ham2 += ham_dis(bin_sign_arr[true_rep], bin_sign_arr[true_algo_reps[0][i-1]])
                        ham_dist_from_true_algo_cluster[singleton] = (ham1 / len(true_algo_reps[0]),
                                                                      0 if ham2 == 0 else ham2 / (len(true_algo_reps[0])-1),
                                                                      true_algo_reps[1], reads_err[singleton][: index_size])
                        log_str += f"singleton cluster id: {singleton_cluster_id} ; "\
                                   f"true_algo_ham_mean: {ham_dist_from_true_algo_cluster[singleton][0]} ; "\
                                   f"origin_ham_mean: {ham_dist_from_true_algo_cluster[singleton][1]} ; "\
                                   f"guessed index: {ham_dist_from_true_algo_cluster[singleton][2]} ; "\
                                   f"singleton index: {ham_dist_from_true_algo_cluster[singleton][3]}\n\n"
                    else: # key is in candidates
                        if key not in true_candidates:
                            count_missed_in_true_candidates += 1
                        if largest_cluster != key:
                            ham = 0
                            true_algo_reps = clusters_reps[key]
                            for i, true_rep in enumerate(true_algo_reps[0]):
                                ham += ham_dis(bin_sign_arr[true_rep], singleton_bin_sign)
                            origin_mean_ham = ham/len(true_algo_reps[0])
                            winner_mean_ham = min(compare)
                            origin_guessed_index = true_algo_reps[1]
                            winner_guessed_index = clusters_reps[largest_cluster][1]
                            singleton_index = reads_err[singleton][:index_size]
                            origin_ed = edit_dis(origin_guessed_index, singleton_index)
                            winner_ed = edit_dis(winner_guessed_index, singleton_index)
                            wrong_right_dict[singleton] = {"origin_mean_ham": origin_mean_ham,
                                                           "min_mean_ham": winner_mean_ham,
                                                           "diff": origin_mean_ham-winner_mean_ham,
                                                           "origin_guessed_index": origin_guessed_index,
                                                           "winner_guessed_index": winner_guessed_index,
                                                           "the singleton index": singleton_index,
                                                           "origin_ed": origin_ed,
                                                           "winner_ed": winner_ed,
                                                           "diff_ed": winner_ed - origin_ed,
                                                           "diff_jacc": jaccard(origin_guessed_index, singleton_index, 2) - jaccard(winner_guessed_index, singleton_index, 2),
                                                           "diff_GPM": GPM_quick_ratio(origin_guessed_index, singleton_index) - GPM_quick_ratio(winner_guessed_index, singleton_index)}

            if largest_cluster != singleton_cluster_id:
                if min(compare) <= 95 \
                   or (edit_dis(clusters_reps[largest_cluster][1], reads_err[singleton][: index_size]) <= 1
                   and min(compare) <= 90):
                    algo_clustering[largest_cluster].extend(algo_clustering[singleton_cluster_id])
                    algo_clustering[singleton_cluster_id] = []
        # print(f'singleton origin cluster in algo clustering = {stat2[origin_id]} ; {singleton_id =}; {candidates =}')
    log_str += f'{ham_dist_candidates_and_origins=}\n'
    log_str += f'{len(ham_dist_candidates_and_origins)=}\n'
    num_of_corrects = sum([int(x[2]) for x in ham_dist_candidates_and_origins.values()])
    sum_of_ham_dists_of_correct_reps = sum([int(x[1]) for x in ham_dist_candidates_and_origins.values() if x[2]])
    max_ham_dists_of_correct_reps = max(ham_dist_candidates_and_origins.values(), key=lambda x: x[1] if x[2] else 0)[1]
    num_of_incorrects = sum([1 for x in ham_dist_candidates_and_origins.values() if not x[2]])
    sum_of_ham_dists_of_incorrect_reps = sum([int(x[1]) for x in ham_dist_candidates_and_origins.values() if not x[2]])
    min_ham_dists_of_incorrect_reps = min(ham_dist_candidates_and_origins.values(), key=lambda x: x[1] if not x[2] else 1000)[1]
    sum_of_origin_ham_dists_of_incorrect_reps = sum([int(x[0]) for x in ham_dist_candidates_and_origins.values() if not x[2]])
    log_str += f'number of corrects: {num_of_corrects}\n'
    log_str += f'mean of corrects ham dist: {sum_of_ham_dists_of_correct_reps/num_of_corrects}\n'
    log_str += f'max ham dist of corrects: {max_ham_dists_of_correct_reps}\n'
    log_str += f'number of incorrects: {num_of_incorrects}\n'
    log_str += f'mean of incorrects ham dist: {sum_of_ham_dists_of_incorrect_reps / num_of_incorrects}\n'
    log_str += f'min ham dist of incorrects: {min_ham_dists_of_incorrect_reps}\n'
    log_str += f'mean of origin incorrects ham dist: {sum_of_origin_ham_dists_of_incorrect_reps / num_of_incorrects}\n'

    for singleton, stat in ham_dist_from_true_algo_cluster.items():
        log_str += f"singleton id: {singleton} ; true_algo_ham_mean: {stat[0]} ; origin_ham_mean: {stat[1]} ; "\
                   f"guessed index: {stat[2]} ; singleton index: {stat[3]}\n\n"
    diff_ed_dict = Counter()
    diff_jac_dict = Counter()
    diff_gpm_dict = Counter()
    for singleton, wrong_dict in wrong_right_dict.items():
        diff_ed_dict.update([wrong_dict['diff_ed']])
        diff_jac_dict.update([wrong_dict['diff_jacc']])
        diff_gpm_dict.update([wrong_dict['diff_GPM']])
        log_str += f"{singleton=}; {wrong_dict}\n\n"
    log_str += f"{sorted(diff_ed_dict.items(), reverse=False, key=lambda x: x[0])=}\n"
    log_str += f"{sorted(diff_jac_dict.items(), reverse=True, key=lambda x: x[0])=}\n"
    log_str += f"{sorted(diff_gpm_dict.items(), reverse=True, key=lambda x: x[0])=}\n"
    log_str += f"{sorted(true_candidates_len_dict.values(), reverse=True)}\n"
    log_str += f"{max(list(true_candidates_len_dict.values()))=}\n"
    log_str += f"{np.mean(list(true_candidates_len_dict.values()))=}\n"
    log_str += f'number of candidates we missed in filter 1: {count_how_many_candidates_we_missed}\n'
    log_str += f"number of candidates we missed in filter 1: {count_missed_in_true_candidates}\n"

    if log:
        return [sorted(x) for x in algo_clustering if x != []], log_str
    return [sorted(x) for x in algo_clustering if x != []]


def handle_singletons_with_index_ver2_5_clean(algo_clustering, orig_cluster_info, bin_sign_arr, index_size, threshold=10):
    algo_clustering_copy = copy.deepcopy(algo_clustering)
    reads_err = orig_cluster_info.reads_err
    true_candidates_len_dict = {}
    # find singletons
    singletons = []
    for cluster_id, cluster in enumerate(algo_clustering):
        if len(cluster) == 1:
            singletons.append((cluster[0], cluster_id))

    # find clusters reps and their index:
    clusters_reps = {}
    for cluster_id, cluster in enumerate(algo_clustering_copy):
        true_index = []
        representatives = random.sample(cluster, min([threshold, len(cluster)]))
        index_mat = np.array([list(reads_err[read][:index_size]) for read in representatives]).T
        for row in index_mat:
            row_list = list(row)
            true_index.append(max(row_list, key=row_list.count))
        str_index = ''.join(true_index)
        clusters_reps[cluster_id] = (representatives, str_index)

    # for each singleton search for his cluster's mates
    for singleton, singleton_cluster_id in singletons:
        largest_cluster = singleton_cluster_id
        singleton_index = reads_err[singleton][:index_size]
        singleton_bin_sign = bin_sign_arr[singleton]
        candidates = []
        for cluster_id, cluster in enumerate(algo_clustering_copy):
            if cluster_id != singleton_cluster_id:
                if cluster[0] != singleton:
                    representatives = clusters_reps[cluster_id][0]
                    str_index = clusters_reps[cluster_id][1]
                    ham_reps = []
                    for rep in representatives:
                        ham_reps.append(ham_dis(bin_sign_arr[rep], singleton_bin_sign))
                    if str_index == reads_err[singleton][: index_size] or \
                            edit_dis(reads_err[singleton][: index_size], str_index) <= 2 or \
                            np.mean(ham_reps) <= 95:
                        candidates.append(cluster_id)
        compare = []
        for cluster_id in candidates:
            representatives = clusters_reps[cluster_id][0]
            ham = 0
            for i, rep in enumerate(representatives):
                ham += ham_dis(bin_sign_arr[rep], singleton_bin_sign)
            compare.append(ham/len(representatives))
        if len(compare) != 0:
            # true candidates:
            min_ham = min(compare)
            min_ham_cluster_id = candidates[np.argmin(compare)]
            largest_cluster = min_ham_cluster_id
            true_candidates = [cluster_id for index, cluster_id in enumerate(candidates) if compare[index] <= min_ham + 6]
            true_candidates_compare = [compare[index] for index, cluster_id in enumerate(candidates) if compare[index] <= min_ham + 6]
            true_candidates_len_dict[singleton] = len(true_candidates)
            # TODO: continue with this
            best_cluster_index = clusters_reps[min_ham_cluster_id][1]
            diffs_dict = {}
            inside_ham_dist = {}
            if len(true_candidates) > 2:
                abs_diffs = {}
                for index, true_can in enumerate(true_candidates):
                    guessed_index = clusters_reps[true_can][1]
                    diff = edit_dis(best_cluster_index, singleton_index) - edit_dis(guessed_index, singleton_index)
                    diffs_dict[true_can] = diff
                    # TODO: check ham_dist inside each candidate to see if it closer to the singleton ham dist
                    ham = 0
                    for i, rep1 in enumerate(clusters_reps[true_can][0]):
                        for j, rep2 in enumerate(clusters_reps[true_can][0]):
                            if i < j:
                                ham += ham_dis(bin_sign_arr[rep1], bin_sign_arr[rep2])
                    inside_ham_dist[true_can] = 0 if ham == 0 else ham/((len(clusters_reps[true_can][0])**2 - len(clusters_reps[true_can][0]))/2)
                    abs_diffs[true_can] = abs(inside_ham_dist[true_can] - true_candidates_compare[index])
                largest_cluster = min(abs_diffs.items(), key=lambda pair: pair[1])[0]
            if largest_cluster != singleton_cluster_id:
                if min(compare) <= 95 \
                   or (edit_dis(clusters_reps[largest_cluster][1], reads_err[singleton][: index_size]) <= 1
                   and min(compare) <= 90):
                    algo_clustering[largest_cluster].extend(algo_clustering[singleton_cluster_id])
                    algo_clustering[singleton_cluster_id] = []
    return [sorted(x) for x in algo_clustering if x != []]


def handle_singletons_with_index_ver3(algo_clustering, orig_cluster_info, index_size, threshold=6):
    reads_err = orig_cluster_info.reads_err
    # find singletons
    singletons = []
    for cluster_id, cluster in enumerate(algo_clustering):
        if len(cluster) == 1:
            singletons.append((cluster[0], cluster_id))
    # for each singleton search for his cluster's mates

    algo_clustering_index_bin_sign = []
    for cluster_id, cluster in enumerate(algo_clustering):
        representatives = random.sample(cluster, min([threshold, len(cluster)]))
        reps_bin_sign = [bin_sig(reads_err[rep][:index_size], 3) for rep in representatives]
        algo_clustering_index_bin_sign.append((reps_bin_sign, cluster_id))

    for singleton, singleton_cluster_id in singletons:
        dist_array = []
        singleton_bin_sign = bin_sig(reads_err[singleton][:index_size], 3)
        for reps_bin_sign, cluster_id in algo_clustering_index_bin_sign:
            if cluster_id == singleton_cluster_id or len(algo_clustering[cluster_id]) == 0:
                continue
            dist_list = []
            for i in range(len(reps_bin_sign)):
                dist_list.append(ham_dis(reps_bin_sign[i], singleton_bin_sign))
            # remove outliers
            if len(dist_list) > threshold/2:
                dist_list.remove(max(dist_list))
            dist_array.append((np.mean(np.array(dist_list)), cluster_id))

        # dist_array_sorted = sorted(dist_array, key=lambda x: x[0])
        # if dist_array_sorted[0][0] == dist_array_sorted[1][0]:
        #     print(singleton)
        if min(dist_array, key=lambda x: x[0])[0] <= 3:
            algo_clustering[min(dist_array, key=lambda x: x[0])[1]].extend(algo_clustering[singleton_cluster_id])
            algo_clustering[singleton_cluster_id] = []
    return [sorted(x) for x in algo_clustering if x != []]


def handle_singletons_with_index_ver4(algo_clustering, orig_cluster_info, bin_sign_arr, index_size, threshold=10):
    algo_clustering_copy = copy.deepcopy(algo_clustering)
    reads_err = orig_cluster_info.reads_err
    stat1, stat2 = find_clusters_stats(algo_clustering, orig_cluster_info)
    # true_candidates_len_dict = {}
    sub_sign_indices = []
    for i in range(5):
        indices = random.sample(list(range(len(bin_sign_arr[0]))), k=20)
        indices_sorted = sorted(indices)
        while indices_sorted in sub_sign_indices:
            indices = random.sample(list(range(len(bin_sign_arr[0]))), k=20)
            indices_sorted = sorted(indices)
        print(indices_sorted)
        sub_sign_indices.append(indices_sorted)

    candidates = {}
    hash_dict = {}
    # find singletons
    singletons = []
    for cluster_id, cluster in enumerate(algo_clustering):
        if len(cluster) == 1:
            singleton_id = cluster[0]
            candidates[singleton_id] = {}
            singletons.append((singleton_id, cluster_id))
            for indices in sub_sign_indices:
                sign_array = np.array(list(bin_sign_arr[singleton_id]))
                str_sub_sign = ''.join(sign_array[indices])

                if str_sub_sign in hash_dict:
                    hash_dict[str_sub_sign]['singletons_ids'].append(singleton_id)
                else:
                    hash_dict[str_sub_sign] = {'cluster_ids': [],
                                               'singletons_ids': [singleton_id]}

                # if singleton_id in candidates[str_sub_sign]:
                #     candidates[str_sub_sign][singleton_id] += 1
                # else:
                #     candidates[str_sub_sign][singleton_id] = 1

    # find clusters reps and their index:
    clusters_reps = {}
    for cluster_id, cluster in enumerate(algo_clustering_copy):
        # candidates[cluster_id] = {}
        true_index = []
        representatives = random.sample(cluster, min([threshold, len(cluster)]))
        index_mat = np.array([list(reads_err[read][:index_size]) for read in representatives]).T
        for row in index_mat:
            row_list = list(row)
            true_index.append(max(row_list, key=row_list.count))
        str_index = ''.join(true_index)
        clusters_reps[cluster_id] = (representatives, str_index)

        for indices in sub_sign_indices:
            temp = []
            for rep in representatives:
                sign_array = np.array(list(bin_sign_arr[rep]))
                temp.append(sign_array[indices])
            mat = np.array(temp).T
            sub_sign = []
            for row in mat:
                row_list = list(row)
                sub_sign.append(max(row_list, key=row_list.count))
            str_sub_sign = ''.join(sub_sign)

            if str_sub_sign in hash_dict:
                hash_dict[str_sub_sign]['cluster_ids'].append(cluster_id)
                for singleton_id in hash_dict[str_sub_sign]['singletons_ids']:
                    if cluster_id in candidates[singleton_id]:
                        candidates[singleton_id][cluster_id] += 1
                    else:
                        candidates[singleton_id][cluster_id] = 1

    # stats:
    # check candidates avg size
    # check number of singletons that we didn't find their true candidates
    count_how_many_candidates_we_missed = 0
    total_number_of_candidats = 0
    for singleton_id, singleton_cluster_id in singletons:
        origin_id = orig_cluster_info.reads_err_original_strand_dict[singleton_id]
        for key, value in stat2[origin_id].items():
            if value[0] >= 0.5 and key != singleton_cluster_id:
                if key not in candidates[singleton_id]:
                    count_how_many_candidates_we_missed += 1
        total_number_of_candidats += len(candidates[singleton_id])

    print(f"number of candidates we missed: {count_how_many_candidates_we_missed}")
    print(f"avg size of candidates: {total_number_of_candidats/len(singletons)}")
    print(f"{hash_dict=}")

    return [sorted(x) for x in algo_clustering if x != []], ''


def handle_singletons_with_index_ver5(algo_clustering, orig_cluster_info, bin_sign_arr, index_size, threshold=10, num_epochs=2):
    algo_clustering_copy = copy.deepcopy(algo_clustering)
    for epoch_id in range(num_epochs):
        print(f"epoch: {epoch_id}")
        reads_err = orig_cluster_info.reads_err
        stat1, stat2 = find_clusters_stats(algo_clustering_copy, orig_cluster_info)
        w = math.ceil(math.log(len(reads_err[0][index_size:]), 4))
        l = math.ceil(math.log(len(reads_err), 4))-3
        hashing_a = []
        for i in range(20):
            a = rand_perm(w)
            hashing_a.append(a)

        candidates = {}
        hash_dict = {}
        # find singletons
        singletons = []
        for cluster_id, cluster in enumerate(algo_clustering_copy):
            if len(cluster) == 1:
                singleton_id = cluster[0]
                candidates[singleton_id] = {}
                singletons.append((singleton_id, cluster_id))
                for a in hashing_a:
                    hash_value = hash_fun(reads_err[singleton_id][index_size:], a, w, l)
                    if hash_value == '':
                        continue
                    if hash_value in hash_dict:
                        hash_dict[hash_value]['singletons_ids'].append(singleton_id)
                    else:
                        hash_dict[hash_value] = {'cluster_ids': [],
                                                 'singletons_ids': [singleton_id]}

                    # if singleton_id in candidates[str_sub_sign]:
                    #     candidates[str_sub_sign][singleton_id] += 1
                    # else:
                    #     candidates[str_sub_sign][singleton_id] = 1
        print(f"there are {len(singletons)} singletons")
        # find clusters reps and their index:
        clusters_reps = {}
        for cluster_id, cluster in enumerate(algo_clustering_copy):
            # candidates[cluster_id] = {}
            true_index = []
            representatives = random.sample(cluster, min([threshold, len(cluster)]))
            index_mat = np.array([list(reads_err[read][:index_size]) for read in representatives]).T
            for row in index_mat:
                row_list = list(row)
                true_index.append(max(row_list, key=row_list.count))
            str_index = ''.join(true_index)
            clusters_reps[cluster_id] = (representatives, str_index)

            for a in hashing_a:
                for rep in representatives:
                    hash_value = hash_fun(reads_err[rep][index_size:], a, w, l)
                    if hash_value in hash_dict:
                        hash_dict[hash_value]['cluster_ids'].append(cluster_id)
                        for singleton_id in hash_dict[hash_value]['singletons_ids']:
                            if cluster_id in candidates[singleton_id]:
                                candidates[singleton_id][cluster_id] += 1
                            else:
                                candidates[singleton_id][cluster_id] = 1

        # stats:
        # check candidates avg size
        # check number of singletons that we didn't find their true candidates
        count_how_many_candidates_we_missed = 0
        total_number_of_candidats = 0
        number_of_zero_candidates = 0
        for singleton_id, singleton_cluster_id in singletons:
            origin_id = orig_cluster_info.reads_err_original_strand_dict[singleton_id]
            for key, value in stat2[origin_id].items():
                if value[0] >= 0.5 and key != singleton_cluster_id and value[0] != 1:
                    if key not in candidates[singleton_id]:
                        count_how_many_candidates_we_missed += 1
                    else:
                        algo_clustering_copy[singleton_cluster_id] = []
                    if len(candidates[singleton_id]) == 0:
                        number_of_zero_candidates += 1
            total_number_of_candidats += len(candidates[singleton_id])
        algo_clustering_copy = [a for a in algo_clustering_copy if a != []]

        print(f"number of candidates we missed: {count_how_many_candidates_we_missed}")
        print(f"avg size of candidates: {total_number_of_candidats/len(singletons)}")
        print(f"number of singletons with zero candidates: {number_of_zero_candidates}\n")
    # print(f"{hash_dict=}")
    exit()
    return [sorted(x) for x in algo_clustering_copy if x != []], ''


def find_best_cluster(c_til, singleton_id, singleton_cluster_id, candidates_ids, reads_err,
                      clusters_reps, bin_sign_arr, index_size=6):
    best_cluster_id = singleton_cluster_id
    # First filter
    candidates_ids = [c for c in candidates_ids if c != singleton_cluster_id]
    compare = []
    true_candidates_sizes = []
    for candidate_id in candidates_ids:
        reps = clusters_reps[candidate_id][0]
        ham = 0
        for i, rep in enumerate(reps):
            ham += ham_dis(bin_sign_arr[rep], bin_sign_arr[singleton_id])
        compare.append(ham / len(reps))
    if len(compare) != 0:
        # true candidates:
        min_ham = min(compare)
        min_ham_cluster_id = candidates_ids[np.argmin(compare)]
        if not (min_ham <= 95
                or (edit_dis(clusters_reps[min_ham_cluster_id][1], reads_err[singleton_id][: index_size]) <= 1
                    and min_ham <= 90)):
            return singleton_cluster_id
        true_candidates = [cluster_id for index, cluster_id in enumerate(candidates_ids) if compare[index] <= min_ham + 6]
        true_candidates_compare = [compare[index] for index, cluster_id in enumerate(candidates_ids) if compare[index] <= min_ham + 6]
        true_candidates_sizes.append(len(true_candidates))
        best_cluster_id = min_ham_cluster_id
        inside_ham_dist = {}
        if len(true_candidates) > 2:
            abs_diffs = {}
            for index, true_can in enumerate(true_candidates):
                ham = 0
                for i, rep1 in enumerate(clusters_reps[true_can][0]):
                    for j, rep2 in enumerate(clusters_reps[true_can][0]):
                        if i < j:
                            ham += ham_dis(bin_sign_arr[rep1], bin_sign_arr[rep2])
                inside_ham_dist[true_can] = 0 if ham == 0 else ham / ((len(clusters_reps[true_can][0]) ** 2 - len(clusters_reps[true_can][0])) / 2)
                abs_diffs[true_can] = abs(inside_ham_dist[true_can] - true_candidates_compare[index])
            best_cluster_id = min(abs_diffs.items(), key=lambda pair: pair[1])[0]

    return best_cluster_id


def handle_singletons_with_index_ver5_5(algo_clustering, orig_cluster_info, bin_sign_arr, index_size,
                                        threshold=10, num_epochs=2, num_of_hashes=20,
                                        return_stats=False):
    algo_clustering_copy = copy.deepcopy(algo_clustering)
    times_for_each_epoch = [0]
    num_of_remaining_singletons = []
    stat1, stat2 = find_clusters_stats(algo_clustering_copy, orig_cluster_info)
    num_of_remaining_singletons.append(len(find_unwanted_singletons(stat1)))
    for epoch_id in range(num_epochs):
        # print(f"epoch: {epoch_id}")
        start = time.perf_counter_ns()
        reads_err = orig_cluster_info.reads_err
        stat1, stat2 = find_clusters_stats(algo_clustering_copy, orig_cluster_info)
        w = math.ceil(math.log(len(reads_err[0][index_size:]), 4))
        l = math.ceil(math.log(len(reads_err), 4))-3
        hashing_a = []
        for i in range(num_of_hashes):
            a = rand_perm(w)
            hashing_a.append(a)

        candidates = {}
        hash_dict = {}
        # find singletons
        singletons = []
        for cluster_id, cluster in enumerate(algo_clustering_copy):
            if len(cluster) == 1:
                singleton_id = cluster[0]
                candidates[singleton_id] = {}
                singletons.append((singleton_id, cluster_id))
                for a in hashing_a:
                    hash_value = hash_fun(reads_err[singleton_id][index_size:], a, w, l)
                    if hash_value == '':
                        continue
                    if hash_value in hash_dict:
                        hash_dict[hash_value]['singletons_ids'].append(singleton_id)
                    else:
                        hash_dict[hash_value] = {'cluster_ids': [],
                                                 'singletons_ids': [singleton_id]}

                    # if singleton_id in candidates[str_sub_sign]:
                    #     candidates[str_sub_sign][singleton_id] += 1
                    # else:
                    #     candidates[str_sub_sign][singleton_id] = 1
        # print(f"there are {len(singletons)} singletons")
        # find clusters reps and their index:
        clusters_reps = {}
        for cluster_id, cluster in enumerate(algo_clustering_copy):
            # candidates[cluster_id] = {}
            true_index = []
            representatives = random.sample(cluster, min([threshold, len(cluster)]))
            index_mat = np.array([list(reads_err[read][:index_size]) for read in representatives]).T
            for row in index_mat:
                row_list = list(row)
                true_index.append(max(row_list, key=row_list.count))
            str_index = ''.join(true_index)
            clusters_reps[cluster_id] = (representatives, str_index)

            for a in hashing_a:
                for rep in representatives:
                    hash_value = hash_fun(reads_err[rep][index_size:], a, w, l)
                    if hash_value in hash_dict:
                        hash_dict[hash_value]['cluster_ids'].append(cluster_id)
                        for singleton_id in hash_dict[hash_value]['singletons_ids']:
                            if cluster_id in candidates[singleton_id]:
                                candidates[singleton_id][cluster_id] += 1
                            else:
                                candidates[singleton_id][cluster_id] = 1

        for singleton_id, singleton_cluster_id in singletons:
            cluster_id = find_best_cluster(algo_clustering_copy, singleton_id, singleton_cluster_id,
                                           list(candidates[singleton_id].keys()), reads_err, clusters_reps,
                                           bin_sign_arr, index_size)
            if cluster_id != singleton_cluster_id:
                algo_clustering_copy[cluster_id].extend(algo_clustering_copy[singleton_cluster_id])
                algo_clustering_copy[singleton_cluster_id] = []

        algo_clustering_copy = [x for x in algo_clustering_copy if x != []]
        end = time.perf_counter_ns()
        elapsed_in_ns = end - start
        times_for_each_epoch.append((elapsed_in_ns * math.pow(10, -9))/60)
        stat1, stat2 = find_clusters_stats(algo_clustering_copy, orig_cluster_info)
        num_of_remaining_singletons.append(len(find_unwanted_singletons(stat1)))
    if return_stats:
        return [sorted(x) for x in algo_clustering_copy if x != []], '',\
               num_of_remaining_singletons, times_for_each_epoch
    return [sorted(x) for x in algo_clustering_copy if x != []], ''


def separate_cluster(cluster, orig_cluster_info, bin_sign_arr, index_size, threshold=100):
    # hoping that in 4 random reps we'll get at least one rep from each cluster
    reads_err = orig_cluster_info.reads_err
    n = min([10, len(cluster)])
    reps = random.sample(cluster, n)
    ham_dist_matrix = np.zeros((n, n))
    # edit_dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            # edit_dist_matrix[i][j] = edit_dis(reads_err[reps[i]][:index_size], reads_err[reps[j]][:index_size])
            ham_dist_matrix[i][j] = ham_dis(bin_sign_arr[reps[i]], bin_sign_arr[reps[j]])
    # edi, edj = np.unravel_index(edit_dist_matrix.argmax(), edit_dist_matrix.shape)
    hami, hamj = np.unravel_index(ham_dist_matrix.argmax(), ham_dist_matrix.shape)
    ham_within_cluster_1 = min([ham_j for j, ham_j in enumerate(ham_dist_matrix[hami]) if hami != j])
    ham_within_cluster_2 = min([ham_i for i, ham_i in enumerate(ham_dist_matrix[hamj]) if hamj != i])
    # mean_ham = np.mean([i for i in ham_dist_matrix.reshape(-1) if i != 0])
    # print(f'{mean_ham=}')
    # print(f'{edit_dist_matrix[edi][edj]=}')
    # print(f'{ham_dist_matrix[hami][hamj]=}')
    # print(f"{(ham_dist_matrix[hami][hamj] - ham_within_cluster_1)=}")
    diff1 = ham_dist_matrix[hami][hamj] - ham_within_cluster_1
    diff2 = ham_dist_matrix[hami][hamj] - ham_within_cluster_2
    if ham_dist_matrix[hami][hamj] < 115 or diff1 < 30 or diff2 < 30:
        return None, None
    # if edit_dist_matrix[edi][edj] >= 3:
    #     cluster1 = [reps[edi]]
    #     cluster2 = [reps[edj]]
    # else:
    cluster1 = [reps[hami]]
    cluster2 = [reps[hamj]]
    for read_id in cluster:
        if read_id != reps[hami] and read_id != reps[hamj]:
            ham1 = ham_dis(bin_sign_arr[cluster1[0]], bin_sign_arr[read_id])
            ham2 = ham_dis(bin_sign_arr[cluster2[0]], bin_sign_arr[read_id])
            if ham1 < ham2:
                cluster1.append(read_id)
            else:
                cluster2.append(read_id)
    # To be sure:
    dists_inside_cluster1 = []
    for i in range(len(cluster1)):
        for j in range(i + 1, len(cluster1)):
            dists_inside_cluster1.append(ham_dis(bin_sign_arr[cluster1[i]], bin_sign_arr[cluster1[j]]))
    dists_inside_cluster2 = []
    for i in range(len(cluster2)):
        for j in range(i + 1, len(cluster2)):
            dists_inside_cluster2.append(ham_dis(bin_sign_arr[cluster2[i]], bin_sign_arr[cluster2[j]]))

    mean_dists_inside_cluster1 = 0 if len(dists_inside_cluster1) == 0 else np.mean(dists_inside_cluster1)
    mean_dists_inside_cluster2 = 0 if len(dists_inside_cluster2) == 0 else np.mean(dists_inside_cluster2)
    if mean_dists_inside_cluster1 > 100 and mean_dists_inside_cluster2 > 100:
        return None, None
    dists_between_clusters = []
    for i in range(len(cluster1)):
        for j in range(len(cluster2)):
            dists_between_clusters.append(ham_dis(bin_sign_arr[cluster1[i]], bin_sign_arr[cluster2[j]]))

    mean_dists_between_clusters = 0 if len(dists_between_clusters) == 0 else np.mean(dists_between_clusters)
    if mean_dists_between_clusters-min([mean_dists_inside_cluster1, mean_dists_inside_cluster2]) <= 25:
        return None, None
    return cluster1, cluster2


def handle_unions(algo_clustering, orig_cluster_info, bin_sign_arr, index_size, threshold=10, log=True):
    avg_cluster_size = sum([len(cluster) for cluster in algo_clustering])/len(algo_clustering)
    count_wrong = 0
    count_right = 0
    sep_dict = {}
    log_str = f'#################################################\n' \
              f'               handle_unions_log                 \n' \
              f'#################################################\n'
    # for checks only
    stats1, stats2 = find_clusters_stats(algo_clustering, orig_cluster_info)
    unwanted_unions = find_unwanted_unions(stats1)
    # print(avg_cluster_size)
    for cluster_id, cluster in enumerate(algo_clustering):
        if len(cluster) >= avg_cluster_size:
            cluster1, cluster2 = separate_cluster(cluster, orig_cluster_info, bin_sign_arr, index_size, threshold=100)
            if cluster1 is not None:
                sep_dict[cluster_id] = (cluster1, cluster2)
                if cluster_id in unwanted_unions:
                    # print(unwanted_unions[cluster_id])
                    count_right += 1
                else:
                    log_str += f"{stats1[cluster_id]}\n"
                    log_str += f"{len(cluster1)=}\n{len(cluster2)=}\n"
                    count_wrong += 1
    log_str += f"{count_wrong=}\n"
    log_str += f"{count_right=}\n"
    for cluster_id, (cluster1, cluster2) in sep_dict.items():
        algo_clustering[cluster_id] = cluster1
        algo_clustering.append(cluster2)
    if log:
        return algo_clustering, log_str
    return algo_clustering


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


def algo_clustering_to_file_aux(input_path, index_size):
    input_file_name = input_path[input_path.rfind("/"):]
    output_path = "files/minion_idt/50000 strands in size 150 with x2 errors and cluster avg of 25/algo_results" + input_file_name.replace(".txt", "_algo_result.txt")
    # output_path = "files/minion_idt/15000 strands in size 150 with x2 errors and cluster avg of 40/algo_results" + input_file_name.replace(".txt", "_algo_result_shuffled.txt")
    clustering_info = ClusteringInfo(file_path=input_path)
    C_til, bin_sig_arr = hash_based_cluster(clustering_info.reads_err, index_size=index_size)
    print(C_til[0])
    with open(output_path, 'w', newline='\n') as f:
        for cluster in C_til:
            for read_id in cluster:
                f.write(f"{read_id}\n")
            f.write("***\n")
        f.write("bin_sign:\n")
        for i in range(len(bin_sig_arr)):
            f.write(f"{bin_sig_arr[i]}\n")


def algo_clustering_to_file(index_size):
    input_path = "files/minion_idt/50000 strands in size 150 with x2 errors and cluster avg of 25/evyat files/evyat0_index.txt"
    for i in range(5):
        # if i != 4:
        if True:
            curr_input_path = input_path.replace("_index.txt", f"{i}_index.txt")
            algo_clustering_to_file_aux(curr_input_path, index_size)


def file_to_algo_clustering(path):
    clustering = []
    bin_sig_arr = []
    cluster_id = 0
    clustering.append([])
    with open(path, 'r') as f:
        # find the clustering
        line = f.readline().strip()
        while line != "bin_sign:":
            if line[0] != '*':
                clustering[cluster_id].append(int(line))
                line = f.readline()
            else:
                line = f.readline().strip()
                if line != "bin_sign:":  # means new cluster
                    clustering.append([])
                    cluster_id += 1

        # find the bin_signs
        line = f.readline().strip()
        while line:
            bin_sig_arr.append(line)
            line = f.readline().strip()
    return clustering, bin_sig_arr


def from_no_index_to_index_via_indices_file(indices_file_path, input_strands_file_path):
    with open(indices_file_path, 'r', newline='\n') as ind_file:
        index_line = ind_file.readline().rstrip()
        index_size = len(index_line)
        with open(input_strands_file_path, 'r', newline='\n') as input_file:
            input_line = input_file.readline().rstrip()
            output_path = input_strands_file_path.replace('.txt', f'_index_{index_size}.txt')
            with open(output_path, 'w', newline='\n') as out:
                while input_line:
                    output_line = index_line + input_line + '\n'
                    out.write(output_line)
                    index_line = ind_file.readline().rstrip()
                    input_line = input_file.readline().rstrip()


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
