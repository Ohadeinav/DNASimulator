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
    output_path = "files/minion_idt/15000 strands in size 150 with x2 errors and cluster avg of 40/algo_results" + input_file_name.replace(".txt", "_algo_result.txt")
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
    input_path = "files/minion_idt/15000 strands in size 150 with x2 errors and cluster avg of 40/evyat files/evyat0_index.txt"
    for i in range(5):
        if i != 4:
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
