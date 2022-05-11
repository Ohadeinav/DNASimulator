import copy
import math

import numpy as np
from clustering import *
# from clustering import ClusteringInfo
from metrics import *
import timeit
import time
import random
from datetime import datetime
# random.seed(datetime.now())


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

    # for index in range(len(algo_clustering)):
    #     if len(algo_clustering[index]) == 0:
    #         algo_clustering.pop(index)
    [cluster for cluster in algo_clustering if len(cluster) != 0]
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


def stats_to_str_dict(stats1):
    str_summery = 'summery:\n'

    # union:
    unwanted_unions = find_unwanted_unions(stats1)
    str_union = ''
    for index, algo_cluster_stat in unwanted_unions.items():
        str_union += f'algo_cluster_index = {index}:\n'
        for orig_cluster_id, stats in algo_cluster_stat.items():
            str_union += (f'orig_cluster_id = {orig_cluster_id} ; size_in_algo_cluster = {stats[1]:0.1f} '
                          f'; size_of_orig_cluster = {stats[2]:0.1f} ; percentage = {stats[0]:0.4f}\n')
        str_union += '\n'
    str_summery_temp = f'The total number of unwanted unions is {len(unwanted_unions.keys())}\n'
    str_union += str_summery_temp
    str_summery += str_summery_temp

    # singletons:
    unwanted_singletons = find_unwanted_singletons(stats1)
    str_singletons = ''
    for index, algo_cluster_stat in unwanted_singletons.items():
        str_singletons += f'algo_cluster_index = {index}:\n'
        for orig_cluster_id, stats in algo_cluster_stat.items():
            str_singletons += f'orig_cluster_id = {orig_cluster_id} ; size_in_algo_cluster = {stats[1]:0.1f} ' \
                              f'; size_of_orig_cluster = {stats[2]:0.1f} ; percentage = {stats[0]:0.4f}\n'
        str_singletons += '\n'
    str_summery_temp = f'The total number of unwanted singletons is {len(unwanted_singletons.keys())}\n'
    str_singletons += str_summery_temp
    str_summery += str_summery_temp

    # rebellious_reads
    sum_of_rebellious_reads = 0
    str_rebellious_reads = ''
    unwanted_rebellious_reads = find_unwanted_rebellious_reads(stats1)
    for index, algo_cluster_stat in unwanted_rebellious_reads.items():
        str_rebellious_reads += f'algo_cluster_index = {index}:\n'
        for orig_cluster_id, stats in algo_cluster_stat.items():
            str_rebellious_reads += f'orig_cluster_id = {orig_cluster_id} ; size_in_algo_cluster = {stats[1]:0.1f} ' \
                                    f'; size_of_orig_cluster = {stats[2]:0.1f} ; percentage = {stats[0]:0.4f}\n'
            sum_of_rebellious_reads += stats[1]
        str_rebellious_reads += '\n'
    num_of_clusters_with_rebellious_reads = len(unwanted_rebellious_reads.keys())
    avg_of_rebellious_reads = sum_of_rebellious_reads / num_of_clusters_with_rebellious_reads
    str_summery_temp = f'The total number of clusters with unwanted rebellious reads is {num_of_clusters_with_rebellious_reads}\n' \
                       f'The avg size of unwanted rebellious reads in a cluster is {avg_of_rebellious_reads:0.4f}\n'
    str_rebellious_reads += str_summery_temp
    str_summery += str_summery_temp
    return {'summery': str_summery, 'unions': str_union, 'singletons': str_singletons,
            'rebellious_reads': str_rebellious_reads}


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


def understanding_singletons():
    file_path = "files/minion_idt/3000 strands in size 150 with x2 errors and cluster avg of 40/evyat00_index.txt"
    result_path = "files/minion_idt/algo_results/evyat00_index_algo_result.txt"
    clustering_info = ClusteringInfo(file_path=file_path)
    C_til = file_to_algo_clustering(result_path)
    stats1, stats2 = find_clusters_stats(C_til, clustering_info)
    unwanted_singletons = find_unwanted_singletons(stats1)
    count = 0
    index_ed_dict = Counter()
    jaccard_dis_dict = Counter()
    gpm_dict = Counter()
    ham_dis_dict = Counter()
    edit_dis_dict = Counter()
    origins_dict = Counter()
    num_of_index_guessed_wrong = 0
    sum_of_representatives_guessed_wrong = 0
    max_cluster_size_of_wrong_guessed_index = 100
    for algo_cluster_id, value in unwanted_singletons.items():
        for origin_cluster_id in value.keys():
            singleton_id = C_til[algo_cluster_id][0]
            singleton = clustering_info.reads_err[singleton_id]
            singleton_index = singleton[:6]
            singleton_origin_id = clustering_info.reads_err_original_strand_dict[singleton_id]
            origins_dict.update([singleton_origin_id])
            singleton_origin = clustering_info.original_strand_dict[singleton_origin_id]
            singleton_origin_index = singleton_origin[:6]
            if singleton_origin_index == singleton_index:
                count += 1
            jaccard_dist = float("{:.3f}".format(jaccard(singleton_origin_index, singleton_index, 2)))
            gpm_dist = float("{:.3f}".format(GPM_quick_ratio(singleton_origin_index, singleton_index)))
            ham_dist = ham_dis(bin_sig(singleton_origin_index, 2), bin_sig(singleton_index, 2))
            edit_dist = edit_dis(singleton_origin_index, singleton_index)
            jaccard_dis_dict.update([jaccard_dist])
            gpm_dict.update([gpm_dist])
            ham_dis_dict.update([ham_dist])
            edit_dis_dict.update([edit_dist])

            singleton_family = stats2[origin_cluster_id]
            algo_family_id = None
            for algo_family, ststs in singleton_family.items():
                if ststs[0] > 0.5:
                    algo_family_id = algo_family
            if algo_family_id is not None:
                representatives = random.sample(C_til[algo_family_id], min([100, len(C_til[algo_family_id])]))
                index_mat = np.array([list(clustering_info.reads_err[read][:6]) for read in representatives]).T
                true_index = []
                for row in index_mat:
                    row_list = list(row)
                    true_index.append(max(row_list, key=row_list.count))
                str_index = ''.join(true_index)
                if str_index != singleton_origin_index:
                    num_of_index_guessed_wrong += 1
                    sum_of_representatives_guessed_wrong += min([100, len(C_til[algo_family_id])])
                    max_cluster_size_of_wrong_guessed_index = len(C_til[algo_family_id])
                    index_ed_dict.update([edit_dis(str_index, singleton_origin_index)])
                print(f'real index = {singleton_origin_index :8} ;   len_algo = {len(C_til[algo_family_id])} guess index = {str_index:8} ;   singleton index = {singleton_index :8}')
            print(f'real index = {singleton_origin_index :8} ;   singleton index = {singleton_index :8}')
    print(f'there are {count} singletons with true index from {len(unwanted_singletons)} singletons')
    print(f'jaccard_dis_dict = {sorted(jaccard_dis_dict.items(), key=lambda pair: pair[0], reverse=True)}')
    print(f'gpm_dict = {sorted(gpm_dict.items(), key=lambda pair: pair[0], reverse=True)}')
    print(f'ham_dis_dict = {sorted(ham_dis_dict.items(), key=lambda pair: pair[0], reverse=False)}')
    print(f'edit_dis_dict = {sorted(edit_dis_dict.items(), key=lambda pair: pair[0], reverse=False)}')
    print(f'origins_dict = {[pair for pair in sorted(origins_dict.items(), key=lambda pair: pair[0], reverse=False) if pair[1]>=2]}')
    print(f'avg cluster size of wrong guessed index {sum_of_representatives_guessed_wrong/num_of_index_guessed_wrong:0.3f}')
    print(f'number of wrong guessed index {num_of_index_guessed_wrong}')
    print(f'max cluster size of wrong guessed index {max_cluster_size_of_wrong_guessed_index}')
    print(f'index_ed_dict = {sorted(index_ed_dict.items(), key=lambda pair: pair[0], reverse=False)}')


def understanding_unions():
    file_path = "files/minion_idt/3000 strands in size 150 with x2 errors and cluster avg of 40/evyat00_index.txt"
    result_path = "files/minion_idt/algo_results/evyat00_index_algo_result.txt"
    clustering_info = ClusteringInfo(file_path=file_path)
    reads_err = clustering_info.reads_err
    index_size = 6
    C_til, bin_sig_arr = file_to_algo_clustering(result_path)
    stats1, stats2 = find_clusters_stats(C_til, clustering_info)
    unwanted_unions = find_unwanted_unions(stats1)
    edit_dist_counter = Counter()
    ham_dist_dict = {}
    reps_dict = {}
    diffs_list = []
    for algo_cluster_id, cluster_stats in unwanted_unions.items():
        reps_dict[algo_cluster_id] = {}
        ham_dist_dict[algo_cluster_id] = {}
        algo_cluster_splited = {}
        guessed_indexes = {}
        for origin_strand_id, stats in cluster_stats.items():
            if stats[0] >= 0.5:
                algo_cluster_splited[origin_strand_id] = []
                for read_id in C_til[algo_cluster_id]:
                    if clustering_info.reads_err_original_strand_dict[read_id] == origin_strand_id:
                        algo_cluster_splited[origin_strand_id].append(read_id)
        for origin_strand_id, cluster in algo_cluster_splited.items():
            true_index = []
            reps = random.sample(cluster, min([10, len(cluster)]))
            reps_dict[algo_cluster_id][origin_strand_id] = reps
            index_mat = np.array([list(reads_err[read_id][:index_size]) for read_id in reps]).T
            for row in index_mat:
                row_list = list(row)
                true_index.append(max(row_list, key=row_list.count))
            str_index = ''.join(true_index)
            ham = 0
            if len(reps) > 1:
                for i in range(len(reps)):
                    for j in range(len(reps)):
                        if i < j:
                            ham += ham_dis(bin_sig_arr[reps[i]], bin_sig_arr[reps[j]])
                ham_dist_dict[algo_cluster_id][origin_strand_id] = ham/((len(reps)**2 - len(reps))/2)
            guessed_indexes[origin_strand_id] = str_index

        reps_list = list(reps_dict[algo_cluster_id].values())
        reps1 = reps_list[0]
        reps2 = reps_list[1]
        ham = 0
        for rep1 in reps1:
            for rep2 in reps2:
                ham += ham_dis(bin_sig_arr[rep1], bin_sig_arr[rep2])
        diff = ham / (len(reps1)*len(reps2)) - max(ham_dist_dict[algo_cluster_id].values())
        ham_dist_dict[algo_cluster_id][-1] = ham / (len(reps1)*len(reps2))
        ham_dist_dict[algo_cluster_id][-2] = diff
        diffs_list.append(diff)
        print(f"{ham_dist_dict[algo_cluster_id]=}")
        print(f"{guessed_indexes}")
        str_indexes = list(guessed_indexes.values())
        edit_dist_counter.update([edit_dis(str_indexes[0], str_indexes[1])])
        print(f"{edit_dis(str_indexes[0], str_indexes[1])}")
        print(len(C_til[algo_cluster_id]))
        # print(f"{ed_matrix}\n")
    print(sorted(edit_dist_counter.items(), reverse=False, key=lambda x: x[0]))
    print(f"{min(diffs_list)=}")

    handle_unions(C_til, clustering_info, bin_sig_arr, index_size=6, threshold=10)


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


#######################################
""" Our version of accuracy """
#######################################


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


def accuracies_cmp():
    origin_path = "files/minion_idt/evyat"
    # for i in range(0, 10):
    #     path = origin_path + str(i) + ".txt"
    #     print(f'file evyat0{i}:')
    #     accuracies_cmp_aux(path)
    accuracies_cmp_aux(origin_path + '.txt')


def accuracies_cmp_aux(file_path):
    clustering_info = ClusteringInfo(file_path=file_path)
    C = file_to_cluster(file_path)
    C_til_new = file_to_cluster(file_path.replace('evyat', 'temp_evyat'))
    C_til_old = hash_based_cluster(clustering_info.reads_err)
    # C_til_old = C_til_new

    for gamma in np.arange(50, 100, 5):
        new_acc = calc_accuracy(C, C_til_new, gamma/100)
        old_acc = calc_acrcy(C_til_old, clustering_info.reads_err, clustering_info.C_dict, clustering_info.C_reps, gamma/100)
        print(f'gamma: {gamma/100}, acc_new: {new_acc}, acc_old: {old_acc}')
