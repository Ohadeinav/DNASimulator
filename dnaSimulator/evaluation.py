import copy
import math

import numpy as np
from clustering import *
import timeit
import time
import random
from datetime import datetime
# random.seed(datetime.now())


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


def test_gen_rand_input():
    res, str0 = gen_rand_input(100, 5, "input/strands_in01.txt")
    for i in range(len(res)):
        print(f'strand {i}:')
        print(res[i])
    print("str:")
    print(str0)


def create_input(file_path, strand_len=150, num_of_strands=1024, index=None):
    strand_list, origin = gen_rand_input(strand_len, num_of_strands)
    index = RegIndex(math.ceil(math.log(num_of_strands, 4)))
    res_str = ''
    if index is not None:
        file_path = file_path.replace('strands_in', f'strands_in_index_{index.len}')
    for i in range(len(strand_list)):
        if index is not None:
            strand_list[i] = index.next() + strand_list[i]
        res_str += strand_list[i] + '\n'
    with open(file_path, 'w', newline='\n') as f:
        f.write(origin)
    if index is not None:
        with open(file_path.replace('.txt', f'_index_{index.len}.txt'), 'w', newline='\n') as f:
            f.write(res_str)


class RegIndex:
    def __init__(self, length=5):
        self.index = 0
        self.len = length
        self.map_dict = {0: "A", 1: "C", 2: "G", 3: "T"}

    def next(self):
        res = decimal_to_dna_str(self.index)
        self.index += 1
        while len(res) != self.len:
            res = 'A' + res
        return res


class IndexStride(RegIndex):
    def __init__(self, stride):
        super.__init__()
        self.stride = stride

    def next(self):
        res = decimal_to_dna_str(self.index)
        self.index += self.stride
        while len(res) != self.len:
            res = 'A' + res
        return res


def test_decimal_to_dna_str():
    for i in range(20):
        print(f'{i} -> {decimal_to_dna_str(i)}')


def test_reg_index():
    index = RegIndex(4)
    for i in range(int(math.pow(4, 4))):
        print(f'{i} -> {index.next()}')


def create_inputs(strand_len=150, num_of_strands=1024):
    origin_file_path = "input/strand_in0.txt"
    for i in range(10):
        file_path = origin_file_path.replace('in0', f'in0{i}')
        create_input(file_path, strand_len, num_of_strands, index=RegIndex())


def test_time_functions():
    s = time.time_ns()
    _s = timeit.default_timer()
    __s = time.perf_counter_ns()
    for i in range(1000):
        pass
    e = time.time_ns()
    _e = timeit.default_timer()
    __e = time.perf_counter_ns()
    print(f'start: {s}, end: {e}, {e - s}')
    print(f'start: {_s}, end: {_e}, {_e - _s}')
    print(f'start: {__s}, end: {__e}, {__e - __s}')


def test_time_and_accuracy_with_index():
    origin_file_path = "files/minion_idt/3000 strands in size 150 with x1.5 errors and cluster avg of 40/evyat0.txt"
    start = time.perf_counter_ns()
    for i in range(0, 10):
        acc_list = []
        acc_index_list = []
        acc_index_regular_list = []
        path_no_index = origin_file_path.replace('.txt', f'{i}.txt')
        path_index = origin_file_path.replace('.txt', f'{i}_index.txt')
        clustering_info_no_index = ClusteringInfo(file_path=path_no_index)
        clustering_info_index = ClusteringInfo(file_path=path_index)
        s1 = time.perf_counter_ns()
        C_til_old_no_index = hash_based_cluster(clustering_info_no_index.reads_err)
        e1 = time.perf_counter_ns()
        elapsed_time_ns1 = e1 - s1
        elapsed_time_sec1 = elapsed_time_ns1 * math.pow(10, -9)
        print(f'file0{i}\nelapsed time: {elapsed_time_sec1:0.4f} sec')
        gammas = [0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
        for gamma in gammas:
            old_acc = calc_acrcy(C_til_old_no_index, clustering_info_no_index.reads_err, clustering_info_no_index.C_dict
                                 , clustering_info_no_index.C_reps, gamma)
            print(f'gamma: {gamma:0.4f}, acc_old: {old_acc:0.6f}')
            acc_list.append(old_acc)

        s2 = time.perf_counter_ns()
        C_til_old_index_regular = hash_based_cluster(clustering_info_index.reads_err)
        e2 = time.perf_counter_ns()
        elapsed_time_ns2 = e2 - s2
        elapsed_time_sec2 = elapsed_time_ns2 * math.pow(10, -9)
        print(f'file0{i}_index_regular\nelapsed time: {elapsed_time_sec2:0.4f} sec')
        for gamma in gammas:
            old_acc = calc_acrcy(C_til_old_index_regular, clustering_info_index.reads_err, clustering_info_index.C_dict
                                 , clustering_info_index.C_reps, gamma)
            acc_index_regular_list.append(old_acc)
            print(f'gamma: {gamma:0.4f}, acc_old: {old_acc:0.6f}')

        s3 = time.perf_counter_ns()
        C_til_old_index = hash_based_cluster(clustering_info_index.reads_err, index_size=5, cond_func=condition4)
        e3 = time.perf_counter_ns()
        elapsed_time_ns3 = e3 - s3
        elapsed_time_sec3 = elapsed_time_ns3 * math.pow(10, -9)
        print(f'file0{i}_index\nelapsed time: {elapsed_time_sec3:0.4f} sec')
        for gamma in gammas:
            old_acc = calc_acrcy(C_til_old_index, clustering_info_index.reads_err, clustering_info_index.C_dict
                                 , clustering_info_index.C_reps, gamma)
            acc_index_list.append(old_acc)
            print(f'gamma: {gamma:0.4f}, acc_old: {old_acc:0.6f}')

        print(f'file0{i}\nsummary:')
        time_improvement_sec_from_no_index = elapsed_time_sec1 - elapsed_time_sec3
        time_improvement_sec_from_index_regular = elapsed_time_sec2 - elapsed_time_sec3
        time_improvement_ns_from_no_index = elapsed_time_ns1 - elapsed_time_ns3
        time_improvement_ns_from_index_regular = elapsed_time_ns2 - elapsed_time_ns3
        improvement_perc1 = (time_improvement_sec_from_no_index / elapsed_time_sec1) * 100
        improvement_perc2 = (time_improvement_sec_from_index_regular / elapsed_time_sec2) * 100
        print(f'elapsed time not indexed     = {elapsed_time_sec1:0.4f} sec\n'
              f'elapsed time indexed regular = {elapsed_time_sec2:0.4f} sec\n'
              f'elapsed time indexed         = {elapsed_time_sec3:0.4f} sec\n'
              f'elapsed time not indexed - elapsed time indexed = {time_improvement_sec_from_no_index:0.4f} sec\n'
              f'elapsed time indexed regular - elapsed time indexed = {time_improvement_sec_from_index_regular:0.4f} sec\n'
              f'improvement percentage with no index = {improvement_perc1:0.4f}\n'
              f'improvement percentage with index regular = {improvement_perc2:0.4f}\n'
              f'gamma not indexed - gamma indexed = {np.array(acc_list) - np.array(acc_index_list)}\n'
              f'gamma indexed regular - gamma indexed = {np.array(acc_index_regular_list) - np.array(acc_index_list)}\n'
              )
    end = time.perf_counter_ns()
    print(f'it took {(end - start)*math.pow(10, -9)} sec to run this')


def algo_clustering_to_file_aux(input_path, index_size):
    input_file_name = input_path[input_path.rfind("/"):]
    output_path = "files/minion_idt/algo_results" + input_file_name.replace(".txt", "_algo_result.txt")
    clustering_info = ClusteringInfo(file_path=input_path)
    C_til = hash_based_cluster(clustering_info.reads_err, index_size=index_size)
    print(C_til[0])
    with open(output_path, 'w', newline='\n') as f:
        for cluster in C_til:
            for read_id in cluster:
                f.write(f"{read_id}\n")
            f.write("***\n")


def algo_clustering_to_file(index_size):
    input_path = "files/minion_idt/3000 strands in size 150 with x2 errors and cluster avg of 40/evyat0_index.txt"
    for i in range(10):
        if i != 4:
            curr_input_path = input_path.replace("_index.txt", f"{i}_index.txt")
            algo_clustering_to_file_aux(curr_input_path, index_size)


def file_to_algo_clustering(path):
    clustering = []
    cluster_id = 0
    clustering.append([])
    with open(path, 'r') as f:
        line = f.readline().strip()
        while line:
            if line[0] != '*':
                clustering[cluster_id].append(int(line))
                line = f.readline()
            else:
                line = f.readline().strip()
                if line:  # means new cluster
                    clustering.append([])
                    cluster_id += 1
    return clustering


def test_file_to_algo_clustering():
    path = "files/minion_idt/algo_results/evyat00_index_algo_result.txt"
    clustering = file_to_algo_clustering(path)
    print(f'{clustering[0]=}')
    print(f'{clustering[-1]=}')


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


def test_stats(unions=False, singletons=False, rebellious_reads=False, summery=True, file_path=None):
    # file_path = "files/minion_idt/3000 strands in size 150 with x2 errors and cluster avg of 40/evyat00_index.txt"
    if file_path is None:
        file_path = "files/minion_idt/3000 strands in size 150 with x2 errors and cluster avg of 40/evyat00_index.txt"
        algo_result_path = "files/minion_idt/algo_results/evyat00_index_algo_result.txt"
    clustering_info = ClusteringInfo(file_path=file_path)
    C_til = file_to_algo_clustering(algo_result_path)
    # C_til1 = handle_singletons_with_index(copy.deepcopy(C_til), clustering_info, index_size=6, threshold=100)
    C_til2 = handle_singletons_with_index_ver2(copy.deepcopy(C_til), clustering_info, index_size=6, threshold=100)
    # C_til3 = handle_singletons_with_index_ver3(copy.deepcopy(C_til), clustering_info, index_size=6, threshold=100)


    stats_ver_0 = stats_to_str_dict(find_clusters_stats(C_til, clustering_info)[0])
    # stats_ver_1 = stats_to_str_dict(find_clusters_stats(C_til1, clustering_info)[0])
    stats_ver_2 = stats_to_str_dict(find_clusters_stats(C_til2, clustering_info)[0])
    # stats_ver_3 = stats_to_str_dict(find_clusters_stats(C_til3, clustering_info)[0])

    if unions:
        print(stats_ver_0['unions'])
    if singletons:
        print(stats_ver_0['singletons'])
    if rebellious_reads:
        print(stats_ver_0['rebellious_reads'])
    if summery:
        acc0 = calc_acrcy(C_til, clustering_info.reads_err, clustering_info.C_dict, clustering_info.C_reps, gamma=0.99)
        # acc1 = calc_acrcy(C_til1, clustering_info.reads_err, clustering_info.C_dict, clustering_info.C_reps, gamma=0.99)
        acc2 = calc_acrcy(C_til2, clustering_info.reads_err, clustering_info.C_dict, clustering_info.C_reps, gamma=0.99)
        # acc3 = calc_acrcy(C_til3, clustering_info.reads_err, clustering_info.C_dict, clustering_info.C_reps, gamma=0.99)
        print(stats_ver_0['summery'])
        print(f'acc_ver_0 = {acc0:0.4f}\n')
        # print(stats_ver_1['summery'])
        # print(f'acc_ver_1 = {acc1:0.4f}\n')
        print(stats_ver_2['summery'])
        print(f'acc_ver_2 = {acc2:0.4f}\n')
        # print(stats_ver_3['summery'])
        # print(f'acc_ver_3 = {acc3:0.4f}\n')


def test_handle_singletons(index_size=6):
    file_path = "files/minion_idt/3000 strands in size 150 with x2 errors and cluster avg of 40/evyat0_index.txt"
    for i in range(10):
        if i == 4:
            continue
        curr_file = file_path.replace("_index.txt", f"{i}_index.txt")
        print(f"file0{i}:")
        test_stats(file_path=curr_file)


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


def test_jaccard():
    s1, _ = gen_rand_input(6, 1)
    s2, _ = gen_rand_input(6, 1)
    # s1 = ['ATAACGAATT']
    # s2 = ['GTAACTAGAT']
    # expected = (2/14) * 100
    res = jaccard(s1[0], s2[0], 2)
    print(f'{s1[0]=}')
    print(f'{s2[0]=}')
    print(f'{res=}')
    # print(f'{expected=}')


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


def main():
    # accuracies_cmp()
    # test_gen_rand_input()
    # test_decimal_to_dna_str()
    # test_reg_index()
    # create_inputs(strand_len=150, num_of_strands=3000)
    # test_time_functions()
    # test_time_and_accuracy_with_index()
    # algo_clustering_to_file(index_size=6)
    # test_file_to_algo_clustering()
    # test_handle_singletons()
    # no_indices_file_path = "input/3000 strands in size 150/strand_in00.txt"
    # indices_file_path = "input/special indices/indices.txt"
    # from_no_index_to_index_via_indices_file(indices_file_path, no_indices_file_path)
    # test_jaccard()
    # test_stats()
    understanding_singletons()


if __name__ == "__main__":
    main()
