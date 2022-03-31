import math

import numpy as np
from clustering import *
import timeit
import time


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
        C_til_old_index = hash_based_cluster(clustering_info_index.reads_err, index_size=5)
        e3 = time.perf_counter_ns()
        elapsed_time_ns3 = e3 - s3
        elapsed_time_sec3 = elapsed_time_ns3 * math.pow(10, -9)
        print(f'file0{i}_index_regular\nelapsed time: {elapsed_time_sec3:0.4f} sec')
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


def main():
    # accuracies_cmp()
    # test_gen_rand_input()
    # test_decimal_to_dna_str()
    # test_reg_index()
    # create_inputs(strand_len=150, num_of_strands=3000)
    # test_time_functions()
    test_time_and_accuracy_with_index()


if __name__ == "__main__":
    main()
