from gen_input import *
from clustering import *
from evaluation import *
import timeit
import time


def test_jaccard():
    """
    test the Jaccard distance function.
    """
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


def test_gen_rand_input():
    """
    test gen_rand_input function
    """
    res, str0 = gen_rand_input(100, 5, "input/strands_in01.txt")
    for i in range(len(res)):
        print(f'strand {i}:')
        print(res[i])
    print("str:")
    print(str0)


def test_decimal_to_dna_str():
    """
    test decimal_to_dna_str function
    """
    for i in range(20):
        print(f'{i} -> {decimal_to_dna_str(i)}')


def test_reg_index():
    """
    test RegIndex class
    """
    index = RegIndex(4)
    for i in range(int(math.pow(4, 4))):
        print(f'{i} -> {index.next()}')


def test_time_functions():
    """
    exploring different types of time libs
    """
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


def time_function(func, **kwargs):
    """
    :return: A dict containing the elapsed time
        it took to run the given function
        in different unit of measurement: ns, sec, minutes.
    :param func: The function which we want to time.
    :param kwargs: Arguments to pass to the func param.
    """
    start = time.perf_counter_ns()
    func_result = func(**kwargs)
    end = time.perf_counter_ns()
    elapsed_in_ns = end - start
    return {'ns': elapsed_in_ns,
            'sec': elapsed_in_ns * math.pow(10, -9),
            'min': (elapsed_in_ns * math.pow(10, -9))/60}, func_result


def test_time_and_accuracy_with_index():
    """
    measures the time and accuracies of some clustering for different values of gamma.
    """
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


def test_file_to_algo_clustering():
    """
    test the file_to_algo_clustering function.
    """
    path = "files/minion_idt/algo_results/evyat00_index_algo_result.txt"
    clustering, bin_sig_arr = file_to_algo_clustering(path)
    print(f'{clustering[0]=}')
    print(f'{bin_sig_arr[-1]=}')


def test_stats(unions=False, singletons=False, rebellious_reads=False, summery=True,
               file_path=None, algo_result_path=None, index_size=0, gamma=0.99):
    """
    get stats and accuracy for algo clustering, after unions handling and after singletons handling.

    :return: str log containing stats.
    :param unions: Whether or not to return the unions full stats. Default False.
    :param singletons: Whether or not to return the singletons full stats. Default False.
    :param rebellious_reads: Whether or not to return the rebellious_reads full stats. Default False.
    :param summery: Whether or not to return the summery stats. Default True.
    :param file_path: The path to an evyat file. Default to None.
    :param algo_result_path: The path to an clustering file. Default to None.
    :param index_size: The number of symbols dedicate for indexing. Default to 0.
    :param gamma: The gamma parameter to pass to the accuracy function. Default to 0.99.
    """
    # file_path = "files/minion_idt/3000 strands in size 150 with x2 errors and cluster avg of 40/evyat00_index.txt"
    if file_path is None:
        file_path = "files/minion_idt/50000 strands in size 150 with x2 errors and cluster avg of 25/evyat files/evyat00_index.txt"
        algo_result_path = "files/minion_idt/50000 strands in size 150 with x2 errors and cluster avg of 25/algo_results/evyat00_index_algo_result.txt"
    clustering_info = ClusteringInfo(file_path=file_path)
    C_til, bin_sig_arr = file_to_algo_clustering(algo_result_path)
    # C_til1 = handle_singletons_with_index(copy.deepcopy(C_til), clustering_info, index_size=index_size, threshold=100)
    C_til1_5, unions_log = handle_unions(copy.deepcopy(C_til), clustering_info, bin_sig_arr, index_size=index_size, threshold=10)
    # C_til2, singletons_log = handle_singletons_with_index_ver2_5(copy.deepcopy(C_til1_5), clustering_info, bin_sig_arr, index_size=index_size, threshold=12)
    # C_til3 = handle_singletons_with_index_ver3(copy.deepcopy(C_til), clustering_info, index_size=index_size, threshold=100)
    C_til4, singletons_log = handle_singletons_with_index_ver5_5(copy.deepcopy(C_til1_5), clustering_info, bin_sig_arr,
                                                                 index_size=index_size, threshold=10, num_epochs=4,
                                                                 converge=True)

    str_log = f'{unions_log}\n{singletons_log}\n'

    stats_ver_0 = stats_to_str_dict(find_clusters_stats(C_til, clustering_info)[0])
    # stats_ver_1 = stats_to_str_dict(find_clusters_stats(C_til1, clustering_info)[0])
    stats_ver_1_5 = stats_to_str_dict(find_clusters_stats(C_til1_5, clustering_info)[0])
    # stats_ver_2 = stats_to_str_dict(find_clusters_stats(C_til2, clustering_info)[0])
    # stats_ver_3 = stats_to_str_dict(find_clusters_stats(C_til3, clustering_info)[0])
    stats_ver_4 = stats_to_str_dict(find_clusters_stats(C_til4, clustering_info)[0])

    if unions:
        str_log += f"{stats_ver_0['unions']}\n"
        # print(stats_ver_0['unions'])
    if singletons:
        str_log += f"{stats_ver_0['singletons']}\n"
        # print(stats_ver_0['singletons'])
    if rebellious_reads:
        str_log += f"{stats_ver_0['rebellious_reads']}\n"
        # print(stats_ver_0['rebellious_reads'])
    if summery:
        acc0 = calc_acrcy(C_til, clustering_info.reads_err, clustering_info.C_dict, clustering_info.C_reps, gamma=gamma)
        # acc1 = calc_acrcy(C_til1, clustering_info.reads_err, clustering_info.C_dict, clustering_info.C_reps, gamma=gamma)
        acc1_5 = calc_acrcy(C_til1_5, clustering_info.reads_err, clustering_info.C_dict, clustering_info.C_reps, gamma=gamma)
        # acc2 = calc_acrcy(C_til2, clustering_info.reads_err, clustering_info.C_dict, clustering_info.C_reps, gamma=gamma)
        # acc3 = calc_acrcy(C_til3, clustering_info.reads_err, clustering_info.C_dict, clustering_info.C_reps, gamma=gamma)
        acc4 = calc_acrcy(C_til4, clustering_info.reads_err, clustering_info.C_dict, clustering_info.C_reps, gamma=gamma)

        str_log += f"{stats_ver_0['summery']}\n"
        str_log += f'acc_ver_0 = {acc0:0.4f}\n\n'
        # str_log += f"{stats_ver_1['summery']}\n"
        # str_log += f'acc_ver_1 = {acc1:0.4f}\n\n'
        str_log += f"{stats_ver_1_5['summery']}\n"
        str_log += f'acc_ver_1_5 = {acc1_5:0.4f}\n\n'
        # str_log += f"{stats_ver_2['summery']}\n"
        # str_log += f'acc_ver_2 = {acc2:0.4f}\n\n'
        # str_log += f"{stats_ver_3['summery']}\n"
        # str_log += f'acc_ver_3 = {acc3:0.4f}\n\n'
        str_log += f"{stats_ver_4['summery']}\n"
        str_log += f'acc_ver_4 = {acc4:0.4f}\n\n'
    return str_log


def test_handle_singletons(index_size=6, log=True, to_print=True):
    """
    running the test_stats function on multiple files.

    :param index_size: The number of symbols dedicate for indexing. Default to 6.
    :param log: Whether or not to save the log in files. Default to True.
    :param to_print: Whether or not to print the log to screen. Default to True.
    """
    file_path = "files/minion_idt/50000 strands in size 150 with x2 errors and cluster avg of 25/evyat files/evyat0_index.txt"
    algo_result_path = "files/minion_idt/50000 strands in size 150 with x2 errors and cluster avg of 25/algo_results/evyat0_index_algo_result.txt"
    log_path = "files/minion_idt/50000 strands in size 150 with x2 errors and cluster avg of 25/stats files/stats0.txt"
    for i in range(5):
        # if i == 4:
        #     continue
        curr_file = file_path.replace("_index.txt", f"{i}_index.txt")
        curr_algo_result_path = algo_result_path.replace("_index_algo_result.txt", f"{i}_index_algo_result.txt")
        # curr_algo_result_path = algo_result_path.replace("_index_algo_result.txt", f"{i}_index_algo_result_shuffled.txt")
        curr_log_path = log_path.replace(".txt", f"{i}.txt")
        print(f"file0{i}:")
        res = test_stats(file_path=curr_file, algo_result_path=curr_algo_result_path, index_size=index_size)
        if log:
            with open(curr_log_path, 'w', newline='\n') as log_file:
                log_file.write(res)
        if to_print:
            print(res)


def test_times(log=True, log_each_file=True, index_size=6):
    """
    Check times for pipe line of functions.

    :return: The time stats string.
    :param log: Whether or not to save all the logs in a single file. Default to True.
    :param log_each_file: Whether or not to save the log for each file in a different file.
        Default to True.
    :param index_size: The number of symbols dedicate for indexing. Default to 6.
    """
    file_path = "files/minion_idt/50000 strands in size 150 with x2 errors and cluster avg of 25/evyat files/evyat0_index.txt"
    algo_result_path = "files/minion_idt/50000 strands in size 150 with x2 errors and cluster avg of 25/algo_results/evyat0_index_algo_result.txt"
    log_path_each_file = "files/minion_idt/50000 strands in size 150 with x2 errors and cluster avg of 25/stats files/times0.txt"
    log_path = "files/minion_idt/50000 strands in size 150 with x2 errors and cluster avg of 25/time stats_ver5_5.txt"

    # functions_to_check = [hash_based_cluster, handle_unions, handle_singletons_with_index_ver2_5_clean]
    functions_to_check = [hash_based_cluster, handle_unions, handle_singletons_with_index_ver5_5]

    time_stats = ''
    for i in range(5):
        # if i != 4 and i != 6:
        if True:
            print(f'file0{i}')
            time_stats_each_file = f'file0{i}\n'
            curr_file_path = file_path.replace("_index.txt", f"{i}_index.txt")
            curr_algo_result_path = algo_result_path.replace("_index_algo_result.txt", f"{i}_index_algo_result.txt")
            curr_log_path = log_path_each_file.replace(".txt", f"{i}.txt")
            clustering_info = ClusteringInfo(file_path=curr_file_path)
            C_til, bin_sig_arr = file_to_algo_clustering(curr_algo_result_path)

            for func in functions_to_check:
                kw = {}
                if func.__name__ == 'hash_based_cluster':
                    kw = {'reads_err': clustering_info.reads_err, 'index_size': index_size}
                elif func.__name__ == 'handle_unions':
                    kw = {'algo_clustering': C_til,
                          'orig_cluster_info': clustering_info,
                          'bin_sign_arr': bin_sig_arr,
                          'index_size': index_size,
                          'threshold': 10, 'log': False}
                elif func.__name__.startswith('handle_singletons_with_index'):
                    kw = {'algo_clustering': C_til,
                          'orig_cluster_info': clustering_info,
                          'bin_sign_arr': bin_sig_arr,
                          'index_size': index_size,
                          'threshold': 10}
                else:
                    print(f"Error: unknown function")

                time_dict, func_result = time_function(func, **kw)
                if func.__name__ != 'hash_based_cluster':
                    C_til = func_result

                time_stats_each_file += f'{func.__name__}:\n' \
                                        f'time in ns: {time_dict["ns"]: 0.3f}\n' \
                                        f'time in sec: {time_dict["sec"]: 0.3f}\n' \
                                        f'time in min: {time_dict["min"]: 0.3f}\n'
                time_stats_each_file += '\n'

            if log_each_file:
                with open(curr_log_path, 'w', newline='\n') as f:
                    f.write(time_stats_each_file)
            time_stats += time_stats_each_file
            time_stats += '\n'

    if log:
        with open(log_path, 'w', newline='\n') as f:
            f.write(time_stats)
    return time_stats


def main():
    # accuracies_cmp()
    # test_gen_rand_input()
    # test_decimal_to_dna_str()
    # test_reg_index()
    # create_inputs(strand_len=150, num_of_strands=9000)
    # test_time_functions()
    # test_time_and_accuracy_with_index()
    algo_clustering_to_file(index_size=9)
    # test_file_to_algo_clustering()
    # no_indices_file_path = "input/3000 strands in size 150/strand_in00.txt"
    # indices_file_path = "input/special indices/indices.txt"
    # from_no_index_to_index_via_indices_file(indices_file_path, no_indices_file_path)
    # test_jaccard()
    # print(test_stats())
    # understanding_singletons()
    # understanding_unions()
    # test_handle_singletons(index_size=8, log=True)
    # test_times(True, True, index_size=8)


if __name__ == "__main__":
    main()
