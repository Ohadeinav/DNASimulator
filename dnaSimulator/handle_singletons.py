from simulator import *
from metrics import *
from evaluation import *
from gen_input import *
from clustering import *
import itertools

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


def get_all_nearby_str_sub_sing(str_sub_sign, dist = 2):
    substr_list = [str_sub_sign]
    str_as_list = list(str_sub_sign)
    for comb in itertools.combinations_with_replacement(range(len(str_sub_sign)), dist):
        nxt_str = str_as_list
        for idx in comb:
            nxt_str[idx] = str(1 - int(str_as_list[idx]))
        next_str = "".join(nxt_str)
        substr_list.append(next_str)

    return substr_list

def handle_singletons_with_index_ver_6(algo_clustering, orig_cluster_info, bin_sign_arr, index_size, threshold=10, num_epochs=2):
    algo_clustering_copy = copy.deepcopy(algo_clustering)
    for epoch_id in range(num_epochs):
        print(f"epoch: {epoch_id}")
        reads_err = orig_cluster_info.reads_err
        stat1, stat2 = find_clusters_stats(algo_clustering_copy, orig_cluster_info)
        #w = math.ceil(math.log(len(reads_err[0][index_size:]), 4))
        #l = math.ceil(math.log(len(reads_err), 4))-3
        HASH_LENGTH = 10
        strt_index = random.sample(list(range(len(bin_sign_arr[0]) - HASH_LENGTH)),k=1)[0]

        candidates = {}
        hash_dict = {}
        # find singletons
        singletons = []
        for cluster_id, cluster in enumerate(algo_clustering_copy):
            if len(cluster) == 1:
                singleton_id = cluster[0]
                candidates[singleton_id] = {}
                singletons.append((singleton_id, cluster_id))
                # gather hash value of singleton
                #sign_array = np.array(list(bin_sign_arr[singleton_id]))
                str_sub_sign = bin_sign_arr[singleton_id][strt_index:strt_index+HASH_LENGTH]

                #find all nearby strings and add them all to hash dict
                str_sub_sign_list = get_all_nearby_str_sub_sing(str_sub_sign, dist = 2)
                for sub_sign in str_sub_sign_list:
                    if sub_sign in hash_dict:
                        hash_dict[sub_sign]['singletons_ids'].append(singleton_id)
                    else:
                        hash_dict[str_sub_sign] = {'cluster_ids': [],
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

            for rep in representatives:
                #rep_sign_array = np.array(list(bin_sign_arr[rep]))
                hash_value = bin_sign_arr[rep][strt_index:strt_index+HASH_LENGTH]
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
    return [sorted(x) for x in algo_clustering_copy if x != []], ''
