import numpy as np
from clustering import *
from timeit import timeit


def accuracies_cmp():
    origin_path = "files/minion_idt/evyat0"
    for i in range(0, 1):
        path = origin_path + str(i) + ".txt"
        print(f'file evyat0{i}:')
        accuracies_cmp_aux(path)


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
    res, str0 = gen_rand_input(100, 5, "input/strand_in01.txt")
    for i in range(len(res)):
        print(f'strand {i}:')
        print(res[i])
    print("str:")
    print(str0)


def main():
    # accuracies_cmp()
    test_gen_rand_input()


if __name__ == "__main__":
    main()
