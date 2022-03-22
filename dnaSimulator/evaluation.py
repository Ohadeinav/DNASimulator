



def comp_clusters(cluster1, cluster2, gamma):
    """
    compute the indicator condition that is a part of the
    accuracy definition in the article section 2 definition 2.1
    Args:
        cluster1: cluster from the algorithm clustering
        cluster2: cluster from the true clustering
        gamma: threshold should be between 0.5 to 1

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
    if gamma < 0.5 or gamma > 1:
        print(f'gamma can only be value between 0.5 to 1. you gave gamma = {gamma}')
        return -1
    accuracy = 0
    # print(len(true_clustering[0]))
    # print(len(algo_clustering[0]))
    # print(true_clustering[1] == algo_clustering[1])
    # for i in range(len(true_clustering.keys())):
    #     if true_clustering[i] != algo_clustering[i]:
    #         print(f'{len(true_clustering[i])}, {len(algo_clustering[i])}')

    for true_cluster in true_clustering.values():
        for algo_cluster in algo_clustering.values():
            if len(algo_cluster) >= 1:
                res = comp_clusters(algo_cluster, true_cluster, gamma)
                accuracy += res
                if res == 1:
                    break
    return accuracy/len(true_clustering)

def file_to_cluster(file_path):
    reads_err = []  # will have all the reads from the sequencing stage
    strand_id = 0
    cluster = {}
    original_strand_dict = {}  # map from orig strand id to the actual strand
    reads_err_original_strand_dict = {}  # map from read_err to it's orig strand id
    # temp_evyat_path = file_path.removesuffix('evyat.txt')
    # temp_evyat_path += 'temp_evyat.txt'
    with open(file_path, 'r') as evyat_f:
        line = "j"
        while line:
            line = evyat_f.readline()
            if not line:
                break
            original_strand_dict.update({strand_id: line})
            cluster[strand_id] = [line.strip()]
            line = evyat_f.readline()  # line == '*****************************\n':
            line = evyat_f.readline()
            while line != '\n':
                reads_err.append(line.strip())
                cluster[strand_id].append(line.strip())
                reads_err_original_strand_dict.update({len(reads_err): strand_id})
                line = evyat_f.readline()
            strand_id = strand_id + 1
            line = evyat_f.readline()
    return cluster

def main():
    origin = "files/minion_idt/evyat.txt"
    algo_output_path = "files/minion_idt/temp_evyat.txt"
    C = file_to_cluster(origin)
    C_til = file_to_cluster(algo_output_path)
    gamma = 0.99999999999999999999999999999999
    print(calc_accuracy(C, C_til, gamma))


if __name__ == "__main__":
    main()
