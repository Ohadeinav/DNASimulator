import argparse
import os
import shutil
import threading

from simulator import *
import time
import math

class HashBasedCluster:

    def __init__(self, chosen_technology):
        self.technology = chosen_technology
        if platform.system() == "Linux":
            self.shuffled_file = '/home_nfs/sgamer8/DNAindex' + str(
                self.index) + '/files/' + self.technology + '/' + 'errors_shuffled.txt'
        elif platform.system() == "Windows":
            self.shuffled_file = 'files/' + self.technology + '/' + 'errors_shuffled.txt'

        if platform.system() == "Linux":
            self.evyat_path = '/home_nfs/sgamer8/DNAindex' + str(self.index) + '/files/' + self.technology + '/' + 'evyat.txt'
        elif platform.system() == "Windows":
            self.evyat_path = 'files/' + self.technology + '/' + 'evyat.txt'

    def hash_fun(self, x, a, w, l):
        """ An implementation of the hash function from the article section 3.2"""
        # assume that x has substring a
        # but what happens if it return -1???
        ind = x.find(a)
        return x[ind:min(len(x), ind + w + l)]

    def ind_st(self, st):
        # All 3-grams: AAA,AAG,AAC,AAT,AGA,AGC,...,TTA,TTG,TTC,TTT
        # A-0, G-1, C-2, T-3
        # index of CAT = Decimal representaion of (203)
        # it's actually 2*4^2 + 0*4^1 + 3*4^0 its quad not decimal
        N_q = {"A": 0, "C": 1, "G": 2, "T": 3}
        dec = 0
        for i in range(0, len(st)):
            dec += N_q[st[i]] * (4 ** (len(st) - i - 1))
        return dec

    def bin_sig(self, x, q):
        bs = [0] * (4 ** q)  # if bs[i] = 1 it means that there is at least 1 substring
                             # of x that in quad base means i
                             # for example: if x is CAT... then CAT is substring in x
                             # and it's value is 203(in base 4) = i(in base 10)
        for i in range(0, len(x) - q + 1):
            st = x[i:i + q]
            bs[self.ind_st(st)] = 1
        bs_str = ''.join(str(e) for e in bs)
        return bs_str

    def ham_dis(self, x, y):
        dis = 0
        for i in range(0, len(x)):
            if x[i] != y[i]:
                dis += 1
        return dis

    def rand_perm(self, w):
        return ''.join(random.choice('ACGT') for _ in range(w))

    def rep_find(self, inp, parent):
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
        while parent[temp] != temp and cnt < 10:
            cnt += 1
            temp = parent[temp]
        return temp

    def edit_dis(self, s1, s2):
        """
        Calculate the edit distance between s1 to s2
        """
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

    def display_parent(self, parent):
        """
        | Args:
        |     parent[i] = the minimum id of all the read in the cluster that contains read i
        |     e.g. if read i is in the cluster {k, i, j} where k < i < j
        |     than parent[k], parent[i], parent[j] = k

        | return:
        |    Clustering that is construct from the parent list.
        |    The output is a dict where the keys are the id of a cluster
        |    and the values are the cluster - list of reads id.
        """
        clstr = {}
        for i in range(0, len(parent)):
            clstr[i] = []
        for i in range(0, len(parent)):
            clstr[self.rep_find(i, parent)].append(i)
        return clstr

    def min_max(self, val1, val2):
        """ return the tuple: (min{val1, val2}, max{val1, val2})"""
        min_val = min(val1, val2)
        max_val = max(val1, val2)
        return min_val, max_val

    def hash_based_cluster(self, number_of_steps, windows_size, similarity_threshold, report_func):
        reads_err = []  # will have all the reads from the sequencing stage
        strand_id = 0
        original_strand_dict = {}  # map from orig strand id to the actual strand
        reads_err_original_strand_dict = {}  # map from read_err to it's orig strand id
        temp_evyat_path = self.evyat_path.removesuffix('evyat.txt')
        temp_evyat_path += 'temp_evyat.txt'
        with open(self.evyat_path, 'r') as evyat_f:
            line = "j"
            while line:
                line = evyat_f.readline()
                if not line:
                    break
                original_strand_dict.update({strand_id: line})
                line = evyat_f.readline()  # line == '*****************************\n':
                line = evyat_f.readline()
                while line != '\n':
                    reads_err.append(line.strip())
                    # TODO: check if it should be {len(reads_err)-1: strand_id}
                    # yes it should !!!!
                    reads_err_original_strand_dict.update({len(reads_err)-1: strand_id})
                    line = evyat_f.readline()
                strand_id = strand_id+1
                line = evyat_f.readline()

        read_err_dict = {}  # map from the index of the read in the read_err list to the real itself
        reads_err_ind = [0] * (len(reads_err))
        parent = [0] * (len(reads_err))
        bin_sig_arr = []

        local_comm = number_of_steps
        tmp_bar = int(local_comm / 10)
        total_bar_size = local_comm + 2 * tmp_bar

        q = windows_size
        # computing the binary signature for all the reads
        for i in range(0, len(reads_err)):
            if i % 100 == 0:  # for the progress bar (GUI)
                report_func(total_bar_size, 1 + int((tmp_bar * i) / len(reads_err)))
            read_err_dict.update({i: reads_err[i]})
            reads_err_ind[i] = (i, reads_err[i])
            parent[i] = i
            bin_sig_arr.append(self.bin_sig(reads_err[i], q))

        C_til = self.display_parent(parent)  # short for C^(~) i.e. C tilda

        dist_arr = []
        # th_low, th_high:
        # finding Hamming distance from the first read
        for i in range(1, len(reads_err)):
            if i % 100 == 0:
                report_func(total_bar_size, tmp_bar + int((tmp_bar * i) / len(reads_err)))
            dist_arr.append(self.ham_dis(bin_sig_arr[i], bin_sig_arr[0]))
        # sort the dist array
        dist_arr.sort()
        # the weird definition of theta_low and theta_high
        # I should probably ask why this specific theta's
        for i in range(0, 1000):
            if dist_arr[i + 1] - dist_arr[i] > 10:
                break
        th_low = min(dist_arr[i] + 5, dist_arr[i + 1] - 8)
        th_high = min(dist_arr[i + 1] - 5, th_low + 10)

        r = similarity_threshold
        # set parameters from the article section 5.1:
        w = math.ceil(math.log(len(reads_err[0]), 4))
        l = math.ceil(math.log(len(reads_err), 4))
        for lcl_step in range(0, local_comm):
            report_func(total_bar_size, (2 * tmp_bar) + lcl_step)
            hash_select = [0] * (len(reads_err))
            # picking random representatives samples from each cluster:
            for i in range(0, len(C_til)):
                if len(C_til[i]) != 0:
                    hash_select[random.choice(C_til[i])] = 1

            # computing the hash value for each representative:
            a = self.rand_perm(w)
            hash_C_til = [0] * (len(reads_err))
            for i in range(0, len(reads_err_ind)):
                if hash_select[i] == 0:
                    hash_C_til[i] = (i, "")
                else:
                    hash_C_til[i] = (i, self.hash_fun(reads_err_ind[i][1], a, w, l))
            # sort the hash values by the strings values
            hash_C_til.sort(key=lambda x: x[1])

            cnt = 0
            # implementation of lines 8-10 in the algorithm in the article:
            for i in range(0, len(hash_C_til) - 1):
                if hash_C_til[i][1] == "":
                    continue
                else:
                    if hash_C_til[i][1] == hash_C_til[i + 1][1]:
                        x = reads_err[hash_C_til[i][0]]
                        y = reads_err[hash_C_til[i + 1][0]]
                        if ((self.ham_dis(bin_sig_arr[hash_C_til[i][0]], bin_sig_arr[hash_C_til[i + 1][0]]) <= th_low) or
                                ((self.ham_dis(bin_sig_arr[hash_C_til[i][0]],
                                               bin_sig_arr[hash_C_til[i + 1][0]]) <= th_high) and
                                 self.edit_dis(x, y) <= r)):
                            cnt += 1
                            min_temp, max_temp = self.min_max(self.rep_find(hash_C_til[i][0], parent),
                                                              self.rep_find(hash_C_til[i + 1][0], parent))
                            # merging x and y clusters
                            C_til[min_temp].extend(C_til[max_temp])
                            C_til[max_temp] = []
                            parent[max_temp] = min_temp

        clusters = [sorted(x) for x in list(C_til.values()) if x != []]
        with open(temp_evyat_path, 'w', newline='\n') as temp_f:
            for cluster in clusters:
                orig_strand_candidates = []
                for cluster_element in cluster:
                    orig_strand_candidates.append(reads_err_original_strand_dict.get(cluster_element))
                    # if cluster_element == 0:
                    #     print('wtf')
                    #     if reads_err_original_strand_dict.get(cluster_element) is None:
                    #         print('wtf2')
                orig_strand_id = max(orig_strand_candidates, key= orig_strand_candidates.count)
                temp_f.write(str(original_strand_dict.get(orig_strand_id)))
                temp_f.write('*****************************\n')
                for cluster_element in cluster:
                    # if cluster_element == 0:
                    #     print('wtf')
                    #     if read_err_dict.get(cluster_element) is None:
                    #         print('wtf2')
                    temp_f.write(read_err_dict.get(cluster_element) + '\n')
                temp_f.write('\n\n')

