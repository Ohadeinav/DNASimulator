import random
import math


def rand_perm(w):
    """ return random DNA strand of length w"""
    return ''.join(random.choice('ACGT') for _ in range(w))


def gen_rand_input(strand_len, num_of_strands, file_path=None):
    """
    Generate num_of_strands random strands each in strand_len length.
    Returns (list of this strands, string of all the strands).
    optionally also write is to a given file.
    """
    strands = []
    res_str = ''
    for i in range(0, num_of_strands):
        strand = rand_perm(strand_len)
        res_str += strand + "\n"
        strands.append(strand)
    if file_path is not None:
        with open(file_path, 'w', newline='\n') as f:
            f.write(res_str)
    return strands, res_str


def decimal_to_dna_str(num):
    """
    return the value of a decimal number in base 4 in which:
            A = 0 | C = 1 | G = 2 | T = 3
    """
    res = ''
    map_dict = {0: "A", 1: "C", 2: "G", 3: "T"}
    i = num
    while i != 0:
        res = map_dict[i % 4] + res
        i = math.floor(i / 4)
    return res


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


def create_inputs(strand_len=150, num_of_strands=1024):
    origin_file_path = "input/strand_in0.txt"
    for i in range(10):
        file_path = origin_file_path.replace('in0', f'in0{i}')
        create_input(file_path, strand_len, num_of_strands, index=RegIndex())