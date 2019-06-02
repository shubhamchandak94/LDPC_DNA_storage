import numpy as np
import random
import argparse
import json


def get_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_file', type=str, default="output.txt")
    parser.add_argument('--sample_file', type=str, default="sample.txt")
    parser.add_argument('--coverage', type=float, default=1.0)
    parser.add_argument('--num_chunks', type=int, default=1000)
    parser.add_argument('--chunk_length', type=int, default=256)
    parser.add_argument('--eps', type=float, default=0.0)
    return parser


def randomly_sample_reads(input_file, sample_file, num_sample_reads, chunk_size, eps):

    f_input = open(input_file, "r")
    input_str = f_input.read().rstrip('\n')
    f_input.close()
    f_sample = open(sample_file, "w")
    num_chunks = len(input_str)//chunk_size
    input_data = np.array([int(c) for c in input_str], dtype=int)
    perm = np.random.permutation(len(input_str))  # random permutation to
    # avoid any strangeness due to code structure and read structure
    # this will be inverted later
    invperm = np.zeros(len(input_str), dtype=int)
    invperm[perm] = np.arange(len(input_str))
    input_data = input_data[perm]
    total_counts = np.zeros(len(input_str))
    # number of zeros received at each pos
    zero_counts = np.zeros(len(input_str))
    unique_reads = set([])
    log_1_minus_eps_by_eps = np.log((1-eps)/eps)
    print(num_sample_reads)
    for i in range(num_sample_reads):
        read_id = np.random.randint(num_chunks)
        unique_reads.add(read_id)
        total_counts[read_id*chunk_size:(read_id+1)*chunk_size] += 1
        zero_counts[read_id*chunk_size:(read_id+1)*chunk_size] += (
            ((input_data[read_id*chunk_size:(read_id+1)*chunk_size]+(np.random.rand(chunk_size) < eps)) % 2) == 0)
    print("Number of unique reads", len(unique_reads))
    llr = (2*zero_counts-total_counts)*log_1_minus_eps_by_eps
    llr = llr[invperm]
    for i in range(len(input_str)):
        f_sample.write(str(llr[i])+' ')
    f_sample.close()


def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    num_samples = int(config.coverage*config.num_chunks)
    randomly_sample_reads(config.output_file, config.sample_file,
                          num_samples, config.chunk_length, config.eps)


if __name__ == '__main__':
    main()
