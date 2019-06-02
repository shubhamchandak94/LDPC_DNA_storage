import numpy as np
import random
import argparse
import json
def get_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_file', type=str, default="output.txt")
    parser.add_argument('--sample_file', type=str, default="sample.txt")
    parser.add_argument('--coverage', type=float, default=1.0)
    parser.add_argument('--num_chunks', type=int, default=100)
    parser.add_argument('--eps', type=float, default=0.0)
    return parser

def randomly_sample_reads(input_file, sample_file,num_sample_reads, eps):
    
	f_input = open(input_file, "r");
	f_sample = open(sample_file, "w");
	input_data = json.load(f_input)
	output_data = dict(input_data)
	output_data['symbols'] = []
	print num_sample_reads
	for i in range(num_sample_reads):
		read = list(random.choice(input_data['symbols']))
		l = [c for c in read[1]]
		for j in range(len(l)):
			if np.random.random() < eps:
				l[j] = '0' if (read[1][j] == '1') else '1'
		read[1] = ''.join(l)
		output_data['symbols'].append(read)
	print("Number of unique reads", len(set([s[0] for s in output_data['symbols']])))
	f_sample.write(json.dumps(output_data, sort_keys = 'False', indent=2, separators=(',', ': ')))

def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    num_samples = int(config.coverage*config.num_chunks)
    randomly_sample_reads( config.output_file, config.sample_file, num_samples, config.eps);

if __name__ == '__main__':
    main()
