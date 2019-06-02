import numpy as np
import random
import argparse
import json
import subprocess
import sys
import bchlib
import binascii
import base64
BCH_POLYNOMIAL = 285


def get_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--recon_file', type=str)
    return parser

def decode(input_file, recon_file):
	decode_raptor_script = "../python-libraptorq/rq --debug decode "    
	f_input = open(input_file, "r");
	data = json.load(f_input)
	f_input.close()
	bch = bchlib.BCH(BCH_POLYNOMIAL, data['BCH_bits'])
	reads_dict = {}
        num_bytes_per_read = len(data['symbols'][0][1])//8 
	# load reads and store according to index
	for read in data['symbols']:
		if read[0] in reads_dict:
			reads_dict[read[0]].append([int(c) for c in read[1]])
		else:
			reads_dict[read[0]] = [[int(c) for c in read[1]]]
	corrected_reads = [] #to store index and corrected bitstring
	bitflips_dict = {} #to store the number of bitflips for each index (for sorting later)
	# convert to numpy arrays, find consensus, do BCH decoding and keep if decoding successful
	for k in reads_dict.keys():
		read_array = np.array(reads_dict[k], dtype = int)
		majority_list = np.array(np.mean(read_array,axis = 0)>0.5,dtype=int)
		majority_str = ''.join([str(c) for c in majority_list])
		majority_str_bytes = binascii.unhexlify(((hex(int(majority_str,2)))[2:-1]).zfill(2*(num_bytes_per_read)))
		maj_data, maj_ecc = majority_str_bytes[:-bch.ecc_bytes], majority_str_bytes[-bch.ecc_bytes:]
		(bitflips, maj_data, maj_ecc) = bch.decode(maj_data, maj_ecc)
		if bitflips >= 0: #success
			corrected_str = base64.urlsafe_b64encode(maj_data)
			corrected_reads.append([k, corrected_str])
			bitflips_dict[k] = bitflips
	corrected_reads.sort(key=lambda x: bitflips_dict[x[0]])
	output_data = dict(data)
	output_data['symbols'] = corrected_reads
	f_intermediate = open("tmpfile", "w");
	f_intermediate.write(json.dumps(output_data, sort_keys = 'False', indent=2, separators=(',', ': ')))
	f_intermediate.close()
	ret = subprocess.call([decode_raptor_script+" tmpfile "+ recon_file], shell=True)
	return ret

def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    ret =  decode( config.input_file, config.recon_file);
    return ret

if __name__ == '__main__': sys.exit(main())
