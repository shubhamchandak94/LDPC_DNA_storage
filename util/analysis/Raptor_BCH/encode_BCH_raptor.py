import numpy as np
import random
import argparse
import json
import subprocess
import sys
import bchlib
import os
import binascii
import base64

BCH_POLYNOMIAL = 285

def get_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--output_file', type=str)
    parser.add_argument('--BCH_bits', type=int)	
    parser.add_argument('--alpha_raptor', type=float)	
    return parser

def encode(input_file, output_file, BCH_bits, alpha_raptor):
	encode_raptor_script = "../python-libraptorq/rq --debug encode "   
	raptor_length = (31-BCH_bits) - (31-BCH_bits)%4 # RaptorQ expects multiples of 4
	arg_string  = " -s " + str(raptor_length)
	arg_string += " -m " + str(1000000)
	arg_string += " --repair-symbols-rate " + str(alpha_raptor)
	arg_string += " "
	arg_string += input_file + " "
	arg_string += "tmpfile"
	encode_raptor_command = encode_raptor_script + arg_string
        print(encode_raptor_command)
	subprocess.call([encode_raptor_command], shell=True)
        assert os.path.isfile("tmpfile"),"The codebook did not get generated"

	bch = bchlib.BCH(BCH_POLYNOMIAL, BCH_bits)
	f_input = open("tmpfile", "r");
	data = json.load(f_input)
	f_input.close()
	for i,s in enumerate(data['symbols']):
                s_byte = base64.urlsafe_b64decode(str(s[1]))
		s_byte_coded = s_byte + bch.encode(s_byte)
		data['symbols'][i][1] = bin(int(binascii.hexlify(s_byte_coded), 16))[2:].zfill(len(s_byte_coded)*8)
	
	data['BCH_bits'] = BCH_bits	
	f_out = open(output_file,'w')
	data = json.dumps(data, sort_keys=True, indent=2, separators=(',', ': '))
	f_out.write(data)
	f_out.close()	

def main():
    parser = get_argument_parser()
    config = parser.parse_args()
    encode( config.input_file, config.output_file, config.BCH_bits, config.alpha_raptor);

if __name__ == '__main__': 
	main()
