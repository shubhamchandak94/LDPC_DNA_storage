import sys
sys.path.append('..')
import dna_storage_genap as dna_storage
import filecmp
import os
import random
from params import *

num_trials = 20 # make sure decoding succeeds for so many random trials

read_file_prefix='/raid/nanopore/shubham/LDPC_DNA_storage_data_genapsys_20201109/Stanford/'
read_file=['S000715_RUN_2020_11_02_11_58_52_sample_005_aligned_'+str(i)+'.reads' for i in range(1,num_files+1)]

num_reads_init = [
7000,
7000,
7000,
]        

eps = [0.04]*3

random.seed(42)

tmpfile_sampled_reads = 'tmpfile_sampled_reads'
tmpfile_barcode_removed = 'tmpfile_barcode_removed'
tmpfile_decoded = 'tmpfile_decoded'

for i in range(num_files):
    index_len = (bin_index_len[i]+5*BCH_bits[i])//2
    sync_pos = (oligo_len-index_len)//2+index_len
    print('sync_pos',sync_pos)
    file_size = os.path.getsize(infile_path_prefix+infile_name[i])
    print('file_size (in bits):', file_size*8)
    num_reads = num_reads_init[i]
    # load all reads into list
    with open(read_file_prefix+read_file[i]) as f:
        all_reads = f.readlines()
    while True:
        print('num_reads:',num_reads)
        print('reading cost (bases/bit):', (num_reads*oligo_len)/(file_size*8))
        success = True
        for t in range(num_trials):
            sampled_reads = random.sample(all_reads,num_reads)
            with open(tmpfile_sampled_reads,'w') as f:
                for sampled_read in sampled_reads:
                    f.write(sampled_read)
            dna_storage.remove_barcodes_flexbar(tmpfile_sampled_reads,start_barcodes[i], end_barcodes[i], tmpfile_barcode_removed) 
            dna_storage.decode_data(tmpfile_barcode_removed,oligo_len,tmpfile_decoded,bin_index_len[i],BCH_bits[i],LDPC_alpha[i],LDPC_path_prefix+LDPC_code[i],file_size,eps=eps[i],sync=sync[i],sync_pos=sync_pos)
            if not filecmp.cmp(infile_path_prefix+infile_name[i],tmpfile_decoded):
                success = False
            os.remove(tmpfile_sampled_reads)
            os.remove(tmpfile_decoded)
            os.remove(tmpfile_barcode_removed)
            if not success:
                break
        if success:
            break
        num_reads += 500
