from params import *
import sys
sys.path.append('..')
import dna_storage
import filecmp
import os

path_to_oligo_files = '/raid/nanopore/shubham/LDPC_DNA_storage_data/oligo_files_1/'

for i in range(num_files):
    index_len = (bin_index_len[i]+6*BCH_bits[i])//2
    sync_pos = (oligo_len-index_len)//2+index_len
    print('sync_pos',sync_pos)
    file_size = os.path.getsize(infile_path_prefix+infile_name[i])
    dna_storage.encode_data(infile_path_prefix+infile_name[i],oligo_len,path_to_oligo_files+'reads.'+str(i),BCH_bits[i],LDPC_alpha[i],LDPC_path_prefix+LDPC_code[i],bin_index_len[i],sync = sync[i], sync_pos = sync_pos)
    # test decoding to see that there is no strange issue
    tmpfile_decoded = 'tmpfile_decoded'
    dna_storage.decode_data(path_to_oligo_files+'reads.'+str(i),oligo_len,tmpfile_decoded,bin_index_len[i],BCH_bits[i],LDPC_alpha[i],LDPC_path_prefix+LDPC_code[i],file_size,0.01,sync=sync[i],sync_pos=sync_pos)
    assert filecmp.cmp(infile_path_prefix+infile_name[i],tmpfile_decoded)
    os.remove(tmpfile_decoded)
    with open(path_to_oligo_files+'reads.'+str(i)) as f_reads, open(path_to_oligo_files+'oligos_'+str(i)+'.fa','w') as f_oligos:
        for j, line in enumerate(f_reads):
            f_oligos.write('>oligos_'+str(i)+'_'+start_barcodes[i]+'_'+end_barcodes[i]+'_'+str(j)+'\n')
            f_oligos.write(start_barcodes[i]+line.rstrip('\n')+end_barcodes[i]+'\n')

# NOTE: the files generated are slightly different due to some random padding, but that doesn't affect the decoding
