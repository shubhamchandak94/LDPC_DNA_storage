import sys
sys.path.append('..')
import dna_storage
import filecmp
import os
from joblib import Parallel, delayed
file_name = '/raid/nanopore/shubham/LDPC_DNA_storage/util/random_data_32KB'
file_size = 32000
bin_index_len = 14
sub_prob = 0.004
del_prob = 0.0085
ins_prob = 0.0005
frac_random_reads = 0.15
oligo_length = 100
sync='AGT'
BCH_bits = 2
index_len = (bin_index_len+6*BCH_bits)//2
sync_pos = (oligo_length-index_len)//2+index_len
LDPC_alpha = 0.5
ldpc_code_prefix='/raid/nanopore/shubham/LDPC_DNA_storage/LDPC-codes/matrices/ldpc_9_new'
tmpfile = 'tmpfile_9'

# encode
dna_storage.encode_data(file_name,oligo_length,tmpfile+'.oligos',BCH_bits,LDPC_alpha,ldpc_code_prefix,bin_index_len,sync=sync,sync_pos=sync_pos)

def run_exp(i):
    tmpfile_prefix = tmpfile+'.'+str(i)
    dna_storage.sample_reads_indel(tmpfile+'.oligos',tmpfile_prefix+'.reads',num_reads,sub_prob=sub_prob, del_prob=del_prob, ins_prob=ins_prob, frac_random_reads=frac_random_reads)
    #decode
    dna_storage.decode_data(tmpfile_prefix+'.reads',oligo_length,tmpfile_prefix+'.decoded',bin_index_len,BCH_bits,LDPC_alpha,ldpc_code_prefix,file_size,0.04,sync=sync,sync_pos=sync_pos)
    if filecmp.cmp(file_name,tmpfile_prefix+'.decoded'):
        correct = 1
    else:
        correct = 0
    os.remove(tmpfile_prefix+'.reads')
    os.remove(tmpfile_prefix+'.decoded')
    return correct

reading_cost = 2.3 
while True:
    num_reads = int(file_size*8*reading_cost/oligo_length)
    print('reading cost',reading_cost)
    print('num_reads',num_reads)
    num_it = 0
    num_correct = 0
    while True:
        correct_arr = Parallel(n_jobs=10)(delayed(run_exp)(i) for i in range(10))
        num_it += len(correct_arr)
        num_correct += sum(correct_arr)
        print('num_it:',num_it)
        print('num_correct:',num_correct)
        if num_it >= 50 and num_it - num_correct >= 10:
            break
    error_rate = 1-num_correct/num_it
    print('error rate',error_rate)
    if error_rate < 1e-4:
        break
    reading_cost += 0.1
