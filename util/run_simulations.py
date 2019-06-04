import sys
sys.path.append('..')
import dna_storage
file_name = '/raid/nanopore/shubham/LDPC_DNA_storage/util/random_data_224KB'
file_size = 224000
bin_index_len = 14
sub_prob = 0.004
del_prob = 0.0085
ins_prob = 0.0005
num_trials = 20
frac_random_reads = 0.15
oligo_length = 100
sync='AGT'
LDPC_alpha = 0.5
ldpc_code_prefix='/raid/nanopore/shubham/LDPC_DNA_storage/LDPC-codes/matrices/ldpc_9_new'

# first vary BCH_bits from 0 to 3
#for BCH_bits in range(3):
#    print('BCH_bits:',BCH_bits)
#    index_len = (bin_index_len+6*BCH_bits)//2
#    sync_pos = (oligo_length-index_len)//2+index_len
#    min_coverage = dna_storage.find_min_coverage(file_name,oligo_length,BCH_bits,LDPC_alpha,ldpc_code_prefix,bin_index_len,file_size,sub_prob,0.04,num_trials,ins_prob = ins_prob, del_prob = del_prob, start_coverage = 2.0, sync = sync, sync_pos=sync_pos,frac_random_reads = frac_random_reads)
#    print('min_coverage',min_coverage)

# for BCH_bits = 3, index len = 16
BCH_bits = 3
bin_index_len = 16
print('BCH_bits:',BCH_bits)
index_len = (bin_index_len+6*BCH_bits)//2
sync_pos = (oligo_length-index_len)//2+index_len
min_coverage = dna_storage.find_min_coverage(file_name,oligo_length,BCH_bits,LDPC_alpha,ldpc_code_prefix,bin_index_len,file_size,sub_prob,0.04,num_trials,ins_prob = ins_prob, del_prob = del_prob, start_coverage = 2.0, sync = sync, sync_pos=sync_pos,frac_random_reads = frac_random_reads)
print('min_coverage',min_coverage)

# without sync
bin_index_len = 14
min_coverage = dna_storage.find_min_coverage(file_name,oligo_length,2,LDPC_alpha,ldpc_code_prefix,bin_index_len,file_size,sub_prob,0.04,num_trials,ins_prob = ins_prob, del_prob = del_prob, start_coverage = 2.0,frac_random_reads = frac_random_reads)
print('min_coverage',min_coverage)

# stress testing 
bin_index_len = 16
index_len = (bin_index_len+6*BCH_bits)//2
sync_pos = (oligo_length-index_len)//2+index_len
sub_prob = 0.02
ins_prob = 0.02
del_prob = 0.02
min_coverage = dna_storage.find_min_coverage(file_name,oligo_length,3,LDPC_alpha,ldpc_code_prefix,bin_index_len,file_size,sub_prob,0.1,num_trials,ins_prob = ins_prob, del_prob = del_prob, start_coverage = 10.0,frac_random_reads = frac_random_reads, sync = sync, sync_pos=sync_pos)
print('min_coverage',min_coverage)
