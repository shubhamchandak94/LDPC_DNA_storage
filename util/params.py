# path containing LDPC matrices
LDPC_path_prefix = '/raid/nanopore/shubham/LDPC_DNA_storage/LDPC-codes/matrices/'
# path containing files to be encoded
infile_path_prefix = '/raid/nanopore/shubham/LDPC_DNA_storage_data/files_encoded/'

num_files = 9
start_barcodes = [
'TCCTGTGCTGCCTGTAATGAGCCAA',
'CAGCTGTCAATGCAATTGAGAATAT',
'GCTACATGTATACTGCGAGACAGAC',
'TCTATCTACTCGTGCTCGCTAGCTG',
'CATCAGCAGTAGAGAGTAGCGCGAT',
'TCACGAGATAGTCACGCTGTCACGT',
'GACGCGTCTAGCTGCATCTCGTACA',
'AGCGATGCTCTCTCAGTAGCATCTA',
'TAGTGTCGTCTGAGCGCAGAGATAT',
]

end_barcodes = [
'AGCATAGAACTGAGACCACGGATTG',
'TTGTGGACATGCGAACATTATGTGT',
'CGATAGTCGCAGTCGCACATCACTC',
'TGAGATCACAGCTACATAGTGAGAG',
'GTGTCACTATATCGCTCTACTGTGA',
'AGCTCGATGCGTGAGTGACTGTGAG',
'GTGTAGTGCGCATATGTCTAGACGT',
'ACACTACGAGCTATAGATCTACTGC',
'TCGCGCTAGACTCTGTATGAGTCGC',
]

infile_name = [
'random_data_160KB',
'random_data_224KB',
'2019-03-28_lossy.jpg.encrypted',        
'2019-03-28_lossy.jpg.encrypted',        
'2019-03-28_lossy.jpg.encrypted',        
'2019-03-28_lossy.jpg.encrypted',        
'2019-03-28_lossy.jpg.encrypted',        
'2019-03-28_lossy.jpg.encrypted',
'288.tar.gz.encrypted',
]

LDPC_code = [
'ldpc_9_old',
'ldpc_33_old',
'ldpc_9_new',
'ldpc_13_new',
'ldpc_33_new',
'ldpc_9_new',
'ldpc_9_new',
'ldpc_9_new',
'ldpc_33_new',
]

LDPC_alpha = [
0.5,
0.1,
0.5,
0.3,
0.1,
0.5,
0.5,
0.5,
0.1
]

BCH_bits = [
2,
2,
2,
2,
2,
3,
1,
2,
1,
]

sync = [
'',
'',
'AGT',
'AGT',
'AGT',
'AGT',
'AGT',
'',
'',
]

bin_index_len = [24]*2+[14]*7

oligo_len = 100
