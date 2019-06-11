# path containing LDPC matrices
LDPC_path_prefix = '/raid/nanopore/shubham/LDPC_DNA_storage/LDPC-codes/matrices/'
# path containing files to be encoded
infile_path_prefix = '/raid/nanopore/shubham/LDPC_DNA_storage/genap/files_encoded/'

num_files = 3
start_barcodes = [
'TCCTGTGCTGCCTGTAATGAGCCAA',
'CAGCTGTCAATGCAATTGAGAATAT',
'GCTACATGTATACTGCGAGACAGAC',
]

end_barcodes = [
'AGCATAGAACTGAGACCACGGATTG',
'TTGTGGACATGCGAACATTATGTGT',
'CGATAGTCGCAGTCGCACATCACTC',
]

infile_name = ['2333main_MM_Image_Feature_19_rs4_resized.jpg.enc']*3

LDPC_code = [
'ldpc_9_genap',
'ldpc_13_genap',
'ldpc_33_genap',
]

LDPC_alpha = [
0.5,
0.3,
0.1,
]

BCH_bits = [
2,
2,
2,
]

sync = [
'',
'',
'',
]

bin_index_len = [14]*3

oligo_len = 50
