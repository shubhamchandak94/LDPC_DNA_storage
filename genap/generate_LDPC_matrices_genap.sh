# The .pchk and .systematic matrices can be downloaded directly from 
# https://github.com/shubhamchandak94/LDPC_DNA_storage_data/tree/master/matrices
# NOTE: make-gen command can take few hours to generate the .gen matrix file.

LDPC_path="../LDPC-codes/"
LDPC_mat_path="../LDPC-codes/matrices"
mkdir $LDPC_mat_path

# matrices generated Jun 10, 2019 with block size 240000

# 50% redundancy (3,9) code
$LDPC_path/make-ldpc $LDPC_mat_path/ldpc_9_genap.pchk 120000 360000 42 evenboth 3 no4cycle
$LDPC_path/make-gen $LDPC_mat_path/ldpc_9_genap.pchk $LDPC_mat_path/ldpc_9_genap.gen dense
$LDPC_path/extract_systematic $LDPC_mat_path/ldpc_9_genap.gen $LDPC_mat_path/ldpc_9_genap.systematic

# 30% redundancy (3,13) code 
$LDPC_path/make-ldpc $LDPC_mat_path/ldpc_13_genap.pchk 72000 312000 42 evenboth 3 no4cycle
$LDPC_path/make-gen $LDPC_mat_path/ldpc_13_genap.pchk $LDPC_mat_path/ldpc_13_genap.gen dense
$LDPC_path/extract_systematic $LDPC_mat_path/ldpc_13_genap.gen $LDPC_mat_path/ldpc_13_genap.systematic

# 10% redundancy (3,33) code
$LDPC_path/make-ldpc $LDPC_mat_path/ldpc_33_genap.pchk 24000 264000 42 evenboth 3 no4cycle
$LDPC_path/make-gen $LDPC_mat_path/ldpc_33_genap.pchk $LDPC_mat_path/ldpc_33_genap.gen dense
$LDPC_path/extract_systematic $LDPC_mat_path/ldpc_33_genap.gen $LDPC_mat_path/ldpc_33_genap.systematic
