LDPC_path="../LDPC-codes/"
LDPC_mat_path="../LDPC-codes/matrices"
mkdir $LDPC_mat_path

# matrices generated earlier (spring 2018)

# 50% redundancy (3,9) code
$LDPC_path/make-ldpc $LDPC_mat_path/ldpc_9_old.pchk 128000 384000 1 evenboth 3
$LDPC_path/make-gen $LDPC_mat_path/ldpc_9_old.pchk $LDPC_mat_path/ldpc_9_old.gen dense
$LDPC_path/extract_systematic $LDPC_mat_path/ldpc_9_old.gen $LDPC_mat_path/ldpc_9_old.systematic

# 30% redundancy (3,13) code
$LDPC_path/make-ldpc $LDPC_mat_path/ldpc_13_old.pchk 76800 332800 1 evenboth 3
$LDPC_path/make-gen $LDPC_mat_path/ldpc_13_old.pchk $LDPC_mat_path/ldpc_13_old.gen dense
$LDPC_path/extract_systematic $LDPC_mat_path/ldpc_9_old.gen $LDPC_mat_path/ldpc_9_old.systematic

# 20% redundancy (3,18) code
$LDPC_path/make-ldpc $LDPC_mat_path/ldpc_18_old.pchk 51200 307200 1 evenboth 3
$LDPC_path/make-gen $LDPC_mat_path/ldpc_18_old.pchk $LDPC_mat_path/ldpc_18_old.gen dense
$LDPC_path/extract_systematic $LDPC_mat_path/ldpc_18_old.gen $LDPC_mat_path/ldpc_18_old.systematic

# 10% redundancy (3,33) code
$LDPC_path/make-ldpc $LDPC_mat_path/ldpc_33_old.pchk 25600 281600 1 evenboth 3
$LDPC_path/make-gen $LDPC_mat_path/ldpc_33_old.pchk $LDPC_mat_path/ldpc_33_old.gen dense
$LDPC_path/extract_systematic $LDPC_mat_path/ldpc_33_old.gen $LDPC_mat_path/ldpc_33_old.systematic

# matrices generated later (winter 2019)

# 50% redundancy (3,9) code
$LDPC_path/make-ldpc $LDPC_mat_path/ldpc_9_new.pchk 128000 384000 42 evenboth 3 no4cycle
$LDPC_path/make-gen $LDPC_mat_path/ldpc_9_new.pchk $LDPC_mat_path/ldpc_9_new.gen dense
$LDPC_path/extract_systematic $LDPC_mat_path/ldpc_9_new.gen $LDPC_mat_path/ldpc_9_new.systematic

# 30% redundancy (3,13) code 
$LDPC_path/make-ldpc $LDPC_mat_path/ldpc_13_new.pchk 76800 332800 42 evenboth 3 no4cycle
$LDPC_path/make-gen $LDPC_mat_path/ldpc_13_new.pchk $LDPC_mat_path/ldpc_13_new.gen dense
$LDPC_path/extract_systematic $LDPC_mat_path/ldpc_13_new.gen $LDPC_mat_path/ldpc_13_new.systematic

# 10% redundancy (3,33) code
$LDPC_path/make-ldpc $LDPC_mat_path/ldpc_33_new.pchk 25600 281600 42 evenboth 3 no4cycle
$LDPC_path/make-gen $LDPC_mat_path/ldpc_33_new.pchk $LDPC_mat_path/ldpc_33_new.gen dense
$LDPC_path/extract_systematic $LDPC_mat_path/ldpc_33_new.gen $LDPC_mat_path/ldpc_33_new.systematic
