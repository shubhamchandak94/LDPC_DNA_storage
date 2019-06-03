import sys
import random
from collections import Counter
import numpy as np
import scipy.stats
#random.seed(42)

if len(sys.argv) != 4:
    print('Incorrect number of arguments')
    sys.exit(1)

with open(str(sys.argv[1])) as f:
    original_ref_list = [l.rstrip('\n') for l in f.readlines()]
num_oligos=int(sys.argv[2])
coverage=float(sys.argv[3])
num_reads = int(num_oligos*coverage)
print('num_oligos',num_oligos)
print('coverage',coverage)
print('num_reads',num_reads)

subsampled_ref_list = random.choices(original_ref_list,k=num_reads)
c = Counter(subsampled_ref_list)

count_list = []
for elem in c:
    count_list.append(c[elem])

count_list += (num_oligos-len(count_list))*[0]
print('mean',np.mean(count_list))
print('coefficient of variation',scipy.stats.variation(count_list))
