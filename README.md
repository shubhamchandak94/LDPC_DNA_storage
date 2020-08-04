# LDPC DNA storage
LDPC codes for Illumina sequencing-based DNA storage. The associated data is available at https://github.com/shubhamchandak94/LDPC_DNA_storage_data.
### Biorxiv preprint: https://www.biorxiv.org/content/10.1101/770032v2
### [Supplmentary material](https://github.com/shubhamchandak94/LDPC_DNA_storage/blob/master/supplementary_material.pdf)

### [Allerton 2019](https://ieeexplore.ieee.org/document/8919890)

Installation instructions (tested on Ubuntu 18.04.1)
```
# clone repository with submodules
git clone --recursive https://github.com/shubhamchandak94/LDPC_DNA_storage
cd LDPC_DNA_storage/
# install flexbar (our code tested with Flexbar 3.0.3, source code: https://github.com/seqan/flexbar)
sudo apt-get install flexbar
# build LDPC codes
cd LDPC-codes/
make
cd ..
# build Kalign MSA
cd kalign2-current/
./configure
make
cd ..
# install BCH codes Python library
cd python-bchlib/
python3 setup.py build
sudo python3 setup.py install

# install joblib for parallelization of kalign consensus
pip3 install --user joblib
```

## Instructions for encoding and decoding data
The code is implemented in Python3.
```
import dna_storage
```
### Encoding
```
dna_storage.encode_data(infile, oligo_length, outfile, BCH_bits, LDPC_alpha, LDPC_prefix, \
bin_index_len, sync='', sync_pos=-1)
```
Parameters:
```
infile:         file to be encoded in DNA (we recommend that the file is compressed and encrypted 
                to randomize the data and avoid issues with DNA synthesis)
oligo_length:   length of oligo
outfile:        file to write oligos to
BCH_bits:       number of bit errors that the BCH code can correct (the number of parity bits used 
                is 6*BCH_bits)
LDPC_alpha:     the redundancy of the LDPC code. The code should map encode 256000 bits into 
                256000(1+alpha) bits.
LDPC_prefix:    the path and prefix of the LDPC matrices. The encoder requires the files 
                LDPC_prefix.gen and LDPC_prefix.systematic.
bin_index_len:  number of bits used for the index (before adding BCH). Should be large enough to
                accomodate the number of oligos.
sync:           (Optional) Synchronization marker string to be added to each oligo (we use 'AGT' in 
                our experiments).
sync_pos:       (Optional) Position to put the synchronization marker. We set it at the center of 
                the payload in the oligo (the part of the oligo after the index).
```
See [`util/encode_files.py`](https://github.com/shubhamchandak94/LDPC_DNA_storage/blob/master/util/encode_files.py) for examples.

### Decoding
Note that the decoder expects a file with only the reads. To extract the reads from a FASTQ file, run (on the shell):
```
sed -n '2~4p' file.fastq > file.reads
```
First we remove the barcodes at the start and end of each reads by running
```
dna_storage.remove_barcodes_flexbar(infile_reads, start_barcode, end_barcode, outfile_reads)
```
This internally uses [Flexbar](https://github.com/seqan/flexbar) to remove the barcodes and also corrects the orientation of any reverse complemented reads. We tested our code with Flexbar 3.0.3.

Then we run the decoder on the trimmed reads:
```
dna_storage.decode_data(infile, oligo_length, outfile, bin_index_len, BCH_bits, LDPC_alpha, \
LDPC_prefix, file_size, eps, sync='', sync_pos=-1, attempt_indel_cor=True)
```
Parameters:
```
infile:             file containing trimmed reads
oligo_length:       same as used in encode_data
outfile:            file to write decoded data to
bin_index_len:      same as used in encode_data
BCH_bits:           same as used in encode_data
LDPC_alpha:         same as used in encode_data
LDPC_prefix:        same as used in encode_data
file_size:          size of the decoded file in bytes (this is required to calculate 
                    the number of LDPC blocks and also for removing any padding)
eps:                eps value used for computing LDPC LLRs (we used 4%, 
                    i.e., 0.04 for most experiments)
sync:               (Optional) same as used in encode_data
sync_pos:           (Optional) same as used in encode_data
attempt_indel_cor:  use indel correction heuristic during BCH decoding
```
See [`util/decode_files.py`](https://github.com/shubhamchandak94/LDPC_DNA_storage/blob/master/util/decode_files.py) for examples.

### Running simulations
To run simulations for testing parameters, we provide a function which encodes the data and finds the minimum reading cost (in bases/bit) needed for successful decoding. The function separates simulations with independent substitution, deletion and insertion errors and also supports adding random reads (to simulate unaligned reads).
```
dna_storage.find_min_coverage(infile_data, oligo_length, BCH_bits, LDPC_alpha, LDPC_prefix, \
bin_index_len, file_size, sub_prob, eps_decode, num_experiments, ins_prob=0.0, del_prob=0.0, \
start_coverage=1.0, sync='', sync_pos=-1, attempt_indel_cor=True, frac_random_reads=0.0)
```
Parameters:
```
infile:             file to be encoded in DNA
oligo_length:       same as used in encode_data
BCH_bits:           see parameters for encode_data
LDPC_alpha:         see parameters for encode_data
LDPC_prefix:        see parameters for encode_data
bin_index_len:      see parameters for encode_data
file_size:          size of infile in bytes
sub_prob:           substitution error probability
eps_decode:         see parameter eps in decode_data
num_experiments:    number of successful trials before declaring success 
                    (e.g., if num_experiments=20, success is declared if 20/20 trials succeed)
ins_prob:           insertion error probability
del_prob:           deletion error probability
start_coverage:     reading cost (bases/bit) to start with (the function then increases it
                    by 0.1 until it succeeds)
sync:               see parameters for encode_data            
sync_pos:           see parameters for encode_data            
attempt_indel_cor:  see parameters for decode_data
frac_random_reads:  fraction of random reads to add while simulating 
```
Returns:
```
minimum reading cost (in bases/bit)
```
See [`util/run_simulations.py`](https://github.com/shubhamchandak94/LDPC_DNA_storage/blob/master/util/run_simulations.py) for examples.
