# LDPC DNA storage
LDPC codes for Illumina sequencing-based DNA storage

Installation instructions (tested on Ubuntu 18.04.1)
```
# clone repository with submodules
git clone --recursive https://github.com/shubhamchandak94/LDPC_DNA_storage
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

