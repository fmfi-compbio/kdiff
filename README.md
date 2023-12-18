# README
```sh
git clone https://github.com/ldenti/mfcnv.git
cd mfcnv
mkdir build ; cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..
```

### How-to
``` sh
cd example

# 1 count kmers from reference
kmc -m16 -t4 -k21 -ci1 -cs65535 -fm reference.fa reference-k21 .

# 2 count kmers from case and control (paired-end fastq)
ls control_1.fq control_2.fq > control.list
kmc -m16 -t4 -k21 -ci2 -cs65535 @control.list control-k21 .

ls case_1.fq case_2.fq > case.list
kmc -m16 -t4 -k21 -ci2 -cs65535 @case.list case-k21 .

# 3 intersect with reference
kmc_tools simple reference-k21 control-k21 intersect control-filter-k21 -ocright
kmc_tools simple reference-k21 case-k21 intersect case-filter-k21 -ocright

# 3 run kdiff
../MFCNV reference.fa control-filter-k21 case-filter-k21 -w 100 > bins.depth
```
