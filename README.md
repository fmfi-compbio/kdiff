# kdiff

kdiff is a tool for the alignment-free detection of differences between sequencing data sets.


### Installation
The following commands will compile kdiff and its dependencies (e.g., KMC3):
```sh
git clone https://github.com/fmfi-compbio/kdiff.git
cd kdiff
mkdir build ; cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..
./kdiff -h
```

### How-to
kdiff requires as input a reference genome and two samples (control/case). Since kdiff is based on k-mers, it is first necessary to compute them (using KMC3). kdiff binary, indeed, requires the reference sequence and three KMC3 databases (reference, control, case).

Here we will explain how to run kdiff using the data provided in the `example` folder. Please refer to `./kdiff -h` for additional options. In this example, we will use the kmc binaries built while compiling kdiff, but you can use your own `kmc` and `kmc_tools` binary.
``` sh
cd example

# 1 count kmers from reference sequence
../build/kmc-prefix/src/kmc/bin/kmc -m16 -t4 -k21 -ci1 -cs65535 -fm reference.fa reference-k21 .

# 2 count kmers from case (assuming paired-end fastq)
ls control_1.fq control_2.fq > control.list
../build/kmc-prefix/src/kmc/bin/kmc -m16 -t4 -k21 -ci2 -cs65535 @control.list control-k21 .

# 3 count kmers from control (assuming paired-end fastq)
ls case_1.fq case_2.fq > case.list
../build/kmc-prefix/src/kmc/bin/kmc -m16 -t4 -k21 -ci2 -cs65535 @case.list case-k21 .

# 3 intersect case/control databases with reference databses to exrtact only those kmers occuring in the reference
../build/kmc-prefix/src/kmc/bin/kmc_tools simple reference-k21 control-k21 intersect control-filter-k21 -ocright
../build/kmc-prefix/src/kmc/bin/kmc_tools simple reference-k21 case-k21 intersect case-filter-k21 -ocright

# 4 run kdiff
../kdiff reference.fa reference-k21 control-filter-k21 case-filter-k21 -w 100 > bins.depth
```

kdiff, by default, partitions the reference genome in windows (of size `-w`, 1000 by default). However, it is also possible to run kdiff on precomputed regions of the genomes. This is done by passing the list list of regions of interest using the `-r` option (BED format, 0-based and half-open).
