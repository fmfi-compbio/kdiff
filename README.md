# README
```sh
git clone --recursive git@github.com:ldenti/mfcnv.git
cd mfcnv
make
```

### How-to
``` sh
# 1 count kmers from reference
kmc -k21 -ci1 -fa reference.fa reference-k21 .

# 2 count kmers from case and control (fasta)

# if you have paired-end, concatenate them
ls case_1.fa case_2.fa > case.fa
ls control_1.fa control_2.fa > control.fa

# if you have fasta, use -fa
kmc -k21 -ci1 -fa case.fa case-k21 .
kmc -k21 -ci1 -fa control.fa control-k21 .

# if you have fastq, use -fq
kmc -k21 -ci1 -fq case.fq case-k21 .
kmc -k21 -ci1 -fq control.fq control-k21 .

# 3 intersect with reference
kmc_tools simple reference-k21 control-k21 intersect control-filter-k21 -ocright
kmc_tools simple reference-k21 case-k21 intersect case-filter-k21 -ocright

# 3 run kdiff
./main reference.fa control-filter-k21 case-filter-k21 300 > bins.depth
```

### Example
```
cd example
kmc -k21 -ci1 -fa reference.fa reference-k21 .

ls control_1.fq control_2.fq > control.lst
ls case_1.fq case_2.fq > case.lst

kmc -k21 -ci1 -fq @control.lst control-k21 .
kmc -k21 -ci1 -fq @case.lst case-k21 .

kmc_tools simple reference-k21 control-k21 intersect control-filter-k21 -ocright
kmc_tools simple reference-k21 case-k21 intersect case-filter-k21 -ocright

../main reference.fa control-filter-k21 case-filter-k21 300 > bins.depth
```
