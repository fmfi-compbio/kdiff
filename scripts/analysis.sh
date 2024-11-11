#!/bin/sh

RD="$(dirname $0)"

fa=$1
region=$2
rk=$3
k1=$4
k2=$5
k3=$6
bw2=$7
bw3=$8
bam1=$9
bam2=${10}
bam3=${11}
wd=${12}
k=${13}

mkdir -p $wd
samtools faidx $fa $region > $wd/region.fa

chrom=$(echo $region | cut -f1 -d':')
s=$(echo $region | cut -f2 -d':' | cut -f1 -d'-')
s=$((s-1))
e=$(echo $region | cut -f2 -d':' | cut -f2 -d'-')
e=$((e-1))

kmc -t2 -k$k -ci1 -cs65535 -fm $wd/region.fa $wd/region-k$k $wd

i=1
for kk in $rk $k1 $k2 $k3
do
  kmc_tools simple $wd/region-k$k $kk intersect $wd/$i -ocright
  kmc_dump $wd/$i $wd/$i.txt
  i=$((i+1))
done

bigWigToBedGraph $bw2 -chrom=$chrom -start=$s -end=$((e+1)) $wd/3.bedgraph
bigWigToBedGraph $bw3 -chrom=$chrom -start=$s -end=$((e+1)) $wd/4.bedgraph

python3 $RD/plot.py $wd/region.fa $wd/1.txt $wd/2.txt $wd/3.txt $wd/4.txt $wd/3.bedgraph $wd/4.bedgraph $bam1 $bam2 $bam3


