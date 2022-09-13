#!/bin/bash
#$ -cwd
#$ -V
#$ -N fetch

#fetch script
srr=$1
root=$2
fetch_path=$3

DUMP_PATH=/slipstream/home/joeboyd/bin/sratoolkit

if [ -z $fetch_path ]; then
fetch_path=$(pwd)
fi

echo srr is $srr
echo root is $root
echo fetch_path is $fetch_path

mkdir -p $fetch_path
cd $fetch_path
 
dl_fq1=${fetch_path}/${srr}"_1.fastq"
dl_fq2=${fetch_path}/${srr}"_2.fastq"
dl_fqse=${fetch_path}/${srr}".fastq"

fq1=${fetch_path}/${root}"_R1_001.fastq"
fq2=${fetch_path}/${root}"_R2_001.fastq"
fqse=${fetch_path}/${root}"_R1_001.fastq"

if [ -f ${fq1}.gz ] && [ -f ${fq2}.gz ]; then
 echo final fastqs are present, will not rerun.;
 exit 0;
else if [ -f ${fqse}.gz ]; then
 echo final fastq is present, will not rerun.;
 exit 0;  
fi

if [ -d ${fetch_path}/${srr} ]; then
echo found prefetch for $srr, will not rerun.;
else
  $DUMP_PATH/prefetch --max-size 500G -O $fetch_path  $srr
fi

if [ -f ${dl_fq1} ] && [ -f ${dl_fq2} ] || [ -f ${dl_fqse} ]; then
  echo found dumped fastqs, will not rerun.
else
  if [ -f ${dl_fqse} ]; then
    echo found dumped fastq, will not rerun.
  else
    $DUMP_PATH/fastq-dump -O $fetch_path --split-e ${fetch_path}/${srr}.sra
  fi
fi

if [ -f ${dl_fq1} ]; then
  mv ${dl_fq1} ${fq1}
  gzip ${fq1}
fi

if [ -f ${dl_fq2} ]; then
  mv ${dl_fq2} ${fq2}
  gzip ${fq2}
fi

if [ -f ${dl_fqse} ]; then
  mv ${dl_fqse} ${fqse}
  gzip ${fqse}
fi
