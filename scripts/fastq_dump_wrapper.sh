#!/bin/bash
#$ -cwd
if [ -z $SRR ]; then
  SRR=$1
fi
if [ -z $FASTQ ]; then
  FASTQ=$2
fi
FASTQ=${FASTQ/.fastq/""}.fastq
echo SRR is $SRR
echo FASTQ is $FASTQ
if [ ! -f $FASTQ ]; then
  /slipstream/usr/local/bin/sratoolkit.2.8.1-ubuntu64/bin//fastq-dump $SRR
  mv ${SRR}.fastq ${FASTQ}
fi
