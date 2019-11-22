#!/bin/bash
#$ -cwd

f=$1
if [ ! -f $f ]; then
  echo cannot find $f. stop.
  exit 1
fi
f_p1=$(basename $f .fastq)"_P1.fastq"
f_p2=$(basename $f .fastq)"_P2.fastq"
lenFull=$2
if [ -z $lenFull ]; then lenFull=202; echo Default lenFull of 202 used; fi
lenHalf=$(( $lenFull / 2 ))
if [ $f == $f_p1 ]; then echo could not derive new filename for $f; stop 1; fi
echo f is $f
echo f_p1 is $f_p1
echo f_p2 is $f_p2
echo lenFull is $lenFull
echo lenHalf is $lenHalf
cat $f | awk -v left_file=$f_p1 -v right_file=$f_p2 -v lenFull=$lenFull -v lenHalf=$lenHalf '{
	sub("length="lenFull, "length="lenHalf); 
	if(NR % 4 == 2 || NR % 4 == 0){
		left_read=substr($0, 1, lenHalf); 
		right_read=substr($0, lenHalf+1, lenFull); 
		print left_read > left_file; 
		print right_read > right_file
	}else{
		print $0 > left_file; 
		print $0 > right_file}
	}'
mkdir -p unsplit_fastq
mv $f unsplit_fastq
