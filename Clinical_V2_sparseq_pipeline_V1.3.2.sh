#/bin/bash

# Dependencies:
#	1. Bowtie 1 (http://bowtie-bio.sourceforge.net/index.shtml)
#	2. HTSeq* (https://htseq.readthedocs.io/en/release_0.11.1/index.html)
#	3. Python3 [os, argparse, re, openpyxl, pandas, datetime]*
#		a. sparseq_analysis_v1.3.X.py
#	*To install pkgs
#		$ pip install --user [package]

run_folder=$1 #./runclX/
resource_folder=$2 #./resources/SPAR_SEQ/

echo "***** Beginning SPAR-Seq Pipeline "$(date)"*****"

echo "	*Unzipping FASTQ files...*"
gunzip $run_folder/R1_files/*.gz

echo "	*Running bowtie alignment and counting reads with HTseq...*"
mkdir -p $run_folder/R1_files/samcounts
for f in $run_folder/R1_files/*.fastq; do bowtie $resource_folder/bowtie_index_V1p2/bw1_sparsq_V1p2_amplicons_index $f $f.sam --best -v 3 -k 1 -m 1 -S; done && for x in $run_folder/R1_files/*.sam; do python3 -m HTSeq.scripts.count -f sam -t CDS $x $resource_folder/bowtie_index_V1p2/sparsq_V1p2_amplicons.gtf > $x.count.txt; done;
mv $run_folder/R1_files/*.count.txt $run_folder/R1_files/samcounts

echo "	*Running analysis...*"
mkdir -p $run_folder/results
python3 $resource_folder/sparseq_analysis_V1.3.2.py $run_folder $resource_folder

echo "	*Zipping FASTQ files...*"
gzip $run_folder/R1_files/*.fastq

echo "	*Remove sam files*" 
rm -v $run_folder/R1_files/*.sam

echo "*****"$(date)" SPAR-Seq Pipeline is now complete. *****"
