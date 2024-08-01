# SPARseq
This is the code for the Wrana Lab's SPARseq analysis pipeline. Originally designed and written by Seda Barutcu [https://github.com/seda-barutcu/MultiplexedPCR-DeepSequence-Analysis], upgraded and streamlined by Kathleen-Rose Zakoor, and maintained and updated by Lauren Caldwell. 

Dependencies
1. Bowtie 0.12.7 (http://bowtie-bio.sourceforge.net/index.shtml)
2. HTSeq 0.13.5 (https://github.com/htseq)
3. Python3 [os, argparse, re, openpyxl, pandas, datetime]
Tested on MacOS 12.7.5


Versions

As covid19 variants emerged, we updated the analysis pipeline. The initial research version was used from January 2021-May 2021. When we began running the analysis daily in June 2021, we updated the pipeline into what is called the V1 pipeline. At the end of June 2021, we replaced the Srbd amplicon with a slightly longer one to have more coverage. 


This is the pipeline for the daily clinical sample processing. 


The run folder and resource folder must be subfolders of the same parent folder:

Run folder structure

	/runclX  
		-/R1_files
			-fastq files

Resource folder structure

	/resources
		-/bowtie_index
				-bowtie index files
				-reference fa file
				-reference GTF file
		-codon_aa_table.csv
		-Clinical_V2_sparseq_analysis_V1.3.2.py
		-Clinical_V2_sparseq_pipeline_V1.3.2.sh

Explanation of files:
codon_aa_table.csv contains nucleotide to amino acid conversion info.
Clinical_V2_sparseq_pipeline_V1.3.2.sh is the wrapper script that is run to launch the rest of the analysis.
Clinical_V2_sparseq_analysis_V1.3.2.py is the analysis script that is run after bowtie and htseq calls to generate output files.


Initial one-time step:
Build bowtie index from within the bowtie_index folder with
bowtie-build clinical_V2_sparsq_V1p2_amplicons.fa clinical_V2_sparsq_V1p2_amplicons

Instructions for one run:
1. Create run folder runclX and subdirectory R1_files; place fastq files into R1_files
2. Run $bash Clinical_V2_sparseq_pipeline_V1.3.2.sh path/to/run/folder/runclXX path/to/resource/folder
		> Wrapper script performs the following: 
			i.   unzips fastqs
			ii.  Creates samcounts directory then runs bowtie and HTseq
			iii. Creates results directory then run Clinical_V2_sparseq_analysis_V1.3.2 which outputs several key results tables 
				-runclxxx_AllVariantDetails.xlsx
				-runclxxx_CountTable.xlsx
				-runclxxx_SelectedVariants.xlsx
				-sparseq_report_runclxxx.xlsx
				-please note variant calling is manually completed by a supervisor based on review of the 4 outputs 
			vi.  Re-zips FASTQ files


