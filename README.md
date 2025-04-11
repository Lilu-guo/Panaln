# Panaln
version 1.3（20240630）

# What is it?
Indexing pan-genome with applications in reads mapping and alignment.

# How to use it?
Panaln consists of three components, data preprocessing, pan-index building, and reads mapping. After preprocessing the data format, you should first build the pan-index with the VCF file (e.g., snp144common.txt from dbSNP dataset) and the reference genome (e.g., GRCh38.fasta), then perform the mapping processing. Since the data is quite large, we put the download link below or you can email us to ask for it.

# Available data:
Illumina reads dataset:   
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/AshkenazimTrio/sequence.index.AJtrio_Illumina_2x250bps_06012016_updated.HG004   
PacBio-CCS reads dataset:   
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_10kb/m54238_180628_014238.Q20.fastq   
Common small variants:   
http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp144Common.txt.gz   
Synthetic reads simulator:   
https://github.com/DaehwanKimLab/hisat2/blob/master/hisat2_simulate_reads.py    
Pipeline of variant calling:   
https://github.com/ksahlin/strobealign/blob/main/evaluation.md   

# Step I. Install
  1. Download (or clone) the source code form https://github.com/Lilu-guo/Panaln
  2. Compile the source code. (Note that you need to compile FM and WFA first)

# Step II. Data preprocessing
  1. Convert vcf format to custom snp format with five columns by running: "./vcf2snp snp144common.txt snp144common.snp". (If your reference genome is an entire sequence, you can split it into independent chromosomes by running: "genom2chr" .)
  2. Embed snp information into the reference genome by running: "./snp_embed". (The default output path is "/home/lab/gll/formatsnp", you can modify it according to your preference.)
  3. Output the snp and indel information by running: "./snp_indel". (The default output path is "/home/lab/gll/formatsnp/panVcf", you can modify it according to your preference.)
  4. Generate the sequence of the linear serialization model by running: "./comb -w 124 genome.fa snp144Comm.fasta snp144Comm_indel.fasta snp144Comm.data"
     
# Step III. Build pan-index
  1. Run the shell command: "./panaln index snp144Comm_indel.fasta".
  2. Get the index file as: "sal.idx3, newD.idx3, lroot.idx3, rroot.idx3, rankl.idx3, and U.idx3"
  
# Step IV. Mapping processing
  1. run the shell command: "./panaln align <pan_index_Name> <fastq_Name> <sam_Name>, where <pan_index_Name> is the index file that has been constructed, <fastq_Name> is the sequencing reads file, and <sam_Name> is the mapping result file.
  
# Feedback
Please report bugs to Email: guolilu@stu.xidian.edu.cn if any questions or suggestions. Your feedback and test results are welcome.