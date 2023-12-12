# Panaln
version 1.0（20230920）

# What is it?
Indexing Pan-genome with Applications in Read Mapping and Alignment.

# How to use it?
Panaln consists of three components, data preprocessing, pan-index building, and read mapping. After preprocessing the data format, you should first build the pan-index with the VCF file (e.g., snp144common.txt from dbSNP dataset) and the reference genome (e.g., GRCh38.fasta), then perform the mapping processing. Since the data is quite large, we put the download link below or you can send me an email to request it.

# Available data:
reference genome: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.4    
vcf file from dbSNP: https://ftp.ncbi.nih.gov/snp/    
sequencing read: https://www.ebi.ac.uk/ena/browser/view/ERR194146     

# Step I. Install
  1. Download (or clone) the source code form https://github.com/Hongweihuo-Lab/Panaln
  2. Compile the source code. (Note that you need to compile FM first)

# Step II. Data preprocessing
  1. Convert vcf format to custom snp format with five columns by running: "./vcf2snp snp144common.txt snp144common.snp". (If your reference genome is an entire sequence, you can split it into independent chromosomes by running: "genom2chr" .)
  2. Embed snp information into the reference genome by running: "./snp_embed". (The default output path is "/home/lab/gll/formatsnp", you can modify it according to your preference.)
  3. Output the snp and indel information by running: "./snp_indel". (The default output path is "/home/lab/gll/formatsnp/panVcf", you can modify it according to your preference.)
  4. Generate the sequence of the linear serialization model by running: "./comb -w 124 genome.fa snp144Comm.fasta snp144Comm_indel.fasta snp144Comm.data"
     
# Step III. Build pan-index
  1. Run the shell command: "./panaln index snp144Comm_indel.fasta".
  2. Get the index file as: "sal.idx.256, newD.idx.256, lroot.idx.256, rroot.idx.256"
  
# Step IV. Mapping processing
  1. run the shell command: "./panaln align -n 3 <processed_fasta_Name> <fastq_Name> <aln_Name>, where processed_fasta_Name is the processed genome name, fastq_Name is the sequencing reads file name, aln_Name is the intermediate alignment file name.
  2. run the shell command: "./panaln aln2sam <processed_fasta_Name> <aln_Name> <sam_Name>, where processed_fasta_Name is the processed genome name, aln_Name is the intermediate alignment file name, sam_Name is the mapping result file.
  
# Feedback
Please report bugs to Email: guolilu@stu.xidian.edu.cn if any questions or suggestions. Your feedback and test results are welcome.
