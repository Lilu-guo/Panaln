# Panaln
 version 2.0（20250410）

# What is it?
 Panaln: Indexing pangenome for read alignment.

# How to use it?
 ```
 Usage:   panaln [combine /index /align] <options>  
 --------------  
 Function: combine  generate pangenome file  
           index    construct a pangenome index  
           align    perform read alignment  
 --------------  
 Feedback Email: <guolilu@stu.xidian.edu.cn> 
 ``` 

# Step I. Install
  1. Download (or clone) the source code form https://github.com/Lilu-guo/Panaln  
  2. Compile the source code. (Note that you need to compile FM and WFA first)  

# Step II. Generate pangenome
 ```
 Usage: panaln combine -s <input.fasta> -v <input.vcf> -b <basename>  
 --------------  
 Specific:  -s STRING [required] reference genome  
            -v STRING [required] vcf file  
            -b STRING [required] basename  
            -c STRING [optional] context size (default:150)  
 Please use ABSOLUTE PATHs when specifying files.  
 --------------  
Feedback Email: <guolilu@stu.xidian.edu.cn>  
```
     
# Step III. Construct index
 ```
 Usage: panaln index -p <input.pan>  
 --------------  
 Specific:  -p STRING [required] pangenome file  
 Please use ABSOLUTE PATHs when specifying files.  
 --------------  
 Feedback Email: <guolilu@stu.xidian.edu.cn>  
 ```
  
# Step IV. Read alignment
 ```
 Usage: panaln align -x <index_basename> -f <input.fastq> -s <output.sam>  
 --------------  
 Specific:  -x STRING [required] basename  
            -f STRING [required] fastq file  
            -s STRING [required] smm file  
 Please use ABSOLUTE PATHs when specifying files.  
 --------------
 Feedback Email: <guolilu@stu.xidian.edu.cn>
 ```
  
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
   
# Feedback
 Please report bugs to Email: guolilu@stu.xidian.edu.cn if any questions or suggestions. 
 Your feedback and test cases are welcome.