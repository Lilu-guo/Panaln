# coding=utf-8
import time
import os
import sys
import subprocess

# read="/home/lab/gll/formatsnp/mason_vcf/mason_genome_10k.fa.2wgsim"                            #模拟，已知VCF
# aln1="/home/lab/gll/bwbble-master/test_data/aln.mason.vcf.10k.n"+sys.argv[1]
# sam1="/home/lab/gll/bwbble-master/test_data/sam.mason.vcf.10k.n"+sys.argv[1]

# read="/home/lab/gll/simPan/mismatch_snp_reads_100M/sim_ms_1.fa.2wgsim"                            #1.ms
# aln1="/home/lab/gll/bwbble-master/test_data/aln.sim_ms-1M.n"+sys.argv[1]
# sam1="/home/lab/gll/bwbble-master/test_data/sam.sim_ms-1M.n"+sys.argv[1]
# samacc="sam-ms-1M"

# # read="/home/lab/gll/simPan/snp_reads_100k/sim_s_1.fa.2wgsim"                                      #2.s
# # aln1="/home/lab/gll/bwbble-master/test_data/aln.sim_s_1-100k.n"+sys.argv[1]
# # sam1="/home/lab/gll/bwbble-master/test_data/sam.sim_s_1-100k.n"+sys.argv[1]
# # samacc="sam-s-100k"

# # read="/home/lab/gll/simPan/mismatch_reads_100k/sim_m_1.fa.2wgsim"                                 #3.m
# # aln1="/home/lab/gll/bwbble-master/test_data/aln.sim_m_1-100k.n"+sys.argv[1]
# # sam1="/home/lab/gll/bwbble-master/test_data/sam.sim_m_1-100k.n"+sys.argv[1]
# # samacc="sam-m-10k"

# # read="/home/lab/gll/simPan/reads_10k/sim_1.fa.2wgsim"                                            #4.exact
# # aln1="/home/lab/gll/bwbble-master/test_data/aln.sim_1-10k.n"+sys.argv[1]
# # sam1="/home/lab/gll/bwbble-master/test_data/sam.sim_1-10k.n"+sys.argv[1]
# # samacc="sam-exact-10k"


# # read="/home/lab/lab/gll/data_aln/bowtie2-master/example-hg19/reads/ERR194146-10000.fastq"      #真实
# # aln1="/home/lab/gll/bwbble-master/test_data/aln.ERR194146-10k.n"+sys.argv[1]
# # sam1="/home/lab/gll/bwbble-master/test_data/sam.ERR194146-10k.n"+sys.argv[1]

# ID="42"
# # ID="111516"
# # ID="1"           #正向read
# # ID="2"
# read="/home/lab/gll/simPan/chr1_ms_2M_101bp_0.2/"+ID+".fastq"   #单个read调试用       

read="/home/lab/gll/simPan/chr1_ms_2M_101bp_0.2/sim_ms_1.fa.2wgsim"                            #BIOINF论文R0用，20250327
aln1="/home/dell198/gll/vg/sam_panaln/2025.4.10.sam"                  #150bp中k用3,250bp中k用7
samacc="sam-101bp-0.2"
# read="/home/lab/gll/simPan/chr1_ms_10k_ccs_10k_Q20/sim_ms_1.fa.2wgsim"                            #BIOINF论文R0用，20250327
# aln1="/home/dell198/gll/vg/sam_panaln/chr1_ms_10k_ccs_10k_Q20.sam"                  #150bp中k用3,250bp中k用7
# samacc="sam-ccs-10k-Q20"

# read="/home/lab/gll/simPan/chr1_ms_2M_101bp_0.1/sim_ms_1.fa.2wgsim"                            #BIOINF论文R0用，20250327
# aln1="/home/dell198/gll/vg/sam_bwbble/chr1_ms_2M_101bp_0.1.aln.n"+sys.argv[1]                  #150bp中k用3,250bp中k用7
# sam1="/home/dell198/gll/vg/sam_bwbble/chr1_ms_2M_101bp_0.1.sam.n"+sys.argv[1]
# samacc="sam-101bp-0.1"
# read="/home/lab/gll/simPan/chr1_ms_2M_250bp_0.2/sim_ms_1.fa.2wgsim"                            #BIOINF论文R0用，20250327
# aln1="/home/dell198/gll/vg/sam_bwbble/chr1_ms_2M_250bp_0.2.aln.n"+sys.argv[1]                  #150bp中k用3,250bp中k用7
# sam1="/home/dell198/gll/vg/sam_bwbble/chr1_ms_2M_250bp_0.2.sam.n"+sys.argv[1]
# samacc="sam-250bp-0.2"
# read="/home/lab/gll/simPan/chr1_ms_2M_250bp_0.1/sim_ms_1.fa.2wgsim"                            #BIOINF论文R0用，20250327
# aln1="/home/dell198/gll/vg/sam_bwbble/chr1_ms_2M_250bp_0.1.aln.n"+sys.argv[1]                  #150bp中k用3,250bp中k用7
# sam1="/home/dell198/gll/vg/sam_bwbble/chr1_ms_2M_250bp_0.1.sam.n"+sys.argv[1]
# samacc="sam-250bp-0.1"

# read="/home/lab/gll/simPan/chr1_ms_10k_ccs_10k_Q20/sim_ms_1.fa.2wgsim"                             # BIOINF补充实验，2025.3.28，CCS_10k_Q20
# aln1="/home/dell198/gll/vg/sam_bwbble/chr1_ms_10k_ccs_10k_Q20.aln.n"+sys.argv[1]
# sam1="/home/dell198/gll/vg/sam_bwbble/chr1_ms_10k_ccs_10k_Q20.sam.n"+sys.argv[1]
# samacc="sam-ccs-10k-Q20"


exe1="/home/lab/gll/panaln/mg-aligner/panaln"
idx1="/home/lab/gll/panaln/mg-aligner/dirPan/guo"                      #2025.4.10 更新版本,使用更方便
# idx1="/home/lab/gll/panaln/mg-ref/1_bubble.fasta"                      #BIOINF论文R0用，20250327
# idx1="/home/lab/gll/bwbble-master/test_data/single_GRCh38.fasta"                               #无VCF


def nnc(str1):
    p=subprocess.Popen(str1,shell=True)
    print(str1)
    p.communicate()

def align():
    # cmd=exe1+" align -P -n "+sys.argv[1]+" "+idx1+" "+read+" "+aln1                         #预计算SA区间加速（2023.9.6）
    # cmd=exe1+" align -n "+sys.argv[1]+" "+idx1+" "+read+" "+aln1                            #Panaln V1
    # cmd=exe1+" align"+" "+idx1+" "+read+" "+aln1                                              #Panaln V2 (未使用PreSA加速, 在Illumina上精度更高)
    cmd=exe1+" align"+" -x "+idx1+" -f "+read+" -s "+aln1   #2025.4.10 更新版本,使用更方便
    start=time.time()
    nnc(cmd)
    end=time.time()
    take1=end-start
    print('running time: ########## %f s' %(take1))
    nnc("/home/lab/gll/protest/mycode/sam4 "+aln1)                                            #1.比对率
    print(' ')
    nnc("/home/lab/gll/bwbble-master/mg-ref/sam_pad "+aln1)                                   #pos从bubble转chr
    print(' ')
    nnc("/home/lab/gll/protest/mycode/"+ samacc +" " +aln1+".pan")                            #2.准确率
    return take1

if __name__=="__main__":
    print('\n')
    print("----------align  "+time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    take1=align()
    # print('\n')
    # print("----------sam  "+time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    # take2=sam()