import sys
import time
import subprocess
import random
def nnc(str1):
    p=subprocess.Popen(str1,shell = True)
    p.communicate()
cntfiles = [
"/home/lab/lab/hzt/testfile/dna",
"/home/lab/lab/hzt/GRCh38.fna",
"/home/lab/lab/hzt/NA12877_10G_Reads",
]
locfiles = [
"/home/lab/lab/hzt/testfile/dna",
]
locG38files = [
"/home/lab/lab/hzt/GRCh38.fna",
"/home/lab/lab/hzt/NA12877_10G_Reads",
]
def cnttest(filename):
	cmd1 = "./fmcount -f {0} -b 256 -r 8 -d 32".format(filename)
	nnc(cmd1)
	cmd1 = "./fmcount -f {0} -b 256 -r 16 -d 32".format(filename)
	nnc(cmd1)
	cmd1 = "./fmcount -f {0} -b 256 -r 32 -d 32".format(filename)
	nnc(cmd1)
	cmd1 = "./fmcount -f {0} -b 256 -r 64 -d 32".format(filename)
	nnc(cmd1)
	return cmd1,
def loctest(filename):
	cmd1 = "./fmloc -f {0} -b 256 -d 32".format(filename)
	nnc(cmd1)
	cmd1 = "./fmloc -f {0} -b 256 -d 64".format(filename)
	nnc(cmd1)
	cmd1 = "./fmloc -f {0} -b 256 -d 128".format(filename)
	nnc(cmd1)
	cmd1 = "./fmloc -f {0} -b 256 -d 256".format(filename)
	nnc(cmd1)
	return cmd1,
def locG38test(filename):
	cmd1 = "./fmlocG38 -f {0} -b 256 -d 32".format(filename)
	nnc(cmd1)
	cmd1 = "./fmlocG38 -f {0} -b 256 -d 64".format(filename)
	nnc(cmd1)
	cmd1 = "./fmlocG38 -f {0} -b 256 -d 128".format(filename)
	nnc(cmd1)
	cmd1 = "./fmlocG38 -f {0} -b 256 -d 256".format(filename)
	nnc(cmd1)
	return cmd1,
if __name__ == "__main__":
	for f in cntfiles:
		cnttest(f)
	for f in locfiles:
		loctest(f)
	for f in locG38files:
		locG38test(f)











