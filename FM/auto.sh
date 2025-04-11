#!/bin/bash
###
 # @Author: your name
 # @Date: 2021-04-02 16:51:32
 # @LastEditTime: 2021-04-28 11:42:32
 # @LastEditors: Please set LastEditors
 # @Description: In User Settings Edit
 # @FilePath: /FM-EF-d/do.sh
### 
echo "auto is carrying!!!"
make clean
make
echo "---------------------------------------------------"
echo "128:"
./fm -f ../testfile/dna -b 128
echo "128 dna is OK!!!"
./fm -f ../testfile/pitches -b 128
echo "128 pitches is OK!!!"
./fm -f ../testfile/proteins -b 128
echo "128 proteins is OK!!!"
./fm -f ../testfile/english -b 128
echo "128 english is OK!!!"
./fm -f ../testfile/sources -b 128
echo "128 sources is OK!!!"
./fm -f ../testfile/dblp.xml -b 128
echo "128 dblp.xml is OK!!!"
./fm -f ../testfile/para -b 128
echo "128 para is OK!!!"
./fm -f ../testfile/influenza -b 128
echo "128 influenza is OK!!!"
./fm -f ../testfile/world_leaders -b 128
echo "128 world_leaders is OK!!!"
./fm -f ../testfile/kernel -b 128
echo "128 kernel is OK!!!"
echo "---------------------------------------------------"
echo "256:"
./fm -f ../testfile/dna -b 256
echo "256 dna is OK!!!"
./fm -f ../testfile/pitches -b 256
echo "256 pitches is OK!!!"
./fm -f ../testfile/proteins -b 256
echo "256 proteins is OK!!!"
./fm -f ../testfile/english -b 256
echo "256 english is OK!!!"
./fm -f ../testfile/sources -b 256
echo "256 sources is OK!!!"
./fm -f ../testfile/dblp.xml -b 256
echo "256 dblp.xml is OK!!!"
./fm -f ../testfile/para -b 256
echo "256 para is OK!!!"
./fm -f ../testfile/influenza -b 256
echo "256 influenza is OK!!!"
./fm -f ../testfile/world_leaders -b 256
echo "256 world_leaders is OK!!!"
./fm -f ../testfile/kernel -b 256
echo "256 kernel is OK!!!"
echo "---------------------------------------------------"
echo "512:"
./fm -f ../testfile/dna -b 512
echo "512 dna is OK!!!"
./fm -f ../testfile/pitches -b 512
echo "512 pitches is OK!!!"
./fm -f ../testfile/proteins -b 512
echo "512 proteins is OK!!!"
./fm -f ../testfile/english -b 512
echo "512 english is OK!!!"
./fm -f ../testfile/sources -b 512
echo "512 sources is OK!!!"
./fm -f ../testfile/dblp.xml -b 512
echo "512 dblp.xml is OK!!!"
./fm -f ../testfile/para -b 512
echo "512 para is OK!!!"
./fm -f ../testfile/influenza -b 512
echo "512 influenza is OK!!!"
./fm -f ../testfile/world_leaders -b 512
echo "512 world_leaders is OK!!!"
./fm -f ../testfile/kernel -b 512
echo "512 kernel is OK!!!"