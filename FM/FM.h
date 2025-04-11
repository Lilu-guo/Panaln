/*
 * @Author: your name
 * @Date: 2020-07-21 10:21:22
 * @LastEditTime: 2021-08-20 21:11:16
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /FM-Adaptive-master/FM.h
 */ 
/*============================================
# Filename: FM.h
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
  External class,you can build a hybrid and adaptive fm-index 
  for a document.and then you can know:
  
  How many times a pattern occs in the document using Counting,
  the reference parmater num will holds the times.
  
  All the positions where the pattern occs,and the times.
  the reference pointer pos will holds all the places,but
  REMEMBER that: prepare and clean the space of pos is your 
  duty.

  Extract a piece of the document using Extracting,the sequence
  will hold T[pos...pos+len-1].T denotes the original document.
  REMEMBER that: prepare and clean space for sequence is your
  duty too.

  Save and Load support serialize a FM object to a file or restore 
  a FM object from a file.

  GetN will tell you the size of the file in byte.
  SizeInByte will tell you the totle space needed for all operations in byte
  SizeInByte_count will tell you sapce only needed for Counting operation
=============================================*/
#ifndef FM_H
#define FM_H
#include"loadkit.h"
#include"savekit.h"
#include"ABS_WT.h"
#include"Huffman_WT.h"
#include"Balance_WT.h"
#include"Hutacker_WT.h"
#include"WT_Handle.h"
class FM
{
	public:
		FM(const char * filename,int blocksize,int r,int D,int shape); //声明（没函数体），没有对wt初始化，在后面实现中有wt初始化
		FM(char * filename); //保存参考序列名，为索引idx3动态命令用，2025.3.23
		~FM(){};
		FM(const FM & h):wt(h.wt){} //有函数体，有对wt的初始化
		FM& operator =(const FM&h){wt=h.wt;return *this;};
		
		void counting(const char *pattern,fm_int &num);
		fm_int * locating(const char *pattern,fm_int & num);
		unsigned char * extracting(fm_int pos,fm_int len);
		
		fm_int FM_Lookup(fm_int pos){return wt.Fm_Lookup(pos);};
		void Fmocc(fm_int l, fm_int r,unsigned char c ,fm_int *Left,fm_int * Right){wt.Fmocc(l,r,c,Left,Right);}; //为bwa，新增的方法
		void Fm16occ(fm_int l, fm_int r,fm_int *Left,fm_int * Right){wt.Fm16occ(l,r,Left,Right);};

		void Occ1p(unsigned char c, fm_int pos1, fm_int pos2, fm_int &occ1, fm_int &occ2){wt.Occ1p(c,pos1,pos2,occ1,occ2);};   //add gll
		void Occ8p(unsigned char c, fm_int pos1, fm_int pos2, fm_int *occ81, fm_int *occ82){wt.Occ8p(c,pos1,pos2,occ81,occ82);}; //add gll
		void Occ16p(fm_int pos1, fm_int pos2, fm_int *occ161, fm_int *occ162){wt.Occ16p(pos1,pos2,occ161,occ162);};               //add gll

		int load(const char * indexfile);
		int loadIdx();                                                                      //add gll
		int save(const char * indexfile);

		fm_int getN();
		fm_int getAlphabetSize();
		double getRuns();
		fm_int sizeInByte();
		fm_int sizeInByteForCount();
		fm_int sizeInByteForLocate();
		fm_int sizeInByteForExtract();
		fm_int sizeInByteForSAL();
		fm_int sizeInByteForRankL();
		fm_int sizeInByteFornewD();
		double compressRatio();
		double compressRatioForCount();
		double compressRatioForLocate();
		double compressRatioForExtract();
		double compressRatioForSal();
		double compressRatioForRankl();
		double compressRatioFornewD();
		fm_int TreePlainBlock(){return wt.TreePlainBlock();};
		fm_int TreeRLGBlock(){return wt.TreeRLGBlock();};
		fm_int TreeALL01Block(){return wt.TreeALL01Block();};
		fm_int TreeEF01Block(){return wt.TreeEF01Block();};
		fm_int TreePlainSize(){return wt.TreePlainSize();};
		fm_int TreeRLGSize(){return wt.TreeRLGSize();};
		fm_int TreeEF01Size(){return wt.TreeEF01Size();};
		fm_int MemSize(){return wt.MemSize();};
		fm_int NodeBlockSize(){return wt.NodeBlockSize();};
		fm_int CodeStyleSize(){return wt.CodeStyleSize();};

	private:
		WT_Handle wt; //属性（代理FM的功能）
};
#endif

