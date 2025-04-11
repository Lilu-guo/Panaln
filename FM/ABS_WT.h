/*============================================
# Filename: ABS_WT.h
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#ifndef ABS_FM_H
#define ABS_FM_H
#define CHAR_SET_SIZE 256
#define CODE_MAX_LEN 256
#include<iomanip>
#include<string.h>
#include"BitMap.h"
#include"InArray.h"
#include"loadkit.h"
#include"savekit.h"
#include"divsufsort64.h"
#include"../io.h"

class ABS_FM
{
	public:
		ABS_FM(const char * filename,int block_size=256,int r=16,int D=32);
		ABS_FM(char * filename){
			this->filename=filename; //保存参考序列名，为索引idx3动态命令用，2025.3.23
		};
		virtual ~ABS_FM();
		void Counting(const char * partten,fm_int &num);
		fm_int * Locating(const char * pattern,fm_int &num);
		fm_int * Locatingnew(const char * pattern,fm_int &num);
		unsigned char* Extracting(fm_int pos,fm_int len);

		void Fmocc(fm_int l, fm_int r,unsigned char c ,fm_int *Left,fm_int* Right);
		// void Fm4occ(fm_int l, fm_int r,fm_int *Left,fm_int* Right);
		void Fm16occ(fm_int l, fm_int r,fm_int *Left,fm_int* Right);

		void Occ1p(unsigned char c, fm_int pos1, fm_int pos2, fm_int &occ1, fm_int &occ2);   //add gll
		void Occ8p(unsigned char c, fm_int pos1, fm_int pos2, fm_int *occ81, fm_int *occ82); //add gll
		void Occ16p(fm_int pos1, fm_int pos2, fm_int *occ161, fm_int *occ162);               //add gll

		fm_int FM_Lookup(fm_int i);

		int Load(loadkit & s);
		int loadIdx();                                                                       //add gll
		int Save(savekit & s);
		int BuildTree(); //WT_Handle初始化时调用
		void buildEFtable(u64 i);

		int GetAlphabetsize(){return alphabetsize;}
		fm_int GetN(){return n;}                                                             //头文件中，直接返回ABS_FM属性n
		fm_int SizeInByte();
		fm_int SizeInByte_count();
		fm_int SizeInByte_locate();
		fm_int SizeInByte_extract();
		fm_int SizeInByte_SAL();
		fm_int SizeInByte_RankL();
		fm_int SizeInByte_newD();
		fm_int TreePlainBlock(BitMap * r);
		fm_int TreeRLGBlock(BitMap * r);
		fm_int TreeALL01Block(BitMap * r);
		fm_int TreeEF01Block(BitMap * r);
		fm_int TreePlainSize(BitMap * r);
		fm_int TreeRLGSize(BitMap * r);
		fm_int TreeEF01Size(BitMap * r);
		fm_int TreeMemSize(BitMap * r);
		fm_int TreeNodeBlockSize(BitMap * r);
		fm_int TreeCodeStyleSize(BitMap * r);
		fm_int PlainBlock(){return TreePlainBlock(root);};
		fm_int RLGBlock(){return TreeRLGBlock(root);};
		fm_int ALL01Block(){return TreeALL01Block(root);};
		fm_int EF01Block(){return TreeEF01Block(root);};
		fm_int PlainSize(){return TreePlainSize(root);};
		fm_int RLGSize(){return TreeRLGSize(root);};
		fm_int EF01Size(){return TreeEF01Size(root);};
		fm_int MemSize(){return TreeMemSize(root);};
		fm_int NodeBlockSize(){return TreeNodeBlockSize(root);};
		fm_int CodeStyleSize(){return TreeCodeStyleSize(root);};
		double GetRuns();

	protected:
		BitMap *root, *lroot, *rroot; //lroot:ACGT, rroot:IUPAC+$
		uchar * Z;
	 	uchar * R;
		uchar * R1;
		uchar * R3;
		Dmap *newD;                  //topSNP bitvector
		Dmap *newd;                  //指采样的bitvector
		unsigned char * T;
		unsigned char * bwt;
		char * filename;
		fm_int	sa0_index;          //BWT $ pos
		fm_int  n;                  //file size in byte
		int block_size;             //the block_size of wavelettree's node.default 256
		int D=32;                   //step of SAL.default 32
		int r;                      //the super_block_size of wavelettree's node,super_block_size=r*block_size,default 16
		InArray * SAL;              //the sampling of SA array
		InArray * RankL;            //the sampling of the Rank Array
		double runs;

		
		bool charMap[256];
		// if charMap[i] is true,char[i] is a member of alphabet
		
		fm_int * C;                 //Cumulative Frequency sum C表
		int *code;                                                                 //???

		fm_int charFreq[CHAR_SET_SIZE];
		int alphabetsize;
		char codeTable[CHAR_SET_SIZE][CODE_MAX_LEN];                               //小波树码表，[字符集大小][码长]
		
		void Inittable();                                                          //初始化解码加速表
		void mycodetable();
		void huffmanTreecode();
		fm_int Occ(fm_int &rank,unsigned char &label,fm_int pos);
		fm_int Occ(unsigned char c,fm_int pos);
		fm_int Occ(unsigned char c, fm_int pos, BitMap *root);                     //add gll
		void Occ4(fm_int pos, fm_int *occ4);                                       //add gll
		void Occ11(fm_int pos, fm_int *occ11);                                     //add gll
		void Occ1pAux(unsigned char c, fm_int pos1, fm_int pos2, fm_int &occ1, fm_int &occ2);           //add gll
		void Occ4p(fm_int pos1, fm_int pos2, fm_int *occ41, fm_int *occ42);        //add gll
		void Occ11p(fm_int pos1, fm_int pos2, fm_int *occ111, fm_int *occ112);     //add gll
		void Occ(unsigned char c,fm_int pos_left,fm_int pos_right,fm_int &rank_left,fm_int &rank_right);
		void Occ4(fm_int pos_left,fm_int pos_right,fm_int *rank_left,fm_int *rank_right);
		fm_int Occ1(fm_int i); //add gll
		fm_int LF(fm_int i);
		fm_int myLF(fm_int i);
		unsigned char myL(fm_int i);
		void DrawBackSearch(const char * pattern,fm_int &Left,fm_int & Right);
		fm_int Lookup(fm_int i);
		void Getpos(fm_int Left,fm_int Right,fm_int *&pos);

		virtual int TreeCode(){return -1;}; //虚函数（由具体树形类处理）

		int BWT(unsigned char * T,saidx64_t * SA,unsigned char * bwt,fm_int n);
		BitMap * CreateWaveletTree(unsigned char * bwt,fm_int n);
		BitMap * FullFillWTNode(unsigned char * bwt,fm_int len,int level);
		int DestroyWaveletTree();
		int blog(fm_int);
		unsigned char * Getfile(const char * filename);

		int SaveNodePosition(BitMap *,int, savekit &);
		int SaveNodeData(BitMap *,savekit &s);
		int SaveWTTree(savekit &s);
		int SaveWTTree(savekit &s, BitMap *root);                      //add gll
		int LoadWTTree(loadkit &s, uchar **tables=NULL);
		int LoadWTTree(loadkit &s, BitMap* &root, uchar **tables=NULL);           //add gll
		int TreeNodeCount(BitMap * root);

		fm_int TreeSizeInByte(BitMap * r);
		
	//	int TreeSizeInByte();
		friend int GammaDecode(u64 * buff,int & index,ABS_FM * t);
		friend int Zeros(u16 x,ABS_FM * t);


};
#endif

