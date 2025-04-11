/*============================================
# Filename: BitMap.h
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#ifndef WT_NODE_H
#define WT_NODE_H
#include<string.h>
#include"loadkit.h"
#include"savekit.h"
#include"InArray.h"
#include"BaseClass.h"
#include"newDmap.h"
#include<math.h>
#include<vector>
#include<iostream>
using namespace std;
class BitMap
{
	public:
		BitMap(unsigned long long int * bitbuff,fm_int bit_len,int level,int block_size=1024,int r=16,unsigned char label='\0',uchar ** tables=NULL); //构建，填充小波树时用 ABS_FM::FullFillWTNode
		//bit_len:0,1串的实际长度，单位bit
		//level:层数
		//block_size:块大小
		//label:当前节点代表的字符.只有叶节点的label域又意义.
		
		BitMap(){}; //从索引，读入小波树时用 ABS_FM::LoadWTTree
		BitMap(uchar ** tables):Z(tables[0]),R(tables[1]),R1(tables[2]),R3(tables[3]){}
		~BitMap();

		fm_int Rank(fm_int pos);
		fm_int Rank(fm_int pos,int &bit);
		void Rank(fm_int pos_left,fm_int pos_right,fm_int &rank_left,fm_int &rank_right);


		void Left(BitMap * left);
		//设置左孩子
		
		BitMap * Left(){return left;};
		//返回左孩子

		void Right(BitMap * right);
		//...

		BitMap * Right(){return right;};
		//...

		unsigned char Label();
		int Load(loadkit & s);
		int Save(savekit & S);
		fm_int SizeInByte();
		fm_int PlainBlock();
		fm_int RLGBlock();
		fm_int ALL01Block();
		fm_int EF01Block();
		fm_int PlainSize();
		fm_int RLGSize();
		fm_int EF01Size();

		fm_int MemSize();
		fm_int NodeBlockSize();
		fm_int CodeStyleSize();
	private:
		uchar* Z;
		uchar *R;
		uchar *R1;
		uchar *R3;
		BitMap(const BitMap &);
		BitMap & operator =(const BitMap& right);
		void Coding();
		void Coding2();
		//得到存储在data中的0,1串中的第index位
		int GetBit(u64 * data,fm_int index);

		u16 Zeros(u16 x){return (Z[x>>8]==8)?Z[x>>8]+Z[(uchar)x]:Z[x>>8];}
		//从buff保存的0.1串中，由index位置开始，返回后续bits位表示的
		//数值.
		u64 GetBits(u64 * buff,fm_int index,int bits);
		void BitCopy1(u64 * temp,i64 &index,u64*data,i64 index1,int len);

		//得到0的runs.
		int GetZerosRuns(u64 * buff,fm_int &index);
		//gamma解
		int GammaDecode(u64 * buff,fm_int &index);
 		
		//得到0,1串中的runs长度，bit标示该runs是针对谁的
		int GetRuns(u64 * data,fm_int &index,int &bit);
		vector<int> ConvertRLtopos(int *runs,int num,int firstbit,int& bit);
		int* ConvertRLtogap(vector<int> pos,int low);
		void GetGapG(vector<int > &pos);
		//index
		void Append_g(u64 * temp,fm_int &index,u64 value);
		void Append_f(u64 * temp,i64 &index,u64 value,u32 maxrl);
		//把u64类型的value拷贝到data串的index处.
		void BitCopy(u64 * temp,fm_int &index,u64 value);

		//返回rl0编码的串中，由index位置开始，长度位bits
		//内的1的个数.
		void RL_Rank(u64 * buff,fm_int &index,int bits_left,int bits_right,fm_int &rank_left,fm_int &rank_right,int rl_type);
		int  RL_Rank(u64 * buff,fm_int &index,int bits_num,int rl_type);
		int  RL_Rank(u64 * buff,fm_int &index,int bits_num,int rl_type,int &bit);
		//GapG编码
		void GapG_Rank(u64 * buff,fm_int &index,int bits_left,int bits_right,fm_int &rank_left,fm_int &rank_right,int rl_type,int cbl);
		int  GapG_Rank(u64 * buff,fm_int &index,int bits_num,int rl_type,int cbl);
		int  GapG_Rank(u64 * buff,fm_int &index,int bits_num,int rl_type,int cbl,int &bit);
		///
		i64 GetCurrentRankvalue(i64 block_anchor,i64 superblock_anchor,i64 offset1,i64 preoffset,i64 radio);

		void EF_Rank(u64 * buff,i64 &index,int bits_left,int bits_right,i64 &rank_left,i64 &rank_right,int ef_type,u64 elemnums,int currentBlockLen);
		u64  EF_Rank(u64 * buff,i64 &index,int bits_num,int ef_type,u64 elemnums,int currentBlockLen);
		u64  EF_Rank(u64 * buff,i64 &index,int bits_num,int ef_type,u64 elemnums,int currentBlockLen,int &bit);
		int EF0_Bit(u64 * buff,i64 &index,int bits,u64 currentblocklen,int elemnums);
		//返回容量编码的串中，由index位置开始，bits位内的1的个数.
		int EF1_Bit(u64 * buff,i64 & index,int bits,u64 currentblocklen,int elemnums);


		int RL0_Rank(u64 * buff,fm_int &index,int bits_num);
		int RL0_Bit(u64 * buff,fm_int &index,int bits);
		int RL0_Rank(u64 * buff,fm_int &index,int bits,int &bit);
		//返回容量编码的串中，由index位置开始，bits位内的1的个数.
		int RL1_Rank(u64 * buff,fm_int &index,int bits);
		int RL1_Bit(u64 * buff,fm_int & index,int bits);
		int RL1_Rank(u64 * buff,fm_int &index,int bits,int &bit);
		//Gapg
		int GapG0_Bit(u64 * buff,fm_int &index,int bits,int cbl);
		int GapG1_Bit(u64 * buff,fm_int & index,int bits,int cbl);
		//buff从index位置开始是直接存储的，从index位置开始，bits
		//位内有几个1.
		void Plain_Rank(u64 *buff,fm_int &index,int bits_left,int bits_right,fm_int &rank_left,fm_int&rank_right);
		int Plain_Rank(u64 * buff,fm_int &index,int bits);
		int Plain_Bit(u64 * buff,fm_int &index,int bits);
		int Plain_Rank(u64 * buff,fm_int &index,int bits,int & bit);

		
		int level;//该串的层数.
		
		unsigned char label;
		//只有叶节点又意义，表示该节点代表的字符

		unsigned long long int * data;
		//0,1串的压缩存储体.

		fm_int bitLen;
		//0,1串的长度，单位bit。

		fm_int memorysize;
		int block_size;
		int r;
		int super_block_size;

		BitMap * left;
		BitMap * right;

		InArray *superblock;//超快偏移量
		InArray *block;//块的偏移量
		InArray *coding_style;//每个块的编码方案.0:plain, 1:RLG0, 2:RLG1;
		SBArray *newblock;
		
		i64 plainblock;
		i64 all01block;
		i64 rlgblock;
		i64 ef01block;
		i64 plainsize;
		i64 rlgsize;
		i64 ef01size;


		//这是个工作变量.
		unsigned long long int * buff;
};

#endif







