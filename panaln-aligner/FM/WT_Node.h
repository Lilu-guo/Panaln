#ifndef WT_NODE_H
#define WT_NODE_H
#include<string.h>
#include"loadkit.h"
#include"savekit.h"
#include"InArray.h"
#include"BaseClass.h"
#include<math.h>
#include<iostream>
using namespace std;
class WT_Node
{
	public:
		WT_Node(unsigned long long int * bitbuff,int bit_len,int level,int block_size=1024,unsigned char label='\0',uchar ** tables=NULL);
		WT_Node(){};
		WT_Node(uchar ** tables):zerostable(tables[0]),R(tables[1]){}
		~WT_Node();
		int Rank(int pos);
		int Rank(int pos,int &bit);
		void Left(WT_Node * left);
		WT_Node * Left(){return left;};
		void Right(WT_Node * right);
		WT_Node * Right(){return right;};
		unsigned char Label();
		int Load(loadkit & s);
		int Save(savekit & S);
		int SizeInByte();
	private:
		uchar* zerostable;
		uchar *R;
		WT_Node(const WT_Node &);
		WT_Node & operator =(const WT_Node& right);
		void Coding();
		int GetBit(u64 * data,int index);
		u64 GetBits(u64 * buff,int &index,int bits);
		int GetZerosRuns(u64 * buff,int &index);
		int FixedDecode(u64 * buff,int &index);
		int GammaDecode(u64 * buff,int &index);
		int GetRuns(u64 * data,int &index,int &bit);
		void Append_g(u64 * temp,int &index,u32 value);
		void BitCopy(u64 * temp,int &index,u64 value);
		int RL0_Rank(u64 * buff,int &index,int bits_num);
		int RL0_Bit(u64 * buff,int &index,int bits);
		int RL0_Rank(u64 * buff,int &index,int bits,int &bit);
		int RL1_Rank(u64 * buff,int &index,int bits);
		int RL1_Bit(u64 * buff,int & index,int bits);
		int RL1_Rank(u64 * buff,int &index,int bits,int &bit);
		int Plain_Rank(u64 * buff,int &index,int bits);
		int Plain_Bit(u64 * buff,int &index,int bits);
		int Plain_Rank(u64 * buff,int &index,int bits,int & bit);		
		int level;
		unsigned char label;
		unsigned long long int * data;
		int bitLen;
		int memorysize;
		int block_size;
		int block_width;
		int super_block_size;
		WT_Node * left;
		WT_Node * right;
		InArray *superblock;
		InArray *block;
		InArray *coding_style;
		unsigned long long int * buff;
};

#endif
