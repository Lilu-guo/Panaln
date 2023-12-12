
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
		BitMap(unsigned long long int * bitbuff,fm_int bit_len,int level,int block_size=1024,int r=16,unsigned char label='\0',uchar ** tables=NULL);
		BitMap(){};
		BitMap(uchar ** tables):Z(tables[0]),R(tables[1]),R1(tables[2]),R3(tables[3]){}
		~BitMap();
		fm_int Rank(fm_int pos);
		fm_int Rank(fm_int pos,int &bit);
		void Rank(fm_int pos_left,fm_int pos_right,fm_int &rank_left,fm_int &rank_right);
		void Left(BitMap * left);
		BitMap * Left(){return left;};
		void Right(BitMap * right);
		BitMap * Right(){return right;};
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
		int GetBit(u64 * data,fm_int index);
		u16 Zeros(u16 x){return (Z[x>>8]==8)?Z[x>>8]+Z[(uchar)x]:Z[x>>8];}
		u64 GetBits(u64 * buff,fm_int index,int bits);
		void BitCopy1(u64 * temp,i64 &index,u64*data,i64 index1,int len);
		int GetZerosRuns(u64 * buff,fm_int &index);
		int GammaDecode(u64 * buff,fm_int &index);
		int GetRuns(u64 * data,fm_int &index,int &bit);
		vector<int> ConvertRLtopos(int *runs,int num,int firstbit,int& bit);
		int* ConvertRLtogap(vector<int> pos,int low);
		void GetGapG(vector<int > &pos);
		void Append_g(u64 * temp,fm_int &index,u64 value);
		void Append_f(u64 * temp,i64 &index,u64 value,u32 maxrl);
		void BitCopy(u64 * temp,fm_int &index,u64 value);
		void RL_Rank(u64 * buff,fm_int &index,int bits_left,int bits_right,fm_int &rank_left,fm_int &rank_right,int rl_type);
		int  RL_Rank(u64 * buff,fm_int &index,int bits_num,int rl_type);
		int  RL_Rank(u64 * buff,fm_int &index,int bits_num,int rl_type,int &bit);
		void GapG_Rank(u64 * buff,fm_int &index,int bits_left,int bits_right,fm_int &rank_left,fm_int &rank_right,int rl_type,int cbl);
		int  GapG_Rank(u64 * buff,fm_int &index,int bits_num,int rl_type,int cbl);
		int  GapG_Rank(u64 * buff,fm_int &index,int bits_num,int rl_type,int cbl,int &bit);
		i64 GetCurrentRankvalue(i64 block_anchor,i64 superblock_anchor,i64 offset1,i64 preoffset,i64 radio);
		void EF_Rank(u64 * buff,i64 &index,int bits_left,int bits_right,i64 &rank_left,i64 &rank_right,int ef_type,u64 elemnums,int currentBlockLen);
		u64  EF_Rank(u64 * buff,i64 &index,int bits_num,int ef_type,u64 elemnums,int currentBlockLen);
		u64  EF_Rank(u64 * buff,i64 &index,int bits_num,int ef_type,u64 elemnums,int currentBlockLen,int &bit);
		int EF0_Bit(u64 * buff,i64 &index,int bits,u64 currentblocklen,int elemnums);
		int EF1_Bit(u64 * buff,i64 & index,int bits,u64 currentblocklen,int elemnums);
		int RL0_Rank(u64 * buff,fm_int &index,int bits_num);
		int RL0_Bit(u64 * buff,fm_int &index,int bits);
		int RL0_Rank(u64 * buff,fm_int &index,int bits,int &bit);
		int RL1_Rank(u64 * buff,fm_int &index,int bits);
		int RL1_Bit(u64 * buff,fm_int & index,int bits);
		int RL1_Rank(u64 * buff,fm_int &index,int bits,int &bit);
		int GapG0_Bit(u64 * buff,fm_int &index,int bits,int cbl);
		int GapG1_Bit(u64 * buff,fm_int & index,int bits,int cbl);
		void Plain_Rank(u64 *buff,fm_int &index,int bits_left,int bits_right,fm_int &rank_left,fm_int&rank_right);
		int Plain_Rank(u64 * buff,fm_int &index,int bits);
		int Plain_Bit(u64 * buff,fm_int &index,int bits);
		int Plain_Rank(u64 * buff,fm_int &index,int bits,int & bit);
		int level;
		unsigned char label;
		unsigned long long int * data;
		fm_int bitLen;
		fm_int memorysize;
		int block_size;
		int r;
		int super_block_size;
		BitMap * left;
		BitMap * right;
		InArray *superblock;
		InArray *block;
		InArray *coding_style;
		SBArray *newblock;
		
		i64 plainblock;
		i64 all01block;
		i64 rlgblock;
		i64 ef01block;
		i64 plainsize;
		i64 rlgsize;
		i64 ef01size;
		unsigned long long int * buff;
};

#endif







