
#include"BitMap.h"
#include<math.h>
#include<iostream>
using namespace std;
#define lookuptable

inline int popcnt(unsigned long long int x)
{
	x = x -((x & 0xAAAAAAAAAAAAAAAA)>>1);
	x = (x & 0x3333333333333333)+((x>>2) & 0x3333333333333333);
	x =((x+(x>>4)) & 0x0F0F0F0F0F0F0F0F);
	return (x*0x0101010101010101)>>56;
}

int blog(fm_int x)
{
	int ans = 0;
	while(x>0)
	{
		ans++;
		x=x>>1;
	}
	return ans;
}

BitMap::BitMap(unsigned long long int * bitbuff,fm_int bit_len,int level,int block_size,int r,unsigned char label,uchar ** tables)
{
	this->data = bitbuff;
	this->bitLen = bit_len;
	this->memorysize = 0;
	this->level = level;
	this->block_size = block_size;
	this->r=r;
	this->super_block_size = r*block_size;
	this->label = label;
	left=NULL;
	right=NULL;
	superblock =NULL;
	block=NULL;
	newblock = NULL;
	coding_style=NULL;
	R=NULL;
	Z=NULL;
	
	if(data!=NULL)
	{
		this->Z = tables[0];
		this->R = tables[1];
		this->R1 = tables[2];
		this->R3 = tables[3];
		Coding();
		buff =NULL;
	}

}

fm_int BitMap::SizeInByte()
{

	fm_int size = 0;
	if(data!=NULL)
	{
		size+= newblock->GetMemorySize();
		size+= memorysize;
	}
	return size;
}

fm_int BitMap::PlainBlock()
{

	fm_int size = 0;
	if(data!=NULL)
		size += plainblock;
	return size;
}

fm_int BitMap::RLGBlock()
{

	fm_int size = 0;
	if(data!=NULL)
		size += rlgblock;
	return size;
}

fm_int BitMap::ALL01Block()
{

	fm_int size = 0;
	if(data!=NULL)
		size += all01block;
	return size;
}

fm_int BitMap::EF01Block()
{

	fm_int size = 0;
	if(data!=NULL)
		size += ef01block;
	return size;
}

fm_int BitMap::PlainSize()
{

	fm_int size = 0;
	if(data!=NULL)
		size += plainsize;
	return size;
}

fm_int BitMap::RLGSize()
{

	fm_int size = 0;
	if(data!=NULL)
		size += rlgsize;
	return size;
}

fm_int BitMap::EF01Size()
{

	fm_int size = 0;
	if(data!=NULL)
		size += ef01size;
	return size;
}


fm_int BitMap::MemSize()
{

	fm_int size = 0;
	if(data!=NULL)
		size += memorysize * 8;
	return size;
}
fm_int BitMap::NodeBlockSize()
{

	fm_int size = 0;
	if(data!=NULL)
	{
//		size += superblock->GetMemorySize()*8;
//		size += block->GetMemorySize()*8;
		size += newblock->GetMemorySize();
	}
	return size;
}
fm_int BitMap::CodeStyleSize()
{

	fm_int size = 0;
	if(data!=NULL)
		size += coding_style->GetMemorySize()*8;
	return size;
}
void BitMap::Coding()
{
	fm_int idxnums = 0, index = 0;
	int step1 = block_size*this->r;
	int step2 = block_size;
	newblock = new SBArray(bitLen/step1+4,blog(bitLen),this->r,blog(step1-step2));
	newblock->SetValue(0,0,1);
	newblock->SetValue(0,0,0);
	fm_int bitnums;
	if(bitLen % 64 == 0)
		bitnums = bitLen / 64;
	else
		bitnums = bitLen / 64 + 1;
	u64 *temp = new u64[bitnums + 1];
	memset(temp,0,(bitnums+1)*8);
	memcpy(temp,data,bitnums * 8);
	delete data;
	data = NULL;
	fm_int rank = 0, pre_rank=0;

	for(;idxnums < bitnums;++idxnums)
	{
		rank += popcnt(temp[idxnums]);
		index += 64;
		if(index >= bitLen) {index -= 64;break;}
		if(index % step1 == 0)
		{
			pre_rank = rank;
			newblock->SetValue(index / step1,pre_rank,1);
//			_sb[sbidx++] = pre_rank;
		}
		if(index % step2 ==0)
		{
			newblock->SetValue(index / step2,rank - pre_rank,0);
//			_b[bidx++] = rank - pre_rank;
		}
	}
	this->memorysize = bitnums * 8;
	data = new u64[bitnums];
	memset(data,0, bitnums * 8);
	memcpy(data,temp,bitnums * 8);
	delete temp;
}

void BitMap::Coding2()
{
	fm_int u64Len =0;
	if(bitLen%64 == 0)
		u64Len = bitLen/64;
	else
		u64Len = bitLen/64+1;
	u64 * temp = new u64[u64Len];
	memset(temp,0,u64Len*8);
	
	fm_int index = 0;
	int step1 = block_size*this->r;
	int step2 = block_size;
	superblock = new InArray(2*(bitLen/step1)+4,blog(bitLen));
	block      = new InArray(2*(bitLen/step2)+4,blog(step1-step2));
	coding_style      = new InArray(bitLen/step2+1,3);
	u64 rank=0;
	u64 space=0;
	u64 palinspace=0;
	u64 rlgspace=0;
	u64 bits =0;
	u64 firstbit;
	int rl_g=0;
	int runs = 0;
	int bit=0;
	int * runs_tmp = new int[block_size];
	u64 k=0;
	i64 index2=0;
	plainblock = 0;
	all01block = 0;
	rlgblock = 0;
	ef01block = 0;
	plainsize = 0;
	rlgsize = 0;
	ef01size = 0;
	u64 pre_rank=0;
	u64 pre_space =0 ;
	superblock->SetValue(0,0);
	superblock->SetValue(1,0);
	block->SetValue(0,0);
	block->SetValue(1,0);

	while(index < bitLen)
	{
		if(index == bitLen)
			break;
		rl_g = 0;
		bits = 0;
		firstbit = 0;
		runs = 0;
		firstbit = GetBit(data,index);
		memset(runs_tmp,0,block_size*4);
		k=0;
 		runs=0;
		while(bits < block_size && index < bitLen)
		{

			runs = GetRuns(data,index,bit);
			bits = bits +runs;
			if(bit ==1)
				rank=rank+runs;
			runs_tmp[k] = runs;
			k++;
		}
		
		if(bits > block_size)
		{
			int step =0;
			index = index -(bits - block_size);
			step = block_size+runs-bits;
			if(bit ==1)
				rank = rank -runs+step;
			runs_tmp[k-1] = step;
		}
		for(int i=0;i<k;i++){
			rl_g = rl_g + 2*blog(runs_tmp[i])-1;
		}
		int codestyle=0;
		int gapglen = block_size;
		vector<int> pos;
		int gapgtype=0;
		int tmplen = 0;
		if(k!=1)
		{
			pos=ConvertRLtopos(runs_tmp,k,firstbit,gapgtype);
			GetGapG(pos);
			for(int i = 0;i < pos.size();++i)
				tmplen = tmplen + 2 * blog(pos[i]) - 1;
			gapglen = tmplen;
		}

		int thred=20;
		int len = min(rl_g,block_size-thred);
		codestyle=gapglen<len?1:0;
		if(len == (block_size-thred) || index == bitLen)
		{
			coding_style->SetValue((index-1)/block_size,2);
			plainblock++;
			i64 index3=0;
			if(index == bitLen)
			{
				space = space + bits;
				plainsize += bits;
				palinspace+=bits;
				index3=index-bits;
				len=bits;
			}
			else
			{
				space = space + block_size;
				plainsize += block_size;
				palinspace+=block_size;
				index3=index-block_size;
				len=block_size;
			}
			BitCopy1(temp,index2,data,index3,len);
		}
		else if(k==1)
		{
			all01block++;
			if(firstbit==0)
				coding_style->SetValue((index-1)/block_size,3);
			else
				coding_style->SetValue((index-1)/block_size,4);
			space = space +0;
		}

		else if(codestyle==0)
		{
			rlgblock++;
			if(firstbit == 0)
				coding_style->SetValue((index-1)/block_size,0);
			else
				coding_style->SetValue((index-1)/block_size,1);
			space =space + rl_g;
			rlgspace+=rl_g;
			rlgsize += rl_g;
			for(int i=0;i<k;i++)
			{
				Append_g(temp,index2,runs_tmp[i]);
			}
		}
		else{
			ef01block++;
			if(gapgtype==0){
				coding_style->SetValue((index-1)/block_size,5);
			}
			else{
				coding_style->SetValue((index-1)/block_size,6);
			}
			space=space+gapglen;
			ef01size += gapglen;
			for(int i=0;i<pos.size();i++)
			{
				Append_g(temp,index2,pos[i]);
			}
		}
		if(index % step1 == 0)
		{
			pre_rank = rank;
			superblock ->SetValue(2*(index/step1),pre_rank);

			pre_space = space;
			superblock->SetValue(2*(index/step1)+1,pre_space);
		}
		if(index % step2 ==0)
		{
			block->SetValue(2*(index/step2),rank - pre_rank);
			block->SetValue(2*(index/step2)+1,space - pre_space);
		}
	}

	delete [] runs_tmp;
	int u64_len_real = 0;
	if(space % 64==0)
		u64_len_real = space /64+1;
	else
		u64_len_real = space /64 +1+1;
	
	this->memorysize = u64_len_real*8;
	delete [] data;
	data = new u64[u64_len_real];

	memset(data,0,u64_len_real*8);
	memcpy(data,temp,(u64_len_real-1)*8);
	delete [] temp;
}

void BitMap::BitCopy1(u64 * temp,i64 &index,u64*data,i64 index1,int len)
{
	while(len>=64){
		u64 value=GetBits(data,index1,64);
		BitCopy(temp,index,value);
		len-=64;
		index1+=64;
	}
	if(len!=0){
		u64 value=GetBits(data,index1,len);
		if(index%64!=0){
			int first = index % 64;
			if(first+len<65){
				int overloop=64-(first+len);
				value=(value<<overloop);
				temp[index/64]=(temp[index/64]|value);
			}
			else{
				first=64-first;
				int second=len-first;
				int third=64-second;
				temp[index/64]=(temp[index/64]|(value>>second));
				temp[index/64+1]=temp[index/64+1]|(value<<third);
			}
		}
		else{
			int overloop=64-len;
			value=(value<<overloop);
			temp[index/64]=value;
		}
		index = index+len;
	}
}

BitMap::~BitMap()
{
	if(left)
		delete left;
	if(right)
		delete right;
	delete [] data;
	delete newblock;
}

vector<int> BitMap::ConvertRLtopos(int *runs,int num,int firstbit,int &bit){
	int k=1;
	int flag=firstbit;
	vector<int> pos_of_ones,pos_of_zeros;
	for(int i=0;i<num;i++){
		for(int j=0;j<runs[i];j++){
			if(flag==1){
				pos_of_ones.push_back(k);
			}else{
				pos_of_zeros.push_back(k);
			}
			k++;
		}
		flag=1-flag;
	}
	if(pos_of_zeros.size()<pos_of_ones.size()){
		bit=0;
		return pos_of_zeros;
	}else{
		bit=1;
		return pos_of_ones;
	}
}

void BitMap::GetGapG(vector<int > &pos)
{
	for(int i = pos.size() - 1;i > 0;--i)
		pos[i] -= pos[i - 1];
}

int* BitMap::ConvertRLtogap(vector<int> pos,int low){
	size_t num=pos.size();
	int *gap=new int[num];
	for(int i=num-1;i>0;i--){
		gap[i]=(pos[i]>>low)-(pos[i-1]>>low);
	}
	gap[0]=pos[0]>>low;
	return gap;
}

fm_int BitMap::Rank(fm_int pos,int & bit)
{
	if(pos < 0 || pos > bitLen)
	{
		cerr<<pos<<"  "<<bitLen<<" BitMap::Rank(i64, int&) error parmater"<<endl;
		exit(0);
	}
	if((pos+1)%block_size!=0)
	{
		i64 block_anchor = (pos+1)/block_size;
		i64 superblock_anchor = (pos+1)/super_block_size;
		i64 rank_base = newblock->GetValue(superblock_anchor,1) + newblock->GetValue(block_anchor,0);
		i64 offset = block_anchor * block_size;
		buff = data +(offset>>6);
		int overloop = (pos+1)%block_size ;
		i64 index = (offset &0x3f);
		i64 rank = 0;
		rank = Plain_Rank(buff,index,overloop,bit);
		return rank_base + rank;
	}
	else
	{
		i64 block_anchor = (pos+1)/block_size;
		i64 superblock_anchor = (pos+1)/super_block_size;
		i64 rank = newblock->GetValue(superblock_anchor,1) + newblock->GetValue(block_anchor,0);
		i64 offset = pos / block_size * block_size;
		i64 index = (offset&0x3f);
		buff = data+(offset>>6);
		i64 overloop = block_size;
		bit=Plain_Bit(buff,index,overloop);
		return rank;
	}
}


i64 BitMap::GetCurrentRankvalue(i64 block_anchor,i64 superblock_anchor,i64 offset1,i64 preoffset,i64 radio)
{
	i64 offset_1 = 0;
	if ((block_anchor + 1) % radio != 0)
		offset_1 = offset1 + block->GetValue(((block_anchor+1) << 1));
	else
		offset_1 = superblock->GetValue(((superblock_anchor+1) << 1) );

	return offset_1 - preoffset;
}

void BitMap::Rank(fm_int pos_left,fm_int pos_right,fm_int &rank_left,fm_int &rank_right)
{
	if(pos_left<0 || pos_right <0 || pos_left > bitLen || pos_right > bitLen)
	{
		cout<<pos_left<<" "<<pos_right<<" "<<pos_left<<" "<<pos_right<<" "<<bitLen<<endl;
		cerr<<"BitMap::Rank(int,int,int&,int&) error parmater"<<endl;
		exit(0);
	}
	i64 block_anchor = (pos_left+1)/block_size;
	if(block_anchor==(pos_right+1)/block_size)
	{
		i64 superblock_anchor=(pos_left+1)/super_block_size;
		i64 rank_base =newblock->GetValue(superblock_anchor,1) + newblock->GetValue(block_anchor,0);
		i64 offset = block_anchor * block_size;
		i64 overloop_left = (pos_left+1)%block_size;
		i64 overloop_right= (pos_right+1)%block_size;
		buff = data+ (offset>>6);
		i64 index = (offset&0x3f);
		rank_left = rank_right =rank_base;
		if(overloop_left!=0)
		{
			Plain_Rank(buff,index,overloop_left,overloop_right,rank_left,rank_right);
			return ;
		}
		rank_right += Plain_Rank(buff,index,overloop_right);
	}
	else
	{
		rank_left= Rank(pos_left);
		rank_right=Rank(pos_right);
	}
}

fm_int BitMap::Rank(fm_int pos)
{
	if (pos<0 || pos > bitLen)
	{
		cerr<<"BitMap::Rank  error paramater"<<endl;
		cerr<<pos<<" "<<bitLen<<endl;
		exit(0);
	}
	i64 block_anchor = (pos+1)/block_size;
	i64 superblock_anchor  = ((pos+1)/super_block_size);
	i64 rank_base = newblock->GetValue(superblock_anchor,1) + newblock->GetValue(block_anchor,0);
	i64 offset = block_anchor * block_size;
	buff = data + (offset>>6);
	int overloop = (pos+1)%block_size;
	i64 index = (offset & 0x3f);
	if(overloop > 0)
		rank_base += Plain_Rank(buff,index,overloop);
	return rank_base;
}

int BitMap::GetBit(u64 * data,fm_int index)
{
	fm_int anchor = index/64;
	int pos = 63-index%64;
	return ((data[anchor] &(0x01ull<<pos))>>pos);
}

int BitMap::GetRuns(u64 * data,fm_int &index,int &bit)
{
	bit = GetBit(data,index);
	index = index +1;
	int totle_runs = 1;
	int runs=0;
	
	while(totle_runs < block_size)
	{
		u16 x= GetBits(data,index,16);
		if(bit==1)
			x=(~x);
		runs = Zeros(x);
		totle_runs +=runs;
		index+=runs;
		if(runs < 16)
			break;
	}
	return totle_runs;
}

void BitMap::Append_g(u64 *temp,fm_int &index,u64 value)
{
	u64 y=value;
	int zerosnum = blog(value)-1;
	index+=zerosnum;
	int onesnum = zerosnum+1;
	if(index%64 + onesnum < 65)
	{
		temp[index/64] = (temp[index/64] | (y<<(64-(index%64 + onesnum))));
	}
	else
	{
		int first = 64 - index%64;
		int second = onesnum - first;
		temp[index/64] = (temp[index/64] | (y>>second));
		temp [index/64 +1] = (temp[index/64+1] | (y<<(64-second)));
	}
	index = index + onesnum;
}

void BitMap::Append_f(u64 *temp,i64 &index,u64 value,u32 len)
{
	if(index%64 + len < 65)
	{
		temp[index/64] = (temp[index/64] | (value<<(64-(index%64 + len))));
	}
	else
	{
		int first = 64 - index%64;
		int second = len - first;
		temp[index/64] = (temp[index/64] | (value>>second));
		temp [index/64 +1] = (temp[index/64+1] | (value<<(64-second)));
	}
	index = index + len;
}

void BitMap::BitCopy(u64 * temp,fm_int & index,u64 value)
{
	if(index%64!=0)
	{
		int first = 64 - index % 64;
		int second = 64 - first;
		temp[index/64] = (temp[index/64] | (value>>second));
		temp[index/64 + 1] = (temp[index/64+1] | (value<<first));
	}
	else
		temp[index/64]  = value;
	index = index +64;
}

int BitMap::RL0_Rank(u64 *buff,fm_int &index,int bits_num)
{
	int bit=0;
	return RL0_Rank(buff,index,bits_num,bit);
}

int BitMap::RL0_Rank(u64 * buff,fm_int & index,int bits_num,int &bit)
{
	int rank = 0;
	int bit_count = 0;
	int bits = 0;
	while(true)
	{
		 bits=GammaDecode(buff,index);
		 bit_count = bit_count + bits;
		 if(bit_count >= bits_num)
		 {
			 bit = 0;
			 return rank;
		 }
		 bits=GammaDecode(buff,index);
		 bit_count = bit_count + bits;
		 rank = rank + bits;
		 if(bit_count >= bits_num)
		 {
			 bit = 1;
			 return rank - (bit_count-bits_num);
		 }
	}
}

int BitMap::RL0_Bit(u64 * buff,fm_int & index,int bits_num)
{
	//int rank = 0;
	int bit_count =0;
	int bits = 0;
	while(true)
	{
		bits=GammaDecode(buff,index);
		bit_count = bit_count + bits;
		if(bit_count >= bits_num)
			return  0;
		bits=GammaDecode(buff,index);
		bit_count = bit_count + bits;
		if(bit_count >= bits_num)
			return 1;
	}
}

int BitMap::RL1_Rank(u64 * buff,fm_int &index,int bits_num)
{
	int bit =0;
	return RL1_Rank(buff,index,bits_num,bit);
}

int BitMap::RL1_Rank(u64 * buff,fm_int &index,int bits_num,int & bit)
{
	int rank = 0;
	int bit_count = 0 ;
	int  bits = 0;
	while(true)
	{
		 bits = GammaDecode(buff,index);
		 bit_count = bit_count + bits;
		 rank = rank + bits;
		 if(bit_count >= bits_num)
		 {
			 bit = 1;
			 return rank - (bit_count-bits_num);
		 }
		 bits = GammaDecode(buff,index);
		 bit_count = bit_count + bits;
		 if(bit_count >= bits_num)
		 {
			 bit = 0;
			 return rank;
		 }
	}
}

int BitMap::RL1_Bit(u64 * buff,fm_int &index,int bits_num)
{
	int bit_count = 0 ;
	int  bits = 0;
	while(true)
	{
		 bits = GammaDecode(buff,index);
		 bit_count = bit_count + bits;
		 if(bit_count >= bits_num)
			 return 1;
		 bits = GammaDecode(buff,index);
		 bit_count = bit_count + bits;
		 if(bit_count >= bits_num)
			 return 0;
	}
}

void BitMap::Plain_Rank(u64 *buff,fm_int &index,int bits_left,int bits_right,fm_int & rank_left,fm_int &rank_right)
{
	int rank = Plain_Rank(buff,index,bits_left);
	rank_left += rank;
	index++;
	buff = buff +(index>>6);
	index = (index &0x3f);
	rank_right+= (Plain_Rank(buff,index,bits_right-bits_left)+rank);
}

int BitMap::Plain_Rank(u64 * buff,fm_int &index,int bits_num,int &bit)
{
	if((index &0x3f) + bits_num < 65)
	{
		u64 temp = (buff[index>>6]<<(index&0x3f))>>(64-bits_num);
		index = index + bits_num -1;
		bit=(buff[index>>6]>>(63-(index&0x3f)))&0x01;
		return popcnt(temp);
	}
	int rank = 0;
	int head = 64 - (index&0x3f);
	u64 temp = (buff[index>>6]<<(index&0x3f));
	rank = rank + popcnt(temp);
	bits_num = bits_num - head;
	
	int times = bits_num>>6;
	int i=0;
	for(i=0;i<times;i++)
	rank = rank + popcnt(buff[i+(index>>6)+1]);
	
	if((bits_num&0x3f)!=0)
		rank = rank + popcnt((buff[i+(index>>6)+1] >> (64-(bits_num&0x3f))));
	
	index = index + head + bits_num - 1;
	bit=(buff[index>>6]>>(63-(index&0x3f)))&0x01;
	return rank;
}

int BitMap::Plain_Bit(u64 * buff,fm_int &index,int bits_num)
{
	index = index + bits_num - 1;
	return (buff[index>>6]>>(63-(index&0x3f)))&0x01;
}

int BitMap::Plain_Rank(u64 * buff,fm_int &index,int bits_num)
{
	int bit=0;
	return Plain_Rank(buff,index,bits_num,bit);
}

int BitMap::GammaDecode(u64 * buff,fm_int & index)
{
	u32 x = GetBits(buff,index,32);
	int runs = Zeros(x>>16);
	int bits = (runs<<1)+1;
	index = index + bits;
	return x>>(32-bits);
}

u64 BitMap::GetBits(u64 * buff,fm_int index,int bits)
{
	if((index & 0x3f) + bits < 65)
		return (buff[index>>6]<<(index &0x3f))>>(64-bits);

	int first = 64 - (index &0x3f);
	int second = bits - first;
	u64 high = (buff[index>>6] & ((0x01ull<<first)-1)) << second;
	return high + (buff[(index>>6)+1]>>(64-second));
}

int BitMap::GetZerosRuns(u64 * buff,fm_int &index)
{
	u32 x = GetBits(buff,index,16);
	int runs = Zeros(x);
	index = index + runs;
	return runs;
}

void BitMap::Left(BitMap * left)
{
     this->left = left;
}

void BitMap::Right(BitMap * right)
{
     this->right = right;
}
 
unsigned char BitMap::Label()
{
    return label;
}

int BitMap::Load(loadkit & s)
{
	s.loadi32(level);
	s.loadu8(label);
	s.loadi64(bitLen);
	s.loadi32(block_size);
	s.loadi32(super_block_size);
	s.loadi64(memorysize);
	this->data=NULL;
	this->superblock=NULL;
	this->block=NULL;
	this->coding_style=NULL;
	this->newblock=NULL;
	if(memorysize!=0)
	{
		this->data = new u64[memorysize/8];
		s.loadu64array(data,memorysize/8);	
		newblock = new SBArray();
		newblock->load(s);
	}
    return 0;
}

int BitMap::Save(savekit & s)
{
	s.writei32(level);
	s.writeu8(label);
	s.writei64(bitLen);
	s.writei32(block_size);
	s.writei32(super_block_size);
	s.writei64(memorysize);
	if(memorysize!=0)
	{
		s.writeu64array(data,memorysize/8);
		newblock->write(s);
	}
	return 0;
}

int BitMap::RL_Rank(u64 * buff,fm_int &index,int bits_num,int rl_type)
{
	int bit=0;
	return RL_Rank(buff,index,bits_num,rl_type,bit);
}

int BitMap::RL_Rank(u64 * buff,fm_int &index,int bits_num,int rl_type,int &bit)
{
	int rank=0;
	int r=0;
	int already = 0;
	u64 x = GetBits(buff,index,64);
	int bits = 0;
	int step = 0;
	int runs =0;
	int runs_num = 0;
	u32 anchor=0;
	rl_type=1-rl_type;
	if(bits_num > 32)
	{
		while(true)
		{
			anchor = (x>>48)<<2;
			runs =R[anchor];
			if(runs >0)
			{
				step =R[1+anchor];
				already = already + step;
				if(already > 64)
				{
					index = index + (already -step);
					x = GetBits(buff,index,64);
					already = 0;
					continue;
				}
				bits = R[2+anchor];
				r=R[3+anchor];
				bits=(bits==0)?256:bits;
				if((runs_num & 0x01) ==rl_type)
					rank = rank + r;
				else
					rank = rank + ( bits-r);
				bits_num = bits_num - bits;
				runs_num = runs_num + runs;
				if(bits_num <=0)
					break;
				x = (x<<step);
			}
			else
			{
				step = 1 + (Zeros(x>>48)<<1);
				already = already + step;
				if(already > 64)
				{
					index = index + (already - step);
					x = GetBits(buff,index,64);
					step = 1 + (Zeros(x>>48)<<1);
					already = step;
				}
				bits_num = bits_num - (x>>( 64 - step));
				if((runs_num &0x01) ==rl_type)
					rank = rank + (x>>( 64 - step));
				if(bits_num <=0)
				{
					if((runs_num &0x01)==rl_type)
					{
						bit=1;
						return rank + bits_num;
					}
					else
					{
						bit=0;
						return rank;
					}
				}
				runs_num++;
				x = (x<<step);
			}
		}
	}
	index = index + (already - step);
	bits_num = bits_num + bits;
	runs_num = runs_num - runs;
	if((runs_num &0x01) ==rl_type)
		rank = rank - r;
	else
		rank = rank - (bits -r );
	already = 0;
	x = GetBits(buff,index,64);
	while(true)
	{
		step = 1+ (Zeros(x>>48)<<1);
		already = already + step;
		if(already > 64)
		{
			index = index + (already - step);
			x =GetBits(buff,index,64);
			step = 1+ (Zeros(x>>48)<<1);
			already = step;
		}
		bits_num = bits_num - (x>>( 64 - step));
		if((runs_num &0x01) ==rl_type)
			rank = rank + (x>>( 64 - step));
		if(bits_num <= 0)
		{
			if((runs_num&0x01)==rl_type)
			{
				bit =1;
				return rank + bits_num;
			}
			else
			{
				bit=0;
				return rank;
			}
		}
		runs_num++;
		x=(x<<step);
	}
}

int BitMap::GapG_Rank(u64 * buff,fm_int &index,int bits_num,int cbl,int rl_type)
{
	int bit=0;
	return GapG_Rank(buff,index,bits_num,rl_type,cbl,bit);
}

int BitMap::GapG_Rank(u64 * buff,fm_int &index,int bits_num,int rl_type,int cbl,int &bit)//
{
	int rank = 0;
	int cur = 0;
	u32 anchor = 0;
	i64 oldindex = index;
	u64 x = GetBits(buff,index,64);
	int runs = 0;
	int bits = 0;
	int step = 0;
	int already = 0;
	while(cbl - cur >= 16)
	{
		anchor = (x >> 48) << 2;
		runs = R[anchor];
		if(runs > 0)
		{
			step = R[1 + anchor];
			already += step;
			cur += step;
			if(already > 64)
			{
				index = index + (already - step);
				cur -= step;
				x = GetBits(buff,index,64);
				already = 0;
				continue;
			}
			bits = R[2+anchor];
			bits=(bits==0)?256:bits;
			if(rl_type == 1)
				rank = rank + runs;
			else
				rank = rank + (bits - runs);
			bits_num = bits_num - bits;
			if(bits_num == 0)
			{
				if(rl_type == 1)
					bit = 1;
				else
					bit = 0;
				return rank;
			}
			else if(bits_num < 0)
			{
				u16 tmp = x >> 48;
				int curtmp = 0;
				if(rl_type == 1)
					rank = rank - runs;
				else
					rank = rank - (bits - runs);
				bits_num += bits;
				while(curtmp < 16)
				{
					step = 1 + (Zeros(tmp) << 1);
					curtmp += step;
					bits_num -= (tmp >> (16 - step));
					if(rl_type == 0)
						rank += (tmp >> (16 - step)) - 1;
					else
						++rank;
					if(bits_num == 0)
					{
						if(rl_type == 1)
							bit = 1;
						else 
							bit = 0;
						return rank;
					}
					else if(bits_num < 0)
					{
						if(rl_type == 1)
						{
							bit = 0;
							return --rank ;
						}
						else
						{
							bit = 1;
							return rank + bits_num + 1;
						}
					}
					tmp = (tmp << step);
				}
			}
			x = (x<<step);
		}
		else
		{
			step = 1 + (Zeros(x>>48) << 1);
			already = already + step;
			cur += step;
			if(already > 64)
			{	
				index = index + (already - step);
				cur -= step;
				already = 0;
				x = GetBits(buff,index,64);	
				continue;
			}
			bits_num -= (x >> (64 - step));
			if(rl_type == 0)
				rank += (x >> (64 - step)) - 1;
			else
				++rank;
			if(bits_num == 0)
			{
				if(rl_type == 1)
					bit = 1;
				else 
					bit = 0;
				return rank;
			}
			else if(bits_num < 0)
			{
				if(rl_type == 1)
				{
					bit = 0;
					return --rank ;
				}
				else
				{
					bit = 1;
					return rank + bits_num + 1;
				}
			}
			x = x << step;
		}
	}
	int tmplen = cbl - cur;
	int curtmp = 0;
	index = oldindex + cur;
	x = GetBits(buff,index,64);	
	while(curtmp < tmplen)
	{
		step = 1 + (Zeros(x>>48) << 1);
		curtmp += step;
		already += step;
		bits_num -= (x >> (64 - step));
		if(rl_type == 0)
			rank += (x >> (64 - step)) - 1;
		else
			++rank;
		if(bits_num == 0)
		{
			if(rl_type == 1)
				bit = 1;
			else 
				bit = 0;
			return rank;
		}
		else if(bits_num < 0)
		{
			if(rl_type == 1)
			{
				bit = 0;
				return --rank ;
			}
			else
			{
				bit = 1;
				return rank + bits_num + 1;
			}
		}
		x = (x << step);
	}
	if(rl_type == 1)
	{
		bit = 0;
		return rank;
	}
	else
	{
		bit = 1;
		return rank + bits_num;
	}
}

void BitMap::EF_Rank(u64 * buff,i64 &index,int bits_left,int bits_right,i64 &rank_left,i64 &rank_right,int ef_type,u64 elemnums,int currentBlockLen){
	u64 rank_left1=EF_Rank(buff,index,bits_left,ef_type,elemnums,currentBlockLen);
	u64 rank_right1=EF_Rank(buff,index,bits_right,ef_type,elemnums,currentBlockLen);
	rank_left+=rank_left1;
	rank_right+=rank_right1;
}
u64  BitMap::EF_Rank(u64 * buff,i64 &index,int bits_num,int ef_type,u64 elemnums,int currentBlockLen){
	int bit=0;
	return EF_Rank(buff,index,bits_num,ef_type,elemnums,currentBlockLen,bit);
}
u64  BitMap::EF_Rank(u64 * buff,i64 &index,int bits_num,int ef_type,u64 elemnums,int currentBlockLen, int &bit){
	bits_num--;
	bit=1-ef_type;
	#ifdef TESTEFRANK
	cout<<"eftype: "<<ef_type<<endl;
	#endif
	int low=blog(block_size/elemnums)-1;
	low=0>low?0:low;
	i64 high_index=index+low*elemnums;
	i64 low_index=index;
	int highbucket=bits_num>>low;
	int curbucket = 0;
	int prebucket = 0;
	int lowbits=bits_num&((1<<low)-1);
	i64 rank=0;
	int highbucketNum = 0;
	int high_currentblocklen=currentBlockLen-low*elemnums;
	i64 upper_index = high_index+high_currentblocklen;

	if(highbucket==0){
		while(high_index<upper_index && GetBit(buff,high_index++)==1){
			int lowbits=GetBits(buff,low_index,low);
			low_index+=low;
			if(lowbits>bits_num){
				break;
			}
			if(lowbits==bits_num){
				bit=ef_type;
			}
			++rank;
		}
	}
	else{
		u16 i = 0;
		while(high_index < upper_index && curbucket<highbucket){
			int getbitslen = 16<upper_index-high_index ? 16:upper_index-high_index;
			i = GetBits(buff,high_index,getbitslen);
			if(getbitslen!=16)
				i = i<<(16-getbitslen);
			prebucket = curbucket;
			curbucket += R1[i];
			rank += R3[i*18]-R1[i];
			high_index += R3[i*18];
		}
		if(high_index==upper_index && curbucket<highbucket){
			return ef_type==1 ? rank:bits_num-rank+1;
		}
		else if(curbucket>=highbucket){
			rank -= R3[i*18]-R1[i];
			//high_index -= R3[i*18];
			for(int j=0;j<highbucket-prebucket;++j){
				rank += R3[i*18+j+1];
			}
			highbucketNum = R3[i*18+1+(highbucket-prebucket)];

			if(curbucket==highbucket && R3[i*18]==16 && high_index<upper_index){
				bool flag = false;
				for(int j=highbucket-prebucket+1;j<=16;++j){
					if(R3[i*18+j+1]>0){
						flag = true;
						break;
					}
				}
				if(!flag){
					int getbitslen = 16<upper_index-high_index ? 16:upper_index-high_index;
					u16 x = GetBits(buff,high_index,getbitslen);
					if(getbitslen!=16)
						x = x<<(16-getbitslen);
					high_index += getbitslen;
					int one_runs = Zeros(~x);
					while(high_index<upper_index && one_runs==16){
						getbitslen = 16<upper_index-high_index ? 16:upper_index-high_index;
						x = GetBits(buff,high_index,getbitslen);
						if(getbitslen!=16)
							x = x<<(16-getbitslen);
						one_runs += Zeros(~x);
						high_index += getbitslen;
					}
					highbucketNum += one_runs;
				}
			}
		}
		int decompressNum = 0;
		low_index += rank*low;
		while(decompressNum<highbucketNum){
			auto lbits = GetBits(buff,low_index,low);
			if(lbits > lowbits){
				break;
			}
			++rank;
			low_index += low;
			++decompressNum;
			if(lbits == lowbits){
				bit = ef_type;
				break;
			}
		}
	}
	return ef_type==1 ? rank:bits_num-rank+1;
}
int BitMap::EF0_Bit(u64 * buff,i64 &index,int bits,u64 currentblocklen,int elemnums)
{
	bits--;
	#ifdef TESTEFRANK
	cout<<"eftype: "<<ef_type<<endl;
	#endif
	int low=blog(block_size/elemnums)-1;
	low=0>low?0:low;
	i64 high_index=index+low*elemnums;
	i64 low_index=index;
	int highbucket=bits>>low;
	int curbucket = 0;
	int prebucket = 0;
	int lowbits=bits&((1<<low)-1);
	i64 rank=0;
	int highbucketNum = 0;
	int high_currentblocklen=currentblocklen-low*elemnums;
	i64 upper_index = high_index+high_currentblocklen;

	if(highbucket==0){
		while(high_index<upper_index && GetBit(buff,high_index++)==1){
			int lowbits=GetBits(buff,low_index,low);
			low_index+=low;
			if(lowbits>bits){
				return 1;
			}
			if(lowbits==bits){
				return 0;
			}
			rank++;
		}
	}
	else{
		u16 i = 0;
		while(high_index < upper_index && curbucket<highbucket){
			int getbitslen = 16<upper_index-high_index ? 16:upper_index-high_index;
			i = GetBits(buff,high_index,getbitslen);
			if(getbitslen!=16)
				i = i<<(16-getbitslen);
			prebucket = curbucket;
			curbucket += R1[i];
			rank += R3[i*18]-R1[i];
			high_index += R3[i*18];
		}
		if(high_index==upper_index && curbucket<highbucket){
			return 1;
		}
		else if(curbucket>=highbucket){
			rank -= R3[i*18]-R1[i];
			for(int j=0;j<highbucket-prebucket;++j){
				rank += R3[i*18+j+1];
			}
			highbucketNum = R3[i*18+1+(highbucket-prebucket)];
			if(curbucket==highbucket && R3[i*18]==16 && high_index<upper_index){
				bool flag = false;
				for(int j=highbucket-prebucket+1;j<=16;++j){
					if(R3[i*18+j+1]>0){
						flag = true;
						break;
					}
				}
				if(!flag){
					int getbitslen = 16<upper_index-high_index ? 16:upper_index-high_index;
					u16 x = GetBits(buff,high_index,getbitslen);
					if(getbitslen!=16)
						x = x<<(16-getbitslen);
					high_index += getbitslen;
					int one_runs = Zeros(~x);
					while(high_index<upper_index && one_runs==16){
						getbitslen = 16<upper_index-high_index ? 16:upper_index-high_index;
						x = GetBits(buff,high_index,getbitslen);
						if(getbitslen!=16)
							x = x<<(16-getbitslen);
						one_runs += Zeros(~x);
						high_index += getbitslen;
					}
					highbucketNum += one_runs;
				}
			}
		}
		int decompressNum = 0;
		int bucketBase = highbucket<<low;
		while(decompressNum<highbucketNum){
			int decompressElem = bucketBase+GetBits(buff,index+rank*low,low);
			if(decompressElem > bits){
				break;
			}
			++rank;
			++decompressNum;
			if(decompressElem==bits){
				return 0;
			}
		}
	}
	return 1;
}
int BitMap::EF1_Bit(u64 * buff,i64 &index,int bits,u64 currentblocklen,int elemnums)
{
	bits--;
	#ifdef TESTEFRANK
	cout<<"eftype: "<<ef_type<<endl;
	#endif
	int low=blog(block_size/elemnums)-1;
	low=0>low?0:low;
	i64 high_index=index+low*elemnums;
	i64 low_index=index;
	int highbucket=bits>>low;
	int curbucket = 0;
	int prebucket = 0;
	int lowbits=bits&((1<<low)-1);
	i64 rank=0;
	int highbucketNum = 0;
	int high_currentblocklen=currentblocklen-low*elemnums;
	i64 upper_index = high_index+high_currentblocklen;

	int highpos=0;
	if(highbucket==0){
		while(high_index<upper_index && GetBit(buff,high_index++)==1){
			int lowbits=GetBits(buff,low_index,low);
			low_index+=low;
			if(lowbits>bits){
				return 0;
			}
			if(lowbits==bits){
				return 1;
			}
			rank++;
			highpos++;
		}
	}
	else{
		u16 i = 0;
		while(high_index < upper_index && curbucket<highbucket){
			int getbitslen = 16<upper_index-high_index ? 16:upper_index-high_index;
			i = GetBits(buff,high_index,getbitslen);
			if(getbitslen!=16)
				i = i<<(16-getbitslen);
			prebucket = curbucket;
			curbucket += R1[i];
			rank += R3[i*18]-R1[i];
			high_index += R3[i*18];
		}
		if(high_index==upper_index && curbucket<highbucket){
			return 0;
		}
		else if(curbucket>=highbucket){
			rank -= R3[i*18]-R1[i];
			for(int j=0;j<highbucket-prebucket;++j){
				rank += R3[i*18+j+1];
			}
			highbucketNum = R3[i*18+1+(highbucket-prebucket)];
			if(curbucket==highbucket && R3[i*18]==16 && high_index<upper_index){
				bool flag = false;
				for(int j=highbucket-prebucket+1;j<=16;++j){
					if(R3[i*18+j+1]>0){
						flag = true;
						break;
					}
				}
				if(!flag){
					int getbitslen = 16<upper_index-high_index ? 16:upper_index-high_index;
					u16 x = GetBits(buff,high_index,getbitslen);
					if(getbitslen!=16)
						x = x<<(16-getbitslen);
					high_index += getbitslen;
					int one_runs = Zeros(~x);
					while(high_index<upper_index && one_runs==16){
						getbitslen = 16<upper_index-high_index ? 16:upper_index-high_index;
						x = GetBits(buff,high_index,getbitslen);
						if(getbitslen!=16)
							x = x<<(16-getbitslen);
						one_runs += Zeros(~x);
						high_index += getbitslen;
					}
					highbucketNum += one_runs;
				}
			}
		}
		int decompressNum = 0;
		int bucketBase = highbucket<<low;
		while(decompressNum<highbucketNum){
			int decompressElem = bucketBase+GetBits(buff,index+rank*low,low);
			if(decompressElem > bits){
				break;
			}
			++rank;
			++decompressNum;
			if(decompressElem==bits){
				return 1;
			}
		}
	}
	return 0;
}


void BitMap::RL_Rank(u64 *buff,fm_int &index,int bits_left,int bits_right,fm_int &rank_left,fm_int &rank_right,int rl_type)
{	
	int rank = 0;
	//int rank_diff = 0;
	int r = 0;
	int already = 0;
	u64 x=GetBits(buff,index,64);
	int bits = 0;
	int step = 0;
	int runs = 0;
	int runs_num = 0;
	u32 anchor = 0;
	bool left=true;
	rl_type = 1-rl_type;

	while(true)
	{
		anchor = (x>>48)<<2;
		runs = R[anchor];
		if(runs>0)
		{
			step = R[anchor+1];
			already = already + step;
			if(already > 64)
			{
				index = index + (already - step);
				x = GetBits(buff,index,64);
				already = 0;
				continue;
			}

			bits = R[anchor+2];
			r = R[anchor+3];
			bits =(bits==0)?256:bits;
			if((runs_num & 0x01)==rl_type)
				rank  = rank + r;
			else
				rank = rank + (bits-r);

			bits_left = bits_left - bits;
			bits_right = bits_right- bits;
			runs_num = runs_num + runs;
			
			if(bits_left <= 0)
				break;
			x = (x<<step);
		}
		else
		{
			step = 1 + (Zeros(x>>48)<<1);
			already = already + step;
			if(already > 64)
			{
				index = index + (already - step);
				x = GetBits(buff,index,64);
				step = 1 + (Zeros(x>>48)<<1);
				already  = step;
			}
			bits = (x>>(64-step));
			bits_left = bits_left - bits;
			bits_right = bits_right - bits;
			r = bits;
			if((runs_num & 0x01)==rl_type)
				rank = rank + bits;
			runs = 1;
			runs_num += runs;
			if(bits_left <=0)
				break;
			x = (x<<step);
		}
	}
	index = index + (already  - step);
	bits_left = bits_left + bits;
	bits_right= bits_right+ bits;
	runs_num = runs_num - runs;
	if((runs_num & 0x01)==rl_type)
		rank = rank - r;
	else
		rank = rank - (bits -r );
	already = 0;
	x = GetBits(buff,index,64);

	while(true)
	{
		step = 1 + (Zeros(x>>48)<<1);
		already  = already + step;
		if(already > 64)
		{
			index = index + (already - step);
			x = GetBits(buff,index,64);
			step = 1 + (Zeros(x>>48)<<1);
			already = step;
		}
		bits  = (x>>(64-step));
		bits_left = bits_left - bits;
		bits_right= bits_right- bits;
		if((runs_num & 0x01)==rl_type)
			rank = rank + bits;
		if(left && bits_left <=0)
		{
			if((runs_num & 0x01)==rl_type)
				rank_left += (rank + bits_left);
			else
				rank_left += rank;
			left = false;
		}

		if(bits_right <=0)
		{
			if((runs_num & 0x01)==rl_type)
				rank_right += (rank + bits_right);
			else
				rank_right += rank;
			return ;
		}
		runs_num++;
		x=(x<<step);
	}
}

void BitMap::GapG_Rank(u64 *buff,fm_int &index,int bits_left,int bits_right,fm_int &rank_left,fm_int &rank_right,int rl_type,int cbl)
{
	
	i64 oldindex = index;
	u64 rank_left1=GapG_Rank(buff,oldindex,bits_left,cbl,rl_type);
	u64 rank_right1=GapG_Rank(buff,index,bits_right,cbl,rl_type);
	rank_left += rank_left1;
	rank_right += rank_right1;
	return;
}

int BitMap::GapG0_Bit(u64 * buff,fm_int & index,int bits_num,int cbl)
{
	int rank = 0;
	int cur = 0;
	u32 anchor = 0;
	i64 oldindex = index;
	u64 x = GetBits(buff,index,64);
	int runs = 0;
	int bits = 0;
	int step = 0;
	int already = 0;
	while(cbl - cur >= 16)
	{
		anchor = (x >> 48) << 2;
		runs = R[anchor];
		if(runs > 0)
		{
			step = R[1 + anchor];
			already += step;
			cur += step;
			if(already > 64)
			{
				index = index + (already - step);
				cur -= step;
				x = GetBits(buff,index,64);
				already = 0;
				continue;
			}
			bits = R[2+anchor];
			bits=(bits==0)?256:bits;
			rank = rank + (bits - runs);
			bits_num = bits_num - bits;
			if(bits_num == 0)
				return 0;
			else if(bits_num < 0)
			{
				u16 tmp = x >> 48;
				int curtmp = 0;
				rank = rank - (bits - runs);
				bits_num += bits;
				while(curtmp < 16)
				{
					step = 1 + (Zeros(tmp) << 1);
					curtmp += step;
					bits_num -= (tmp >> (16 - step));
					rank += (tmp >> (16 - step)) - 1;
					if(bits_num == 0)
						return 0;
					else if(bits_num < 0)
						return 1;
					tmp = (tmp << step);
				}
			}
			x = (x<<step);
		}
		else
		{
			step = 1 + (Zeros(x>>48) << 1);
			already = already + step;
			cur += step;
			if(already > 64)
			{	
				index = index + (already - step);
				cur -= step;
				already = 0;
				x = GetBits(buff,index,64);	
				continue;
			}
			bits_num -= (x >> (64 - step));
			rank += (x >> (64 - step)) - 1;
			if(bits_num == 0)
				return 0;
			else if(bits_num < 0)
				return 1;
			x = x << step;
		}
	}
	int tmplen = cbl - cur;
	int curtmp = 0;
	x = GetBits(buff,oldindex + cur,64);	
	while(curtmp < tmplen)
	{
		step = 1 + (Zeros(x>>48) << 1);
		curtmp += step;
		already += step;
		if(already > 64)
		{
			index = index + (already - step);
			curtmp -= step;
			x = GetBits(buff,index,64);
			already = 0;
			continue;
		}
		bits_num -= (x >> (64 - step));
		rank += (x >> (64 - step)) - 1;
		if(bits_num == 0)
			return 0;
		else if(bits_num < 0)
			return 1;
		x = (x << step);
	}
	return 1;
}
int BitMap::GapG1_Bit(u64 * buff,fm_int & index,int bits_num,int cbl)
{
	int rank = 0;
	int cur = 0;
	u32 anchor = 0;
	i64 oldindex = index;
	u64 x = GetBits(buff,index,64);
	int runs = 0;
	int bits = 0;
	int step = 0;
	int already = 0;
	while(cbl - cur >= 16)
	{
		anchor = (x >> 48) << 2;
		runs = R[anchor];
		if(runs > 0)
		{
			step = R[1 + anchor];
			already += step;
			cur += step;
			if(already > 64)
			{
				index = index + (already - step);
				cur -= step;
				x = GetBits(buff,index,64);
				already = 0;
				continue;
			}
			bits = R[2+anchor];
			bits=(bits==0)?256:bits;
			rank = rank + (bits - runs);
			bits_num = bits_num - bits;
			if(bits_num == 0)
				return 1;
			else if(bits_num < 0)
			{
				u16 tmp = x >> 48;
				int curtmp = 0;
				rank = rank - (bits - runs);
				bits_num += bits;
				while(curtmp < 16)
				{
					step = 1 + (Zeros(tmp) << 1);
					curtmp += step;
					bits_num -= (tmp >> (16 - step));
					rank += (tmp >> (16 - step)) - 1;
					if(bits_num == 0)
						return 1;
					else if(bits_num < 0)
						return 0;
					tmp = (tmp << step);
				}
			}
			x = (x<<step);
		}
		else
		{
			step = 1 + (Zeros(x>>48) << 1);
			already = already + step;
			cur += step;
			if(already > 64)
			{	
				index = index + (already - step);
				cur -= step;
				already = 0;
				x = GetBits(buff,index,64);	
				continue;
			}
			bits_num -= (x >> (64 - step));
			rank += (x >> (64 - step)) - 1;
			if(bits_num == 0)
				return 1;
			else if(bits_num < 0)
				return 0;
			x = x << step;
		}
	}
	int tmplen = cbl - cur;
	int curtmp = 0;
	x = GetBits(buff,oldindex + cur,64);	
	while(curtmp < tmplen)
	{
		step = 1 + (Zeros(x>>48) << 1);
		curtmp += step;
		already += step;
		if(already > 64)
		{
			index = index + (already - step);
			curtmp -= step;
			x = GetBits(buff,index,64);
			already = 0;
			continue;
		}
		bits_num -= (x >> (64 - step));
		rank += (x >> (64 - step)) - 1;
		if(bits_num == 0)
			return 1;
		else if(bits_num < 0)
			return 0;
		x = (x << step);
	}
	return 0;
}