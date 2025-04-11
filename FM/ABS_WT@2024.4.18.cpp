/*============================================
# Filename: ABS_WT.cpp
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
//2024.4.18 当前版本SA为下标采样，且不支持extract操作，索引也没有RankL结构和D结构
#include"ABS_WT.h"
#include<string.h>
#include<map>
u64 GetBits(u64 * buff,int &index,int bits)
{
	if((index & 0x3f) + bits < 65)
		return (buff[index>>6] >>( 64 -((index&0x3f) + bits))) & ((0x01ull<<bits)- 1);
	int first = 64 - (index &0x3f);
	int second = bits - first;
	u64 high = (buff[index>>6] & ((0x01ull<<first)-1)) << second;
	return high + (buff[(index>>6)+1]>>(64-second));
}

int Zeros(u16 x,ABS_FM *t)
{
	if(t->Z[x>>8]==8)
		return t->Z[x>>8]+t->Z[(uchar)x];
	else
		return t->Z[x>>8];
}

int GammaDecode(u64 * buff,int & index,ABS_FM * t)
{
	u32 x = GetBits(buff,index,32);
	int runs = Zeros(x>>16,t);
	int bits = (runs<<1)+1;
	index = index + bits;
	return x>>(32-bits);
}

ABS_FM::ABS_FM(const char * filename,int block_size,int r,int D)
{
	this->block_size = block_size;
	this->r=r;
	this->D =D;
	this->T=NULL;
	T = Getfile(filename);                                                                        //扫描本文，构建C表
	Inittable();
}

ABS_FM::~ABS_FM()
{
	DestroyWaveletTree();
	if(T)
		delete [] T;
	if(bwt)
		delete [] bwt;
	if(SAL)
		delete SAL;
	if(RankL)
		delete RankL;
	if(C)
		delete [] C;
	if(code)
		delete [] code;
	if(Z)
		delete [] Z;
	if(R)
		delete [] R;
}

fm_int ABS_FM::SizeInByte()
{
	return TreeSizeInByte(root) + SAL->GetMemorySize() + RankL->GetMemorySize()+newD->GetMemorySize();
}

fm_int ABS_FM::SizeInByte_count()
{
	return TreeSizeInByte(root);
}

fm_int ABS_FM::SizeInByte_locate()
{
	return TreeSizeInByte(root) + SAL->GetMemorySize() + RankL->GetMemorySize() +newD->GetMemorySize();
}

fm_int ABS_FM::SizeInByte_extract()
{
	return TreeSizeInByte(root) + SAL->GetMemorySize() + RankL->GetMemorySize()+newD->GetMemorySize();
}
fm_int ABS_FM::SizeInByte_SAL()
{
	return SAL->GetMemorySize();
}
fm_int ABS_FM::SizeInByte_RankL()
{
	return RankL->GetMemorySize();
}
fm_int ABS_FM::SizeInByte_newD()
{
	return newD->GetMemorySize();
}

double ABS_FM::GetRuns(){
	return runs;
}

int ABS_FM::TreeNodeCount(BitMap * r)
{
	if(r==NULL)
		return 0;
	return TreeNodeCount(r->Left()) + TreeNodeCount(r->Right()) + 1; //中间节点 + 叶子节点，2N-1（N为字符表大小）
}

fm_int ABS_FM::TreeSizeInByte(BitMap * r)
{
	fm_int size = 0;
	if(r->Left())
		size += TreeSizeInByte(r->Left());
	if(r->Right())
		size+=TreeSizeInByte(r->Right());
	size = size + r->SizeInByte();
	return size;
}

fm_int ABS_FM::TreePlainBlock(BitMap * r)
{
	fm_int size = 0;
	if(r->Left())
		size += TreePlainBlock(r->Left());
	if(r->Right())
		size+=TreePlainBlock(r->Right());
	size = size + r->PlainBlock();
	return size;
}

fm_int ABS_FM::TreeRLGBlock(BitMap * r)
{
	fm_int size = 0;
	if(r->Left())
		size += TreeRLGBlock(r->Left());
	if(r->Right())
		size+=TreeRLGBlock(r->Right());
	size = size + r->RLGBlock();
	return size;
}

fm_int ABS_FM::TreeALL01Block(BitMap * r)
{
	fm_int size = 0;
	if(r->Left())
		size += TreeALL01Block(r->Left());
	if(r->Right())
		size+=TreeALL01Block(r->Right());
	size = size + r->ALL01Block();
	return size;
}

fm_int ABS_FM::TreeEF01Block(BitMap * r)
{
	fm_int size = 0;
	if(r->Left())
		size += TreeEF01Block(r->Left());
	if(r->Right())
		size+=TreeEF01Block(r->Right());
	size = size + r->EF01Block();
	return size;
}

fm_int ABS_FM::TreePlainSize(BitMap * r)
{
	fm_int size = 0;
	if(r->Left())
		size += TreePlainSize(r->Left());
	if(r->Right())
		size+=TreePlainSize(r->Right());
	size = size + r->PlainSize();
	return size;
}

fm_int ABS_FM::TreeRLGSize(BitMap * r)
{
	fm_int size = 0;
	if(r->Left())
		size += TreeRLGSize(r->Left());
	if(r->Right())
		size+=TreeRLGSize(r->Right());
	size = size + r->RLGSize();
	return size;
}

fm_int ABS_FM::TreeEF01Size(BitMap * r)
{
	fm_int size = 0;
	if(r->Left())
		size += TreeEF01Size(r->Left());
	if(r->Right())
		size+=TreeEF01Size(r->Right());
	size = size + r->EF01Size();
	return size;
}

fm_int ABS_FM::TreeMemSize(BitMap * r)
{
	fm_int size = 0;
	if(r->Left())
		size += TreeMemSize(r->Left());
	if(r->Right())
		size+= TreeMemSize(r->Right());
	size = size + r->MemSize();
	return size;
}
fm_int ABS_FM::TreeNodeBlockSize(BitMap * r)
{
	fm_int size = 0;
	if(r->Left())
		size += TreeNodeBlockSize(r->Left());
	if(r->Right())
		size+= TreeNodeBlockSize(r->Right());
	size = size + r->NodeBlockSize();
	return size;
}
fm_int ABS_FM::TreeCodeStyleSize(BitMap * r)
{
	fm_int size = 0;
	if(r->Left())
		size += TreeCodeStyleSize(r->Left());
	if(r->Right())
		size+= TreeCodeStyleSize(r->Right());
	size = size + r->CodeStyleSize();
	return size;
}














void ABS_FM::DrawBackSearch(const char * pattern,fm_int & Left,fm_int &Right)
{
	int len = strlen(pattern);
	fm_int occ_left=0;
	fm_int occ_right=0;
	if(len <=0)
	{
		Left =1;
		Right = 0;
		return;
	}
	int i = len -1;
	unsigned char c = pattern[i];
	int coding = code[c];
	if (coding==-1)
	{
		Left = 1;
		Right = 0;
		return ;
	}
	Left = C[coding];
	Right = C[coding+1]-1;
	i=i-1;
	while ((Left <= Right) && (i>=0))
	{
		c = pattern[i];
		coding = code[c];
		if(coding == -1)
		{
			Left = 1;
			Right = 0;
			return;
		}
		Occ(c,Left-1,Right,occ_left,occ_right);
		Left = C[coding]+occ_left;
		Right = C[coding]+occ_right-1;
		i=i-1;
	}
	if(Right < Left)
	{
		Left = 1;
		Right = 0;
		return ;
	}
	return;
}


void ABS_FM::Counting(const char * pattern,fm_int & num)
{
	fm_int Left=1;
	fm_int Right =0;
	DrawBackSearch(pattern,Left,Right);
	num = Right -Left +1;
}


fm_int * ABS_FM::Locating(const char * pattern,fm_int &num)
{
	fm_int Left=1;
	fm_int Right = 0;
	DrawBackSearch(pattern,Left,Right);
	if(Right < Left)
		return NULL;
	num = Right - Left + 1;
	fm_int *pos =new fm_int[num];
	for(fm_int i=0;i<num;i++){
		pos[i]=Lookup(Left + i);
//		cout<<pos[i]<<endl;
	}
	return pos;
}

fm_int * ABS_FM::Locatingnew(const char * pattern,fm_int &num)
{
	fm_int Left=1;
	fm_int Right = 0;
	DrawBackSearch(pattern,Left,Right);
	if(Right < Left)
		return NULL;
	num = Right - Left + 1;
	fm_int *pos =new fm_int[num];
	Getpos(Left,Right,pos);
	return pos;
}

void ABS_FM::Getpos(fm_int left,fm_int right,fm_int *&pos){
	fm_int num=right-left+1;
	fm_int *dis=new fm_int[num];//距离数组，标记当前查询的位置与上一次查询到当前位置的跳步次数
	fm_int *pred=new fm_int[num];//记录跳步到当前查询位置的上一个位置
	for(fm_int i=0;i<num;i++){
		pred[i]=-1;
		pos[i]=-1;
	}
	fm_int f=0;
	fm_int i=0;
	fm_int step=0;
	fm_int q=0;
	fm_int s=0;
	for (fm_int j = right; j >= left; j--){
		f = 0;
		i = j;
		step = 0;
		// cout<<"start####"<<j<<endl;
		while (newD->getD(i)!=1){
			i = LF(i);
			step++;
			//当跳步达到的位置在左右边界内，则记录跳步次数，及再次跳步到查询区间的查询位置
			if (left <= i&&i <= right){
				dis[j - left] = step;
				pred[i - left] = j;
				f = 1;
				break;
			}
			
		}
		//标志位为0，表示查询位置跳步到采样殿过程中没有再次进入查询区间，由采样点和跳步次数获取出现位置
		if (f == 0){
			fm_int r=newD->rank(i);
			fm_int sal=0;
			//若当前文本长度不为采样步长的整数倍，则表示在采样时，标记数组对于sa[0]=n-1的位置多记录了一次，因此要r-2,否则r-1
			if((n-1)%D==0){
				sal=SAL->GetValue(r-1);
			}else{
				sal=SAL->GetValue(r-2);
			}
			pos[j - left] = (sal*D+step)%n;
		}
	}
	for(fm_int j=left;j<=right;j++){
		if(pos[j-left]!=-1){
			q=j;
			while(pred[q-left]!=-1){ 
				s=pos[q-left];
				i=pred[q-left];
				step=dis[i-left];
				pos[i-left]=s+step;
				pred[q-left]=-1;
				q=i;
			}
		}
	}
	delete [] pred;
	delete [] dis;
}

void ABS_FM::Fmocc(fm_int l, fm_int r,unsigned char c ,fm_int *L,fm_int * R)               //给[l,r]，LF更新出[L,R]，指定c
{
	int coding = code[c];                                                                  //给c，转为小波树上的码值
	// cout<<"Fmocc... "<<"  l:"<<l<<"  r:"<<r<<"  c:"<<(int)c;
	/* for(int i=0;i<16;i++){
		cout<<"C["<<i<<"]..."<<(char)iupacChar[i]<<"..."<<C[code[i]]<<"..."<<code[i]<<endl;                //打印C表
	}
	getchar(); */
	// if (coding==-1)
	// {
	// 	*L = 1;
	// 	*R = 0;
	// 	return ;
	// }
	// if(l == 0 && r == n - 1)                                                              //这里要求，给的[L,R]初值必须是[0,n-1]
	if(l == -1 && r == n)
	{
		if(code[c]<0){                                       //不存在的字符（C表上映射到下一个）
			*L = 1 +C[code[c+1]];
			*R = 0 +C[code[c+1]];
			return;
		}
		*L = C[code[c]]+1;                              //按bwt的习惯（20220910）
		*R = C[code[c] + 1];
		return;
	}
//	long long Left = *L, Right = *R;
	// Occ(c,l - 1,r,*L,*R);                                                                 //bwt上的Rank（小优化，L>R，区间为空，提前终止）
	if(code[c]<0){                                       //不存在的字符（C表上映射到下一个）
		*L = 1 +C[code[c+1]];
		*R = 0 +C[code[c+1]];
		return;
	}
	// Occ(c,l-1,r-1,*L,*R);                                                                 //bwt上的Rank（小优化，L>R，区间为空，提前终止）
	Occ(c,l,r,*L,*R);                                                                 
	if(*L==-1){
		*L=0;
	}
	// if(*L!=*R){
	// 	cout<<"l:"<<l<<"  r:"<<r<<"  *L:"<<*L<<"  *R:"<<*R;
	// }
	// *L += C[coding];
	*L += C[coding] + 1;                                                                      //再加上C[]
	*R += C[coding];
	// if(*L!=*R+1){
		// cout<<"  char:"<<(char)c<<"  C[]:"<<C[coding]<<"  L:"<<*L<<"  R:"<<*R<<"  num:"<<*R-*L+1<<endl;
		// cout<<"  char:"<<(char)iupacChar[c]<<"  C[]:"<<C[coding]<<"  L:"<<*L<<"  R:"<<*R<<"  num:"<<*R-*L+1<<endl;
		// cout<<"  L:"<<*L<<"  R:"<<*R<<"  num:"<<*R-*L+1<<endl;
	// }
	// if(*R < *L)
	// {
	// 	*L = 1;
	// 	*R = 0;
	// 	return ;
	// }
	
	return;
}

// void ABS_FM::Fm4occ(fm_int l, fm_int r,fm_int *L,fm_int * R) //同上，但不指定c，同时更新ACGT
// {
// 	int coding[4];
// 	for(int i = 0;i < 4;++i)
// 		coding[i] = code[i];                                 //转为小波树上的码值
// 		// coding[i] = code[20 + i];                            //转字符编码，hzt
// 	if(l == 0 && r == n - 1)                                 //处理初始的区间
// 	{
// 		for(int i = 0;i < 4;++i)
// 		{
// 			L[i] = C[coding[i]];
// 			R[i] = C[coding[i] + 1] - 1;
// 		}
// 		return ;
// 	}
// 	Occ4(l - 1,r,L,R);                                       //Rank部分
// //	Occ(20,l - 1,r,L[0],R[0]);
// 	for(int i = 0;i < 4;++i)
// 	{
// 		L[i] += C[coding[i]];                                //加上C[]
// 		R[i] += C[coding[i]] - 1;
// 		if(R[i] < L[i])
// 		{
// 			L[i] = 1;
// 			R[i] = 0;
// 		}
// 	}	
// 	return;
// }

void ABS_FM::Fm16occ(fm_int l, fm_int r,fm_int *L,fm_int * R) //同上，但不指定c，同时更新ACGT
{	
	/* cout<<"print...code...l:"<<l<<"  r:"<<r<<"  n:"<<n<<endl;                            //打印code
	for(int i=0;i<16;i++){
		cout<<i<<":"<<code[i]<<endl;
	}
	getchar(); */
	if(l == -1 && r == n)                                  //处理初始的区间（主要是加速用）
	{
		for(int i=1; i<16; ++i)
		{
			// cout<<"fm16occ...speed"<<endl;
			// L[i] = C[code[i]];
			if(code[i]<0){                                       //不存在的字符（C表上映射到下一个）
				L[i] = 1 +C[code[i+1]];
				R[i] = 0 +C[code[i+1]];
				continue;
			}
			L[i] = C[code[i]]+1;                              //按bwt的习惯（20220910）
			R[i] = C[code[i] + 1];
		}
		return ;
	}
	// Occ4(l - 1,r,L,R);                                    //Rank部分
    // Occ(20,l - 1,r,L[0],R[0]);
	for(int i=1; i<16; ++i)                                  //这里为12
	{
		if(code[i]<0){                                       //不存在的字符（C表上映射到下一个）
			L[i] = 1 +C[code[i+1]];
			R[i] = 0 +C[code[i+1]];
			continue;
		}
		// Occ(i,l-1,r-1,L[i],R[i]);                          //调用单个的Occ（fm的12，更新到bwt的16）
		Occ(i,l,r,L[i],R[i]);                          //调用单个的Occ（fm的12，更新到bwt的16）    
		L[i] += C[code[i]] + 1;                                  //加上C[]
		R[i] += C[code[i]];
		// if(R[i] < L[i])
		// {
		// 	L[i] = 1;
		// 	R[i] = 0;
		// }
	}	
	return;
}

// void ABS_FM::Occ4(fm_int pos_left,fm_int pos_right,fm_int *rank_left,fm_int *rank_right)    //限定死两层，不循环找叶子节点
// {
// 	BitMap *r = root;                                                                       //获取小波树根节点（第一层）
// 	fm_int r1_l0_left, r1_l0_right, r0_l0_left, r0_l0_right;
// 	int l_flag = 1, r_flag = 1;
// 	if(pos_left > -1 && pos_right > -1) 
// 	{
// 		fm_int tmp_left, tmp_right;	
// 		r->Rank(pos_left, pos_right, tmp_left, tmp_right);                                 //输入左右边界，输出左右边界前1的个数
// 		r1_l0_left = tmp_left - 1;                                                         //左边界前，1的个数
// 		r1_l0_right = tmp_right - 1;                                                       //右边界前，1的个数
// 		r0_l0_left = pos_left - tmp_left;                                                  //左边界前，0的个数
// 		r0_l0_right = pos_right - tmp_right;                                               //右边界前，0的个数
// 	}
// 	else if(pos_right > -1)//pos_left = -1
// 	{
// 		r1_l0_right = r->Rank(pos_right)-1;                                          
// 		r0_l0_right = pos_right - r1_l0_right - 1;                                         //同样传出这四个变量
// 		r1_l0_left = pos_left;
// 		r0_l0_left = pos_left;
// 	}
// 	else
// 	{	
// 		for(int i = 0;i < 4;++i)
// 		{
// 			rank_left[i] = 0;
// 			rank_right[i] = 0;
// 		}
// 		return ;                                                                          //直接结束，返回rank_left 和right数组 都为0
// 	}

// 	if(r1_l0_left == r1_l0_right)//early stop                                             //左右边界，1的个数相等
// 	{
// 		r_flag = 0;                                                                       //标记全1，不需要后续rank
// 		for(int i = 2;i < 4;++i)
// 		{
// 			rank_left[i] = r1_l0_left + 1;                                                //更新rank_left 和right数组
// 			rank_right[i] = r1_l0_left + 1;
// 			//return ;
// 		}
// 	}
// 	if(r0_l0_left == r0_l0_right)//early stop                                             //左右边界，0的个数相等
// 	{
// 		l_flag = 0;                                                                       //标记全0，不需要后续rank
// 		for(int i = 0;i < 2;++i)
// 		{
// 			rank_left[i] = r0_l0_left + 1;                                                //更新rank_left 和right数组
// 			rank_right[i] = r0_l0_left + 1;
// 			//return ;
// 		}
// 	}

// 	if(r_flag)                                                                            //非全1（同开头的过程）
// 	{
// 		BitMap *rr = r->Right();                                                          //向右（第二层）
// 		fm_int r1_lr_left, r1_lr_right, r0_lr_left, r0_lr_right;
// 		fm_int tmp_left, tmp_right;	
// 		if(r1_l0_left > -1 && r1_l0_right > -1) 
// 		{
// 			int dollar_l = 0, dollar_r = 0;
// 			rr->Rank(r1_l0_left, r1_l0_right, tmp_left, tmp_right, dollar_l, dollar_r);
// 			r1_lr_left = tmp_left -1;
// 			r1_lr_right = tmp_right -1;
// 			r0_lr_left = r1_l0_left - tmp_left - dollar_l;
// 			r0_lr_right = r1_l0_right - tmp_right - dollar_r;
// 		}
// 		else if(r1_l0_right > -1)
// 		{
// 			int dollar = 0;
// 			r1_lr_right = rr->Rank(r1_l0_right,dollar,1) - 1;
// 			r0_lr_right = r1_l0_right - r1_lr_right - 1 - dollar;
// 			r1_lr_left = r1_l0_left;
// 			r0_lr_left = r1_l0_left;
// 		}
// 		else
// 		{
// 			for(int i = 2;i < 4;++i)
// 			{
// 				rank_left[i] = 0;
// 				rank_right[i] = 0;
// 			}
// 		} 
// 		rank_left[2] = r0_lr_left + 1;
// 		rank_right[2] = r0_lr_right + 1;
// 		rank_left[3] = r1_lr_left + 1;
// 		rank_right[3] = r1_lr_right + 1;
// 	}
// 	if(l_flag)
// 	{
// 		BitMap *rl = r->Left();                                                           //向左（第二层）
// 		fm_int r1_ll_left, r1_ll_right, r0_ll_left, r0_ll_right;
// 		fm_int tmp_left, tmp_right;	
// 		if(r0_l0_left > -1 && r0_l0_right > -1) 
// 		{
// 			rl->Rank(r0_l0_left, r0_l0_right, tmp_left, tmp_right);
// 			r1_ll_left = tmp_left -1;
// 			r1_ll_right = tmp_right -1;
// 			r0_ll_left = r0_l0_left - tmp_left;
// 			r0_ll_right = r0_l0_right - tmp_right;
// 		}
// 		else if(r0_l0_right > -1)
// 		{
// 			r1_ll_right = rl->Rank(r0_l0_right)-1;
// 			r0_ll_right = r0_l0_right - r1_ll_right - 1;
// 			r1_ll_left = r0_l0_left;
// 			r0_ll_left = r0_l0_left;
// 		}
// 		else
// 		{
// 			for(int i = 0;i < 2;++i)
// 			{
// 				rank_left[i] = 0;
// 				rank_right[i] = 0;
// 			}
// 		}
// 		rank_left[0] = r0_ll_left + 1;
// 		rank_right[0] = r0_ll_right + 1;
// 		rank_left[1] = r1_ll_left + 1;
// 		rank_right[1] = r1_ll_right + 1;
// 	}
// 	return ;
// }

unsigned char* ABS_FM::Extracting(fm_int pos,fm_int len)
{
	if(pos + len > n-1 || pos <0)
	{
		cout<<pos<<"  "<<len<<endl;
		cout<<pos+len<<" "<<n-1<<" "<<pos<<endl;
		cout<<"ABS_FM::Extracting  error parmaters"<<endl;
		return NULL;
	}

	unsigned char * sequence=new unsigned char[len+1];
	sequence[len]='\0';
	fm_int end = pos + len;
	fm_int anchor = 0;
	int overloop = 0;
	int step = this->D;
	overloop = end%step; //距离采样点多远
	anchor = end/step; //上取整
	if(overloop!=0){
		anchor++;
		if(anchor!=RankL->GetNum()-1){
			overloop=step-overloop;
		}else{
			overloop=(n-1)%step-overloop;
		}
	}
	fm_int i= RankL->GetValue(anchor);
	i=newD->select(i); //采样点处的排位
	for(int j=0;j<overloop;j++) //跳到pos+len处的排位
		i = LF(i);

	for(int j=0;j<len;j++)
	{
		sequence[len-1-j]=L(i); //索引上取字符
		i = LF(i);
	}
	return sequence;
}


fm_int ABS_FM::Lookup(fm_int i) //值采样，未使用，2023-4-12
{
	int step = 0;
	int D = this->D;
	while(newD->getD(i)!=1)
	{
		i=LF(i);
		step =step +1;
	}
	fm_int r=newD->rank(i);
	fm_int sal=0;
	if((n-1)%D==0){
		sal=SAL->GetValue(r-1);
	}else{
		sal=SAL->GetValue(r-2);
	}
	return (sal*D+step)%n;
}

// fm_int ABS_FM::FM_Lookup(fm_int i)         //值采样（仍有bug 20220913）
// {
// 	cout<<"locate...i:"<<i<<endl;
// 	int step = 0;
// 	int D = this->D;
// 	while(newD->getD(i)!=1)                //值采样
// 	{
// 		printf("jump:%d  i:%d \n",step,i);
// 		// cout<<"jump:"<<step<<"  i:"<<i<<endl;
// 		i=LF(i);                           //LF映射
// 		step =step +1;
// 	}
// 	fm_int r=newD->rank(i);
// 	fm_int sal=0;
// 	if((n-1)%D==0){                        //处理边界值
// 		sal=SAL->GetValue(r-1);
// 	}else{
// 		sal=SAL->GetValue(r-2);
// 	}
// 	cout<<"step:"<<step;
// 	int j=5;
// 	for(int i=-j; i<j; i++){
// 		cout<<" sa"<<i<<":"<<SAL->GetValue(r+i)*D;
// 	}
// 	cout<<endl;
// 	return (sal*D+step)%n;
// }

fm_int ABS_FM::FM_Lookup(fm_int i)          //下标采样
{
	/* ///检查C表
	for(int i=0; i<=16; i++){
		cout<<i<<"... "<<this->C[i]<<endl;
	}
	cout<<"sa0_index:"<<sa0_index<<endl;
	// getchar();
	/// */
	int step = 0;
	int D = this->D;                       //SA采样率
	// cout<<"locate...i:"<<i<<"  D:"<<D<<endl;
	// while(newD->getD(i)!=1)                //值采样
	while(i%D !=0)
	{
		// cout<<newD->getD(i)<<"  ";
		i =myLF(i);                           //LF映射
		step =step +1;
		// cout<<"jump:"<<step<<"  i:"<<i<<endl; //1在右树，0在左树，方便定位错误
	}
	// fm_int r=newD->rank(i);
	// fm_int sal=0;
	// if((n-1)%D==0){                        //处理边界值
	// 	sal=SAL->GetValue(r-1);
	// }else{
	// 	sal=SAL->GetValue(r-2);
	// }
	fm_int sal =SAL->GetValue(i/D); //这里报错 336672  169750249 InArray:GetValue: index out of boundary
	// cout<<"sal:"<<sal<<endl;
	/* int j=5;
	for(int k=-j; k<j; k++){ //打印附件的sal，看是否有正确的
		cout<<" sa"<<k<<":"<<SAL->GetValue(i/D +k);
	}
	cout<<endl; */
	return (sal+step)%n;
}

//返回L串中c字符在位置pos_left 和pos_right之前出现的次数，结果由rank_left 和rank_right带回. （count用）
void ABS_FM::Occ(unsigned char c,fm_int pos_left,fm_int pos_right,fm_int &rank_left,fm_int &rank_right)
{
	BitMap *r = root;
	int level=0;
	char code = '0';
	while(r->Left())                                                        //循环，直到叶子节点
	{
		code = codeTable[c][level];
		/* cout<<"print codeTable"<<endl;
		for(int i=0; i<16; i++){                                            //打印编码表codeTable
			cout<<"i..."<<i<<"...";
			for(int j=0; j<10; j++){
				int tmp=(int)codeTable[i][j]-48;
				if(tmp==0 || tmp==1){
					cout<<tmp;                                              //0放"0"，对应asc码48
				}
			}
			cout<<endl;
		}
		getchar(); */
		if(code == '1')                                                     //向左
		{
			if(pos_left>-1 && pos_right >-1) 
			{
				r->Rank(pos_left,pos_right,rank_left,rank_right);           //传出 rank_left和 rank_right
				pos_left = rank_left -1;
				pos_right = rank_right -1;                                  //用传出结果，更新pos_left，作为下轮输入
			}
			else if(pos_right > -1)
			{
				pos_right=r->Rank(pos_right)-1;
			}
			else
			{
				break;
			}
			//-----------kkzone-----------
			if(pos_left==pos_right)
			{
				rank_left=pos_left+1;
				rank_right=pos_right+1;
				return;
			}
			//-----------kkzone-----------
			r= r->Right();
		}
		else                                                                //向右
		{
			if(pos_left>-1 && pos_right >-1)
			{
				r->Rank(pos_left,pos_right,rank_left,rank_right);
				pos_left = (pos_left+1) - rank_left-1;
				pos_right= (pos_right+1)- rank_right-1;
			}
			else if(pos_right > -1)
			{
				pos_right = (pos_right+1)-r->Rank(pos_right)-1;
			}
			else
			{
				break;
			}
			//-----------kkzone-----------
			if(pos_left==pos_right)
			{
				rank_left=pos_left+1;
				rank_right=pos_right+1;
				return;
			}
			//-----------kkzone-----------
			r=r->Left();
		}
		level++;
	}
	rank_left = pos_left+1;
	rank_right= pos_right+1;
	return ;
}

fm_int ABS_FM::Occ(unsigned char c,fm_int pos)                                    //LF跳步用，更新排位
{
	BitMap * r = root;
	int level = 0;
	char code ='0';
	while(r->Left() && pos > -1)
	{
		code = codeTable[c][level];                                              //[c]哪个字符，[level]第几层
		if(code == '1')
		{
			pos = r->Rank(pos) - 1;
			r = r->Right();
		}
		else
		{
			pos = (pos+1) - r->Rank(pos)-1;
			r = r->Left();
		}
		level = level +1;
	}
	return pos+1;
}

fm_int ABS_FM::Occ(unsigned char c,fm_int pos,BitMap* root)
{
	BitMap * r = root;
	int level = 0;
	char code ='0';
	while(r->Left() && pos > -1)
	{
		code = codeTable[c][level];
		if(code == '1')
		{
			pos = r->Rank(pos) - 1;
			r = r->Right();
		}
		else
		{
			pos = (pos+1) - r->Rank(pos)-1;
			r = r->Left();
		}
		level = level +1;
	}
	return pos+1;
}

//辅助Occ8p，本身缺少newD->rank2(pos)
void ABS_FM::Occ1pAux(unsigned char c, fm_int pos1, fm_int pos2, fm_int &occ1, fm_int &occ2) //传入c为gray码
{
	BitMap *r=lroot;
	int level=0;
	char code = '0';
	while(r->Left())                                             //循环，直到叶子节点
	{
		code = codeTable[c][level];
		if(code == '1')                                          //向左
		{
			if(pos1 >-1 && pos2 >-1) {
				r->Rank(pos1, pos2, occ1, occ2);                 //传出 rank_left和 rank_right
				pos1 = occ1 -1;
				pos2 = occ2 -1;                                  //用传出结果，更新pos_left，作为下轮输入
			}else if(pos2 > -1) {
				pos2 = r->Rank(pos2) - 1;
			}else{
				break;
			}
			//-----------kkzone-----------
			if(pos1 == pos2){
				occ1 = pos1 + 1;
				occ2 = pos2 + 1;
				return;
			}
			//-----------kkzone-----------
			r = r->Right();
		}
		else                                                     //向右
		{
			if(pos1 > -1 && pos2 > -1) {
				r->Rank(pos1, pos2, occ1, occ2);
				pos1 = (pos1 + 1) - occ1 - 1;
				pos2 = (pos2 + 1) - occ2 - 1;
			}else if(pos2 > -1) {
				pos2 = (pos2 + 1) - r->Rank(pos2) - 1;
			}else {
				break;
			}
			//-----------kkzone-----------
			if(pos1 == pos2) {
				occ1 = pos1 + 1;
				occ2 = pos2 + 1;
				return;
			}
			//-----------kkzone-----------
			r = r->Left();
		}
		level++;
	}
	occ1 = pos1 + 1;
	occ2 = pos2 + 1;
	return ;
}

//myLF中做locate定位用，2023-4-12
fm_int ABS_FM::Occ1(fm_int i)
{
	fm_int ii =i;
	BitMap *r;
	int bit =0;
	fm_int rank =0;
	unsigned char c =0; //传出i处的字符，用于在C表取值
	if(newD->getD(i)){  //1.去右树
		r = rroot;
		if(i>0) i = newD->rank2(i) - 1; //默认
	}else{              //2.去左树
		r = lroot;
		if(i>0) i = i - newD->rank2(i);
	}
	while(r->Left())    //一直走到叶子节点
	{
		rank = r->Rank(i,bit);
		if(bit==1)
		{
			i = rank -1;
			r =r ->Right();
		}else
		{
			i = (i+1) - rank -1;
			r = r->Left();
		}
	}
	c = r->Label();
	if(c==0 && ii>sa0_index) i--;    //说明bwt中还是放入结束符0的，2023-4-13，调试正确
	i =i + ((c<10)?C[c]:C[c-1]) + 1; //加上C表（由于N随机化，编码不存在，所以大于10向前顺移1）
	return i;
}

void ABS_FM::Occ1p(unsigned char c, fm_int pos1, fm_int pos2, fm_int &occ1, fm_int &occ2) //传入c为gray码
{
	if(pos1 == -1 && pos2 == n){                                 //初始字符，查C表
		occ1 = ((c<10)?C[c]:C[c-1]) + 1;                           //由于N随机化，编码不存在，所以大于10向前顺移1
		occ2 = (c<10)?C[c+1]:C[c];
		return;
	}

	BitMap *r;
	// cout<<"...1"<<endl;
	if(is_snp_$[c]){                                             //去左树，还是右树
		r = rroot;
		// cout<<"...2"<<endl;
		if(pos1>0) pos1 = newD->rank2(pos1) - 1;
		if(pos2>0) pos2 = newD->rank2(pos2) - 1;
		//-----------kkzone-----------
			if(pos1 == pos2){ //必须要判断相等，并提前返回，要不然会报错（20230104）
				occ1 = pos1 + 1;
				occ2 = pos2 + 1;
				occ1+=((c<10)?C[c]:C[c-1]) + 1; //提前中止，把C表加上，L要多加1
				occ2+= (c<10)?C[c]:C[c-1];
				return;
			}
		//-----------kkzone-----------
	}else{
		r = lroot;
		// cout<<"...3"<<" pos1:"<<pos1<<" pos2:"<<pos2<<endl;
		if(pos1>0) pos1 = pos1 - newD->rank2(pos1);
		if(pos2>0) pos2 = pos2 - newD->rank2(pos2);
		//-----------kkzone-----------
			if(pos1 == pos2){
				occ1 = pos1 + 1;
				occ2 = pos2 + 1;
				occ1+=((c<10)?C[c]:C[c-1]) + 1; //提前中止，把C表加上，L要多加1
				occ2+= (c<10)?C[c]:C[c-1];
				return;
			}
		//-----------kkzone-----------
	}

	// cout<<"pos1:"<<pos1<<" pos2:"<<pos2<<endl;
	int level=0;
	char code = '0';
	while(r->Left())                                             //循环，直到叶子节点
	{
		code = codeTable[c][level];
		if(code == '1')                                          //向左
		{
			// cout<<"...5"<<endl;
			if(pos1 >-1 && pos2 >-1) {
				r->Rank(pos1, pos2, occ1, occ2);                 //传出 rank_left和 rank_right
				pos1 = occ1 -1;
				pos2 = occ2 -1;                                  //用传出结果，更新pos_left，作为下轮输入
			}else if(pos2 > -1) {
				pos2 = r->Rank(pos2) - 1;
			}else{
				break;
			}
			//-----------kkzone-----------
			if(pos1 == pos2){
				occ1 = pos1 + 1;
				occ2 = pos2 + 1;
				occ1+=((c<10)?C[c]:C[c-1]) + 1; //提前中止，把C表加上，L要多加1
				occ2+= (c<10)?C[c]:C[c-1];
				return;
			}
			//-----------kkzone-----------
			r = r->Right();
		}
		else                                                     //向右
		{
			// cout<<"...6"<<endl;
			if(pos1 > -1 && pos2 > -1) {
				r->Rank(pos1, pos2, occ1, occ2);
				pos1 = (pos1 + 1) - occ1 - 1;
				pos2 = (pos2 + 1) - occ2 - 1;
			}else if(pos2 > -1) {
				pos2 = (pos2 + 1) - r->Rank(pos2) - 1;
			}else {
				break;
			}
			//-----------kkzone-----------
			if(pos1 == pos2) {
				occ1 = pos1 + 1;
				occ2 = pos2 + 1;
				occ1+=((c<10)?C[c]:C[c-1]) + 1; //提前中止，把C表加上，L要多加1
				occ2+= (c<10)?C[c]:C[c-1];
				return;
			}
			//-----------kkzone-----------
			r = r->Left();
		}
		level++;
	}
	// cout<<"...7"<<endl;
	occ1 = pos1 + 1; //pos1和pos2都在不断的更新
	occ2 = pos2 + 1;
	// cout<<"...L:"<<occ1<<" R:"<<occ2;
	occ1+=((c<10)?C[c]:C[c-1]) + 1; //加上C表（由于N随机化，编码不存在，所以大于10向前顺移1）
	occ2+= (c<10)?C[c]:C[c-1];
	// cout<<" C[]:"<<((c<10)?C[c]:C[c-1])<<" c:"<<(int)c<<endl;
	return ;
}

void ABS_FM::Occ8p(unsigned char c, fm_int pos1, fm_int pos2, fm_int* occ81, fm_int* occ82) //exact用（传入c为gray码）
{ 
	if(pos1 == -1 && pos2 == n){ //首个字符，查C表
		switch (c) //七个元素，N已经跳过
		{
		case 0: /*A*/ //对应read上取字符，0:A，1:G，2:C，3:T
			occ81[0]=C[8]+1; occ81[1]=C[9]+1; occ81[2]=C[10]+1; occ81[3]=C[11]+1; occ81[4]=C[12]+1; occ81[5]=C[13]+1; occ81[6]=C[14]+1; //MHVRDWA
			occ82[0]=C[9];   occ82[1]=C[10];  occ82[2]=C[11];   occ82[3]=C[12];   occ82[4]=C[13];   occ82[5]=C[14];   occ82[6]=C[15];
			break;
		case 2: /*C*/
			occ81[0]=C[4]+1; occ81[1]=C[5]+1; occ81[2]=C[6]+1; occ81[3]=C[7]+1; occ81[4]=C[8]+1; occ81[5]=C[9]+1; occ81[6]+=C[10]+1; //SBYCMHV
			occ82[0]=C[5];   occ82[1]=C[6];   occ82[2]=C[7];   occ82[3]=C[8];   occ82[4]=C[9];   occ82[5]=C[10];  occ82[6]+=C[11];
			break;
		case 1: /*G*/
			occ81[0]=C[2]+1; occ81[1]=C[3]+1; occ81[2]=C[4]+1; occ81[3]=C[5]+1; occ81[4]=C[10]+1; occ81[5]=C[11]+1; occ81[6]=C[12]+1; //KGSBVRD
			occ82[0]=C[3];   occ82[1]=C[4];   occ82[2]=C[5];   occ82[3]=C[6];   occ82[4]=C[11];   occ82[5]=C[12];   occ82[6]=C[13];
			break;
		case 3: /*T*/
			occ81[0]=C[1]+1; occ81[1]=C[2]+1; occ81[2]+=C[5]+1; occ81[3]+=C[6]+1; occ81[4]+=C[9]+1; occ81[5]+=C[12]+1; occ81[6]+=C[13]+1; //TKBYHDW
			occ82[0]=C[2];   occ82[1]=C[3];   occ82[2]+=C[6];   occ82[3]+=C[7];   occ82[4]+=C[10];  occ82[5]+=C[13];   occ82[6]+=C[14];
			break;
		default:
			break;
		}
		return;
	}

	long long p1R=0, p1L=0, p2R=0, p2L=0, occ1=0, occ2=0, occ111[11]={0}, occ112[11]={0};
	p1R=newD->rank2(pos1); p1L=pos1-p1R; p1R--;
	p2R=newD->rank2(pos2); p2L=pos2-p2R; p2R--;
	if(p1L !=p2L) Occ1pAux(nt4_gray[c], p1L, p2L, occ1, occ2); //传入A0 G1 C2 T3，转为A0 C1 G2 T3
	if(p1R !=p2R) Occ11p(p1R, p2R, occ111, occ112);
	switch (c)
	{
	case 0: /*A*/ //对应read上取字符，0:A，1:G，2:C，3:T
		occ81[0]=occ111[5]; occ81[1]=occ111[6]; occ81[2]=occ111[7]; occ81[3]=occ111[8]; occ81[4]=occ111[9]; occ81[5]=occ111[10]; occ81[6]=occ1; //MHVRDWA
		occ82[0]=occ112[5]; occ82[1]=occ112[6]; occ82[2]=occ112[7]; occ82[3]=occ112[8]; occ82[4]=occ112[9]; occ82[5]=occ112[10]; occ82[6]=occ2;
		occ81[0]+=C[8]+1; occ81[1]+=C[9]+1; occ81[2]+=C[10]+1; occ81[3]+=C[11]+1; occ81[4]+=C[12]+1; occ81[5]+=C[13]+1; occ81[6]+=C[14]+1;      //加上C表
		occ82[0]+=C[8];   occ82[1]+=C[9];   occ82[2]+=C[10];   occ82[3]+=C[11];   occ82[4]+=C[12];   occ82[5]+=C[13];   occ82[6]+=C[14];
		break;
	case 2:  /*C*/
		occ81[0]=occ111[2]; occ81[1]=occ111[3]; occ81[2]=occ111[4]; occ81[3]=occ1; occ81[4]=occ111[5]; occ81[5]=occ111[6]; occ81[6]=occ111[7];  //SBYCMHV
		occ82[0]=occ112[2]; occ82[1]=occ112[3]; occ82[2]=occ112[4]; occ82[3]=occ2; occ82[4]=occ112[5]; occ82[5]=occ112[6]; occ82[6]=occ112[7];
		// cout<<"0 ...L:"<<occ81[0]<<" R:"<<occ82[0]<<" C[]:"<<C[4]<<" c:"<<4<<endl;
		// cout<<"1 ...L:"<<occ81[1]<<" R:"<<occ82[1]<<" C[]:"<<C[5]<<" c:"<<5<<endl;
		// cout<<"2 ...L:"<<occ81[2]<<" R:"<<occ82[2]<<" C[]:"<<C[6]<<" c:"<<6<<endl;
		// cout<<"3 ...L:"<<occ81[3]<<" R:"<<occ82[3]<<" C[]:"<<C[7]<<" c:"<<7<<endl;
		// cout<<"4 ...L:"<<occ81[4]<<" R:"<<occ82[4]<<" C[]:"<<C[8]<<" c:"<<8<<endl;
		// cout<<"5 ...L:"<<occ81[5]<<" R:"<<occ82[5]<<" C[]:"<<C[9]<<" c:"<<9<<endl;
		// cout<<"6 ...L:"<<occ81[6]<<" R:"<<occ82[6]<<" C[]:"<<C[10]<<" c:"<<10<<endl;
		occ81[0]+=C[4]+1; occ81[1]+=C[5]+1; occ81[2]+=C[6]+1; occ81[3]+=C[7]+1; occ81[4]+=C[8]+1; occ81[5]+=C[9]+1; occ81[6]+=C[10]+1;          //加上C表
		occ82[0]+=C[4];   occ82[1]+=C[5];   occ82[2]+=C[6];   occ82[3]+=C[7];   occ82[4]+=C[8];   occ82[5]+=C[9];   occ82[6]+=C[10];
		break;
	case 1:  /*G*/
		occ81[0]=occ111[1]; occ81[1]=occ1; occ81[2]=occ111[2]; occ81[3]=occ111[3]; occ81[4]=occ111[7]; occ81[5]=occ111[8]; occ81[6]=occ111[9];  //KGSBVRD
		occ82[0]=occ112[1]; occ82[1]=occ2; occ82[2]=occ112[2]; occ82[3]=occ112[3]; occ82[4]=occ112[7]; occ82[5]=occ112[8]; occ82[6]=occ112[9];
		// cout<<"0 ...L:"<<occ81[0]<<" R:"<<occ82[0]<<" C[]:"<<C[2]<<" c:"<<2<<endl;
		// cout<<"1 ...L:"<<occ81[1]<<" R:"<<occ82[1]<<" C[]:"<<C[3]<<" c:"<<3<<endl;
		// cout<<"2 ...L:"<<occ81[2]<<" R:"<<occ82[2]<<" C[]:"<<C[4]<<" c:"<<4<<endl;
		// cout<<"3 ...L:"<<occ81[3]<<" R:"<<occ82[3]<<" C[]:"<<C[5]<<" c:"<<5<<endl;
		// cout<<"4 ...L:"<<occ81[4]<<" R:"<<occ82[4]<<" C[]:"<<C[10]<<" c:"<<10<<endl;
		// cout<<"5 ...L:"<<occ81[5]<<" R:"<<occ82[5]<<" C[]:"<<C[11]<<" c:"<<11<<endl;
		// cout<<"6 ...L:"<<occ81[6]<<" R:"<<occ82[6]<<" C[]:"<<C[12]<<" c:"<<12<<endl;
		occ81[0]+=C[2]+1; occ81[1]+=C[3]+1; occ81[2]+=C[4]+1; occ81[3]+=C[5]+1; occ81[4]+=C[10]+1; occ81[5]+=C[11]+1; occ81[6]+=C[12]+1;        //加上C表
		occ82[0]+=C[2];   occ82[1]+=C[3];   occ82[2]+=C[4];   occ82[3]+=C[5];   occ82[4]+=C[10];   occ82[5]+=C[11];   occ82[6]+=C[12];
		break;
	case 3: /*T*/
		occ81[0]=occ1; occ81[1]=occ111[1]; occ81[2]=occ111[3]; occ81[3]=occ111[4]; occ81[4]=occ111[6]; occ81[5]=occ111[9]; occ81[6]=occ111[10]; //TKBYHDW
		occ82[0]=occ2; occ82[1]=occ112[1]; occ82[2]=occ112[3]; occ82[3]=occ112[4]; occ82[4]=occ112[6]; occ82[5]=occ112[9]; occ82[6]=occ112[10];
		occ81[0]+=C[1]+1; occ81[1]+=C[2]+1; occ81[2]+=C[5]+1; occ81[3]+=C[6]+1; occ81[4]+=C[9]+1; occ81[5]+=C[12]+1; occ81[6]+=C[13]+1;         //加上C表
		occ82[0]+=C[1];   occ82[1]+=C[2];   occ82[2]+=C[5];   occ82[3]+=C[6];   occ82[4]+=C[9];   occ82[5]+=C[12];   occ82[6]+=C[13];
		break;
	default:
		break;
	}
}

void ABS_FM::Occ16p(fm_int pos1, fm_int pos2, fm_int* occ161, fm_int* occ162) //$TKGSBYCMHVRDWA（gray序），inexact用
{
	if(pos1 == -1 && pos2 == n){
		for(int i=0; i<16; i++){ //原本从1开始，这样导致0的值丢失，已纠正 2023-4-3
			if(i==10){
				occ161[10]=0; //19
				occ162[10]=-1; //18
				continue; //N直接跳过，初始值为0
			}
			occ161[i] = ((i<10)?C[i]:C[i-1]) + 1;                           //由于N随机化，编码不存在，所以大于10向前顺移1
			occ162[i] = (i<10)?C[i+1]:C[i];
		}
		return;
	}

	long long p1R=0, p1L=0, p2R=0, p2L=0, occ41[4]={0}, occ42[4]={0}, occ111[11]={0}, occ112[11]={0};
	// cout<<"newD...rank2...pos1: "<<pos1<<"  pos2:"<<pos2<<endl; //下面容易报错，rank2导致的20230201
	p1R=newD->rank2(pos1); p1L=pos1-p1R; p1R--; //左树的pos1、pos2
	p2R=newD->rank2(pos2); p2L=pos2-p2R; p2R--;
	// cout<<"Occ4p... p1L:"<<setw(10)<<left<<p1L<<"p2L:"<<setw(10)<<left<<p2L<<endl;
	// cout<<"plL:"<<p1L<<" p2L:"<<p2L<<" p1R:"<<p1R<<" p2R:"<<p2R<<endl;
	// if(p1L == p2L){ //ACGT
	// 	p1L = p1L + 1;
	// 	p2L = p2L + 1;
	// 	for(int i=1; i<16; i++){
	// 		occ161[i] =p1L + ((i<10)?C[i]:C[i-1]) + 1;
	// 		occ162[i] =p2L +  (i<10)?C[i]:C[i-1];
	// 	}
	// 	return;
	// }
	// if(p1R == p2R){
	// 	p1R = p1R + 1;
	// 	p2R = p2R + 1;
	// 	for(int i=1; i<16; i++){
	// 		occ161[i] =p1R + ((i<10)?C[i]:C[i-1]) + 1;
	// 		occ162[i] =p2R +  (i<10)?C[i]:C[i-1];
	// 	}
	// 	return;
	// }

	if(p1L !=p2L) Occ4p(p1L, p2L, occ41, occ42);            //ACGT
	// cout<<"Occ11p...p1R:"<<setw(10)<<left<<p1R<<"p2R:"<<setw(10)<<left<<p2R<<endl;
	if(p1R !=p2R) Occ11p(p1R, p2R, occ111, occ112);         //$KSBYMHVRDW（gray序）
	// cout<<"ok"<<endl;
	occ161[0]=occ111[0]; occ161[1]=occ41[3];  occ161[2]=occ111[1]; occ161[3]=occ41[2];   occ161[4]=occ111[2];  occ161[5]=occ111[3];  occ161[6]=occ111[4];   occ161[7]=occ41[1]; //10返回L:19 R:18（无符号型，无法区分0和-1）
	occ161[8]=occ111[5]; occ161[9]=occ111[6]; occ161[10]=19;       occ161[11]=occ111[7]; occ161[12]=occ111[8]; occ161[13]=occ111[9]; occ161[14]=occ111[10]; occ161[15]=occ41[0]; //$TKGSBYCMHVRDWA（gray序）
	occ162[0]=occ112[0]; occ162[1]=occ42[3];  occ162[2]=occ112[1]; occ162[3]=occ42[2];   occ162[4]=occ112[2];  occ162[5]=occ112[3];  occ162[6]=occ112[4];   occ162[7]=occ42[1];
	occ162[8]=occ112[5]; occ162[9]=occ112[6]; occ162[10]=18;       occ162[11]=occ112[7]; occ162[12]=occ112[8]; occ162[13]=occ112[9]; occ162[14]=occ112[10]; occ162[15]=occ42[0];

	occ161[0]+=C[0]+1; occ161[1]+=C[1]+1; occ161[2]+=C[2]+1; occ161[3]+=C[3]+1;   occ161[4]+=C[4]+1;   occ161[5]+=C[5]+1;   occ161[6]+=C[6]+1;   occ161[7]+=C[7]+1; //加上C表
	occ161[8]+=C[8]+1; occ161[9]+=C[9]+1;                    occ161[11]+=C[10]+1; occ161[12]+=C[11]+1; occ161[13]+=C[12]+1; occ161[14]+=C[13]+1; occ161[15]+=C[14]+1;
	occ162[0]+=C[0]; occ162[1]+=C[1]; occ162[2]+=C[2]; occ162[3]+=C[3];   occ162[4]+=C[4];   occ162[5]+=C[5];   occ162[6]+=C[6];   occ162[7]+=C[7]; 
	occ162[8]+=C[8]; occ162[9]+=C[9];                  occ162[11]+=C[10]; occ162[12]+=C[11]; occ162[13]+=C[12]; occ162[14]+=C[13]; occ162[15]+=C[14];    
}

void ABS_FM::Occ4p(fm_int pos1, fm_int pos2, fm_int *occ41, fm_int *occ42) //ACGT（为其它服务，左子树）
{   
	BitMap *r = lroot;
	int f_k1_h0=1, f_k0_h0=1;
	fm_int p1_k0_h0, p1_k1_h0, p2_k0_h0, p2_k1_h0;
	fm_int p1_k0_h1_1, p1_k1_h1_1, p2_k0_h1_1/*G*/, p2_k1_h1_1/*T*/;
	fm_int p1_k0_h1_0, p1_k1_h1_0, p2_k0_h1_0/*A*/, p2_k1_h1_0/*C*/;
	
	//第0层
	if(pos1 > -1 && pos2 > -1){
		r->Rank(pos1, pos2, p1_k1_h0, p2_k1_h0);
		p1_k0_h0 = pos1 - p1_k1_h0; p2_k0_h0 = pos2 - p2_k1_h0;
	}else if(pos2 > -1){
		p2_k1_h0 = r->Rank(pos2); p2_k0_h0 = pos2 - p2_k1_h0;
		p1_k1_h0 = 0; p1_k0_h0 = -1;
	}else{
		memset(occ41, 0, 4*sizeof(fm_int));
		memset(occ42, 0, 4*sizeof(fm_int));
		return;
	}
	p1_k1_h0--; p2_k1_h0--;
	if(p1_k1_h0 == p2_k1_h0){                                                      //早停（两pos间的occ，不需要最终排位）
		f_k1_h0 = 0;
		p1_k1_h1_1 = p2_k1_h1_1 = p1_k0_h1_1 = p2_k0_h1_1 = -1;                    //整个右支的TG都置-1
	}
	if(p1_k0_h0 == p2_k0_h0){
		f_k0_h0 = 0;
		p1_k1_h1_0 = p2_k1_h1_0 = p1_k0_h1_0 = p2_k0_h1_0 = -1;                    //整个左支的AC都置-1
	}

	//第1层
	if(f_k1_h0){
		BitMap *rRoot = r->Right();
		if(p1_k1_h0 > -1 && p2_k1_h0 > -1){
			rRoot->Rank(p1_k1_h0, p2_k1_h0, p1_k1_h1_1, p2_k1_h1_1);                   //11 T
			p1_k0_h1_1 = p1_k1_h0 - p1_k1_h1_1; p2_k0_h1_1 = p2_k1_h0 - p2_k1_h1_1;    //10 G
		}else if(p2_k1_h0 > -1){
			p2_k1_h1_1 = rRoot->Rank(p2_k1_h0); p2_k0_h1_1 = p2_k1_h0 - p2_k1_h1_1;
			p1_k1_h1_1 = 0; p1_k0_h1_1 = -1;
		}else{
			p1_k1_h1_1 = 0; p1_k0_h1_1 = -1;
			p2_k1_h1_1 = 0; p2_k0_h1_1 = -1;
		}
		p1_k1_h1_1--; p2_k1_h1_1--;
	}
	if(f_k0_h0){
		BitMap *lRoot = r->Left();
		if(p1_k0_h0 > -1 && p2_k0_h0 > -1){
			lRoot->Rank(p1_k0_h0, p2_k0_h0, p1_k1_h1_0, p2_k1_h1_0);                   //01 C
			p1_k0_h1_0 = p1_k0_h0 - p1_k1_h1_0; p2_k0_h1_0 = p2_k0_h0 - p2_k1_h1_0;    //00 A
		}else if(p2_k0_h0 > -1){
			p2_k1_h1_0 = lRoot->Rank(p2_k0_h0); p2_k0_h1_0 = p2_k0_h0 - p2_k1_h1_0;
			p1_k1_h1_0 = 0; p1_k0_h1_0 = -1;
		}else{
			p1_k1_h1_0 = 0; p1_k0_h1_0 = -1;
			p2_k1_h1_0 = 0; p2_k0_h1_0 = -1;
		}
		p1_k1_h1_0--; p2_k1_h1_0--;
	}

	occ41[0] = p1_k0_h1_0; occ41[1] = p1_k1_h1_0; occ41[2] = p1_k0_h1_1; occ41[3] = p1_k1_h1_1; //ACGT
	occ42[0] = p2_k0_h1_0; occ42[1] = p2_k1_h1_0; occ42[2] = p2_k0_h1_1; occ42[3] = p2_k1_h1_1;

	for(int i=0; i<4; i++){
		occ41[i]++; occ42[i]++;
	}
}

void ABS_FM::Occ11p(fm_int pos1, fm_int pos2, fm_int *occ111, fm_int *occ112) //$KSBYMHVRDW（gray序），为其它服务，右子树
{	
	BitMap *r = rroot, *rRoot, *lRoot, *rRoot0, *lRoot0, *rRoot01, *rRoot00, *lRoot001, *rRoot0010, *lRoot0010;
	int f_k1_h0=1, f_k0_h0=1, f_k1_h1_0=1, f_k0_h1_0=1, f_k1_h2_01=1, f_k1_h2_00=1, f_k0_h3_001=1, f_k1_h4_0010=1, f_k0_h4_0010=1; //9个开关
	fm_int p1_k0_h0, p1_k1_h0, p2_k0_h0, p2_k1_h0;
	fm_int p1_k0_h1_1, p1_k1_h1_1, p2_k0_h1_1/*R*/, p2_k1_h1_1/*Y*/;
	fm_int p1_k0_h1_0, p1_k1_h1_0, p2_k0_h1_0, p2_k1_h1_0;
	fm_int p1_k0_h2_01, p1_k1_h2_01, p2_k0_h2_01/*$*/, p2_k1_h2_01;
	fm_int p1_k0_h2_00, p1_k1_h2_00, p2_k0_h2_00/*S*/, p2_k1_h2_00;
	fm_int p1_k0_h3_011, p1_k1_h3_011, p2_k0_h3_011/*K*/, p2_k1_h3_011/*M*/;
	fm_int p1_k0_h3_001, p1_k1_h3_001, p2_k0_h3_001, p2_k1_h3_001/*W*/;
	fm_int p1_k0_h4_0010, p1_k1_h4_0010, p2_k0_h4_0010, p2_k1_h4_0010;
	fm_int p1_k0_h5_00101, p1_k1_h5_00101, p2_k0_h5_00101/*B*/, p2_k1_h5_00101/*V*/;
	fm_int p1_k0_h5_00100, p1_k1_h5_00100, p2_k0_h5_00100/*D*/, p2_k1_h5_00100/*H*/; 

	//第0层
	if(pos1 > -1 && pos2 > -1){
		r->Rank(pos1, pos2, p1_k1_h0, p2_k1_h0);
		p1_k0_h0 = pos1 - p1_k1_h0; p2_k0_h0 = pos2 - p2_k1_h0;
	}else if(pos2 > -1){
		p2_k1_h0 = r->Rank(pos2); p2_k0_h0 = pos2 - p2_k1_h0;
		p1_k1_h0 = 0; p1_k0_h0 = -1;
	}else{
		memset(occ111, 0, 11*sizeof(fm_int));
		memset(occ112, 0, 11*sizeof(fm_int));
		return;
	}
	p1_k1_h0--; p2_k1_h0--;
	if(p1_k1_h0 == p2_k1_h0){                                                                                      //早停（两pos间的occ，不需要最终排位）
		f_k1_h0 = 0;
		p1_k1_h1_1 = p2_k1_h1_1 = p1_k0_h1_1 = p2_k0_h1_1 = -1; //RY
	}
	if(p1_k0_h0 == p2_k0_h0){
		f_k0_h0 = f_k1_h1_0 = f_k0_h1_0 = f_k1_h2_01 = f_k1_h2_00 = f_k0_h3_001 = f_k1_h4_0010 = f_k0_h4_0010 = 0; //关8个
		p1_k1_h3_011 = p2_k1_h3_011 = p1_k0_h3_011 = p2_k0_h3_011 = p1_k0_h2_01 = p2_k0_h2_01 = p1_k1_h3_001 = p2_k1_h3_001 = p1_k1_h5_00101 = p2_k1_h5_00101 \
		= p1_k0_h5_00101 = p2_k0_h5_00101 = p1_k1_h5_00100 = p2_k1_h5_00100 = p1_k0_h5_00100 = p2_k0_h5_00100 = p1_k0_h2_00 = p2_k0_h2_00 = -1; //MK$WVBHDS
	}

	//第1层
	if(f_k1_h0){
		rRoot = r->Right();
		if(p1_k1_h0 > -1 && p2_k1_h0 > -1){
			rRoot->Rank(p1_k1_h0, p2_k1_h0, p1_k1_h1_1, p2_k1_h1_1);                                                   //11 Y
			p1_k0_h1_1 = p1_k1_h0 - p1_k1_h1_1; p2_k0_h1_1 = p2_k1_h0 - p2_k1_h1_1;                                    //10 R
		}else if(p2_k1_h0 > -1){
			p2_k1_h1_1 = rRoot->Rank(p2_k1_h0); p2_k0_h1_1 = p2_k1_h0 - p2_k1_h1_1;
			p1_k1_h1_1 = 0; p1_k0_h1_1 = -1;
		}else{
			p1_k1_h1_1 = 0; p1_k0_h1_1 = -1;
			p2_k1_h1_1 = 0; p2_k0_h1_1 = -1;
		}
		p1_k1_h1_1--; p2_k1_h1_1--;
	}
	if(f_k0_h0){
		lRoot = r->Left();
		if(p1_k0_h0 > -1 && p2_k0_h0 > -1){
			lRoot->Rank(p1_k0_h0, p2_k0_h0, p1_k1_h1_0, p2_k1_h1_0);
			p1_k0_h1_0 = p1_k0_h0 - p1_k1_h1_0; p2_k0_h1_0 = p2_k0_h0 - p2_k1_h1_0;
		}else if(p2_k0_h0 > -1){
			p2_k1_h1_0 = lRoot->Rank(p2_k0_h0); p2_k0_h1_0 = p2_k0_h0 - p2_k1_h1_0;
			p1_k1_h1_0 = 0; p1_k0_h1_0 = -1;
		}else{
			p1_k1_h1_0 = 0; p1_k0_h1_0 = -1;
			p2_k1_h1_0 = 0; p2_k0_h1_0 = -1; //由于左右相等，后续分支开关关闭
		}
		p1_k1_h1_0--; p2_k1_h1_0--;
		if(p1_k1_h1_0 == p2_k1_h1_0){
			f_k1_h1_0 = f_k1_h2_01 = 0;
			p1_k1_h3_011 = p2_k1_h3_011 = p1_k0_h3_011 = p2_k0_h3_011 = p1_k0_h2_01 = p2_k0_h2_01 = -1; //MK$
		}
		if(p1_k0_h1_0 == p2_k0_h1_0){
			f_k0_h1_0 = f_k1_h2_00 = f_k0_h3_001 = f_k1_h4_0010 = f_k0_h4_0010 = 0;
			p1_k1_h3_001 = p2_k1_h3_001 = p1_k1_h5_00101 = p2_k1_h5_00101 = p1_k0_h5_00101 = p2_k0_h5_00101 = p1_k1_h5_00100 = p2_k1_h5_00100 \
			= p1_k0_h5_00100 = p2_k0_h5_00100 = p1_k0_h2_00 = p2_k0_h2_00 = -1; //WVBHDS
		}
	}

	//第2层
	if(f_k1_h1_0){
		rRoot0 = lRoot->Right();
		if(p1_k1_h1_0 > -1 && p2_k1_h1_0 > -1){
			rRoot0->Rank(p1_k1_h1_0, p2_k1_h1_0, p1_k1_h2_01, p2_k1_h2_01);
			p1_k0_h2_01 = p1_k1_h1_0 - p1_k1_h2_01; p2_k0_h2_01 = p2_k1_h1_0 - p2_k1_h2_01;                //010 $
		}else if(p2_k1_h1_0 > -1){
			p2_k1_h2_01 = rRoot0->Rank(p2_k1_h1_0); p2_k0_h2_01 = p2_k1_h1_0 - p2_k1_h2_01;
			p1_k1_h2_01 = 0; p1_k0_h2_01 = -1;
		}else{
			p1_k1_h2_01 = 0; p1_k0_h2_01 = -1;
			p2_k1_h2_01 = 0; p2_k0_h2_01 = -1;
		}
		p1_k1_h2_01--; p2_k1_h2_01--;
		if(p1_k1_h2_01 == p2_k1_h2_01){
			f_k1_h2_01 = 0;
			p1_k1_h3_011 = p2_k1_h3_011 = p1_k0_h3_011 = p2_k0_h3_011 = -1; //MK
		}
	}
	if(f_k0_h1_0){
		lRoot0 = lRoot->Left();
		if(p1_k0_h1_0 > -1 && p2_k0_h1_0 > -1){
			lRoot0->Rank(p1_k0_h1_0, p2_k0_h1_0, p1_k1_h2_00, p2_k1_h2_00);
			p1_k0_h2_00 = p1_k0_h1_0 - p1_k1_h2_00; p2_k0_h2_00 = p2_k0_h1_0 - p2_k1_h2_00;                //000 S
		}else if(p2_k0_h1_0 > -1){
			p2_k1_h2_00 = lRoot0->Rank(p2_k0_h1_0); p2_k0_h2_00 = p2_k0_h1_0 - p2_k1_h2_00;
			p1_k1_h2_00 = 0; p1_k0_h2_00 = -1;
		}else{
			p1_k1_h2_00 = 0; p1_k0_h2_00 = -1;
			p2_k1_h2_00 = 0; p2_k0_h2_00 = -1;
		}
		p1_k1_h2_00--; p2_k1_h2_00--;
		if(p1_k1_h2_00 == p2_k1_h2_00){
			f_k1_h2_00 = f_k0_h3_001 = f_k1_h4_0010 = f_k0_h4_0010 = 0;
			p1_k1_h3_001 = p2_k1_h3_001 = p1_k1_h5_00101 = p2_k1_h5_00101 = p1_k0_h5_00101 = p2_k0_h5_00101 = p1_k1_h5_00100 = p2_k1_h5_00100 = p1_k0_h5_00100 = p2_k0_h5_00100 = -1; //WVBHD
		}
	}

	//第3层
	if(f_k1_h2_01){
		rRoot01 = rRoot0->Right();
		if(p1_k1_h2_01 > -1 && p2_k1_h2_01 > -1){
			rRoot01->Rank(p1_k1_h2_01, p2_k1_h2_01, p1_k1_h3_011, p2_k1_h3_011);                           //0111 M
			p1_k0_h3_011 = p1_k1_h2_01 - p1_k1_h3_011; p2_k0_h3_011 = p2_k1_h2_01 - p2_k1_h3_011;          //0110 K
		}else if(p2_k1_h2_01 > -1){
			p2_k1_h3_011 = rRoot01->Rank(p2_k1_h2_01); p2_k0_h3_011 = p2_k1_h2_01 - p2_k1_h3_011;
			p1_k1_h3_011 = 0; p1_k0_h3_011 = -1;
		}else{
			p1_k1_h3_011 = 0; p1_k0_h3_011 = -1;
			p2_k1_h3_011 = 0; p2_k0_h3_011 = -1;
		}
		p1_k1_h3_011--; p2_k1_h3_011--;
	}
	if(f_k1_h2_00){
		rRoot00 = lRoot0->Right();
		if(p1_k1_h2_00 > -1 && p2_k1_h2_00 > -1){
			rRoot00->Rank(p1_k1_h2_00, p2_k1_h2_00, p1_k1_h3_001, p2_k1_h3_001);                           //0011 W
			p1_k0_h3_001 = p1_k1_h2_00 - p1_k1_h3_001; p2_k0_h3_001 = p2_k1_h2_00 - p2_k1_h3_001;
		}else if(p2_k1_h2_00 > -1){
			p2_k1_h3_001 = rRoot00->Rank(p2_k1_h2_00); p2_k0_h3_001 = p2_k1_h2_00 - p2_k1_h3_001;
			p1_k1_h3_001 = 0; p1_k0_h3_001 = -1;
		}else{
			p1_k1_h3_001 = 0; p1_k0_h3_001 = -1;
			p2_k1_h3_001 = 0; p2_k0_h3_001 = -1;
		}
		p1_k1_h3_001--; p2_k1_h3_001--;
		if(p1_k0_h3_001 == p2_k0_h3_001){
			f_k0_h3_001 = f_k1_h4_0010 = f_k0_h4_0010 = 0;
			p1_k1_h5_00101 = p2_k1_h5_00101 = p1_k0_h5_00101 = p2_k0_h5_00101 = p1_k1_h5_00100 = p2_k1_h5_00100 = p1_k0_h5_00100 = p2_k0_h5_00100 = -1; //VBHD
		}
	}

	//第4层
	if(f_k0_h3_001){
		lRoot001 = rRoot00->Left();
		if(p1_k0_h3_001 > -1 && p2_k0_h3_001 > -1){
			lRoot001->Rank(p1_k0_h3_001, p2_k0_h3_001, p1_k1_h4_0010, p2_k1_h4_0010);
			p1_k0_h4_0010 = p1_k0_h3_001 - p1_k1_h4_0010; p2_k0_h4_0010 = p2_k0_h3_001 - p2_k1_h4_0010;
		}else if(p2_k0_h3_001 > -1){
			p2_k1_h4_0010 = lRoot001->Rank(p2_k0_h3_001); p2_k0_h4_0010 = p2_k0_h3_001 - p2_k1_h4_0010;
			p1_k1_h4_0010 = 0; p1_k0_h4_0010 = -1;
		}else{
			p1_k1_h4_0010 = 0; p1_k0_h4_0010 = -1;
			p2_k1_h4_0010 = 0; p2_k0_h4_0010 = -1;
		}
		p1_k1_h4_0010--; p2_k1_h4_0010--;
		if(p1_k1_h4_0010 == p2_k1_h4_0010){
			f_k1_h4_0010 = 0;
			p1_k1_h5_00101 = p2_k1_h5_00101 = p1_k0_h5_00101 = p2_k0_h5_00101 = -1; //VB
		}
		if(p1_k0_h4_0010 == p2_k0_h4_0010){
			f_k0_h4_0010 = 0;
			p1_k1_h5_00100 = p2_k1_h5_00100 = p1_k0_h5_00100 = p2_k0_h5_00100 = -1; //HD
		}
	}

	//第5层
	if(f_k1_h4_0010){
		rRoot0010 = lRoot001->Right();
		if(p1_k1_h4_0010 > -1 && p2_k1_h4_0010 > -1){
			rRoot0010->Rank(p1_k1_h4_0010, p2_k1_h4_0010, p1_k1_h5_00101, p2_k1_h5_00101);                       //001011 V
			p1_k0_h5_00101 = p1_k1_h4_0010 - p1_k1_h5_00101; p2_k0_h5_00101 = p2_k1_h4_0010 - p2_k1_h5_00101;    //001010 B
		}else if(p2_k1_h4_0010 > -1){
			p2_k1_h5_00101 = rRoot0010->Rank(p2_k1_h4_0010); p2_k0_h5_00101 = p2_k1_h4_0010 - p2_k1_h5_00101;
			p1_k1_h5_00101 = 0; p1_k0_h5_00101 = -1;
		}else{
			p1_k1_h5_00101 = 0; p1_k0_h5_00101 = -1;
			p2_k1_h5_00101 = 0; p2_k0_h5_00101 = -1;
		}
		p1_k1_h5_00101--; p2_k1_h5_00101--;
	}
	if(f_k0_h4_0010){
		lRoot0010 = lRoot001->Left();
		if(p1_k0_h4_0010 > -1 && p2_k0_h4_0010 > -1){
			lRoot0010->Rank(p1_k0_h4_0010, p2_k0_h4_0010, p1_k1_h5_00100, p2_k1_h5_00100);                       //001001 H
			p1_k0_h5_00100 = p1_k0_h4_0010 - p1_k1_h5_00100; p2_k0_h5_00100 = p2_k0_h4_0010 - p2_k1_h5_00100;    //001000 D
		}else if(p2_k0_h4_0010 > -1){
			p2_k1_h5_00100 = lRoot0010->Rank(p2_k0_h4_0010); p2_k0_h5_00100 = p2_k0_h4_0010 - p2_k1_h5_00100;
			p1_k1_h5_00100 = 0; p1_k0_h5_00100 = -1;
		}else{
			p1_k1_h5_00100 = 0; p1_k0_h5_00100 = -1;
			p2_k1_h5_00100 = 0; p2_k0_h5_00100 = -1;
		}
		p1_k1_h5_00100--; p2_k1_h5_00100--;
	}

	occ111[0] = p1_k0_h2_01; occ111[1] = p1_k0_h3_011; occ111[2] = p1_k0_h2_00; occ111[3] = p1_k0_h5_00101; occ111[4] = p1_k1_h1_1; occ111[5] = p1_k1_h3_011;
	occ111[6] = p1_k1_h5_00100; occ111[7] = p1_k1_h5_00101; occ111[8] = p1_k0_h1_1; occ111[9] = p1_k0_h5_00100; occ111[10] = p1_k1_h3_001;
	occ112[0] = p2_k0_h2_01; occ112[1] = p2_k0_h3_011; occ112[2] = p2_k0_h2_00; occ112[3] = p2_k0_h5_00101; occ112[4] = p2_k1_h1_1; occ112[5] = p2_k1_h3_011;
	occ112[6] = p2_k1_h5_00100; occ112[7] = p2_k1_h5_00101; occ112[8] = p2_k0_h1_1; occ112[9] = p2_k0_h5_00100; occ112[10] = p2_k1_h3_001; //$KSBYMHVRDW

	for(int i=0; i<11; i++){
		occ111[i]++; occ112[i]++;
	}
}

void ABS_FM::Occ4(fm_int pos, fm_int* occ4)
{   //第0层
	BitMap *r = lroot;
	fm_int k0_h0, k1_h0;
	k1_h0 = r->Rank(pos) - 1;
	k0_h0 = pos - k1_h0 - 1;
	
	//第1层
	BitMap *rRoot = r->Right();
	fm_int k0_h1_1, k1_h1_1;
	k1_h1_1 = rRoot->Rank(k1_h0) - 1;           //11 T
	k0_h1_1 = k1_h0 - k1_h1_1 - 1;              //10 G

	BitMap *lRoot = r->Left();
	fm_int k0_h1_0, k1_h1_0;
	k1_h1_0 = lRoot->Rank(k0_h0) - 1;           //01 C
	k0_h1_0 = k0_h0 - k1_h1_0 - 1;              //00 A

	occ4[0] = k0_h1_0; occ4[1] = k1_h1_0; occ4[2] = k0_h1_1; occ4[3] = k1_h1_1; //ACGT
}

void ABS_FM::Occ11(fm_int pos, fm_int* occ11)
{	//第0层
	BitMap *r = rroot;
	fm_int k0_h0, k1_h0;
	k1_h0 = r->Rank(pos) - 1;
	k0_h0 = pos - k1_h0 - 1;

	//第1层
	BitMap *rRoot = r->Right();
	fm_int k0_h1_1, k1_h1_1;
	k1_h1_1 = rRoot->Rank(k1_h0) - 1;               //11 Y
	k0_h1_1 = k1_h0 - k1_h1_1 - 1;                  //10 R

	BitMap *lRoot = r->Left();
	fm_int k0_h1_0, k1_h1_0;
	k1_h1_0 = lRoot->Rank(k0_h0) - 1;
	k0_h1_0 = k0_h0 - k1_h1_0 - 1;

	//第2层
	BitMap *rRoot0 = lRoot->Right();
	fm_int k0_h2_01, k1_h2_01;
	k1_h2_01 = rRoot0->Rank(k1_h1_0) - 1;
	k0_h2_01 = k1_h1_0 - k1_h2_01 - 1;             //010 $

	BitMap *lRoot0 = lRoot->Left();
	fm_int k0_h2_00, k1_h2_00;
	k1_h2_00 = lRoot0->Rank(k0_h1_0) - 1;
	k0_h2_00 = k0_h1_0 - k1_h2_00 - 1;             //000 S

	//第3层
	BitMap *rRoot01 = rRoot0->Right();
	fm_int k0_h3_011, k1_h3_011;
	k1_h3_011 = rRoot01->Rank(k1_h2_01) - 1;       //0111 M
	k0_h3_011 = k1_h2_01 - k1_h3_011 - 1;          //0110 K

	BitMap *rRoot00 = lRoot0->Right();
	fm_int k0_h3_001, k1_h3_001;
	k1_h3_001 = rRoot00->Rank(k1_h2_00) - 1;       //0011 W
	k0_h3_001 = k1_h2_00 - k1_h3_001 - 1;

	//第4层
	BitMap *lRoot001 = rRoot00->Left();
	fm_int k0_h4_0010, k1_h4_0010;
	k1_h4_0010 = lRoot001->Rank(k0_h3_001) - 1;
	k0_h4_0010 = k0_h3_001 - k1_h4_0010 - 1;

	//第5层
	BitMap *rRoot0010 = lRoot001->Right();
	fm_int k0_h5_00101, k1_h5_00101;
	k1_h5_00101 = rRoot0010->Rank(k1_h4_0010) - 1; //001011 V
	k0_h5_00101 = k1_h4_0010 - k1_h5_00101 - 1;    //001010 B

	BitMap *lRoot0010 = lRoot001->Left();
	fm_int k0_h5_00100, k1_h5_00100;
	k1_h5_00100 = lRoot0010->Rank(k0_h4_0010) - 1; //001001 H
	k0_h5_00100 = k0_h4_0010 - k1_h5_00100 - 1;    //001000 D


	occ11[0] = k0_h2_01; occ11[1] = k0_h3_011; occ11[2] = k0_h2_00; occ11[3] = k0_h5_00101; occ11[4] = k1_h1_1; occ11[5] = k1_h3_011;
	occ11[6] = k1_h5_00100; occ11[7] = k1_h5_00101; occ11[8] = k0_h1_1; occ11[9] = k0_h5_00100; occ11[10] = k1_h3_001; //$KSBYMHVRDW
}

fm_int ABS_FM::myLF(fm_int i)  //跳步，更新排位，使用新的小波树，2023-4-11
{	
	if(i==sa0_index){
		return 0;
	}
	return Occ1(i);
}

fm_int ABS_FM::LF(fm_int i)                                                   //跳步，更新排位
{	
	if(i==sa0_index){
		return 0;
	}
	fm_int occ =0;
	unsigned char label =0;
	Occ(occ,label,i);                                                         //原始
	// Occ(occ,label,i-1);                                                         //调用Occ（传出i处的label）
	int coding = code[label];                                                 //由label得到小波树上的编码
	/* printf("............c:%c  C[]:%lld  rank:%lld  label:%d  sa0_index:%lld \n",iupacChar[label],C[coding],occ,(int)label,sa0_index); */
	// return occ + C[coding] - 1;                                               //原始
	return occ + C[coding];                                                //更新i
}


unsigned char ABS_FM::L(fm_int i)
{
	BitMap * r = root;
	int bit =0;
	fm_int rank = 0;
	
	while(r->Left())
	{
		rank = r->Rank(i,bit);
		if(bit ==1)
		{
			i = rank -1;
			r=r->Right();
		}
		else
		{
			i = (i+1) - rank -1;
			r=r->Left();
		}
	}
	return r->Label();

}

fm_int ABS_FM::Occ(fm_int & occ , unsigned char & label,fm_int pos)          //重载3（传出occ 和 label）
{
	BitMap * r = root;
	int bit =0;
	fm_int rank =0;
	while(r->Left())                                                         //一直走到叶子节点
	{
		rank = r->Rank(pos,bit);
		
		if(bit==1)
		{
			pos = rank -1;
			r =r ->Right();
		}
		else
		{
			pos = (pos+1) - rank -1;
			r = r->Left();
		}
	}
	occ = pos +1;
	label = r->Label();                                                     //取叶子节点上的Label
	return 0;
}


unsigned char * ABS_FM::Getfile(const char *filename)
{
	// this->filename = filename;
	this->filename = (char*) malloc(strlen(filename)); 
	strcpy(this->filename,filename);

	FILE * fp = fopen(filename,"r+");
	if(fp==NULL)
	{
		cout<<"Be sure the file is available "<<filename<<endl;
		exit(0);
	}
	fseek(fp,0,SEEK_END);
	// this->n = ftell(fp)+1;
	this->n = ftell(fp);                          //去结尾的$符
	cout<<"ftell...n:"<<ftell(fp)<<endl;
	// unsigned char * T = new unsigned char[n];     //读文本
	unsigned char * T = new unsigned char[n];     //读文本
	fseek(fp,0,SEEK_SET);

	fm_int e=0;
	fm_int num=0;
	// while((e=fread(T+num,sizeof(uchar),n-1-num,fp))!=0){     //逐字符读
	while((e=fread(T+num,sizeof(uchar),n-num,fp))!=0){          //逐字符读
		// if(T[num]==0){
		// 	cout<<"text...include...$"<<endl;                   //原文中并没有$
		// }
		num = num +e;
	}
	// if(num!=n-1)
	if(num!=n)
	{
		cout<<num<<"  "<<n<<"  Read source file failed"<<endl;
		exit(0);
	}
	// T[n-1]=0;                                               //去结尾的$符
	fclose(fp);

	memset(charFreq,0,256*sizeof(long long));                  //bug已改 gll 20220831
	memset(charMap,0,256*sizeof(bool));
	for(fm_int i=0;i<n;i++){                                   //统计字符频率
		charFreq[T[i]]++;                                      //按asc码值，映射到256中（后面Huffman的nodeinit函数中用）
		// cout<<"i:"<<i<<"  "<<(int)T[i]<<"  "<<T[i]<<endl;
		// if(T[i]==0){
		// 	cout<<"pos...$:"<<i<<endl;         //看$符的位置
		// }	
	}
	this->alphabetsize = 0;
	for(int i=0;i<256;i++)
		if(charFreq[i]!=0)
		{
			this->alphabetsize++;                              //频率表非空，则字符表加1
			this->charMap[i]=true;                             //bool标记该ASC码，有映射
		}
	this->code = new int[256];//kkzone-bai debug
	this->C = new fm_int[alphabetsize+1];                      //C表（字符个数+1）
	memset(C,0,(alphabetsize+1)*4);
	this->C[alphabetsize] = n;                                 //0至alphabetsize，最后一个放文本总长度n
	this->C[0] = 0;                                            //第一个放0
	int k=1;
	fm_int pre =0;
	for(int i=0;i<256;i++)                                     //i对应asc码，k对应C表下标
	{
		if(charFreq[i]!=0)
		{
			code[i]=k-1;                                       //code[asc码]=C表下标（字典序）
			C[k]=pre + charFreq[i];
			pre = C[k];
			k++;
		}
		else
			code[i]=-1;                                        //字符为空的，code[asc码]赋值-1
	}

	cout<<"alphaSize: "<< alphabetsize <<endl;
	for(int i=0; i<alphabetsize+1; i++){
		cout<<"i:"<< i <<"   C[]:"<< C[i] <<endl;                //打印C表
	}
	// getchar();
	return T;
}


int ABS_FM::BWT(unsigned char *T,saidx64_t * SA,unsigned char * bwt,fm_int len)
{
	// fm_int i=0;
	// fm_int index=0;
	// for(i=0;i<len;i++)
	// {
	// 	index = (SA[i]-1+len)%len;
	// 	bwt[i]=T[index];                                    //bwt=T[SA-1]
	// }
	for (bwtint_t i = 0; i <= len; i++) {
		if (SA[i] == 0) {                                                         //SA为0，即原文本序列行，也即bwt为$的行
			sa0_index = i;                                                        //记录行号
			bwt[i] = nt16_table['$'];                                              //将$添加到bwt串中
		}
		else {
			bwt[i] = T[SA[i] - 1];                                                 //生成bwt串，复用SA结构
		}
	}
	return 0;
}


/* int ABS_FM::BuildTree()                                  //下标采样（20220913，还有bug）
{
	saidx64_t *SA = new saidx64_t[n];
	divsufsort64(T,SA,n);		

	int step =this->D;
	SAL=new InArray(n/step+1,blog(n/step)+1);
	RankL=new InArray(n/step+2,blog(n/step)+1);
	newD=new Dmap(n);
	fm_int i=0;
	fm_int j=0;
	for(i=0;i<n;i++){
		if(SA[i]%step==0){
			SAL->SetValue(j++,SA[i]/step);                //SA采样
			newD->SetValue(i,1);                          //记录D
		}else if(i==0){
			newD->SetValue(i,1);
		}else{
			newD->SetValue(i,0);
		}
	}
	newD->constructrank(256,16);
	if((n-1)%step!=0){
		RankL->SetValue((n-1)/step+1,newD->rank(0));
	}
	for(i=0;i<n;i++)
	{
		if(SA[i]%step==0){
			RankL->SetValue(SA[i]/step,newD->rank(i));   //逆SA采样
		}
	}
	
	bwt = new unsigned char[n];                          //用于存放bwt串
	BWT(T,SA,bwt,n);                                     //生成bwt
	
	cout<<"write bwt..."<<n<<endl;
	FILE* bwtFile = (FILE*) fopen("/home/lab/gll/panaln/test_data/fasta0RevComp-bwt", "wt");             //写bwt串
	for(int i=0; i<n; i++){putc(bwt[i],bwtFile);}                                               
	fclose(bwtFile);
	cout<<"finish..."<<endl;
	getchar();


	runs=0.0;
	for(fm_int i=0;i<n-1;i++)
		if(bwt[i]!=bwt[i+1])
			runs++;
	runs=n/runs;
	
	TreeCode();
	root=CreateWaveletTree(bwt,n);
	delete [] T;
	T=NULL;
	delete [] SA;
	SA=NULL;
	delete[] bwt;
	bwt=NULL;

	return 0;
} */

void ABS_FM::mycodetable(){ //codeTable给的是gray码
	memset(codeTable, 0, CHAR_SET_SIZE * CHAR_SET_SIZE);
	// codeTable['A'][0] = '0'; codeTable['A'][1] = '0';
	// codeTable['C'][0] = '0'; codeTable['C'][1] = '1';
	// codeTable['G'][0] = '1'; codeTable['G'][1] = '0';
	// codeTable['T'][0] = '1'; codeTable['T'][1] = '1';

	// codeTable['R'][0] = '1'; codeTable['R'][1] = '0';
	// codeTable['Y'][0] = '1'; codeTable['Y'][1] = '1';
	// codeTable['$'][0] = '0'; codeTable['$'][1] = '1'; codeTable['$'][2] = '0';
	// codeTable['S'][0] = '0'; codeTable['S'][1] = '0'; codeTable['S'][2] = '0';
	// codeTable['K'][0] = '0'; codeTable['K'][1] = '1'; codeTable['K'][2] = '1'; codeTable['K'][3] = '0';
	// codeTable['M'][0] = '0'; codeTable['M'][1] = '1'; codeTable['M'][2] = '1'; codeTable['M'][3] = '1';
	// codeTable['W'][0] = '0'; codeTable['W'][1] = '0'; codeTable['W'][2] = '1'; codeTable['W'][3] = '1';
	// codeTable['B'][0] = '0'; codeTable['B'][1] = '0'; codeTable['B'][2] = '1'; codeTable['B'][3] = '0'; codeTable['B'][4] = '1'; codeTable['B'][5] = '0';
	// codeTable['V'][0] = '0'; codeTable['V'][1] = '0'; codeTable['V'][2] = '1'; codeTable['V'][3] = '0'; codeTable['V'][4] = '1'; codeTable['V'][5] = '1';
	// codeTable['D'][0] = '0'; codeTable['D'][1] = '0'; codeTable['D'][2] = '1'; codeTable['D'][3] = '0'; codeTable['D'][4] = '0'; codeTable['D'][5] = '0';
	// codeTable['H'][0] = '0'; codeTable['H'][1] = '0'; codeTable['H'][2] = '1'; codeTable['H'][3] = '0'; codeTable['H'][4] = '0'; codeTable['H'][5] = '1';

	codeTable[nt16_table['A']][0] = '0'; codeTable[nt16_table['A']][1] = '0'; //A:0 C:1 G:2 T:3
	codeTable[nt16_table['C']][0] = '0'; codeTable[nt16_table['C']][1] = '1';
	codeTable[nt16_table['G']][0] = '1'; codeTable[nt16_table['G']][1] = '0'; //G与R冲突，不影响（分别在两个树中）
	codeTable[nt16_table['T']][0] = '1'; codeTable[nt16_table['T']][1] = '1'; //T与Y冲突，不影响

	codeTable[nt16_table['R']][0] = '1'; codeTable[nt16_table['R']][1] = '0';
	codeTable[nt16_table['Y']][0] = '1'; codeTable[nt16_table['Y']][1] = '1';
	codeTable[nt16_table['$']][0] = '0'; codeTable[nt16_table['$']][1] = '1'; codeTable[nt16_table['$']][2] = '0';
	codeTable[nt16_table['S']][0] = '0'; codeTable[nt16_table['S']][1] = '0'; codeTable[nt16_table['S']][2] = '0'; //由于'N'的码，默认全0，故与S相同
	codeTable[nt16_table['K']][0] = '0'; codeTable[nt16_table['K']][1] = '1'; codeTable[nt16_table['K']][2] = '1'; codeTable[nt16_table['K']][3] = '0';
	codeTable[nt16_table['M']][0] = '0'; codeTable[nt16_table['M']][1] = '1'; codeTable[nt16_table['M']][2] = '1'; codeTable[nt16_table['M']][3] = '1';
	codeTable[nt16_table['W']][0] = '0'; codeTable[nt16_table['W']][1] = '0'; codeTable[nt16_table['W']][2] = '1'; codeTable[nt16_table['W']][3] = '1';
	codeTable[nt16_table['B']][0] = '0'; codeTable[nt16_table['B']][1] = '0'; codeTable[nt16_table['B']][2] = '1'; codeTable[nt16_table['B']][3] = '0'; codeTable[nt16_table['B']][4] = '1'; codeTable[nt16_table['B']][5] = '0';
	codeTable[nt16_table['V']][0] = '0'; codeTable[nt16_table['V']][1] = '0'; codeTable[nt16_table['V']][2] = '1'; codeTable[nt16_table['V']][3] = '0'; codeTable[nt16_table['V']][4] = '1'; codeTable[nt16_table['V']][5] = '1';
	codeTable[nt16_table['D']][0] = '0'; codeTable[nt16_table['D']][1] = '0'; codeTable[nt16_table['D']][2] = '1'; codeTable[nt16_table['D']][3] = '0'; codeTable[nt16_table['D']][4] = '0'; codeTable[nt16_table['D']][5] = '0';
	codeTable[nt16_table['H']][0] = '0'; codeTable[nt16_table['H']][1] = '0'; codeTable[nt16_table['H']][2] = '1'; codeTable[nt16_table['H']][3] = '0'; codeTable[nt16_table['H']][4] = '0'; codeTable[nt16_table['H']][5] = '1';

	memset(C, 0, 16 * 4);
	// C[0]=0;          C[1]=3077225;    C[2]=2048660993;  C[3]=2051143964;  C[4]=3474995542;  C[5]=3477492930;  C[6]=3477621711;  C[7]=3487542170; //C表 gray序（bwt串含结尾处的$）
	// C[8]=4911393748; C[9]=4913876719; C[10]=4913977057; C[11]=4914105838; C[12]=4924026297; C[13]=4924126635; C[14]=4926282491; C[15]=6971866259;
	C[0]=0;          C[1]=3077224;    C[2]=2048660992;  C[3]=2051143963;  C[4]=3474995541;  C[5]=3477492929;  C[6]=3477621710;  C[7]=3487542169; //C表 gray序（不含结尾处的$）
	C[8]=4911393747; C[9]=4913876718; C[10]=4913977056; C[11]=4914105837; C[12]=4924026296; C[13]=4924126634; C[14]=4926282490; C[15]=6971866258;
}

void ABS_FM::huffmanTreecode(){
	memset(codeTable, 0, CHAR_SET_SIZE * CHAR_SET_SIZE);

	// codeTable[nt16_table['A']][0] = '1'; codeTable[nt16_table['A']][1] = '0';
	// codeTable[nt16_table['C']][0] = '0'; codeTable[nt16_table['C']][1] = '1'; codeTable[nt16_table['C']][2] = '1';
	// codeTable[nt16_table['G']][0] = '0'; codeTable[nt16_table['G']][1] = '0';
	// codeTable[nt16_table['T']][0] = '1'; codeTable[nt16_table['T']][1] = '1';
	// codeTable[nt16_table['R']][0] = '0'; codeTable[nt16_table['R']][1] = '1'; codeTable[nt16_table['R']][2] = '0'; codeTable[nt16_table['R']][3] = '1'; codeTable[nt16_table['R']][4] = '0';
	// codeTable[nt16_table['Y']][0] = '0'; codeTable[nt16_table['Y']][1] = '1'; codeTable[nt16_table['Y']][2] = '0'; codeTable[nt16_table['Y']][3] = '1'; codeTable[nt16_table['Y']][4] = '1';
	// codeTable[nt16_table['$']][0] = '0'; codeTable[nt16_table['$']][1] = '1'; codeTable[nt16_table['$']][2] = '0'; codeTable[nt16_table['$']][3] = '0'; codeTable[nt16_table['$']][4] = '1'; codeTable[nt16_table['$']][5] = '0';
	// codeTable[nt16_table['S']][0] = '0'; codeTable[nt16_table['S']][1] = '1'; codeTable[nt16_table['S']][2] = '0'; codeTable[nt16_table['S']][3] = '0'; codeTable[nt16_table['S']][4] = '0'; codeTable[nt16_table['S']][5] = '0';
	// codeTable[nt16_table['K']][0] = '0'; codeTable[nt16_table['K']][1] = '1'; codeTable[nt16_table['K']][2] = '0'; codeTable[nt16_table['K']][3] = '0'; codeTable[nt16_table['K']][4] = '1'; codeTable[nt16_table['K']][5] = '1'; codeTable[nt16_table['K']][6] = '0';
	// codeTable[nt16_table['M']][0] = '0'; codeTable[nt16_table['M']][1] = '1'; codeTable[nt16_table['M']][2] = '0'; codeTable[nt16_table['M']][3] = '0'; codeTable[nt16_table['M']][4] = '1'; codeTable[nt16_table['M']][5] = '1'; codeTable[nt16_table['M']][6] = '1';
	// codeTable[nt16_table['W']][0] = '0'; codeTable[nt16_table['W']][1] = '1'; codeTable[nt16_table['W']][2] = '0'; codeTable[nt16_table['W']][3] = '0'; codeTable[nt16_table['W']][2] = '0'; codeTable[nt16_table['W']][5] = '1'; codeTable[nt16_table['W']][6] = '1';
	// codeTable[nt16_table['B']][0] = '0'; codeTable[nt16_table['B']][1] = '1'; codeTable[nt16_table['B']][2] = '0'; codeTable[nt16_table['B']][3] = '0'; codeTable[nt16_table['B']][4] = '0'; codeTable[nt16_table['B']][5] = '1'; codeTable[nt16_table['B']][6] = '0'; codeTable[nt16_table['B']][7] = '1'; codeTable[nt16_table['B']][8] = '0';
	// codeTable[nt16_table['V']][0] = '0'; codeTable[nt16_table['V']][1] = '1'; codeTable[nt16_table['V']][2] = '0'; codeTable[nt16_table['V']][3] = '0'; codeTable[nt16_table['V']][4] = '0'; codeTable[nt16_table['V']][5] = '1'; codeTable[nt16_table['V']][6] = '0'; codeTable[nt16_table['V']][7] = '1'; codeTable[nt16_table['V']][8] = '1';
	// codeTable[nt16_table['D']][0] = '0'; codeTable[nt16_table['D']][1] = '1'; codeTable[nt16_table['D']][2] = '0'; codeTable[nt16_table['D']][3] = '0'; codeTable[nt16_table['D']][4] = '0'; codeTable[nt16_table['D']][5] = '1'; codeTable[nt16_table['D']][6] = '0'; codeTable[nt16_table['D']][7] = '0'; codeTable[nt16_table['D']][8] = '0';
	// codeTable[nt16_table['H']][0] = '0'; codeTable[nt16_table['H']][1] = '1'; codeTable[nt16_table['H']][2] = '0'; codeTable[nt16_table['H']][3] = '0'; codeTable[nt16_table['H']][4] = '0'; codeTable[nt16_table['H']][5] = '1'; codeTable[nt16_table['H']][6] = '0'; codeTable[nt16_table['H']][7] = '0'; codeTable[nt16_table['H']][8] = '1';

	//TA、GC、YR、HD 交换
	codeTable[nt16_table['T']][0] = '1'; codeTable[nt16_table['T']][1] = '0';
	codeTable[nt16_table['G']][0] = '0'; codeTable[nt16_table['G']][1] = '1'; codeTable[nt16_table['G']][2] = '1';
	codeTable[nt16_table['C']][0] = '0'; codeTable[nt16_table['C']][1] = '0';
	codeTable[nt16_table['A']][0] = '1'; codeTable[nt16_table['A']][1] = '1';
	codeTable[nt16_table['Y']][0] = '0'; codeTable[nt16_table['Y']][1] = '1'; codeTable[nt16_table['Y']][2] = '0'; codeTable[nt16_table['Y']][3] = '1'; codeTable[nt16_table['Y']][4] = '0';
	codeTable[nt16_table['R']][0] = '0'; codeTable[nt16_table['R']][1] = '1'; codeTable[nt16_table['R']][2] = '0'; codeTable[nt16_table['R']][3] = '1'; codeTable[nt16_table['R']][4] = '1';
	codeTable[nt16_table['$']][0] = '0'; codeTable[nt16_table['$']][1] = '1'; codeTable[nt16_table['$']][2] = '0'; codeTable[nt16_table['$']][3] = '0'; codeTable[nt16_table['$']][4] = '1'; codeTable[nt16_table['$']][5] = '0';
	codeTable[nt16_table['S']][0] = '0'; codeTable[nt16_table['S']][1] = '1'; codeTable[nt16_table['S']][2] = '0'; codeTable[nt16_table['S']][3] = '0'; codeTable[nt16_table['S']][4] = '0'; codeTable[nt16_table['S']][5] = '0';
	codeTable[nt16_table['K']][0] = '0'; codeTable[nt16_table['K']][1] = '1'; codeTable[nt16_table['K']][2] = '0'; codeTable[nt16_table['K']][3] = '0'; codeTable[nt16_table['K']][4] = '1'; codeTable[nt16_table['K']][5] = '1'; codeTable[nt16_table['K']][6] = '0';
	codeTable[nt16_table['M']][0] = '0'; codeTable[nt16_table['M']][1] = '1'; codeTable[nt16_table['M']][2] = '0'; codeTable[nt16_table['M']][3] = '0'; codeTable[nt16_table['M']][4] = '1'; codeTable[nt16_table['M']][5] = '1'; codeTable[nt16_table['M']][6] = '1';
	codeTable[nt16_table['W']][0] = '0'; codeTable[nt16_table['W']][1] = '1'; codeTable[nt16_table['W']][2] = '0'; codeTable[nt16_table['W']][3] = '0'; codeTable[nt16_table['W']][2] = '0'; codeTable[nt16_table['W']][5] = '1'; codeTable[nt16_table['W']][6] = '1';
	codeTable[nt16_table['B']][0] = '0'; codeTable[nt16_table['B']][1] = '1'; codeTable[nt16_table['B']][2] = '0'; codeTable[nt16_table['B']][3] = '0'; codeTable[nt16_table['B']][4] = '0'; codeTable[nt16_table['B']][5] = '1'; codeTable[nt16_table['B']][6] = '0'; codeTable[nt16_table['B']][7] = '1'; codeTable[nt16_table['B']][8] = '0';
	codeTable[nt16_table['V']][0] = '0'; codeTable[nt16_table['V']][1] = '1'; codeTable[nt16_table['V']][2] = '0'; codeTable[nt16_table['V']][3] = '0'; codeTable[nt16_table['V']][4] = '0'; codeTable[nt16_table['V']][5] = '1'; codeTable[nt16_table['V']][6] = '0'; codeTable[nt16_table['V']][7] = '1'; codeTable[nt16_table['V']][8] = '1';
	codeTable[nt16_table['H']][0] = '0'; codeTable[nt16_table['H']][1] = '1'; codeTable[nt16_table['H']][2] = '0'; codeTable[nt16_table['H']][3] = '0'; codeTable[nt16_table['H']][4] = '0'; codeTable[nt16_table['H']][5] = '1'; codeTable[nt16_table['H']][6] = '0'; codeTable[nt16_table['H']][7] = '0'; codeTable[nt16_table['H']][8] = '0';
	codeTable[nt16_table['D']][0] = '0'; codeTable[nt16_table['D']][1] = '1'; codeTable[nt16_table['D']][2] = '0'; codeTable[nt16_table['D']][3] = '0'; codeTable[nt16_table['D']][4] = '0'; codeTable[nt16_table['D']][5] = '1'; codeTable[nt16_table['D']][6] = '0'; codeTable[nt16_table['D']][7] = '0'; codeTable[nt16_table['D']][8] = '1';
}

int ABS_FM::BuildTree()                                    //值采样，测试索引构建时间，2023-4-10
{	
	saidx64_t *SA = new saidx64_t[n+1];
	SA[0] = n;
	cout<<"build tree n:"<<n<<"  D:"<<this->D<<endl;       //这里的D为32，没问题
	// getchar();

	divsufsort64(T,SA+1,n);		
	int step =this->D;
	SAL=new InArray((n+1)/step+2,blog(n+1)+1);                //多了$的一个
	fm_int i=0, j=0;
	// SAL->SetValue(j++,n);                                  //由于没存$，所以这里补上
	for(i=0; i<=n; i++){
		if(i%step ==0){
			SAL->SetValue(j++,SA[i]);                
		}
	}
	cout<<"construct SAL ..."<<endl;
	string fsal = "sal.idx"; //写索引1
	savekit s4(fsal.c_str());
	SAL->write(s4); //保存SA采样
	cout<<"SAL space: "<<SAL->GetMemorySize()/1024/1024<<" MB"<<endl;
	SAL->~InArray(); //用完析构
	SAL=NULL;
	s4.close();

	bwt = new unsigned char[n+1];                                                                                  //用于存放bwt串
	BWT(T,SA,bwt,n);                                                                                               //生成bwt（传入n，但bwt的长度为n+1）
	delete [] T; //释放T
	T=NULL;
	delete [] SA; //释放SA
	SA=NULL;

	///////////////////////////////////////////
	newD = new Dmap(n+1);	//为了与直接读.bwtSeq结果一致，2023-4-11
	unsigned char* lptr = new unsigned char[n+1];                        //左序列
	unsigned char* rptr = new unsigned char[n+1];                        //右序列
	memset(lptr, 0, n+1);  memset(rptr, 0, n+1);
	long rightLen = 0, leftLen = 0;
	unsigned char c/* , cZ */;
	for(long i=0; i<n+1; i++){
		/* cZ = bwt[i];
		c = nt16_table[cZ]; */ //IUPAC，ref中已经是gray
		c =bwt[i];
		if(is_snp_$[c]){
			newD->SetValue(i,1);
			rptr[rightLen++] = c;                                      //右子树SNP
		}else{
			lptr[leftLen++] = c;                                       //左子树ACGT
		}
	}
	// cout<<"read sequence finish"<<endl;
	delete[] bwt; //释放bwt
	bwt=NULL;

	cout<<"construct newD ..."<<endl;
	newD->constructrank(1024,16);
	cout<<"Dmap space: "<<newD->GetMemorySize()/1024/1024<<" MB"<<endl;
	string fnewD = "newD.idx"; //写索引2
	savekit s(fnewD.c_str());
	newD->write(s);
	newD->~Dmap(); //析构Dmap
	s.close();

	mycodetable(); //重编码
	cout<<"construct Ltree ..."<<endl;
	lroot = CreateWaveletTree(lptr, leftLen);
	cout<<"Ltree space: "<<TreeSizeInByte(lroot)/1024/1024<<" MB"<<endl;
	string flroot = "lroot.idx"; //写索引3
	savekit s2(flroot.c_str());
	SaveWTTree(s2, lroot);
	delete[] lptr; //释放lptr
	lptr=NULL;
	s2.close();

	mycodetable(); //重编码
	cout<<"construct Rtree ..."<<endl;
	rroot = CreateWaveletTree(rptr, rightLen);
	cout<<"Rtree space: "<<TreeSizeInByte(rroot)/1024/1024<<" MB"<<endl<<endl;
	string frroot = "rroot.idx"; //写索引4
	savekit s3(frroot.c_str());
	SaveWTTree(s3, rroot);
	delete[] rptr; //释放rptr
	rptr=NULL;
	s.close();

	cout<<"finish index save, pause"<<endl;
	// getchar();

	return 0;
}

// 这个是目前用的，2023-4-10
// int ABS_FM::BuildTree0(){ //WT_Handle初始化时调用，这个是正确的，直接读bwt串，2023-4-10
// 	FILE* bwtFile = (FILE*) fopen("/home/lab/gll/panaln/mg-ref/multigenome_bubble_GRCh38.fasta.ref.bwtSeq","r+"); //ASC码明文
// 	if(bwtFile == NULL) cout<<"file dost not exit"<<endl;
// 	fseek(bwtFile, 0, SEEK_END);
// 	n = ftell(bwtFile);
// 	cout<<"n:"<<n<<endl;
// 	fseek(bwtFile, 0, SEEK_SET);
// 	newD = new Dmap(n);	
// 	memset(charFreq, 0, 256*sizeof(long long));
// 	unsigned char* seq = new unsigned char[n];                         //整序列
// 	unsigned char* lptr = new unsigned char[n];                        //左序列
// 	unsigned char* rptr = new unsigned char[n];                        //右序列
// 	memset(lptr, 0, n);  memset(rptr, 0, n); memset(seq, 0, n);
// 	long rightLen = 0, leftLen = 0;
	
// 	unsigned char c, cZ;
// 	for(long i=0; i<n; i++){
// 		cZ = fgetc(bwtFile);
// 		c = nt16_table[cZ]; //IUPAC
// 		seq[i] = c; charFreq[c]++;		
// 		if(is_snp_$[c]){
// 			newD->SetValue(i,1);
// 			rptr[rightLen++] = c;                                      //右子树SNP
// 		}else{
// 			lptr[leftLen++] = c;                                       //左子树ACGT
// 		}
// 	}
// 	cout<<"read sequence finish"<<endl;
// 	alphabetsize = 0;
// 	for(int i=0; i<256; i++){                                          //统计字符表大小
// 		if(charFreq[i] != 0){
// 			alphabetsize++;
// 		}
// 	}
// 	cout<<"Alphabet size: "<<alphabetsize<<endl<<endl;
// 	/* C = new long long[alphabetsize + 1];                               //构建C表
// 	memset(C, 0, (alphabetsize + 1) * sizeof(long long));
// 	C[0] = 0;
// 	int k = 1;
// 	for(int i=0; i<256; i++){
// 		if(charFreq[i] != 0){
// 			code[i] = k - 1;
// 			C[k] = pre + charFreq[i];
// 			pre = C[k];
// 			k++;
// 		}
// 		else
// 			code[i] = -1;
// 	} */

// 	// TreeCode(); //为了测试索引构建时间，暂时去掉root的部分，2023-4-10
// 	// cout<<"construct tree ..."<<endl;
// 	// root = CreateWaveletTree(seq, n);
// 	// cout<<"Wtree space: "<<TreeSizeInByte(root)/1024/1024<<" MB"<<endl;
// 	// string froot = "root.idx";
// 	// savekit s1(froot.c_str());
// 	// SaveWTTree(s1, root);
// 	// s1.close();
// 	/* {
// 	cout<<"load tree ..."<<endl;
// 	string froot2 = "root.idx";
// 	loadkit s1(froot2.c_str());
// 	LoadWTTree(s1, root);
// 	s1.close();
// 	} */

// 	cout<<"construct newD ..."<<endl;
// 	newD->constructrank(1024,16);
// 	cout<<"Dmap space: "<<newD->GetMemorySize()/1024/1024<<" MB"<<endl;
// 	string fnewD = "newD.idx";
// 	savekit s(fnewD.c_str());
// 	newD->write(s);
// 	s.close();
// 	/* {
// 	cout<<"load newD ..."<<endl;
// 	newD = new Dmap(6971866259);
// 	string fnewD2 = "newD.idx";
// 	loadkit s2(fnewD2.c_str());
// 	newD->load(s2);
// 	s2.close();
// 	} */

// 	mycodetable(); //重编码
// 	cout<<"construct Ltree ..."<<endl;
// 	lroot = CreateWaveletTree(lptr, leftLen);
// 	cout<<"Ltree space: "<<TreeSizeInByte(lroot)/1024/1024<<" MB"<<endl;
// 	string flroot = "lroot.idx";
// 	savekit s2(flroot.c_str());
// 	SaveWTTree(s2, lroot);
// 	s2.close();
// 	/* {
// 	cout<<"load Ltree ..."<<endl;
// 	string flroot2 = "lroot.idx";
// 	loadkit s3(flroot2.c_str());
// 	LoadWTTree(s3, lroot);
// 	s3.close();
// 	} */

// 	mycodetable(); //重编码
// 	cout<<"construct Rtree ..."<<endl;
// 	rroot = CreateWaveletTree(rptr, rightLen);
// 	cout<<"Rtree space: "<<TreeSizeInByte(rroot)/1024/1024<<" MB"<<endl<<endl;
// 	// cout<<"rroot addr:"<<rroot<<endl;
// 	string frroot = "rroot.idx";
// 	savekit s3(frroot.c_str());
// 	SaveWTTree(s3, rroot);
// 	s.close();
// 	/* {
// 	cout<<"load Rtree ..."<<endl;
// 	string frroot2 = "rroot.idx";
// 	loadkit s4(frroot2.c_str());
// 	LoadWTTree(s4, rroot);
// 	s4.close();
// 	} */

// 	cout<<"finish index save, pause"<<endl;
// 	getchar();

// 	// long EPOCH = 1000 * 1000 * 100;
// 	long EPOCH = 1;
// 	srand48(9);
// 	cout<<"############ 1. test new"<<endl;
// 	mycodetable();
// 	long long occ161[16]={0}, occ162[16]={0}, occ1=0, occ2=0, occ81[8]={0}, occ82[8]={0};
// 	long plus=888;
// 	long pos1=500000+plus, pos2=600000+plus;
// 	// long pos1=0, pos2=0, Max=6971866259;
// 	clock_t t1 = clock();
// 	for(long e=0; e<EPOCH; e++){
// 		// pos1 = (long) (drand48() * Max); pos2 = ((pos1 + 10000) > n) ? n : (pos1 + 10000); 

// 		Occ16p(pos1, pos2, occ161, occ162);

// 		// Occ8p(nt16_table['T'], pos1, pos2, occ81, occ82);     //2. 测试Occ8p正确性

// 		// if(e==10) continue;
// 		// Occ1p(e, pos1, pos2, occ1, occ2, 1);                  //1. 测试Occ1p正确性
// 		// cout<<setw(5)<<left<<e<<"occ1: "<<setw(12)<<left<<occ1<<"occ2: "<<setw(12)<<left<<occ2<<endl;
// 	}
// 	for(int i=0; i<16; i++){
// 		if(i==10) continue;
// 		cout<<setw(5)<<left<<i<<"occ1: "<<setw(10)<<left<<occ161[i]<<"occ2: "<<setw(10)<<left<<occ162[i]<<endl;
// 	}
// 	// cout<<"1. cost time: " << setprecision(6) << (float)(clock() - t1) / CLOCKS_PER_SEC<<endl;
	

// 	cout<<endl<<endl;
// 	cout<<"############ 2. test old"<<endl;
// 	huffmanTreecode();
// 	// unsigned char a = nt16_table['A'];
// 	// Occ(a, pos1, pos2, rank1, rank2);
// 	// cout<<"rank1: "<<setw(10)<<left<<rank1<<"rank2: "<<setw(10)<<left<<rank2<<endl;
// 	clock_t t2 = clock();
// 	for(long e=0; e<EPOCH; e++){
// 		// pos1 = (long) (drand48() * Max); pos2 = ((pos1 + 10000) > n) ? n : (pos1 + 10000);
// 		for(int i=0; i<16; i++){
// 			// Occ(nucl_bases_table[3][i], pos1, pos2, occ1, occ2); //AGCT顺序
// 			if(i==10) continue;
// 			Occ(i, pos1, pos2, occ1, occ2);
// 			cout<<setw(5)<<left<<i<<"occ1: "<<setw(10)<<left<<occ1<<"occ2: "<<setw(10)<<left<<occ2<<endl;
// 		}

// 		// if(e==10) continue;
// 		// Occ(e, pos1, pos2, rank1, rank2);
// 		// cout<<setw(5)<<left<<e<<"occ1: "<<setw(12)<<left<<rank1<<"occ2: "<<setw(12)<<left<<rank2<<endl;
// 	}
// 	// cout<<"2. cost time: "<< setprecision(6) << (float)(clock() - t2) / CLOCKS_PER_SEC<<endl;
// 	// cout<<endl<<"############ test occ11"<<endl;
// 	// long long occ111[11]={0}, occ112[11]={0};
// 	// long pos3=2648, pos4=3170;
// 	// Occ11p(pos3, pos4, occ111, occ112);
// 	// for(int i=0; i<11; i++){
// 	// 	cout<<setw(10)<<left<<occ111[i];
// 	// }
// 	// cout<<endl;
// 	// for(int i=0; i<11; i++){
// 	// 	cout<<setw(10)<<left<<occ112[i];
// 	// }

// 	getchar();
// 	// fclose(bwtFile);
// 	// delete[] seq;
// 	// delete[] lptr;
// 	// delete[] rptr;
// }



// int ABS_FM::BuildTree0()                                    //值采样
// {	
// 	if(0)
// 	{	//测试正确20221106（直接读bwt文件）
// 		long pos = 500000;
// 		cout<<"#########################"<<endl;
// 		// FILE* bwtFile = (FILE*) fopen("/home/lab/gll/panaln/mg-ref/1k.bwtSeq","r+");
// 		FILE* bwtFile = (FILE*) fopen("/home/lab/gll/panaln/mg-ref/multigenome_bubble_GRCh38.fasta.ref.bwtSeq","r+");
// 		if(bwtFile == NULL) cout<<"file dost not exit"<<endl;
// 		fseek(bwtFile, 0, SEEK_END);
// 		n = ftell(bwtFile);
// 		fseek(bwtFile, 0, SEEK_SET);
		
// 		newD = new Dmap(n);
// 		unsigned char* lptr = new unsigned char[n];
// 	    unsigned char* rptr = new unsigned char[n];
// 		memset(lptr, 0, n);  memset(rptr, 0, n);
// 		long rightLen = 0, leftLen = 0;

// 		unsigned char c;
// 		for(long i=0; i<n; i++){
// 			c = fgetc(bwtFile);
// 			if(is_snp_$[nt16_table[c]]){
// 				newD->SetValue(i,1);
// 				rptr[rightLen++] = c;                                      //右子树SNP
// 			}else{
// 				lptr[leftLen++] = c;                                       //左子树ACGT
// 			}
// 		}
// 		fclose(bwtFile);

// 		cout<<"construct newD ..."<<endl;
// 		newD->constructrank(256,16); cout<<endl;
// 		mycodetable();
// 		cout<<"construct Treel ..."<<endl;
// 		BitMap* lroot = CreateWaveletTree(lptr, leftLen);
// 		cout<<"construct Treer ..."<<endl;
// 		BitMap* rroot = CreateWaveletTree(rptr, rightLen);


// 		long topr = newD->rank2(pos);
// 		cout<<endl<<"test occ ..."<<endl;
// 		long topl = pos - topr;
// 		long occA = Occ('A', topl, lroot); long occC = Occ('C', topl, lroot); long occG = Occ('G', topl, lroot); long occT = Occ('T', topl, lroot);
// 		long occR = Occ('R', topr, rroot); long occY = Occ('Y', topr, rroot); long occ$ = Occ('$', topr, rroot); long occS = Occ('S', topr, rroot); long occK = Occ('K', topr, rroot);
// 		long occM = Occ('M', topr, rroot); long occW = Occ('W', topr, rroot); long occB = Occ('B', topr, rroot); long occV = Occ('V', topr, rroot); long occD = Occ('D', topr, rroot); long occH = Occ('H', topr, rroot);
// 		cout<<"occA:"<<setw(15)<<left<<occA<<"occC:"<<setw(15)<<left<<occC<<"occG:"<<setw(15)<<left<<occG<<"occT:"<<setw(15)<<left<<occT<<endl; 
// 		// cout<<"topl:"<<topl<<endl<<endl;
// 		cout<<"occY:"<<setw(15)<<left<<occY<<"occR:"<<setw(15)<<left<<occR<<"occ$:"<<setw(15)<<left<<occ$<<"occS:"<<setw(15)<<left<<occS<<"occM:"<<setw(15)<<left<<occM<<endl;
// 		cout<<"occK:"<<setw(15)<<left<<occK<<"occW:"<<setw(15)<<left<<occW<<"occV:"<<setw(15)<<left<<occV<<"occB:"<<setw(15)<<left<<occB<<"occH:"<<setw(15)<<left<<occH<<"occD:"<<setw(15)<<left<<occD<<endl;
// 		// cout<<"topr:"<<topr<<endl<<endl;
// 		if((topl == occA+occC+occG+occT) && (topr == occR+occY+occ$+occS+occK+occM+occW+occB+occV+occD+occH)){
// 			cout<<"pass!"<<endl;
// 		}else{
// 			cout<<"acgt: "<<setw(15)<<left<<occA+occC+occG+occT<<"topl:"<<setw(15)<<left<<topl<<endl;
// 			cout<<"iupac:"<<setw(15)<<left<<occR+occY+occ$+occS+occK+occM+occW+occB+occV+occD+occH<<"topr:"<<setw(15)<<left<<topr<<endl;
// 		}

// 		cout<<endl<<"test occ15..."<<endl;
// 		long long occ4[4], occ11[11];
// 		Occ4(topl, lroot, occ4); Occ11(topr, rroot, occ11);
// 		cout<<"occ[A]:"<<setw(15)<<left<<occ4[0]<<"occ[C]:"<<setw(15)<<left<<occ4[1]<<"occ[G]:"<<setw(15)<<left<<occ4[2]<<"occ[T]:"<<setw(15)<<left<<occ4[3]<<endl;
// 		cout<<"occ[Y]:"<<setw(15)<<left<<occ11[0]<<"occ[R]:"<<setw(15)<<left<<occ11[1]<<"occ[$]:"<<setw(15)<<left<<occ11[2]<<"occ[S]:"<<setw(15)<<left<<occ11[3]<<"occ[M]:"<<setw(15)<<left<<occ11[4]<<endl
// 			<<"occ[K]:"<<setw(15)<<left<<occ11[5]<<"occ[W]:"<<setw(15)<<left<<occ11[6]<<"occ[V]:"<<setw(15)<<left<<occ11[7]<<"occ[B]:"<<setw(15)<<left<<occ11[8]<<"occ[H]:"<<setw(15)<<left<<occ11[9]<<"occ[D]:"<<setw(15)<<left<<occ11[10]<<endl;

// 		getchar();
// 		delete newD;
// 		//清理小波树
// 	}

// 	/* {      //这个版本的plan_rank报错
// 		cout<<endl<<"construct D... bitmap"<<endl;
// 		n = 6971866258;  n = n + 1;                                        //后面不再加1
// 		// n = 996;
// 		long u64Len = 0;
// 		if(n%64 == 0)
// 			u64Len = n/64 + 1;
// 		else
// 			u64Len = n/64 + 2;
// 		unsigned long long int *bitBuff = new unsigned long long int[u64Len]; //长64的桶
// 		memset(bitBuff, 0, u64Len*8);
// 		long bytePos = 0;
// 		int bitOffset = 0;
// 		FILE* bwtFile = (FILE*) fopen("/home/lab/gll/panaln/mg-ref/multigenome_bubble_GRCh38.fasta.ref.bwtSeq","r+");
// 		// FILE* bwtFile = (FILE*) fopen("/home/lab/gll/panaln/mg-ref/1k.bwtSeq","r+");
// 		if(bwtFile == NULL) cout<<"file dost not exit"<<endl;
// 		unsigned char c;
// 		for(long i=0; i<n; i++){                           //构建0/1串
// 			c = fgetc(bwtFile);
// 			if(is_snp_$[nt16_table[c]]){
// 				bitBuff[bytePos] |= (0x01ull<<(63-bitOffset));
// 			}
// 			bitOffset++;
// 			if(bitOffset == 64){
// 				bytePos++;
// 				bitOffset = 0;
// 			}
// 		}

// 		uchar* tables[4] = {this->Z, this->R, this->R1, this->R3};
// 		BitMap* node = new BitMap(bitBuff, n, 0, 256, 16, 0U, tables);         //栈内存会崩
// 		delete[] bitBuff;
// 		cout<<"finish..."<<endl;
// 		int bit = 1;
// 		cout<<"Bitmap space: "<<node->SizeInByte()/1024/1024 <<" MB"<<endl;
// 		cout<<"test...rank: "<<node->Rank(pos, bit)<<endl;
// 		delete node;
// 		getchar();
// 	} */

// 	saidx64_t *SA = new saidx64_t[n+1];
// 	SA[0] = n;
// 	divsufsort64(T,SA+1,n);		

// 	/* cout<<"write sa..."<<n+1<<endl;
// 	FILE* saFile = (FILE*) fopen("/home/lab/gll/panaln/test_data/fasta0RevComp-sa-panaln", "wt");             //写SA（只是用于检查）
// 	for(long i=0; i<n+1; i++){printf("%lld ",SA[i]);}                                              
// 	fclose(saFile);
// 	cout<<"finish..."<<endl;
// 	// getchar(); */

// 	int step =this->D;
// 	SAL=new InArray((n+1)/step+2,blog(n+1)+1);                //多了$的一个
// 	// SAL=new InArray(n/step+1,blog(n)+1);                   //下标采样要扩容
// 	// SAL=new InArray(n/step+1,blog(n/step)+1);
// 	/* RankL=new InArray(n/step+2,blog(n/step)+1);
// 	newD=new Dmap(n); */
// 	fm_int i=0, j=0;
// 	// SAL->SetValue(j++,n);                                  //由于没存$，所以这里补上
// 	for(i=0;i<=n;i++){
// 		if(i%step==0){
// 			/* SAL->SetValue(j++,SA[i]/step);                //SA采样 */
// 			// SAL->SetValue(j++,SA[i-1]);                //SA采样
// 			SAL->SetValue(j++,SA[i]);                
// 			// printf("generate...sal \n");
// 			// getchar();
// 			// for(int i=0;i<50;i++){
// 			// 	printf("sa...i:%d  %d \n",i,SA[i-1]);
// 			// }
// 			// printf("sa...i:%d  %d \n",i,SA[i]);
// 			/* newD->SetValue(i,1);                          //记录D */
// 		}/* else if(i==0){
// 			newD->SetValue(i,1);
// 		}else{
// 			newD->SetValue(i,0);
// 		} */
// 	}
// 	/* newD->constructrank(256,16);
// 	if((n-1)%step!=0){
// 		RankL->SetValue((n-1)/step+1,newD->rank(0));
// 	}
// 	for(i=0;i<n;i++)
// 	{
// 		if(SA[i]%step==0){
// 			RankL->SetValue(SA[i]/step,newD->rank(i));   //逆SA采样
// 		}
// 	} */
	
// 	bwt = new unsigned char[n+1];                                                                                  //用于存放bwt串
// 	BWT(T,SA,bwt,n);                                                                                               //生成bwt（传入n，但bwt的长度为n+1）

// 	cout<<"write bwt..." << n+1 << endl;
// 	char *fbwt = (char*) malloc(strlen(this->filename) + 7); strcpy(fbwt,this->filename); char *bwtSuf=".bwtSeq"; strcat(fbwt,bwtSuf);
// 	cout<<"bwt file: "<< fbwt << endl;
// 	FILE* bwtFile = (FILE*) fopen(fbwt, "wt");                                                                     //写bwt串
// 	// for(long i=0; i<n+1; i++){putc(bwt[i],bwtFile);}                                                            //gray码（0~15）     
// 	for(long i=0; i<n+1; i++){
// 		putc(iupacChar[bwt[i]], bwtFile);                                                                          //ASC码
// 		if(is_snp_$[bwt[i]]){
// 			newD->SetValue(i,1);                                                                                   //记录到D中
// 		}
// 	}
// 	fclose(bwtFile);
// 	cout << "finish... bwt"<< endl;
// 	// getchar();

// 	runs=0.0;
// 	for(fm_int i=0;i<n;i++)
// 		if(bwt[i]!=bwt[i+1])
// 			runs++;
// 	runs=n/runs;
	
// 	TreeCode();                                            //处理树形相关，产生codeTable
// 	root=CreateWaveletTree(bwt,n+1);                       //给bwt序列，创建小波树
// 	delete [] T;
// 	T=NULL;
// 	delete [] SA;
// 	SA=NULL;
// 	delete[] bwt;
// 	bwt=NULL;

// 	return 0;
// }

BitMap * ABS_FM::CreateWaveletTree(unsigned char * bwt,fm_int n)
{
	BitMap * root = NULL;                                  //每个Huffman节点，对应一个BitMap结构
	root = FullFillWTNode(bwt,n,0);                        //填充各个节点的0/1串             
	if(!root)
	{
		cout<<"FullfillWTNode failed"<<endl;               //失败
		exit(0);
	}
	return root;                                           //返回根节点
}


BitMap * ABS_FM::FullFillWTNode(unsigned char * buff,fm_int len,int level)            //填充各个节点的0/1串
{
//	cout<<level<<endl;
	// cout<<"in...FullFillWTNode...block_size:"<<block_size<<endl;
	int CurrentLevel = level;                                                         //当前层
	long long CurrentBitLen = len;                                                    //当前位长
	unsigned char CurrentLabel = '\0';                                                //当前字符
	unsigned long long int *CurrentBitBuff = NULL;                                    //当前位串
	
	if ((int)strlen((const char*)codeTable[buff[0]])==level)                                                            //该层串的第一个字符，编码长度=层数?
	{
		CurrentLabel = buff[0];                                                                                         //bwt串上取字符
		CurrentBitBuff = NULL;
		//uchar * tables[5]={this->zerostable,this->R1,this->R2,this->R3,this->R4};
		uchar * tables[4] ={this->Z,this->R, this->R1, this->R3};                                                       //解码加速表
		BitMap * node = new BitMap(CurrentBitBuff,CurrentBitLen,CurrentLevel,block_size,r,CurrentLabel,tables);         //新建节点（含这些结构）
		node->Left(NULL);                                                                                               //左孩子
		node->Right(NULL);                                                                                              //右孩子
		return node;                                                                                                    //特殊情况，直接返回
	}
	
	fm_int u64Len=0;
	if(len%64==0)
		u64Len = len/64+1;
	else
		u64Len = len/64+2;
	CurrentBitBuff = new unsigned long long int[u64Len];                                                                //64bit的桶，数组（一个字符对应一位）
	memset(CurrentBitBuff,0,u64Len*8);                                                                                  //置0
	unsigned char * lptr=NULL;                                                                                          //放分到左边的字符（左右可以不等）
	unsigned char * rptr=NULL;                                                                                          //放分到右边的
	fm_int leftLen=0;                                                                                                   //分到左边的字符个数
	fm_int rightLen=0;                                                                                                  //分到右边的

	lptr = new unsigned char[len];
	rptr = new unsigned char[len];
	memset(lptr,0,len);
	memset(rptr,0,len);

	//computer bitvect;

	fm_int i=0;
	fm_int bytePos=0;
	int bitOffset=0;
	u64 last = 0;
	for(i=0;i<len;i++)                                                                                                  //逐个处理，该层串上字符
	{
		/* if(charFreq[buff[i]]==0){                                                                                       //bwt串第i个字符的频率为0，报错（检查用，GetFile中得）
			cout<<"FullFillWTNode error"<<endl;
			exit(0);
		} */
		if(codeTable[buff[i]][level]=='1')                                                                              //bwt字符的第i层为1
		{
			CurrentBitBuff[bytePos] |= (0x01ull<<(63-bitOffset));                                                       //或，向64的桶中放（把1从左向右放，0为本身不管）
			rptr[rightLen++]=buff[i];                                                                                   //右，该层串上字符（给下一层用）
			last = 0;
		}
		else                                                                                                            //bwt字符的第i层为0
		{
			lptr[leftLen++]=buff[i];                                                                                    //左，该层串上字符（给下一层用）
			last = 1;
		}
		bitOffset++;                                                                                                    //桶中位偏移
		if(bitOffset == 64)
		{
			bytePos++;                                                                                                  //64的桶满
			bitOffset = 0;
		}
	}
	CurrentBitBuff[bytePos] |= (last<<(63-bitOffset));                                                                 //这里last有什么作用？
	uchar * tables[4] = {this->Z,this->R, this->R1, this->R3};
	BitMap * node = new BitMap(CurrentBitBuff,CurrentBitLen,CurrentLevel,block_size,r,CurrentLabel,tables);            //同上
	if(leftLen !=0)
	{
		BitMap * left =FullFillWTNode(lptr,leftLen,level+1);                                                           //递归构建，对上层level位为0的所有组合字符串
		node->Left(left);                                                                                              //挂为左孩子
		delete [] lptr;
		lptr=NULL;
	}
	if(rightLen!=0)
	{
		BitMap * right = FullFillWTNode(rptr,rightLen,level+1);                                                        //递归构建，对上层level位为1的所有组合字符串
		node->Right(right);                                                                                            //挂为右孩子
		delete [] rptr;
		rptr=NULL;
	}
	return node;                                                                                                       //返回根节点（应该是从上到下构建？）
}


int ABS_FM::DestroyWaveletTree()
{
	delete root ;
	root=NULL;
	return 0;
}


int ABS_FM::blog(fm_int x)
{
	int ans=0;
	while(x>0)
	{
		ans++;
		x=(x>>1);
	}
	return ans;
}


void ABS_FM::Inittable()
{
	this -> Z = new uchar[1<<8];
	int tablewidth = 8;
	for(int i=0;i<tablewidth;i++)
		for(int j=(1<<i);j<(2<<i);j++)
			Z[j] = tablewidth-1-i;
	Z[0] = tablewidth;
	
	u64 tablesize = (1<<16);
	R  = new uchar[tablesize<<2];
	R1 = new uchar[tablesize]();
	R3 = new uchar[tablesize*18]();
	//查找表的初始化：在16bits的0,1串上模拟gamma解码的过程，得到
	//这些表
	u64 B[2]={0xffffffffffffffffull,0xffffffffffffffffull};
	int sum =0;//gamma编码的和,含义为原串被编码的bits数目。
	int step=0;//16bits可以完整解码的bits数,该值不会大于16.
	int rank = 0;//16bits表示的几个完整的runs,假设第一个runs是1-runs,这几个runs的rank值。
	int runs = 0 ;//runs 个数.
	
	int x = 0;//工作变量，保存本次的gamma解码值.
	int prestep = 0;//前一次正确解码的bits数(累加),<=16.
	for(u64 i=0;i<tablesize;i++)
	{
		B[0] = (i<<48);
		sum  =0 ;
		step = 0;
		prestep=0;
		rank = 0;
		runs = 0;
		while(1)
		{
			x = GammaDecode(B,step,this);//step会联动.
			if(step > 16)
				break;
			sum = sum + x;
			prestep = step;
			runs ++;
			if(runs%2==1)
				rank = rank + x;
		}
		R[i<<2] = runs;//r4
		R[(i<<2)+1] = prestep;//r2
		R[(i<<2)+2] = sum; //r1;
		R[(i<<2)+3] = rank;//r3

		if(i==0){
			R1[0] = 16;
			R3[0] = 16;
			continue;
		}
		buildEFtable(i);
	}
}

inline int popcnt1(unsigned long long int x) //默认的pop64，软实现 2023-2-17
{
	x = x - ((x & 0xAAAAAAAAAAAAAAAA) >> 1);
	x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
	x = ((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F);
	return (x * 0x0101010101010101) >> 56;
}

// inline int popcnt1(unsigned long long int x) //改进的pop64，硬实现 2023-2-17
// {
// 	return __builtin_popcountll(x); //硬实现
// }

void ABS_FM::buildEFtable(u64 i){
	int decodelen = 0;
	u16 x = i;
	while(1){
		if(x&1) break;
		x = x>>1;
		decodelen++;
	}
	R3[i*18] = 16-decodelen;

	int onesnum = popcnt1(i);
	R1[i] = R3[i*18]-onesnum;
	//cout<<(int)R1[i]<<" "<<(int)R2[i]<<" "<<(int)R3[i*18]<<endl;

	x = i;
	int zeronums = 0;
	int alreadyBits = 0;
	// int index = 0;
	while(alreadyBits < R3[i*18]){
		int num_0 = Zeros(x,this);
		zeronums += num_0;
		alreadyBits += num_0;
		x = x<<num_0;
		int num_1 = Zeros(~x,this);
		R3[i*18+zeronums+1] = num_1;//桶中元素个数
		//[i*35+18+zeronums] = alreadyBits;//桶的起始位置
		alreadyBits += num_1;
		x = x<<num_1;
	}
}

//递归保存节点的编号信息
int ABS_FM::SaveNodePosition(BitMap *r, int position, savekit &s)
{
	if(!r)
		return 1;
	s.writei32(position);                            //写当前编号: cur（先序遍历，编号不连续）
	SaveNodePosition(r->Left(),  2*position, s);     //左孩子编号: 2*cur
	SaveNodePosition(r->Right(), 2*position+1, s);   //右孩子编号: 2*cur + 1
	return 0;
}

//递归保存节点的数据信息
int ABS_FM::SaveNodeData(BitMap *r,savekit &s)
{
	if(!r)
		return 1 ;
	r->Save(s);                                     //写BitMap结构（先序遍历，用编号对应顺序）
	SaveNodeData(r->Left(),s);                      //左孩子
	SaveNodeData(r->Right(),s);                     //右孩子
	return 0;
}

int ABS_FM::SaveWTTree(savekit &s, BitMap *root) //gll add
{
	//保存编号信息
	int nodecount = TreeNodeCount(root);
	s.writei32(nodecount);
	cout<<"nodecount:"<<nodecount<<endl;
	SaveNodePosition(root,1,s);                  //根节点编号为1， 左孩子为2*cur，右孩子为2*cur + 1

	//保存节点数据信息
	SaveNodeData(root,s);                        
	return 0;
}

int ABS_FM::SaveWTTree(savekit &s)
{
	//保存编号信息
	int nodecount = TreeNodeCount(root); //2N-1，N为字符表大小
	s.writei32(nodecount);
	SaveNodePosition(root,1,s);

	//保存节点数据信息
	SaveNodeData(root,s);
	return 0;
}

int ABS_FM::LoadWTTree(loadkit &s, BitMap* &root, uchar **tables) //gll add（原样读进来，跟原始的blockSize一样，不用管 2023-2-11）
{
	//读取数据，map的int域对应该节点的位置
	int nodecount = 0;
	s.loadi32(nodecount);
	// cout<<"load nodecount:"<<nodecount<<endl;
	// int nodecount = 2*alphabetsize -1;
	// cout<<"alphabetsize:"<<alphabetsize<<endl;
	int *p = new int[nodecount];
	s.loadi32array(p, nodecount);                        //读入节点位置数组（cur，2*cur，2*cur + 1）
	map<int, BitMap *> pmap;                             //map容器，键为int，值为BitMap指针
	BitMap * r = NULL;                                   //临时根节点
	for(int i=0;i<nodecount;i++)
	{
		if(tables)
			r = new BitMap(tables);                      //table为4列的加速表
		else
			r = new BitMap();
		r->Load(s);                                      //逐个读入节点结构
		pmap[p[i]] = r;                                  //将节点结构，映射到节点编号上（先序遍历，写p的顺序与bitmap顺序一致）
	}
	//挂链
	map<int ,BitMap *>::iterator iter;                   //map迭代器
	map<int ,BitMap *>::iterator f_iter;                 //临时
	for(iter=pmap.begin(); iter!=pmap.end(); iter++)
	{
		// cout<<"each size:"<<iter->second->SizeInByte()<<endl;
		f_iter = pmap.find(2*iter->first);               //定位数据，返回出现位置的迭代器（若无，则返回end函数返回的迭代器）
		if(f_iter != pmap.end())
			iter->second->Left(f_iter->second);          //2*cur，左孩子挂链
		else
			iter->second->Left(NULL);                    //左孩子为空
		
		f_iter = pmap.find(2*iter->first +1);            //2*cur，右孩子挂链
		if(f_iter!=pmap.end())
			iter->second->Right(f_iter->second);
		else
			iter->second->Right(NULL);
	}
	f_iter = pmap.find(1);                               //cur=1，为根节点
	if(f_iter !=pmap.end()){
		root = f_iter->second;                           //把根节点的结构，赋给ABS_FM的root属性
	}else{
		cerr<<"Load WTTree error"<<endl;
		root = NULL;
		exit(0);                                         //失败，直接退出了
	}
	return 0;
}

int ABS_FM::LoadWTTree(loadkit &s, uchar **tables)
{
	//读取数据，map的int域对应该节点的位置
	int nodecount = 0;
	s.loadi32(nodecount);
	// int nodecount = 2*alphabetsize -1;
	// cout<<"alphabetsize:"<<alphabetsize<<endl;
	int *p = new int[nodecount];
	s.loadi32array(p, nodecount);
	map<int, BitMap *> pmap;                             //map容器，键为int，值为BitMap指针
	BitMap * r=NULL;
	for(int i=0;i<nodecount;i++)
	{
		if(tables)
			r = new BitMap(tables);                      //table为4列的加速表
		else
			r = new BitMap(); //读入各个bitmap，但是并没有指定blocksize 2023-2-11
		r->Load(s);
		pmap[p[i]] = r;
	}
	//挂链
	map<int ,BitMap *>::iterator iter;
	map<int ,BitMap *>::iterator f_iter;
	for(iter = pmap.begin();iter!=pmap.end();iter++)
	{
		f_iter = pmap.find(2*iter->first);
		if(f_iter != pmap.end())
			iter->second->Left(f_iter->second);
		else
			iter->second->Left(NULL);
		
		f_iter = pmap.find(2*iter->first +1);
		if(f_iter!=pmap.end())
			iter->second->Right(f_iter->second);
		else
			iter->second->Right(NULL);
	}
	f_iter = pmap.find(1);
	if(f_iter !=pmap.end())
		this->root = f_iter->second;
	else
	{
		cerr<<"Load WTTree error"<<endl;
		this->root = NULL;
		exit(0);
	}
	return 0;
}

int ABS_FM::loadIdx()
{
	cout<<"loadIdx..."<<endl;

	// cout<<"load tree ..."<<endl;
	// string froot2 = "/home/lab/gll/panaln/mg-aligner/FM/root.idx"; //测试正确性时，载入的单棵小波树，实际运行时不需要20230202
	// loadkit s1(froot2.c_str());
	// LoadWTTree(s1, root); //载入后，用root指针指向（ABS_FM属性）
	// s1.close();

	cout<<"load newD ..."<<endl;
	newD = new Dmap(6971866259);
	string fnewD2 = "/home/lab/gll/panaln/mg-aligner/FM/newD.idx.256"; //blockSize大小256,128,64,32
	loadkit s2(fnewD2.c_str());
	newD->load(s2); //载入后，用newD指针指向（ABS_FM属性）
	s2.close();

	cout<<"load Ltree ..."<<endl;
	string flroot2 = "/home/lab/gll/panaln/mg-aligner/FM/lroot.idx.256";
	loadkit s3(flroot2.c_str());
	LoadWTTree(s3, lroot); //载入后，用lroot指针指向（ABS_FM属性）
	s3.close();

	cout<<"load Rtree ..."<<endl;
	string frroot2 = "/home/lab/gll/panaln/mg-aligner/FM/rroot.idx.256";
	loadkit s4(frroot2.c_str());
	LoadWTTree(s4, rroot); //载入后，用rroot指针指向（ABS_FM属性）
	s4.close();

	cout<<"load SAL ..."<<endl;
	string fsal2 = "/home/lab/gll/panaln/mg-aligner/FM/sal.idx.256";
	loadkit s5(fsal2.c_str());
	this->SAL = new InArray();
	SAL->load(s5); //载入后，用SAL指针指向（ABS_FM属性）
	s5.close();

	this->sa0_index =2932456607;
	this->D =32; //sa采样
	this->n =6971866258; //给ABS_FM属性赋值
	this->alphabetsize =15; //除去N
	this->C =new fm_int[alphabetsize+1]; //多1是，C表的尾部为总数
	mycodetable(); //给C表，编码表赋值
}

int ABS_FM::Load(loadkit &s)
{
	s.loadi64(this->n);
	s.loadi32(this->alphabetsize);
	s.loadi32(this->D);
	s.loadi32(this->r);
	s.loaddouble(this->runs);
	s.loadi64(sa0_index);
	
	//for C
	this->C = new fm_int[alphabetsize+1];
	s.loadi64array(this->C,alphabetsize+1);
	//for code
	this->code = new int[256];//kkzone-bai debug
	s.loadi32array(this->code,256);
	//for codeTable;
	memset(codeTable,0,sizeof(codeTable));
	for(int i=0;i<256;i++)
	{
		uchar len=0;
		s.loadu8(len);
		if(len!=0)
		{
			int bytes = len%8?len/8+1:len/8;
			uchar * bits = new uchar[bytes];
			s.loadu8array(bits,bytes);
			int in_index =0;
			int off_index =0;
			for(int j=0;j<len;j++)
			{
				if(bits[off_index] & (0x01<<(7-in_index)))
					codeTable[i][j] = '1';
				else
					codeTable[i][j] = '0';
				in_index++;
				if(in_index==8)
				{
					in_index =0;
					off_index ++;
				}
			}
		}
	}
	
	//for SAL
	this->SAL = new InArray();
	this->SAL->load(s);
	/* //for Rankl                             //值采样，这部分去掉（对应Save）
	this->RankL = new InArray();           
	this->RankL->load(s);
	//for newD
	this->newD=new Dmap();
	newD->load(s); */
	Inittable();
	uchar * par[4]={Z,R,R1,R3};
	LoadWTTree(s,par);                         //载入小波树
	T=NULL;
	bwt=NULL;
	return 0;
}

int ABS_FM::Save(savekit &s)
{
	s.writei64(n);
	s.writei32(alphabetsize);
	s.writei32(D);//SA的采样率
	s.writei32(r);
	s.writedouble(runs);
	s.writei64(sa0_index);
	
	//C表
	s.writei64array(C,alphabetsize+1);
	//code表
	s.writei32array(code,256);//kkzone-bai debug
	
	//codeTable
	for(int i=0;i<256;i++)
	{
		uchar len = strlen(codeTable[i]);
		s.writeu8(len);
		if(0!=len)
		{
			int bytes = len%8?len/8+1:len/8;
			uchar *bits = new uchar[bytes];
			memset(bits,0,bytes);
			int off_index=0;
			int in_index =0;
			for(int j=0;j<len;j++)
			{
				if(codeTable[i][j]=='1')
					bits[off_index] = bits[off_index]|(0x01<<(7-in_index));
				in_index++;
				if(8==in_index)
				{
					in_index=0;
					off_index++;
				}
			}
			s.writeu8array(bits,bytes);
		}
	}
	SAL->write(s);
	// RankL->write(s);                                     //值采样，去掉这些部分
	// newD->write(s);
	SaveWTTree(s);                                          //保存小波树
	return 0;
}
