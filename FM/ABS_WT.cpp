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
	T = Getfile(filename);                                                                        
	this->filename=(char *)filename; 
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
	return TreeNodeCount(r->Left()) + TreeNodeCount(r->Right()) + 1; 
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
	fm_int *dis=new fm_int[num];
	fm_int *pred=new fm_int[num];
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
		while (newD->getD(i)!=1){
			i = LF(i);
			step++;
			if (left <= i&&i <= right){
				dis[j - left] = step;
				pred[i - left] = j;
				f = 1;
				break;
			}
		}
		if (f == 0){
			fm_int r=newD->rank(i);
			fm_int sal=0;
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
void ABS_FM::Fmocc(fm_int l, fm_int r,unsigned char c ,fm_int *L,fm_int * R)               
{
	int coding = code[c];                                                                  
	if(l == -1 && r == n)
	{
		if(code[c]<0){                                       
			*L = 1 +C[code[c+1]];
			*R = 0 +C[code[c+1]];
			return;
		}
		*L = C[code[c]]+1;                              
		*R = C[code[c] + 1];
		return;
	}
	if(code[c]<0){                                       
		*L = 1 +C[code[c+1]];
		*R = 0 +C[code[c+1]];
		return;
	}
	Occ(c,l,r,*L,*R);                                                                 
	if(*L==-1){
		*L=0;
	}
	*L += C[coding] + 1;                                                                      
	*R += C[coding];
	return;
}
void ABS_FM::Fm16occ(fm_int l, fm_int r,fm_int *L,fm_int * R) 
{	
	if(l == -1 && r == n)                                  
	{
		for(int i=1; i<16; ++i)
		{
			if(code[i]<0){                                       
				L[i] = 1 +C[code[i+1]];
				R[i] = 0 +C[code[i+1]];
				continue;
			}
			L[i] = C[code[i]]+1;                              
			R[i] = C[code[i] + 1];
		}
		return ;
	}
	for(int i=1; i<16; ++i)                                  
	{
		if(code[i]<0){                                       
			L[i] = 1 +C[code[i+1]];
			R[i] = 0 +C[code[i+1]];
			continue;
		}
		Occ(i,l,r,L[i],R[i]);                          
		L[i] += C[code[i]] + 1;                                  
		R[i] += C[code[i]];
	}	
	return;
}
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
	overloop = end%step; 
	anchor = end/step; 
	if(overloop!=0){
		anchor++;
		if(anchor!=RankL->GetNum()-1){
			overloop=step-overloop;
		}else{
			overloop=(n-1)%step-overloop;
		}
	}
	fm_int i= RankL->GetValue(anchor);
	i=newd->select(i); 
	for(int j=0;j<overloop;j++) 
		i = myLF(i);
	for(int j=0;j<len;j++)
	{
		sequence[len-1-j]=myL(i); 
		i = myLF(i);
	}
	return sequence;
}
fm_int ABS_FM::Lookup(fm_int i) 
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
fm_int ABS_FM::FM_Lookup(fm_int i)          
{
	int step = 0;
	int D = this->D;                       
	while(newd->getD(i)!=1)                
	{
		i =myLF(i);                           
		step =step +1;
	}
	fm_int r=newd->rank2(i);
	fm_int sal=0;
	if((n-1)%D==0){                        
		sal=SAL->GetValue(r-1);
	}else{
		sal=SAL->GetValue(r-2);
	}
	return (sal*D+step)%n;
}
void ABS_FM::Occ(unsigned char c,fm_int pos_left,fm_int pos_right,fm_int &rank_left,fm_int &rank_right)
{
	BitMap *r = root;
	int level=0;
	char code = '0';
	while(r->Left())                                                        
	{
		code = codeTable[c][level];
		if(code == '1')                                                     
		{
			if(pos_left>-1 && pos_right >-1) 
			{
				r->Rank(pos_left,pos_right,rank_left,rank_right);           
				pos_left = rank_left -1;
				pos_right = rank_right -1;                                  
			}
			else if(pos_right > -1)
			{
				pos_right=r->Rank(pos_right)-1;
			}
			else
			{
				break;
			}
			if(pos_left==pos_right)
			{
				rank_left=pos_left+1;
				rank_right=pos_right+1;
				return;
			}
			r= r->Right();
		}
		else                                                                
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
			if(pos_left==pos_right)
			{
				rank_left=pos_left+1;
				rank_right=pos_right+1;
				return;
			}
			r=r->Left();
		}
		level++;
	}
	rank_left = pos_left+1;
	rank_right= pos_right+1;
	return ;
}
fm_int ABS_FM::Occ(unsigned char c,fm_int pos)                                    
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
void ABS_FM::Occ1pAux(unsigned char c, fm_int pos1, fm_int pos2, fm_int &occ1, fm_int &occ2) 
{
	BitMap *r=lroot;
	int level=0;
	char code = '0';
	while(r->Left())                                             
	{
		code = codeTable[c][level];
		if(code == '1')                                          
		{
			if(pos1 >-1 && pos2 >-1) {
				r->Rank(pos1, pos2, occ1, occ2);                 
				pos1 = occ1 -1;
				pos2 = occ2 -1;                                  
			}else if(pos2 > -1) {
				pos2 = r->Rank(pos2) - 1;
			}else{
				break;
			}
			if(pos1 == pos2){
				occ1 = pos1 + 1;
				occ2 = pos2 + 1;
				return;
			}
			r = r->Right();
		}
		else                                                     
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
			if(pos1 == pos2) {
				occ1 = pos1 + 1;
				occ2 = pos2 + 1;
				return;
			}
			r = r->Left();
		}
		level++;
	}
	occ1 = pos1 + 1;
	occ2 = pos2 + 1;
	return ;
}
fm_int ABS_FM::Occ1(fm_int i)
{
	fm_int ii =i;
	BitMap *r;
	int bit =0;
	fm_int rank =0;
	unsigned char c =0; 
	if(newD->getD(i)){  
		r = rroot;
		if(i>0) i = newD->rank2(i) - 1; 
	}else{              
		r = lroot;
		if(i>0) i = i - newD->rank2(i);
	}
	while(r->Left())    
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
	if(c==0 && ii>sa0_index) i--;    
	i =i + ((c<10)?C[c]:C[c-1]) + 1; 
	return i;
}
void ABS_FM::Occ1p(unsigned char c, fm_int pos1, fm_int pos2, fm_int &occ1, fm_int &occ2) 
{
	if(pos1 == -1 && pos2 == n){                                 
		occ1 = ((c<10)?C[c]:C[c-1]) + 1;                           
		occ2 = (c<10)?C[c+1]:C[c];
		return;
	}
	BitMap *r;
	if(is_snp_$[c]){                                             
		r = rroot;
		if(pos1>0) pos1 = newD->rank2(pos1) - 1;
		if(pos2>0) pos2 = newD->rank2(pos2) - 1;
			if(pos1 == pos2){ 
				occ1 = pos1 + 1;
				occ2 = pos2 + 1;
				occ1+=((c<10)?C[c]:C[c-1]) + 1; 
				occ2+= (c<10)?C[c]:C[c-1];
				return;
			}
	}else{
		r = lroot;
		if(pos1>0) pos1 = pos1 - newD->rank2(pos1);
		if(pos2>0) pos2 = pos2 - newD->rank2(pos2);
			if(pos1 == pos2){
				occ1 = pos1 + 1;
				occ2 = pos2 + 1;
				occ1+=((c<10)?C[c]:C[c-1]) + 1; 
				occ2+= (c<10)?C[c]:C[c-1];
				return;
			}
	}
	int level=0;
	char code = '0';
	while(r->Left())                                             
	{
		code = codeTable[c][level];
		if(code == '1')                                          
		{
			if(pos1 >-1 && pos2 >-1) {
				r->Rank(pos1, pos2, occ1, occ2);                 
				pos1 = occ1 -1;
				pos2 = occ2 -1;                                  
			}else if(pos2 > -1) {
				pos2 = r->Rank(pos2) - 1;
			}else{
				break;
			}
			if(pos1 == pos2){
				occ1 = pos1 + 1;
				occ2 = pos2 + 1;
				occ1+=((c<10)?C[c]:C[c-1]) + 1; 
				occ2+= (c<10)?C[c]:C[c-1];
				return;
			}
			r = r->Right();
		}
		else                                                     
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
			if(pos1 == pos2) {
				occ1 = pos1 + 1;
				occ2 = pos2 + 1;
				occ1+=((c<10)?C[c]:C[c-1]) + 1; 
				occ2+= (c<10)?C[c]:C[c-1];
				return;
			}
			r = r->Left();
		}
		level++;
	}
	occ1 = pos1 + 1; 
	occ2 = pos2 + 1;
	occ1+=((c<10)?C[c]:C[c-1]) + 1; 
	occ2+= (c<10)?C[c]:C[c-1];
	return ;
}
void ABS_FM::Occ8p(unsigned char c, fm_int pos1, fm_int pos2, fm_int* occ81, fm_int* occ82) 
{ 
	if(pos1 == -1 && pos2 == n){ 
		switch (c) 
		{
		case 0:  
			occ81[0]=C[8]+1; occ81[1]=C[9]+1; occ81[2]=C[10]+1; occ81[3]=C[11]+1; occ81[4]=C[12]+1; occ81[5]=C[13]+1; occ81[6]=C[14]+1; 
			occ82[0]=C[9];   occ82[1]=C[10];  occ82[2]=C[11];   occ82[3]=C[12];   occ82[4]=C[13];   occ82[5]=C[14];   occ82[6]=C[15];
			break;
		case 2: 
			occ81[0]=C[4]+1; occ81[1]=C[5]+1; occ81[2]=C[6]+1; occ81[3]=C[7]+1; occ81[4]=C[8]+1; occ81[5]=C[9]+1; occ81[6]+=C[10]+1; 
			occ82[0]=C[5];   occ82[1]=C[6];   occ82[2]=C[7];   occ82[3]=C[8];   occ82[4]=C[9];   occ82[5]=C[10];  occ82[6]+=C[11];
			break;
		case 1: 
			occ81[0]=C[2]+1; occ81[1]=C[3]+1; occ81[2]=C[4]+1; occ81[3]=C[5]+1; occ81[4]=C[10]+1; occ81[5]=C[11]+1; occ81[6]=C[12]+1; 
			occ82[0]=C[3];   occ82[1]=C[4];   occ82[2]=C[5];   occ82[3]=C[6];   occ82[4]=C[11];   occ82[5]=C[12];   occ82[6]=C[13];
			break;
		case 3: 
			occ81[0]=C[1]+1; occ81[1]=C[2]+1; occ81[2]+=C[5]+1; occ81[3]+=C[6]+1; occ81[4]+=C[9]+1; occ81[5]+=C[12]+1; occ81[6]+=C[13]+1; 
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
	if(p1L !=p2L) Occ1pAux(nt4_gray[c], p1L, p2L, occ1, occ2); 
	if(p1R !=p2R) Occ11p(p1R, p2R, occ111, occ112);
	switch (c)
	{
	case 0:  
		occ81[0]=occ111[5]; occ81[1]=occ111[6]; occ81[2]=occ111[7]; occ81[3]=occ111[8]; occ81[4]=occ111[9]; occ81[5]=occ111[10]; occ81[6]=occ1; 
		occ82[0]=occ112[5]; occ82[1]=occ112[6]; occ82[2]=occ112[7]; occ82[3]=occ112[8]; occ82[4]=occ112[9]; occ82[5]=occ112[10]; occ82[6]=occ2;
		occ81[0]+=C[8]+1; occ81[1]+=C[9]+1; occ81[2]+=C[10]+1; occ81[3]+=C[11]+1; occ81[4]+=C[12]+1; occ81[5]+=C[13]+1; occ81[6]+=C[14]+1;      
		occ82[0]+=C[8];   occ82[1]+=C[9];   occ82[2]+=C[10];   occ82[3]+=C[11];   occ82[4]+=C[12];   occ82[5]+=C[13];   occ82[6]+=C[14];
		break;
	case 2:  
		occ81[0]=occ111[2]; occ81[1]=occ111[3]; occ81[2]=occ111[4]; occ81[3]=occ1; occ81[4]=occ111[5]; occ81[5]=occ111[6]; occ81[6]=occ111[7];  
		occ82[0]=occ112[2]; occ82[1]=occ112[3]; occ82[2]=occ112[4]; occ82[3]=occ2; occ82[4]=occ112[5]; occ82[5]=occ112[6]; occ82[6]=occ112[7];
		occ81[0]+=C[4]+1; occ81[1]+=C[5]+1; occ81[2]+=C[6]+1; occ81[3]+=C[7]+1; occ81[4]+=C[8]+1; occ81[5]+=C[9]+1; occ81[6]+=C[10]+1;          
		occ82[0]+=C[4];   occ82[1]+=C[5];   occ82[2]+=C[6];   occ82[3]+=C[7];   occ82[4]+=C[8];   occ82[5]+=C[9];   occ82[6]+=C[10];
		break;
	case 1:  
		occ81[0]=occ111[1]; occ81[1]=occ1; occ81[2]=occ111[2]; occ81[3]=occ111[3]; occ81[4]=occ111[7]; occ81[5]=occ111[8]; occ81[6]=occ111[9];  
		occ82[0]=occ112[1]; occ82[1]=occ2; occ82[2]=occ112[2]; occ82[3]=occ112[3]; occ82[4]=occ112[7]; occ82[5]=occ112[8]; occ82[6]=occ112[9];
		occ81[0]+=C[2]+1; occ81[1]+=C[3]+1; occ81[2]+=C[4]+1; occ81[3]+=C[5]+1; occ81[4]+=C[10]+1; occ81[5]+=C[11]+1; occ81[6]+=C[12]+1;        
		occ82[0]+=C[2];   occ82[1]+=C[3];   occ82[2]+=C[4];   occ82[3]+=C[5];   occ82[4]+=C[10];   occ82[5]+=C[11];   occ82[6]+=C[12];
		break;
	case 3: 
		occ81[0]=occ1; occ81[1]=occ111[1]; occ81[2]=occ111[3]; occ81[3]=occ111[4]; occ81[4]=occ111[6]; occ81[5]=occ111[9]; occ81[6]=occ111[10]; 
		occ82[0]=occ2; occ82[1]=occ112[1]; occ82[2]=occ112[3]; occ82[3]=occ112[4]; occ82[4]=occ112[6]; occ82[5]=occ112[9]; occ82[6]=occ112[10];
		occ81[0]+=C[1]+1; occ81[1]+=C[2]+1; occ81[2]+=C[5]+1; occ81[3]+=C[6]+1; occ81[4]+=C[9]+1; occ81[5]+=C[12]+1; occ81[6]+=C[13]+1;         
		occ82[0]+=C[1];   occ82[1]+=C[2];   occ82[2]+=C[5];   occ82[3]+=C[6];   occ82[4]+=C[9];   occ82[5]+=C[12];   occ82[6]+=C[13];
		break;
	default:
		break;
	}
}
void ABS_FM::Occ16p(fm_int pos1, fm_int pos2, fm_int* occ161, fm_int* occ162) 
{
	if(pos1 == -1 && pos2 == n){
		for(int i=0; i<16; i++){ 
			if(i==10){
				occ161[10]=0; 
				occ162[10]=-1; 
				continue; 
			}
			occ161[i] = ((i<10)?C[i]:C[i-1]) + 1;                           
			occ162[i] = (i<10)?C[i+1]:C[i];
		}
		return;
	}
	long long p1R=0, p1L=0, p2R=0, p2L=0, occ41[4]={0}, occ42[4]={0}, occ111[11]={0}, occ112[11]={0};
	p1R=newD->rank2(pos1); p1L=pos1-p1R; p1R--; 
	p2R=newD->rank2(pos2); p2L=pos2-p2R; p2R--;
	if(p1L !=p2L) Occ4p(p1L, p2L, occ41, occ42);            
	if(p1R !=p2R) Occ11p(p1R, p2R, occ111, occ112);         
	occ161[0]=occ111[0]; occ161[1]=occ41[3];  occ161[2]=occ111[1]; occ161[3]=occ41[2];   occ161[4]=occ111[2];  occ161[5]=occ111[3];  occ161[6]=occ111[4];   occ161[7]=occ41[1]; 
	occ161[8]=occ111[5]; occ161[9]=occ111[6]; occ161[10]=19;       occ161[11]=occ111[7]; occ161[12]=occ111[8]; occ161[13]=occ111[9]; occ161[14]=occ111[10]; occ161[15]=occ41[0]; 
	occ162[0]=occ112[0]; occ162[1]=occ42[3];  occ162[2]=occ112[1]; occ162[3]=occ42[2];   occ162[4]=occ112[2];  occ162[5]=occ112[3];  occ162[6]=occ112[4];   occ162[7]=occ42[1];
	occ162[8]=occ112[5]; occ162[9]=occ112[6]; occ162[10]=18;       occ162[11]=occ112[7]; occ162[12]=occ112[8]; occ162[13]=occ112[9]; occ162[14]=occ112[10]; occ162[15]=occ42[0];
	occ161[0]+=C[0]+1; occ161[1]+=C[1]+1; occ161[2]+=C[2]+1; occ161[3]+=C[3]+1;   occ161[4]+=C[4]+1;   occ161[5]+=C[5]+1;   occ161[6]+=C[6]+1;   occ161[7]+=C[7]+1; 
	occ161[8]+=C[8]+1; occ161[9]+=C[9]+1;                    occ161[11]+=C[10]+1; occ161[12]+=C[11]+1; occ161[13]+=C[12]+1; occ161[14]+=C[13]+1; occ161[15]+=C[14]+1;
	occ162[0]+=C[0]; occ162[1]+=C[1]; occ162[2]+=C[2]; occ162[3]+=C[3];   occ162[4]+=C[4];   occ162[5]+=C[5];   occ162[6]+=C[6];   occ162[7]+=C[7]; 
	occ162[8]+=C[8]; occ162[9]+=C[9];                  occ162[11]+=C[10]; occ162[12]+=C[11]; occ162[13]+=C[12]; occ162[14]+=C[13]; occ162[15]+=C[14];    
}
void ABS_FM::Occ4p(fm_int pos1, fm_int pos2, fm_int *occ41, fm_int *occ42) 
{   
	BitMap *r = lroot;
	int f_k1_h0=1, f_k0_h0=1;
	fm_int p1_k0_h0, p1_k1_h0, p2_k0_h0, p2_k1_h0;
	fm_int p1_k0_h1_1, p1_k1_h1_1, p2_k0_h1_1, p2_k1_h1_1;
	fm_int p1_k0_h1_0, p1_k1_h1_0, p2_k0_h1_0, p2_k1_h1_0;
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
	if(p1_k1_h0 == p2_k1_h0){                                                      
		f_k1_h0 = 0;
		p1_k1_h1_1 = p2_k1_h1_1 = p1_k0_h1_1 = p2_k0_h1_1 = -1;                    
	}
	if(p1_k0_h0 == p2_k0_h0){
		f_k0_h0 = 0;
		p1_k1_h1_0 = p2_k1_h1_0 = p1_k0_h1_0 = p2_k0_h1_0 = -1;                    
	}
	if(f_k1_h0){
		BitMap *rRoot = r->Right();
		if(p1_k1_h0 > -1 && p2_k1_h0 > -1){
			rRoot->Rank(p1_k1_h0, p2_k1_h0, p1_k1_h1_1, p2_k1_h1_1);                   
			p1_k0_h1_1 = p1_k1_h0 - p1_k1_h1_1; p2_k0_h1_1 = p2_k1_h0 - p2_k1_h1_1;    
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
			lRoot->Rank(p1_k0_h0, p2_k0_h0, p1_k1_h1_0, p2_k1_h1_0);                   
			p1_k0_h1_0 = p1_k0_h0 - p1_k1_h1_0; p2_k0_h1_0 = p2_k0_h0 - p2_k1_h1_0;    
		}else if(p2_k0_h0 > -1){
			p2_k1_h1_0 = lRoot->Rank(p2_k0_h0); p2_k0_h1_0 = p2_k0_h0 - p2_k1_h1_0;
			p1_k1_h1_0 = 0; p1_k0_h1_0 = -1;
		}else{
			p1_k1_h1_0 = 0; p1_k0_h1_0 = -1;
			p2_k1_h1_0 = 0; p2_k0_h1_0 = -1;
		}
		p1_k1_h1_0--; p2_k1_h1_0--;
	}
	occ41[0] = p1_k0_h1_0; occ41[1] = p1_k1_h1_0; occ41[2] = p1_k0_h1_1; occ41[3] = p1_k1_h1_1; 
	occ42[0] = p2_k0_h1_0; occ42[1] = p2_k1_h1_0; occ42[2] = p2_k0_h1_1; occ42[3] = p2_k1_h1_1;
	for(int i=0; i<4; i++){
		occ41[i]++; occ42[i]++;
	}
}
void ABS_FM::Occ11p(fm_int pos1, fm_int pos2, fm_int *occ111, fm_int *occ112) 
{	
	BitMap *r = rroot, *rRoot, *lRoot, *rRoot0, *lRoot0, *rRoot01, *rRoot00, *lRoot001, *rRoot0010, *lRoot0010;
	int f_k1_h0=1, f_k0_h0=1, f_k1_h1_0=1, f_k0_h1_0=1, f_k1_h2_01=1, f_k1_h2_00=1, f_k0_h3_001=1, f_k1_h4_0010=1, f_k0_h4_0010=1; 
	fm_int p1_k0_h0, p1_k1_h0, p2_k0_h0, p2_k1_h0;
	fm_int p1_k0_h1_1, p1_k1_h1_1, p2_k0_h1_1, p2_k1_h1_1;
	fm_int p1_k0_h1_0, p1_k1_h1_0, p2_k0_h1_0, p2_k1_h1_0;
	fm_int p1_k0_h2_01, p1_k1_h2_01, p2_k0_h2_01, p2_k1_h2_01;
	fm_int p1_k0_h2_00, p1_k1_h2_00, p2_k0_h2_00, p2_k1_h2_00;
	fm_int p1_k0_h3_011, p1_k1_h3_011, p2_k0_h3_011, p2_k1_h3_011;
	fm_int p1_k0_h3_001, p1_k1_h3_001, p2_k0_h3_001, p2_k1_h3_001;
	fm_int p1_k0_h4_0010, p1_k1_h4_0010, p2_k0_h4_0010, p2_k1_h4_0010;
	fm_int p1_k0_h5_00101, p1_k1_h5_00101, p2_k0_h5_00101, p2_k1_h5_00101;
	fm_int p1_k0_h5_00100, p1_k1_h5_00100, p2_k0_h5_00100, p2_k1_h5_00100; 
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
	if(p1_k1_h0 == p2_k1_h0){                                                                                      
		f_k1_h0 = 0;
		p1_k1_h1_1 = p2_k1_h1_1 = p1_k0_h1_1 = p2_k0_h1_1 = -1; 
	}
	if(p1_k0_h0 == p2_k0_h0){
		f_k0_h0 = f_k1_h1_0 = f_k0_h1_0 = f_k1_h2_01 = f_k1_h2_00 = f_k0_h3_001 = f_k1_h4_0010 = f_k0_h4_0010 = 0; 
		p1_k1_h3_011 = p2_k1_h3_011 = p1_k0_h3_011 = p2_k0_h3_011 = p1_k0_h2_01 = p2_k0_h2_01 = p1_k1_h3_001 = p2_k1_h3_001 = p1_k1_h5_00101 = p2_k1_h5_00101 \
		= p1_k0_h5_00101 = p2_k0_h5_00101 = p1_k1_h5_00100 = p2_k1_h5_00100 = p1_k0_h5_00100 = p2_k0_h5_00100 = p1_k0_h2_00 = p2_k0_h2_00 = -1; 
	}
	if(f_k1_h0){
		rRoot = r->Right();
		if(p1_k1_h0 > -1 && p2_k1_h0 > -1){
			rRoot->Rank(p1_k1_h0, p2_k1_h0, p1_k1_h1_1, p2_k1_h1_1);                                                   
			p1_k0_h1_1 = p1_k1_h0 - p1_k1_h1_1; p2_k0_h1_1 = p2_k1_h0 - p2_k1_h1_1;                                    
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
			p2_k1_h1_0 = 0; p2_k0_h1_0 = -1; 
		}
		p1_k1_h1_0--; p2_k1_h1_0--;
		if(p1_k1_h1_0 == p2_k1_h1_0){
			f_k1_h1_0 = f_k1_h2_01 = 0;
			p1_k1_h3_011 = p2_k1_h3_011 = p1_k0_h3_011 = p2_k0_h3_011 = p1_k0_h2_01 = p2_k0_h2_01 = -1; 
		}
		if(p1_k0_h1_0 == p2_k0_h1_0){
			f_k0_h1_0 = f_k1_h2_00 = f_k0_h3_001 = f_k1_h4_0010 = f_k0_h4_0010 = 0;
			p1_k1_h3_001 = p2_k1_h3_001 = p1_k1_h5_00101 = p2_k1_h5_00101 = p1_k0_h5_00101 = p2_k0_h5_00101 = p1_k1_h5_00100 = p2_k1_h5_00100 \
			= p1_k0_h5_00100 = p2_k0_h5_00100 = p1_k0_h2_00 = p2_k0_h2_00 = -1; 
		}
	}
	if(f_k1_h1_0){
		rRoot0 = lRoot->Right();
		if(p1_k1_h1_0 > -1 && p2_k1_h1_0 > -1){
			rRoot0->Rank(p1_k1_h1_0, p2_k1_h1_0, p1_k1_h2_01, p2_k1_h2_01);
			p1_k0_h2_01 = p1_k1_h1_0 - p1_k1_h2_01; p2_k0_h2_01 = p2_k1_h1_0 - p2_k1_h2_01;                
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
			p1_k1_h3_011 = p2_k1_h3_011 = p1_k0_h3_011 = p2_k0_h3_011 = -1; 
		}
	}
	if(f_k0_h1_0){
		lRoot0 = lRoot->Left();
		if(p1_k0_h1_0 > -1 && p2_k0_h1_0 > -1){
			lRoot0->Rank(p1_k0_h1_0, p2_k0_h1_0, p1_k1_h2_00, p2_k1_h2_00);
			p1_k0_h2_00 = p1_k0_h1_0 - p1_k1_h2_00; p2_k0_h2_00 = p2_k0_h1_0 - p2_k1_h2_00;                
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
			p1_k1_h3_001 = p2_k1_h3_001 = p1_k1_h5_00101 = p2_k1_h5_00101 = p1_k0_h5_00101 = p2_k0_h5_00101 = p1_k1_h5_00100 = p2_k1_h5_00100 = p1_k0_h5_00100 = p2_k0_h5_00100 = -1; 
		}
	}
	if(f_k1_h2_01){
		rRoot01 = rRoot0->Right();
		if(p1_k1_h2_01 > -1 && p2_k1_h2_01 > -1){
			rRoot01->Rank(p1_k1_h2_01, p2_k1_h2_01, p1_k1_h3_011, p2_k1_h3_011);                           
			p1_k0_h3_011 = p1_k1_h2_01 - p1_k1_h3_011; p2_k0_h3_011 = p2_k1_h2_01 - p2_k1_h3_011;          
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
			rRoot00->Rank(p1_k1_h2_00, p2_k1_h2_00, p1_k1_h3_001, p2_k1_h3_001);                           
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
			p1_k1_h5_00101 = p2_k1_h5_00101 = p1_k0_h5_00101 = p2_k0_h5_00101 = p1_k1_h5_00100 = p2_k1_h5_00100 = p1_k0_h5_00100 = p2_k0_h5_00100 = -1; 
		}
	}
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
			p1_k1_h5_00101 = p2_k1_h5_00101 = p1_k0_h5_00101 = p2_k0_h5_00101 = -1; 
		}
		if(p1_k0_h4_0010 == p2_k0_h4_0010){
			f_k0_h4_0010 = 0;
			p1_k1_h5_00100 = p2_k1_h5_00100 = p1_k0_h5_00100 = p2_k0_h5_00100 = -1; 
		}
	}
	if(f_k1_h4_0010){
		rRoot0010 = lRoot001->Right();
		if(p1_k1_h4_0010 > -1 && p2_k1_h4_0010 > -1){
			rRoot0010->Rank(p1_k1_h4_0010, p2_k1_h4_0010, p1_k1_h5_00101, p2_k1_h5_00101);                       
			p1_k0_h5_00101 = p1_k1_h4_0010 - p1_k1_h5_00101; p2_k0_h5_00101 = p2_k1_h4_0010 - p2_k1_h5_00101;    
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
			lRoot0010->Rank(p1_k0_h4_0010, p2_k0_h4_0010, p1_k1_h5_00100, p2_k1_h5_00100);                       
			p1_k0_h5_00100 = p1_k0_h4_0010 - p1_k1_h5_00100; p2_k0_h5_00100 = p2_k0_h4_0010 - p2_k1_h5_00100;    
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
	occ112[6] = p2_k1_h5_00100; occ112[7] = p2_k1_h5_00101; occ112[8] = p2_k0_h1_1; occ112[9] = p2_k0_h5_00100; occ112[10] = p2_k1_h3_001; 
	for(int i=0; i<11; i++){
		occ111[i]++; occ112[i]++;
	}
}
void ABS_FM::Occ4(fm_int pos, fm_int* occ4)
{   
	BitMap *r = lroot;
	fm_int k0_h0, k1_h0;
	k1_h0 = r->Rank(pos) - 1;
	k0_h0 = pos - k1_h0 - 1;
	BitMap *rRoot = r->Right();
	fm_int k0_h1_1, k1_h1_1;
	k1_h1_1 = rRoot->Rank(k1_h0) - 1;           
	k0_h1_1 = k1_h0 - k1_h1_1 - 1;              
	BitMap *lRoot = r->Left();
	fm_int k0_h1_0, k1_h1_0;
	k1_h1_0 = lRoot->Rank(k0_h0) - 1;           
	k0_h1_0 = k0_h0 - k1_h1_0 - 1;              
	occ4[0] = k0_h1_0; occ4[1] = k1_h1_0; occ4[2] = k0_h1_1; occ4[3] = k1_h1_1; 
}
void ABS_FM::Occ11(fm_int pos, fm_int* occ11)
{	
	BitMap *r = rroot;
	fm_int k0_h0, k1_h0;
	k1_h0 = r->Rank(pos) - 1;
	k0_h0 = pos - k1_h0 - 1;
	BitMap *rRoot = r->Right();
	fm_int k0_h1_1, k1_h1_1;
	k1_h1_1 = rRoot->Rank(k1_h0) - 1;               
	k0_h1_1 = k1_h0 - k1_h1_1 - 1;                  
	BitMap *lRoot = r->Left();
	fm_int k0_h1_0, k1_h1_0;
	k1_h1_0 = lRoot->Rank(k0_h0) - 1;
	k0_h1_0 = k0_h0 - k1_h1_0 - 1;
	BitMap *rRoot0 = lRoot->Right();
	fm_int k0_h2_01, k1_h2_01;
	k1_h2_01 = rRoot0->Rank(k1_h1_0) - 1;
	k0_h2_01 = k1_h1_0 - k1_h2_01 - 1;             
	BitMap *lRoot0 = lRoot->Left();
	fm_int k0_h2_00, k1_h2_00;
	k1_h2_00 = lRoot0->Rank(k0_h1_0) - 1;
	k0_h2_00 = k0_h1_0 - k1_h2_00 - 1;             
	BitMap *rRoot01 = rRoot0->Right();
	fm_int k0_h3_011, k1_h3_011;
	k1_h3_011 = rRoot01->Rank(k1_h2_01) - 1;       
	k0_h3_011 = k1_h2_01 - k1_h3_011 - 1;          
	BitMap *rRoot00 = lRoot0->Right();
	fm_int k0_h3_001, k1_h3_001;
	k1_h3_001 = rRoot00->Rank(k1_h2_00) - 1;       
	k0_h3_001 = k1_h2_00 - k1_h3_001 - 1;
	BitMap *lRoot001 = rRoot00->Left();
	fm_int k0_h4_0010, k1_h4_0010;
	k1_h4_0010 = lRoot001->Rank(k0_h3_001) - 1;
	k0_h4_0010 = k0_h3_001 - k1_h4_0010 - 1;
	BitMap *rRoot0010 = lRoot001->Right();
	fm_int k0_h5_00101, k1_h5_00101;
	k1_h5_00101 = rRoot0010->Rank(k1_h4_0010) - 1; 
	k0_h5_00101 = k1_h4_0010 - k1_h5_00101 - 1;    
	BitMap *lRoot0010 = lRoot001->Left();
	fm_int k0_h5_00100, k1_h5_00100;
	k1_h5_00100 = lRoot0010->Rank(k0_h4_0010) - 1; 
	k0_h5_00100 = k0_h4_0010 - k1_h5_00100 - 1;    
	occ11[0] = k0_h2_01; occ11[1] = k0_h3_011; occ11[2] = k0_h2_00; occ11[3] = k0_h5_00101; occ11[4] = k1_h1_1; occ11[5] = k1_h3_011;
	occ11[6] = k1_h5_00100; occ11[7] = k1_h5_00101; occ11[8] = k0_h1_1; occ11[9] = k0_h5_00100; occ11[10] = k1_h3_001; 
}
fm_int ABS_FM::myLF(fm_int i)  
{	
	if(i==sa0_index){
		return 0;
	}
	return Occ1(i);
}
fm_int ABS_FM::LF(fm_int i)                                                   
{	
	if(i==sa0_index){
		return 0;
	}
	fm_int occ =0;
	unsigned char label =0;
	Occ(occ,label,i);                                                         
	int coding = code[label];                                                 
	return occ + C[coding];                                                
}
unsigned char ABS_FM::myL(fm_int i) 
{
	fm_int ii =i;
	BitMap *r;
	int bit =0;
	fm_int rank =0;
	unsigned char c =0; 
	if(newD->getD(i)){  
		r = rroot;
		if(i>0) i = newD->rank2(i) - 1; 
	}else{              
		r = lroot;
		if(i>0) i = i - newD->rank2(i);
	}
	while(r->Left())    
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
	return c;
}
fm_int ABS_FM::Occ(fm_int & occ , unsigned char & label,fm_int pos)          
{
	BitMap * r = root;
	int bit =0;
	fm_int rank =0;
	while(r->Left())                                                         
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
	label = r->Label();                                                     
	return 0;
}
unsigned char * ABS_FM::Getfile(const char *filename)
{
	this->filename = (char*) malloc(strlen(filename)+1); 
	strcpy(this->filename,filename);
	FILE * fp = fopen(this->filename,"r+");
	if(fp==NULL)
	{
		cout<<"Be sure the file is available "<<this->filename<<endl;
		exit(0);
	}
	fseek(fp,0,SEEK_END);
	this->n = ftell(fp);                    
	unsigned char * T = new unsigned char[n];     
	fseek(fp,0,SEEK_SET);
	fm_int e=0;
	fm_int num=0;
	while((e=fread(T+num,sizeof(uchar),n-num,fp))!=0){          
		num = num +e;
	}
	if(num!=n)
	{
		cout<<num<<"  "<<n<<"  Read source file failed"<<endl;
		exit(0);
	}
	fclose(fp);
	memset(charFreq,0,256*sizeof(long long));                  
	memset(charMap,0,256*sizeof(bool));
	for(fm_int i=0;i<n;i++){                                   
		charFreq[T[i]]++;                                      
	}
	this->alphabetsize = 0;
	for(int i=0;i<256;i++)
		if(charFreq[i]!=0)
		{
			this->alphabetsize++;                              
			this->charMap[i]=true;                             
		}
	this->code = new int[256];
	this->C = new fm_int[alphabetsize+1];                      
	memset(C,0,(alphabetsize+1)*4);
	this->C[alphabetsize] = n;                                 
	this->C[0] = 0;                                            
	int k=1;
	fm_int pre =0;
	for(int i=0;i<256;i++)                                     
	{
		if(charFreq[i]!=0)
		{
			code[i]=k-1;                                       
			C[k]=pre + charFreq[i];
			pre = C[k];
			k++;
		}
		else
			code[i]=-1;                                        
	}
	return T;
}
int ABS_FM::BWT(unsigned char *T,saidx64_t * SA,unsigned char * bwt,fm_int len)
{
	for (bwtint_t i = 0; i <= len; i++) {
		if (SA[i] == 0) {                                                         
			sa0_index = i;                                                        
			bwt[i] = nt16_table['$'];                                              
		}
		else {
			bwt[i] = T[SA[i] - 1];                                                 
		}
	}
	return 0;
}
void ABS_FM::mycodetable(){ 
	memset(codeTable, 0, CHAR_SET_SIZE * CHAR_SET_SIZE);
	codeTable[nt16_table['A']][0] = '0'; codeTable[nt16_table['A']][1] = '0'; 
	codeTable[nt16_table['C']][0] = '0'; codeTable[nt16_table['C']][1] = '1';
	codeTable[nt16_table['G']][0] = '1'; codeTable[nt16_table['G']][1] = '0'; 
	codeTable[nt16_table['T']][0] = '1'; codeTable[nt16_table['T']][1] = '1'; 
	codeTable[nt16_table['R']][0] = '1'; codeTable[nt16_table['R']][1] = '0';
	codeTable[nt16_table['Y']][0] = '1'; codeTable[nt16_table['Y']][1] = '1';
	codeTable[nt16_table['$']][0] = '0'; codeTable[nt16_table['$']][1] = '1'; codeTable[nt16_table['$']][2] = '0';
	codeTable[nt16_table['S']][0] = '0'; codeTable[nt16_table['S']][1] = '0'; codeTable[nt16_table['S']][2] = '0'; 
	codeTable[nt16_table['K']][0] = '0'; codeTable[nt16_table['K']][1] = '1'; codeTable[nt16_table['K']][2] = '1'; codeTable[nt16_table['K']][3] = '0';
	codeTable[nt16_table['M']][0] = '0'; codeTable[nt16_table['M']][1] = '1'; codeTable[nt16_table['M']][2] = '1'; codeTable[nt16_table['M']][3] = '1';
	codeTable[nt16_table['W']][0] = '0'; codeTable[nt16_table['W']][1] = '0'; codeTable[nt16_table['W']][2] = '1'; codeTable[nt16_table['W']][3] = '1';
	codeTable[nt16_table['B']][0] = '0'; codeTable[nt16_table['B']][1] = '0'; codeTable[nt16_table['B']][2] = '1'; codeTable[nt16_table['B']][3] = '0'; codeTable[nt16_table['B']][4] = '1'; codeTable[nt16_table['B']][5] = '0';
	codeTable[nt16_table['V']][0] = '0'; codeTable[nt16_table['V']][1] = '0'; codeTable[nt16_table['V']][2] = '1'; codeTable[nt16_table['V']][3] = '0'; codeTable[nt16_table['V']][4] = '1'; codeTable[nt16_table['V']][5] = '1';
	codeTable[nt16_table['D']][0] = '0'; codeTable[nt16_table['D']][1] = '0'; codeTable[nt16_table['D']][2] = '1'; codeTable[nt16_table['D']][3] = '0'; codeTable[nt16_table['D']][4] = '0'; codeTable[nt16_table['D']][5] = '0';
	codeTable[nt16_table['H']][0] = '0'; codeTable[nt16_table['H']][1] = '0'; codeTable[nt16_table['H']][2] = '1'; codeTable[nt16_table['H']][3] = '0'; codeTable[nt16_table['H']][4] = '0'; codeTable[nt16_table['H']][5] = '1';
}
void ABS_FM::huffmanTreecode(){
	memset(codeTable, 0, CHAR_SET_SIZE * CHAR_SET_SIZE);
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
int ABS_FM::BuildTree()                                    
{	
	saidx64_t *SA = new saidx64_t[n+1];
	SA[0] = n;
	divsufsort64(T,SA+1,n);		
	int step =this->D; 
	SAL=new InArray((n+1)/step+2,blog((n+1)/step)+1);
	RankL=new InArray((n+1)/step+2,blog((n+1)/step)+1);
	newd=new Dmap(n+1);
	fm_int i=0, j=0;
	for(i=0; i<=n; i++){
		if(SA[i]%step ==0){
			SAL->SetValue(j++,SA[i]/step);                
			newd->SetValue(i,1);                          
		}else if(i==0){
			newd->SetValue(i,1);
		}else{
			newd->SetValue(i,0);
		}
	}
	newd->constructrank(256,16); 
	if((n)%step!=0){ 
		RankL->SetValue((n)/step+1,newd->rank2(0)); 
	}
	for(i=0; i<=n; i++)
	{
		if(SA[i]%step==0){
			RankL->SetValue(SA[i]/step,newd->rank2(i));   
		}
	}
	string prefix(this->filename);
	prefix = prefix.substr(0, prefix.find('.')); 
	cout<<"Construct Sal ..."<<endl;
	string fsal = prefix+".1.idx"; 
	savekit s4(fsal.c_str());
	SAL->write(s4); 
	SAL->~InArray();
	SAL=NULL;
	s4.close();
	cout<<"Construct Rankl ..."<<endl;
	string frankl = prefix+".2.idx";
	savekit s5(frankl.c_str());
	RankL->write(s5);
	RankL->~InArray();
	RankL=NULL;
	s5.close();
	cout<<"Construct D ..."<<endl;
	string fnewd = prefix+".3.idx";
	savekit s6(fnewd.c_str());
	newd->write(s6);
	newd->~Dmap();
	newd=NULL;
	s6.close();
	bwt = new unsigned char[n+1];                                                                                  
	BWT(T,SA,bwt,n);
	string finfo = prefix+".7.idx"; 
	ofstream inf(finfo);
	inf << this->sa0_index <<endl;
	inf << this->D <<endl;
	inf << this->n <<endl;
	inf << this->alphabetsize <<endl;
	for(int i=0; i<alphabetsize+1; i++){
		inf<< this->C[i] <<endl;                
	}
	inf.close();
	delete [] T; 
	T=NULL;
	delete [] SA; 
	SA=NULL;
	newD = new Dmap(n+1);	
	unsigned char* lptr = new unsigned char[n+1];                        
	unsigned char* rptr = new unsigned char[n+1];                        
	memset(lptr, 0, n+1);  memset(rptr, 0, n+1);
	long rightLen = 0, leftLen = 0;
	unsigned char c;
	for(long i=0; i<n+1; i++){
		c =bwt[i];
		if(is_snp_$[c]){
			newD->SetValue(i,1);
			rptr[rightLen++] = c;                                      
		}else{
			lptr[leftLen++] = c;                                       
		}
	}
	delete[] bwt; 
	bwt=NULL;
	cout<<"Construct U ..."<<endl;
	newD->constructrank(1024,16); 
	string fnewD = prefix+".4.idx";
	savekit s(fnewD.c_str());
	newD->write(s);
	newD->~Dmap(); 
	s.close();
	mycodetable(); 
	cout<<"Construct Ltree ..."<<endl;
	lroot = CreateWaveletTree(lptr, leftLen);
	string flroot = prefix+".5.idx";
	savekit s2(flroot.c_str());
	SaveWTTree(s2, lroot);
	delete[] lptr; 
	lptr=NULL;
	s2.close();
	mycodetable(); 
	cout<<"Construct Rtree ..."<<endl;
	rroot = CreateWaveletTree(rptr, rightLen);
	string frroot = prefix+".6.idx";
	savekit s3(frroot.c_str());
	SaveWTTree(s3, rroot);
	delete[] rptr; 
	rptr=NULL;
	s.close();
	std::cerr << "The generated file is in: "<< prefix.substr(0, prefix.find_last_of('/'))<<"/"<< std::endl;
	return 0;
}
BitMap * ABS_FM::CreateWaveletTree(unsigned char * bwt,fm_int n)
{
	BitMap * root = NULL;                                  
	root = FullFillWTNode(bwt,n,0);                        
	if(!root)
	{
		cout<<"FullfillWTNode failed"<<endl;               
		exit(0);
	}
	return root;                                           
}
BitMap * ABS_FM::FullFillWTNode(unsigned char * buff,fm_int len,int level)            
{
	int CurrentLevel = level;                                                         
	long long CurrentBitLen = len;                                                    
	unsigned char CurrentLabel = '\0';                                                
	unsigned long long int *CurrentBitBuff = NULL;                                    
	if ((int)strlen((const char*)codeTable[buff[0]])==level)                                                            
	{
		CurrentLabel = buff[0];                                                                                         
		CurrentBitBuff = NULL;
		uchar * tables[4] ={this->Z,this->R, this->R1, this->R3};                                                       
		BitMap * node = new BitMap(CurrentBitBuff,CurrentBitLen,CurrentLevel,block_size,r,CurrentLabel,tables);         
		node->Left(NULL);                                                                                               
		node->Right(NULL);                                                                                              
		return node;                                                                                                    
	}
	fm_int u64Len=0;
	if(len%64==0)
		u64Len = len/64+1;
	else
		u64Len = len/64+2;
	CurrentBitBuff = new unsigned long long int[u64Len];                                                                
	memset(CurrentBitBuff,0,u64Len*8);                                                                                  
	unsigned char * lptr=NULL;                                                                                          
	unsigned char * rptr=NULL;                                                                                          
	fm_int leftLen=0;                                                                                                   
	fm_int rightLen=0;                                                                                                  
	lptr = new unsigned char[len];
	rptr = new unsigned char[len];
	memset(lptr,0,len);
	memset(rptr,0,len);
	fm_int i=0;
	fm_int bytePos=0;
	int bitOffset=0;
	u64 last = 0;
	for(i=0;i<len;i++)                                                                                                  
	{
		if(codeTable[buff[i]][level]=='1')                                                                              
		{
			CurrentBitBuff[bytePos] |= (0x01ull<<(63-bitOffset));                                                       
			rptr[rightLen++]=buff[i];                                                                                   
			last = 0;
		}
		else                                                                                                            
		{
			lptr[leftLen++]=buff[i];                                                                                    
			last = 1;
		}
		bitOffset++;                                                                                                    
		if(bitOffset == 64)
		{
			bytePos++;                                                                                                  
			bitOffset = 0;
		}
	}
	CurrentBitBuff[bytePos] |= (last<<(63-bitOffset));                                                                 
	uchar * tables[4] = {this->Z,this->R, this->R1, this->R3};
	BitMap * node = new BitMap(CurrentBitBuff,CurrentBitLen,CurrentLevel,block_size,r,CurrentLabel,tables);            
	if(leftLen !=0)
	{
		BitMap * left =FullFillWTNode(lptr,leftLen,level+1);                                                           
		node->Left(left);                                                                                              
		delete [] lptr;
		lptr=NULL;
	}
	if(rightLen!=0)
	{
		BitMap * right = FullFillWTNode(rptr,rightLen,level+1);                                                        
		node->Right(right);                                                                                            
		delete [] rptr;
		rptr=NULL;
	}
	return node;                                                                                                       
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
	u64 B[2]={0xffffffffffffffffull,0xffffffffffffffffull};
	int sum =0;
	int step=0;
	int rank = 0;
	int runs = 0 ;
	int x = 0;
	int prestep = 0;
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
			x = GammaDecode(B,step,this);
			if(step > 16)
				break;
			sum = sum + x;
			prestep = step;
			runs ++;
			if(runs%2==1)
				rank = rank + x;
		}
		R[i<<2] = runs;
		R[(i<<2)+1] = prestep;
		R[(i<<2)+2] = sum; 
		R[(i<<2)+3] = rank;
		if(i==0){
			R1[0] = 16;
			R3[0] = 16;
			continue;
		}
		buildEFtable(i);
	}
}
inline int popcnt1(unsigned long long int x) 
{
	x = x - ((x & 0xAAAAAAAAAAAAAAAA) >> 1);
	x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
	x = ((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F);
	return (x * 0x0101010101010101) >> 56;
}
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
	x = i;
	int zeronums = 0;
	int alreadyBits = 0;
	while(alreadyBits < R3[i*18]){
		int num_0 = Zeros(x,this);
		zeronums += num_0;
		alreadyBits += num_0;
		x = x<<num_0;
		int num_1 = Zeros(~x,this);
		R3[i*18+zeronums+1] = num_1;
		alreadyBits += num_1;
		x = x<<num_1;
	}
}
int ABS_FM::SaveNodePosition(BitMap *r, int position, savekit &s)
{
	if(!r)
		return 1;
	s.writei32(position);                            
	SaveNodePosition(r->Left(),  2*position, s);     
	SaveNodePosition(r->Right(), 2*position+1, s);   
	return 0;
}
int ABS_FM::SaveNodeData(BitMap *r,savekit &s)
{
	if(!r)
		return 1 ;
	r->Save(s);                                     
	SaveNodeData(r->Left(),s);                      
	SaveNodeData(r->Right(),s);                     
	return 0;
}
int ABS_FM::SaveWTTree(savekit &s, BitMap *root) 
{
	int nodecount = TreeNodeCount(root);
	s.writei32(nodecount);
	SaveNodePosition(root,1,s);                  
	SaveNodeData(root,s);                        
	return 0;
}
int ABS_FM::SaveWTTree(savekit &s)
{
	int nodecount = TreeNodeCount(root); 
	s.writei32(nodecount);
	SaveNodePosition(root,1,s);
	SaveNodeData(root,s);
	return 0;
}
int ABS_FM::LoadWTTree(loadkit &s, BitMap* &root, uchar **tables) 
{
	int nodecount = 0;
	s.loadi32(nodecount);
	int *p = new int[nodecount];
	s.loadi32array(p, nodecount);                        
	map<int, BitMap *> pmap;                             
	BitMap * r = NULL;                                   
	for(int i=0;i<nodecount;i++)
	{
		if(tables)
			r = new BitMap(tables);                      
		else
			r = new BitMap();
		r->Load(s);                                      
		pmap[p[i]] = r;                                  
	}
	map<int ,BitMap *>::iterator iter;                   
	map<int ,BitMap *>::iterator f_iter;                 
	for(iter=pmap.begin(); iter!=pmap.end(); iter++)
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
	if(f_iter !=pmap.end()){
		root = f_iter->second;                           
	}else{
		cerr<<"Load WTTree error"<<endl;
		root = NULL;
		exit(0);                                         
	}
	return 0;
}
int ABS_FM::LoadWTTree(loadkit &s, uchar **tables)
{
	int nodecount = 0;
	s.loadi32(nodecount);
	int *p = new int[nodecount];
	s.loadi32array(p, nodecount);
	map<int, BitMap *> pmap;                             
	BitMap * r=NULL;
	for(int i=0;i<nodecount;i++)
	{
		if(tables)
			r = new BitMap(tables);                      
		else
			r = new BitMap(); 
		r->Load(s);
		pmap[p[i]] = r;
	}
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
	string prefix(this->filename);
	string fsal2 = prefix+".1.idx";
	loadkit s1(fsal2.c_str());
	this->SAL = new InArray();
	SAL->load(s1); 
	s1.close();
	string rankl2 = prefix+".2.idx";
	loadkit s2(rankl2.c_str());
	this->RankL =new InArray();
	RankL->load(s2);
	s2.close();
	this->newd = new Dmap(); 
	string fnewD2 = prefix+".3.idx";
	loadkit s3(fnewD2.c_str());
	newd->load(s3); 
	s3.close();
	string fnewd2 = prefix+".4.idx";
	loadkit s4(fnewd2.c_str());
	this->newD = new Dmap(); 
	newD->load(s4);
	s4.close();
	string flroot2 = prefix+".5.idx";
	loadkit s5(flroot2.c_str());
	LoadWTTree(s5, lroot); 
	s5.close();
	string frroot2 = prefix+".6.idx";
	loadkit s6(frroot2.c_str());
	LoadWTTree(s6, rroot); 
	s6.close();
	string finfo = prefix+".7.idx"; 
    ifstream inf(finfo);
    inf >> this->sa0_index;
    inf >> this->D;
    inf >> this->n;
    inf >> this->alphabetsize;
	this->C =new fm_int[this->alphabetsize +1]; 
    for(int i=0; i<alphabetsize+1; i++){
        inf >> this->C[i];                
    }
    inf.close();
	mycodetable(); 
	return 0;
}
int ABS_FM::Load(loadkit &s)
{
	s.loadi64(this->n);
	s.loadi32(this->alphabetsize);
	s.loadi32(this->D);
	s.loadi32(this->r);
	s.loaddouble(this->runs);
	s.loadi64(sa0_index);
	this->C = new fm_int[alphabetsize+1];
	s.loadi64array(this->C,alphabetsize+1);
	this->code = new int[256];
	s.loadi32array(this->code,256);
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
	this->SAL = new InArray();
	this->SAL->load(s);
	Inittable();
	uchar * par[4]={Z,R,R1,R3};
	LoadWTTree(s,par);                         
	T=NULL;
	bwt=NULL;
	return 0;
}
int ABS_FM::Save(savekit &s)
{
	s.writei64(n);
	s.writei32(alphabetsize);
	s.writei32(D);
	s.writei32(r);
	s.writedouble(runs);
	s.writei64(sa0_index);
	s.writei64array(C,alphabetsize+1);
	s.writei32array(code,256);
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
	SaveWTTree(s);                                          
	return 0;
}
