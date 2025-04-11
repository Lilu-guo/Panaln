/*============================================
# Filename: Dmap.cpp
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#include"newDmap.h"

fm_int Dmap:: GetMemorySize()
{	
	if(Md==NULL)
		return 0;
	else
		return Md->GetMemorySize() + Sd->GetMemorySize() + SBd->GetMemorySize() + Bd->GetMemorySize(); //四个结构，字节为单位
};

int ablog(fm_int x)
{
	int ans = 0;
	while (x > 0)
	{
		ans++;
		x = x >> 1;
	}
	return ans;
}

void Dmap:: write(savekit &s){     //保存
	s.writei64(datanum);
	s.writei32(rankflag);
	s.writei32(Blength);
	s.writei32(SBlength);
	Md->write(s);
	Sd->write(s);
	SBd->write(s);
	Bd->write(s);	
};

void Dmap::load(loadkit &s){       //读入
	s.loadi64(this->datanum);
	s.loadi32(this->rankflag);
	s.loadi32(this->Blength);
	s.loadi32(this->SBlength);
	this->Md  = new InArray(); //D1flg，块中是否有1
	this->Sd  = new InArray(); //D1pos，1在块中位置
	this->SBd = new InArray(); //超块rank
	this->Bd  = new InArray(); //小块rank
	this->Md->load(s);
	this->Sd->load(s);
	this->SBd->load(s);
	this->Bd->load(s);
};

fm_int Dmap::rank(fm_int i)                            //i一定要是标记1的pos（是通过查D1pos与i，判断停止的）
{
	fm_int Mdnum = i >> Blength;                       //定位小块号
	fm_int SBnum = Mdnum >> SBlength;                  //定位超块号
	fm_int Sbdbefore = SBd->GetValue(SBnum);           //超块rank
	fm_int Sdnowi = Sbdbefore + Bd->GetValue(Mdnum);   //超块rank + 小块rank
	while (1)
	{
		fm_int Sdnow = Sd->GetValue(Sdnowi);           //查D1pos（从超块rank + 小块rank，对应的pos下标处开始）
		if (((Mdnum << Blength) + Sdnow) == i)         //小块号*小块长 + 块内相对位置（递增）
			break;                                     //碰到i，返回rank（sdnowi）
		Sdnowi++;
	}
	return Sdnowi + 1;
}

fm_int Dmap::rank2(fm_int i)
{
    u64 Mdnum = i >> Blength;                         //块号
    u64 low = i & ((1 << Blength) - 1);               //留低位，小块中要扫的位置个数（含0和1的位置）
    u64 SBnum = Mdnum >> SBlength;                    //超块号
	// cout<<"i:"<<i<<"  SBnum:"<<SBnum<<"  Mdnum:"<<Mdnum<<endl; //错误调试用20230201
	// cout<<"rank2...1"<<endl;
    u64 Sbdbefore = SBd->GetValue(SBnum); //这里出问题？？？
	// cout<<"rank2...2"<<endl;
    u64 Sbdafter = SBd->GetValue(SBnum+1);            //下一个超快rank
	// cout<<"rank2...3"<<endl;
    u64 Sdnowi = Sbdbefore + Bd->GetValue(Mdnum);     //超块rank + 小块rank

	// cout<<"rank2...4"<<endl;
    if(Md->GetValue(Mdnum) == 0)                      //1. 小块中无1，返回（超快rank + 小块rank）
        return Sdnowi;

    if((Mdnum+1)%16 == 0)                             //为大块中的最后一个小块
    {
        while(1) //块内，对每个1的位置，逐个判断，直至碰到i或跳过i.
        {
            u64 Sdnow = Sd->GetValue(Sdnowi);         //块内首个1的相对位置（用Sdnowi更新取下一个）
            if(low == Sdnow)                          //2.1. 碰到i（等于）
            {
                Sdnowi++;                             //记录rank数
                break;
            }
            else if(low > Sdnow)                      //2.2. 还没到i（小于），low为待查长度
            {
                Sdnowi++;
                if(Sdnowi >= Sbdafter)                //超出了小块
                    break;
            }
            else                                      //2.3. 跳过了i（大于）
                break;
        }
    }
	else{
        u64 Sdagap = Bd->GetValue(Mdnum+1) - Bd->GetValue(Mdnum);    //小块中有几个1（通过小块rank计算）
        i64 k = 0;

        while(1)
        {
            u64 Sdnow = Sd->GetValue(Sdnowi);   //从D1pos中，依次取出1的位置(块内相对位置)
            if(low == Sdnow)
            {
                Sdnowi++;
                break;
            }
            else if(low > Sdnow)                                   //2. 还没到i
            {
                Sdnowi++;
                k++;
                if(k >= Sdagap)                                    //超出了小块
                    break;
            }
            else
                break;
        }
    }
    return Sdnowi;
}

u64 Dmap::GetValue(fm_int index)                                   //直接从位串取0/1
{
	if (index > datanum - 1 || index < 0)
	{
		cerr << "InArray:GetValue: index out of boundary" << endl;
		exit(0);
	}
	fm_int anchor = index >> 5;
	int overloop=index % 32;
	return (data[anchor]>>(31-overloop))&1;
}

int Dmap::getD(fm_int i)
{
	if(i<0||i>=datanum){
		cout<<"Dmap: Index out of boundary"<<i<<"  "<<datanum<<endl;
		exit(0);
	}
	fm_int Mdnum = i >> Blength;                                          //定位小块号
	if (Md->GetValue(Mdnum) == 0)                                         //查D1flg，看块中是否有1
		return 0;
	else
	{
		fm_int SBnum = Mdnum >> SBlength;                                 //定位超块号
		fm_int Sbdbefore = SBd->GetValue(SBnum);                          //得超块记录的rank
		fm_int Sdnowi = Sbdbefore + Bd->GetValue(Mdnum);                  //超块rank + 小块rank
		fm_int Sdnow = Sd->GetValue(Sdnowi);                              //查D1pos，看块内相对位置
		if (((Mdnum << Blength) + Sdnow) > i)                             //1. 块内的首个1，出现在i之后
			return 0;
		else if (((Mdnum << Blength) + Sdnow) == i)                       //2. 块内的首个1，刚好在i处
			return 1;
		else                                                              //3. 块内的首个1，出现在i之前
		{
			fm_int temnum1 = 0;                                           //需要在D1pos中扫的位置个数
			if (Mdnum != Md->GetNum() - 1)                                //没碰到最后一个小块
			{
				fm_int nextBd = Bd->GetValue(Mdnum + 1);                  //取下个小块的rank
				if (nextBd == 0)
					temnum1 = SBd->GetValue(SBnum + 1) - Sbdbefore - Bd->GetValue(Mdnum); //下个超块的rank - 超块的rank - 小块的rank
				else
					temnum1 = nextBd - Bd->GetValue(Mdnum);
			}
			else
			{
				temnum1 = Sd->GetNum() - Sdnowi;
			}
			for (fm_int j = 1; j < temnum1; j++)                           //逐个扫，直到碰到或大于i
			{
				Sdnow = Sd->GetValue(Sdnowi + j);                          //从小块起始rank对应的pos，向后
				if (((Mdnum << Blength) + Sdnow) > i)                      //超过i
					return 0;
				else if (((Mdnum << Blength) + Sdnow) == i)                //碰到i
					return 1;
			}
			return 0;
		}
	}
}

Dmap::~Dmap(void)
{
	delete[] data;
	if (rankflag)
	{
		delete Md;
		delete Sd;
		delete SBd;
		delete Bd;
	}
}

Dmap::Dmap()
{
	Md=NULL;
	Sd=NULL;
	SBd=NULL;
	Bd=NULL;
}

Dmap::Dmap(fm_int len)
{
	// cout<<ablog(4)<<endl;
	if (len < 0)
		cout << "Dmap initialization variable error" << endl;
	else
	{
		this->datanum = len;
		u64 totlesize = datanum;
		if (totlesize % 32 == 0)
			totlesize = totlesize / 32 + 1;
		else
			totlesize = (totlesize / 32) + 2;
		this->data = new u32[totlesize];
		memset(data, 0, 4 * totlesize);
	}
}

void Dmap::constructrank(int Blength,int SBlength)
{
	rankflag=1;
	this->Blength=ablog(Blength)-1;
	this->SBlength=ablog(SBlength)-1;
	i64 totlesize=datanum;                                          //0/1串长度
	if(totlesize%32==0)
		totlesize=totlesize/32;
	else
		totlesize=(totlesize/32)+1;
	fm_int sum1=0;
	for(int i=0;i<totlesize;i++)
		sum1+=(u32)__builtin_popcount(data[i]);                     //统计1的个数（popc每次32位）
	fm_int Mdm;
	if(datanum%Blength==0)Mdm=datanum/Blength+1;
	else Mdm=datanum/Blength+2;
	Md=new InArray(Mdm,1);                                          //D1flg，块中是否有1（bitvector位向量）
	Sd=new InArray(sum1,ablog(Blength));                            //D1pos，1在块中位置（个数，位宽定长）
	SBd=new InArray((datanum/(Blength*SBlength)+2),ablog(sum1));    //超块rank
	Bd=new InArray(Mdm,ablog(Blength*SBlength+1));                  //小块rank
	fm_int Sdcount=0;
	int tem;
	fm_int now=-1;
	for(u64 i=0;i<Mdm;i++)
	{
		int flag=0;
		for(int j=0;j<Blength;j++)
		{
			now++;
			if(now>=datanum) break;
			tem=GetValue(now);
			if(tem==1)
			{
				flag=1;
				Sd->SetValue(Sdcount,j);
				Sdcount++;
			}
		}
		if(flag!=0)Md->SetValue(i,1);
		else Md->SetValue(i,0);
	}
	fm_int SBnum;
	if(Mdm%SBlength==0)SBnum=Mdm/SBlength+1;
	else SBnum=Mdm/SBlength+2;
	fm_int SBsum1=0;
    int flag=0;
	now=-1;
	for(fm_int i=0;i<SBnum;i++)
	{
		SBd->SetValue(i,SBsum1);
		fm_int Bdsum1=0;
        if(flag==1){
            break;
        }
		for(int j=0;j<SBlength;j++)
		{
			Bd->SetValue(i*SBlength+j,Bdsum1);
            if(flag==1){
                break;
            }
			for(int k=0;k<Blength;k++)
			{
				now++;
				if(now>=datanum){
                    flag=1;
                    break;
                }
				else{
					int tem=GetValue(now);
					if(tem==1)
					{
						SBsum1++;
						Bdsum1++;
					}
				}
			}
		}

	}
	// cout<<"Md:      "<<Md->GetMemorySize()/1024/1024<<" MB"<<endl;
	// cout<<"Sd:      "<<Sd->GetMemorySize()/1024/1024<<" MB"<<endl;
	// cout<<"SBd:     "<<SBd->GetMemorySize()/1024/1024<<" MB"<<endl;
	// cout<<"Bd:      "<<Bd->GetMemorySize()/1024/1024<<" MB"<<endl;
	// cout<<"Sd...num: "<<Sd->GetNum()<<"  width: "<<Sd->GetDataWidth()<<endl;
	// cout<<"rank test..."<<endl;
	// test();   //rank测试
}

void Dmap::SetValue(fm_int index, u64 v)
{
	if (index > datanum - 1 || index < 0)
	{
		cerr << "InArray:index out of boundary" << endl;
		exit(0);
	}
	else
	{
		fm_int anchor=index/32;
		int overloop=index%32;
		data[anchor]=data[anchor]|(v<<(31-overloop));
	}
}

fm_int Dmap::select(fm_int i)//查找第i个1所在的位置
{
	if(i<=0||i>Sd->GetNum()){
		cout<<"select out of boundary"<<endl;
		return 0;
	}
	int b=1<<Blength;
	int sb=1<<SBlength;
	fm_int Mdnum;
	if(datanum%b==0)Mdnum=datanum/b+1;
	else Mdnum=datanum/b+2;
	fm_int SBnum;
	if(datanum%(sb*b)==0)SBnum=datanum/(sb*b)+1;
	else SBnum=datanum/(sb*b)+2;
	fm_int Sbpos=BS(SBd,0,SBnum-1,i);
	// cout<<"Sbpos: "<<Sbpos<<endl;
	fm_int l=Sbpos*sb;
	fm_int r=(Sbpos+1)*sb-1;
	if(r>=Mdnum)
		r=Mdnum-1;
	fm_int k=i-SBd->GetValue(Sbpos);
	fm_int Bpos=BS(Bd,l,r,k);
	// cout<<"bpos: "<<Bpos<<endl;
	fm_int finalpos=Bpos*b+Sd->GetValue(i-1);
	return finalpos;
}
fm_int Dmap::BS(InArray *B,fm_int l,fm_int r,fm_int i)//二分查询第i个1所在的块
{
	if(B->GetValue(l)==i) return l-1;
	if(B->GetValue(r)==i) {
		while(B->GetValue(r)==i){
			r--;
		}
		return r;
	}
	if(B->GetValue(r)<i) return r;
	while(l<=r){
		auto midpos=l+((r-l)>>1);
		// cout<<l<<" ## "<<midpos<<" ## "<<r<<endl;
		if(midpos==l){
			if(B->GetValue(l)>=i){
				return l-1;
			}
			else{
				return l;
			}
		}
		auto mid=B->GetValue(midpos);
		if(mid<i) l=midpos+1;
		else if(mid>=i) r=midpos;
	}
}

void Dmap::test()
{
	fm_int miss=0;
    fm_int num=0;
    for(fm_int i=0;i<datanum;i++){
        if(getD(i)==1){
            num++;                                                           //若是1，记录个数
            // if(rank(i)!=num){                                                //与rank函数比较
			if(rank2(i)!=num){
                // cout<<"rank error:"<<num<<"  "<<rank(i)<<"  "<<i<<endl;
				miss++;
            }
        }
    }
	cout<<"rank miss "<<miss<<endl;
    cout<<"#####################rank test2 finished################"<<endl;

	/* fm_int n=Sd->GetNum();                                                   //1总的出现次数，查D1pos
	miss=0;
    if(n!=num){                                                                 //检查D1pos是否出错（只看个数）
        cout<<"one num error"<<endl;
    }
	fm_int t=0;
	cout<<"number of 1:  "<<n<<endl;
	for(fm_int i=1;i<=n;i++){
		t=select(i);
		if(rank(t)!=i){                                                         //用rank检查select
			// cout<<"select error"<<i<<"  "<<rank(t)<<endl;
			miss++;
		}
	}

	cout<<"error:"<<miss<<endl;
    cout<<"####################select test finished###############"<<endl; */

	// miss=0;
	// for(fm_int i=0;i<datanum;i++){
	// 	int a=getD(i);                                                         //利用D结构，取0/1
	// 	int b=GetValue(i);                                                     //直接从位串，拿0/1，比较
	// 	if(a!=b){
	// 		cout<<a<<"  "<<b<<"  "<<i<<endl;
	// 		miss++;
	// 	}
	// }
	// cout<<"getD error:"<<miss<<endl;

}