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
#include"Dmap.h"
fm_int Dmap:: GetMemorySize()
{
	if(Md==NULL)
		return 0;
	else
		return Md->GetMemorySize() + Sd->GetMemorySize() + SBd->GetMemorySize() + Bd->GetMemorySize();
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
void Dmap:: write(savekit &s){
	s.writei64(datanum);
	s.writei32(rankflag);
	s.writei32(Blength);
	s.writei32(SBlength);
	Md->write(s);
	Sd->write(s);
	SBd->write(s);
	Bd->write(s);	
};
void Dmap::load(loadkit &s){
	s.loadi64(this->datanum);
	s.loadi32(this->rankflag);
	s.loadi32(this->Blength);
	s.loadi32(this->SBlength);
	this->Md  = new InArray();
	this->Sd  = new InArray();
	this->SBd = new InArray();
	this->Bd  = new InArray();
	this->Md->load(s);
	this->Sd->load(s);
	this->SBd->load(s);
	this->Bd->load(s);
};
fm_int Dmap::rank(fm_int i)
{
	fm_int Mdnum = i >> Blength;
	fm_int SBnum = Mdnum >> SBlength;
	fm_int Sbdbefore = SBd->GetValue(SBnum);
	fm_int Sdnowi = Sbdbefore + Bd->GetValue(Mdnum);
	while (1)
	{
		fm_int Sdnow = Sd->GetValue(Sdnowi);
		if (((Mdnum << Blength) + Sdnow) == i)
			break;
		Sdnowi++;
	}
	return Sdnowi + 1;
}
u64 Dmap::GetValue(fm_int index)
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
	fm_int Mdnum = i >> Blength;
	if (Md->GetValue(Mdnum) == 0)
		return 0;
	else
	{
		fm_int SBnum = Mdnum >> SBlength;
		fm_int Sbdbefore = SBd->GetValue(SBnum);
		fm_int Sdnowi = Sbdbefore + Bd->GetValue(Mdnum);
		fm_int Sdnow = Sd->GetValue(Sdnowi);
		if (((Mdnum << Blength) + Sdnow) > i)
			return 0;
		else if (((Mdnum << Blength) + Sdnow) == i)
			return 1;
		else
		{
			fm_int temnum1 = 0;
			if (Mdnum != Md->GetNum() - 1)
			{
				fm_int nextBd = Bd->GetValue(Mdnum + 1);
				if (nextBd == 0)
					temnum1 = SBd->GetValue(SBnum + 1) - Sbdbefore - Bd->GetValue(Mdnum);
				else
					temnum1 = nextBd - Bd->GetValue(Mdnum);
			}
			else
			{
				temnum1 = Sd->GetNum() - Sdnowi;
			}
			for (fm_int j = 1; j < temnum1; j++)
			{
				Sdnow = Sd->GetValue(Sdnowi + j);
				if (((Mdnum << Blength) + Sdnow) > i)
					return 0;
				else if (((Mdnum << Blength) + Sdnow) == i)
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
	cout<<ablog(4)<<endl;
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
	i64 totlesize=datanum;
	if(totlesize%32==0)
		totlesize=totlesize/32;
	else
		totlesize=(totlesize/32)+1;
	fm_int sum1=0;
	for(int i=0;i<totlesize;i++)
		sum1+=(u32)__builtin_popcount(data[i]);
	fm_int Mdm;
	if(datanum%Blength==0)Mdm=datanum/Blength+1;
	else Mdm=datanum/Blength+2;
	Md=new InArray(Mdm,1);
	Sd=new InArray(sum1,ablog(Blength));
	SBd=new InArray((datanum/(Blength*SBlength)+2),ablog(sum1));
	Bd=new InArray(Mdm,ablog(Blength*SBlength+1));
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
	// test();
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
            num++;
            if(rank(i)!=num){
                // cout<<"rank error:"<<num<<"  "<<rank(i)<<"  "<<i<<endl;
				miss++;
            }
        }
    }
	cout<<"rank miss "<<miss<<endl;
    cout<<"#####################rank test finished################"<<endl;
	fm_int n=Sd->GetNum();
	miss=0;
    if(n!=num){
        cout<<"one num error"<<endl;
    }
	fm_int t=0;
	cout<<"number of 1:  "<<n<<endl;
	for(fm_int i=1;i<=n;i++){
		t=select(i);
		if(rank(t)!=i){
			// cout<<"select error"<<i<<"  "<<rank(t)<<endl;
			miss++;
		}
	}

	cout<<"error:"<<miss<<endl;
    cout<<"####################select test finished###############"<<endl;
	// miss=0;
	// for(fm_int i=0;i<datanum;i++){
	// 	int a=getD(i);
	// 	int b=GetValue(i);
	// 	if(a!=b){
	// 		cout<<a<<"  "<<b<<"  "<<i<<endl;
	// 		miss++;
	// 	}
	// }
	// cout<<"getD error:"<<miss<<endl;

}