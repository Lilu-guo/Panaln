/*============================================
# Filename: InArray.cpp
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#include "InArray.h"
#include<string.h>
#include<stdlib.h>
#include<iostream>
#include<math.h>
using namespace std;
InArray::~InArray(void)
{
	delete [] data;
}
InArray::InArray()
{
	this->datanum = 0;
	this->datawidth = 0;
}
InArray::InArray(fm_int data_num, int data_width)
{
	if(data_num<=0||data_width<=0){
		cout<<data_num<<"  "<<data_width<<"InArray build error: data_num<=0||data_width<=0"<<endl;
		exit(0);
	}
	else
	{
		this->datanum =data_num;
		this->datawidth =data_width;
	    fm_int totlesize=datanum*datawidth;
		if(totlesize%32==0)
			totlesize=totlesize/32+1;
		else
			totlesize=(totlesize/32)+2;
		this->data =new u32[totlesize];//valgrand warns
		memset(data,0,4*totlesize);
		mask = (u64)pow(2,datawidth)-1;
	}
}

void InArray::SetValue (fm_int index, u64 v)
{

	if(index>datanum-1|| index<0)
	{
		cerr<<"InArray:index out of boundary:"<<datanum<<"  "<<index<<endl;
		exit(0);
	}
	else if(v>mask)
	{
		cerr<<"InArray:value:"<<v<<" is out of boundary"<<endl;
		exit(0);
	}
	else
	{
		u64 value=v;
		fm_int anchor=(index*datawidth)>>5;
		u64 temp1=data[anchor];
		u32 temp2=data[anchor+1];
		temp1=(temp1<<32)+temp2;
		i32 overloop=((anchor+2)<<5)-(index+1)*datawidth;
		if(overloop<0)
		{
			value=(value>>(-overloop));//35
			temp1=temp1+value;
			data[anchor+1]=(temp1&(0xffffffff));
			data[anchor]=(temp1>>32)&(0xffffffff);
			data[anchor+2]=(v<<(32+overloop))&(0xffffffff);
		}
		else
		{
			value=(value<<overloop);
			temp1=temp1+value;
			data[anchor+1]=(temp1&(0xffffffff));
			data[anchor]=(temp1>>32)&(0xffffffff);
		}
	}
}

fm_int InArray::GetNum () //获取元素个数
{
	return datanum;
}
fm_int InArray::GetMemorySize()
{
	return (datanum*datawidth)/8; //字节为单位
}


int InArray::GetDataWidth()
{
	return datawidth;
}

u64 InArray::GetValue(fm_int index)
{
	if(index>datanum-1||index<0)
	{
		cerr<<datanum-1<<"  "<<index<<" InArray:GetValue: index out of boundary"<<endl; //InArray中，输出两个值的报错
		exit(0); //直接退出了
	}

	fm_int anchor=(index*datawidth)>>5;
	u64 temp1=data[anchor];
	u32 temp2=data[anchor+1];
	u32 temp3=data[anchor+2];
	i32 overloop=((anchor+2)<<5)-(index+1)*datawidth;
	temp1=(temp1<<32)+temp2;
	if(overloop<0)
	{
	   temp1 = (temp1<<(-overloop))+(temp3>>(32+overloop));
		return temp1&mask;
	}
	else
	{
		return (temp1>>overloop)&mask;
	}
}
i64 InArray::write(savekit & s)
{

	s.writei64(datanum);
	s.writei32(datawidth);
	u64 len=(datanum*datawidth);
	if(len%32==0)
		len=len/32+1;
	else
		len=len/32+2;
	s.writeu64(len);
	s.writeu32array(data,len);
	return 1;
}
i64 InArray::load(loadkit & s)
{
	if(s.loadi64(datanum)!=1){
		cout<<"InArray datanum read error"<<endl;
		exit(0);
	}
	if(s.loadi32(datawidth)!=1){
		cout<<"InArray datawidth read error"<<endl;
		exit(0);		
	}
	u64 len=0;
	if(s.loadu64(len)!=1){
		cout<<"InArray len read error"<<endl;
		exit(0);
	}
	data=new u32[len];
	i32 f=s.loadu32array(data,len);
	if(f!=len){
		cout<<len<<" "<<f<<"InArray data read error"<<endl;
		exit(0);
	}
	
	mask=(u64)pow(2,datawidth)-1;
	return 1;
}

SBArray::~SBArray(void)
{
	delete [] data;
}
SBArray::SBArray()
{
	this->data = NULL;
	this->sbnum = 0;
	this->sbwidth = 0;
	this->radio = 0;
	this->bwidth = 0;
	this->sbmask = 0;
	this->bmask = 0;
}

SBArray::SBArray(fm_int sb_num,int sb_width,int _radio,int b_width)
{
	if(sb_num<=0||sb_width<=0||_radio<=0||b_width<=0){
		cout<<sb_num<<"  "<<sb_width<<"  "<<_radio<<"  "<<b_width<<"SBArray build error: data_num<=0||data_width<=0"<<endl;
		exit(0);
	}
	else
	{
		this->sbnum =sb_num;
		this->sbwidth =sb_width;
		this->radio = _radio;
		this->bwidth = b_width; 
	    fm_int totlesize=sbnum*(sbwidth+radio*bwidth);
		if(totlesize%64==0)
			totlesize=totlesize/64+1;
		else
			totlesize=(totlesize/64)+2;
		len = totlesize;
		this->data =new u64[totlesize];//valgrand warns
		memset(data,0,8*totlesize);
		sbmask = (u64)pow(2,sbwidth)-1;
		bmask = (u64)pow(2,bwidth)-1;
	}
}

void SBArray::SetValue(fm_int index, u64 v,int issb)
{
	if(issb)
	{
		if(index > sbnum - 1 || index<0)
		{
			cerr<<"SBArray:SBindex out of boundary:"<<sbnum<<"  "<<index<<endl;
			exit(0);
		}
		else if(v > sbmask)
		{
			cerr<<"SBArray:SBvalue:"<<v<<" is out of boundary"<<endl;
			exit(0);
		}
		else 
		{
			u64 value=v;
			fm_int anchor=(index*(sbwidth + radio * bwidth))>>6;
			int overloop = (index*(sbwidth + radio * bwidth)) % 64;
			if(overloop + sbwidth <= 64)
				data[anchor] |= (value << (64 - overloop - sbwidth)); 
			else
			{
				int right = overloop + sbwidth - 64;
				data[anchor] |= (value >> right);
				data[anchor + 1] |= (value << (64 - right));
			}
				
		}
	}
	else
	{
		if(index > sbnum * radio - 1 || index<0)
		{
			cerr<<"BArray:Bindex out of boundary:"<<sbnum * radio<<"  "<<index<<endl;
			exit(0);
		}
		else if(v > bmask)
		{
			cerr<<"BArray:Bvalue:"<<v<<" is out of boundary"<<endl;
			exit(0);
		}
		else 
		{
			u64 value=v;
			fm_int tmp = index / radio * (sbwidth + radio * bwidth) + sbwidth + (index % radio * bwidth);
			fm_int anchor = tmp >> 6;
			int overloop = tmp % 64;
			if(overloop + bwidth <= 64)
				data[anchor] |= (value << (64 - overloop - bwidth)); 
			else
			{
				int right = overloop + bwidth - 64;
				data[anchor] |= (value >> right);
				data[anchor + 1] |= (value << (64 - right));
			}			
		}
	}
}
u64 SBArray::GetValue(fm_int index,int issb)
{
	if(issb)
	{
		if(index>sbnum-1||index<0)
		{
			cerr<<sbnum-1<<"  "<<index<<" SBArray:SB GetValue: index out of boundary"<<endl;
			exit(0);
		}
		u64 value = 0;
		fm_int anchor=(index*(sbwidth + radio * bwidth))>>6;
		int overloop = (index*(sbwidth + radio * bwidth)) % 64;
		if(overloop + sbwidth <= 64)
		{
			u64 tmp = data[anchor] >> (64 - overloop - sbwidth);
			value |= (tmp & ((1ull << sbwidth) - 1));
			return value;
		}
		else
		{
			int right = overloop + sbwidth - 64;
			value |= ((data[anchor] & ((1ull << (sbwidth - right)) - 1)) << right);
			value |= (data[anchor + 1] >> (64 - right));
			return value;
		}
	}
	else
	{
		if(index>sbnum * radio-1||index<0)
		{
			cerr<<sbnum * radio-1<<"  "<<index<<" SBArray:B GetValue: index out of boundary"<<endl;
			exit(0);
		}
		u64 value = 0;
		fm_int ttmp = index / radio * (sbwidth + radio * bwidth) + sbwidth + (index % radio * bwidth);
		fm_int anchor = ttmp >> 6;
		int overloop = ttmp % 64;
		if(overloop + bwidth <= 64)
		{
			u64 tmp = data[anchor] >> (64 - overloop - bwidth);
			value |= (tmp & ((1ull << bwidth) - 1));
			return value;
		}
		else
		{
			int right = overloop + bwidth - 64;
			value |= ((data[anchor] & ((1ull << (bwidth - right)) - 1)) << right);
			value |= (data[anchor + 1] >> (64 - right));
			return value;
		}
	}
}
fm_int SBArray::GetSBNum ()
{
	return sbnum;
}
int SBArray::GetSBWidth()
{
	return sbwidth;
}
fm_int SBArray::GetRadio()
{
	return radio;
}
int SBArray::GetBWidth()
{
	return bwidth;
}
fm_int SBArray:: GetMemorySize()
{
	return len * sizeof(u64);
}

i64 SBArray::write(savekit & s)
{

	s.writei64(sbnum);
	s.writei32(sbwidth);
	s.writei32(radio);
	s.writei32(bwidth);
	s.writei64(len);
	s.writeu64array(data,len);
	return 1;
}
i64 SBArray::load(loadkit & s)
{
	if(s.loadi64(sbnum)!=1){
		cout<<"SBArray sbnum read error"<<endl;
		exit(0);
	}
	if(s.loadi32(sbwidth)!=1){
		cout<<"SBArray sbwidth read error"<<endl;
		exit(0);		
	}
	if(s.loadi32(radio)!=1){
		cout<<"SBArray radio read error"<<endl;
		exit(0);		
	}
	if(s.loadi32(bwidth)!=1){
		cout<<"BArray bwidth read error"<<endl;
		exit(0);		
	}
	if(s.loadi64(len)!=1){
		cout<<"SBArray len read error"<<endl;
		exit(0);		
	}
	data=new u64[len];
	i32 f=s.loadu64array(data,len);
	if(f!=len){
		cout<<len<<" "<<f<<"InArray data read error"<<endl;
		exit(0);
	}
	
	sbmask=(u64)pow(2,sbwidth)-1;
	bmask=(u64)pow(2,bwidth)-1;
	return 1;
}
