/*============================================
# Filename: WT_Handle.cpp
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#include"WT_Handle.h"
#include"Huffman_WT.h"
#include"Balance_WT.h"
#include"Hutacker_WT.h"

WT_Handle::WT_Handle(char * filename):fm(new ABS_FM(filename)),u(){} //读索引，比对时用

WT_Handle::WT_Handle(const char * filename,int block_size,int r,int D,int shape) //构建索引用
{
	if(block_size<=0 || shape<0 || shape >2)
	{
		cout<<"WT_Handle::WT_handle error parmater"<<endl;
		exit(0);
	}

	switch(shape)
	{
		case 0: fm =new Hutacker_FM(filename,block_size,r,D);break;
		case 1: fm =new Huffman_FM(filename,block_size,r,D);break; //间接通过ABS_FM的getfile来，读入fastq
		case 2: fm =new Balance_FM(filename,block_size,r,D);break;
		default: fm=new Hutacker_FM(filename,block_size,r,D);break;
	}
	fm->BuildTree(); //建小波树（由ABS_FM去做）；读文件getFile、表初始化initTable在ABS_FM构造方法中完成
}

WT_Handle::WT_Handle(const WT_Handle &h):fm(h.fm),u(h.u){}

WT_Handle & WT_Handle::WT_Handle:: operator =(const WT_Handle & h)
{
	if(u.reattach(h.u))
	{
		delete fm;
	}
	fm = h.fm;
	return * this;
}

WT_Handle::~WT_Handle()
{
	if(u.only())
		delete fm;
}


