/*============================================
# Filename: loadkit.cpp
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#include"loadkit.h"
void loadkit::close()
{
	if(r!=NULL)
		// fclose(r);
		r->close();
	r=NULL;
}

loadkit::~loadkit()
{
	if(r!=NULL)
		// fclose(r);
		r->close();
}

loadkit::loadkit(const char * file)
{
	// this->r=fopen(file,"rb");
	this->r=new ifstream(file,ifstream::binary);
	if(r==NULL)
	{
		cout<<"Fopen error"<<endl;
		exit(0);
	}
}

i32 loadkit::loadi64(i64 & value)
{
	// i32 num=fread(&value,sizeof(i64),1,r);
	r->read((char*)&value,sizeof(i64));
	return r->gcount()/sizeof(i64);
}
i32 loadkit::loadu64(u64 & value)
{
	// i32 num=fread(&value,sizeof(u64),1,r);
	r->read((char*)&value,sizeof(u64));
	return r->gcount()/sizeof(u64);
}
i32 loadkit::loadi32(i32 & value)
{
	// i32 num=fread(&value,sizeof(i32),1,r);
	// return num;
	r->read((char*)&value,sizeof(i32));
	return r->gcount()/sizeof(i32);
}
i32 loadkit::loadu32(u32 & value)
{
	// i32 num=fread(&value,sizeof(u32),1,r);
	// return num;
	r->read((char*)&value,sizeof(u32));
	return r->gcount()/sizeof(u32);
}
i32 loadkit::loadi16(i16 & value)
{
	// i32 num=fread(&value,sizeof(i16),1,r);
	// return num;
	r->read((char*)&value,sizeof(i16));
	return r->gcount()/sizeof(i16);
}
i32 loadkit::loadu16(u16 & value)
{
	// i32 num=fread(&value,sizeof(u16),1,r);
	// return num;
	r->read((char*)&value,sizeof(u16));
	return r->gcount()/sizeof(u16);
}
i32 loadkit::loadu8(u8 & value)
{
	// i32 num = fread(&value,sizeof(u8),1,r);
	// return num;
	r->read((char*)&value,sizeof(u8));
	return r->gcount()/sizeof(u8);
}
i32 loadkit::loaddouble(double & value){
	r->read((char*)&value,sizeof(double));
	return r->gcount()/sizeof(double);
}
i32 loadkit::loadi64array(i64 * value,i32 len)
{
	// i32 num=fread(value,sizeof(i64),len,r);
	// return num;
	r->read((char*)value,len*sizeof(i64));
	return r->gcount()/sizeof(i64);
}
i32 loadkit::loadu64array(u64 * value,i32 len)
{
	// i32 num=fread(value,sizeof(u64),len,r);
	// return num;
	r->read((char*)value,len*sizeof(u64));
	return r->gcount()/sizeof(u64);
}
i32 loadkit::loadi32array(i32 * value,i32 len)
{
	// i32 num=fread(value,sizeof(i32),len,r);
	// return num;
	r->read((char*)value,len*sizeof(i32));
	return r->gcount()/sizeof(i32);
}
i32 loadkit::loadu32array(u32 * value,i32 len)
{
	// i32 num=fread(value,sizeof(u32),len,r);
	// return num;
	r->read((char*)value,len*sizeof(u32));
	return r->gcount()/sizeof(u32);
}
i32 loadkit::loadi16array(i16 * value,i32 len)
{
	// i32 num=fread(value,sizeof(i16),len,r);
	// return num;
	r->read((char*)value,len*sizeof(i16));
	return r->gcount()/sizeof(i16);
}
i32 loadkit::loadu16array(u16 * value,i32 len)
{
	// i32 num=fread(value,sizeof(u16),len,r);
	// return num;
	r->read((char*)value,len*sizeof(u16));
	return r->gcount()/sizeof(u16);
}

i32 loadkit::loadu8array(u8 * value,i32 len)
{
	// i32 num = fread(value,sizeof(u8),len,r);
	// return num;
	r->read((char*)value,len*sizeof(u8));
	return r->gcount()/sizeof(u8);
}












