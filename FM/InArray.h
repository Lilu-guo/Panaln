/*============================================
# Filename: InArray.h
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#ifndef _Inarray
#define _Inarray
#include"BaseClass.h"
#include"savekit.h"
#include"loadkit.h"

class InArray
{
public:
	InArray();
	InArray(fm_int data_num,int data_width);
	~InArray(void);
	u64 GetValue(fm_int index);
	// u64 operator[](u64 index){return GetValue(index);}
	void SetValue(fm_int index, u64 value);
	fm_int GetNum();
	int GetDataWidth();
	fm_int GetMemorySize();
	i64 write(savekit & s);
	i64 load(loadkit & s);

private:
	u32 * data;
    fm_int datanum;
	int datawidth;
	u64 mask;
};

class SBArray
{
public:
	SBArray();
	SBArray(fm_int sb_num,int sb_width,int _radio,int b_width);
	~SBArray(void);
	u64 GetValue(fm_int index,int issb);
	// u64 operator[](u64 index){return GetValue(index);}
	void SetValue(fm_int index, u64 value, int issb);
	fm_int GetSBNum();
	int GetSBWidth();
	fm_int GetRadio();
	int GetBWidth();
	fm_int GetMemorySize();
	i64 write(savekit & s);
	i64 load(loadkit & s);

private:
	u64 * data;
    fm_int sbnum;
	int sbwidth;
	int radio;
	int bwidth;
	u64 sbmask;
	u64 bmask;
	i64 len;
};


#endif