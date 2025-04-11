/*============================================
# Filename: Dmap.h
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#include<string.h>
#include"loadkit.h"
#include"savekit.h"
#include"InArray.h"
#include"BaseClass.h"
#include<math.h>
#include<iostream>
#include <bitset>
using namespace std;
//typedef long long fm_int;
class Dmap
{
  public:
	Dmap();
	Dmap(fm_int data_num);
    fm_int GetMemorySize();
	void write(savekit &s);
	void load(loadkit &s);
	~Dmap(void);
	void SetValue(fm_int index, u64 v);
	void constructrank(int Blength, int SBlength);
	fm_int select(fm_int i);
	fm_int BS(InArray *B,fm_int l,fm_int r,fm_int i);
	void test();
	u64 GetValue(fm_int index);
	int getD(fm_int i);
	fm_int rank(fm_int i);
    fm_int rank2(fm_int i);
  private:
	u32 *data;
	fm_int datanum;
	int rankflag;
	InArray *Md;//标记当前小块中是否有1
	InArray *Sd;//标记每个1在当前块中的位置
	InArray *SBd;//当前超块前有多少个1
	InArray *Bd;//所在超块中，当前小块前的小块中有多少个1
	int Blength;
	int SBlength;
};