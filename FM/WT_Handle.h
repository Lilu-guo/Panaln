 /*============================================
# Filename: WT_Handle.h
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
  A handle class for pointer fm(ABS_FM *),ABS_FM will
  actual points to a specific subclass.
=============================================*/
#ifndef WTHANDLE_H
#define WTHANDLE_H
#include"UseCount.h"
#include"ABS_WT.h"
class WT_Handle
{
	private:
		ABS_FM * fm;                                                               //FM -> WT_Handle -> ABS_FM（都通过列表初始化，调用构造函数）
		UseCount u;
	public:
		WT_Handle(char * filename); //保存参考序列名，为索引idx3动态命令用，2025.3.23
		WT_Handle(const char * filename,int block_size=256,int r=16,int D=32,int shape=1);
		WT_Handle(const WT_Handle &);
		WT_Handle & operator = (const WT_Handle & );
		~WT_Handle();
		
		void Counting(const char * pattern,fm_int &num) { fm->Counting(pattern,num); };
		fm_int * Locating(const char * pattern,fm_int &num){ return fm->Locatingnew(pattern,num); };
		unsigned char *Extracting(fm_int pos,fm_int len){ return fm->Extracting(pos,len);};
		int Load(loadkit & s) { return fm->Load(s);};
		int loadIdx(){return fm->loadIdx();};                                                //add gll
		int Save(savekit & s){ return fm->Save(s);};

		fm_int Fm_Lookup(fm_int pos){return fm->FM_Lookup(pos);};
		void Fmocc(fm_int l, fm_int r,unsigned char c,fm_int *Left,fm_int * Right){fm->Fmocc(l,r,c,Left,Right);};    //对应fm.h，为bwa新增的方法
		void Fm16occ(fm_int l, fm_int r,fm_int *Left,fm_int * Right){fm->Fm16occ(l,r,Left,Right);};

		void Occ1p(unsigned char c, fm_int pos1, fm_int pos2, fm_int &occ1, fm_int &occ2){fm->Occ1p(c,pos1,pos2,occ1,occ2);};      //add gll
		void Occ8p(unsigned char c, fm_int pos1, fm_int pos2, fm_int *occ81, fm_int *occ82){fm->Occ8p(c,pos1,pos2,occ81,occ82);};  //add gll
		void Occ16p(fm_int pos1, fm_int pos2, fm_int *occ161, fm_int *occ162){fm->Occ16p(pos1,pos2,occ161,occ162);};               //add gll
		
		fm_int GetN(){ return fm->GetN();}
		int GetAlphabetsize(){return fm->GetAlphabetsize();}
		fm_int SizeInByte(){ return fm->SizeInByte();};
		fm_int SizeInByte_count() { return fm->SizeInByte_count();};
		fm_int SizeInByte_locate() { return fm->SizeInByte_locate();};
		fm_int SizeInByte_extract() { return fm->SizeInByte_extract();};
		fm_int SizeinByte_sal(){return fm->SizeInByte_SAL();};
		fm_int SizeinByte_rankl(){return fm->SizeInByte_RankL();};
		fm_int SizeinByte_newD(){return fm->SizeInByte_newD();};
		fm_int TreePlainBlock(){return fm->PlainBlock();};
		fm_int TreeRLGBlock(){return fm->RLGBlock();};
		fm_int TreeALL01Block(){return fm->ALL01Block();};
		fm_int TreeEF01Block(){return fm->EF01Block();};
		fm_int TreePlainSize(){return fm->PlainSize();};
		fm_int TreeRLGSize(){return fm->RLGSize();};
		fm_int TreeEF01Size(){return fm->EF01Size();};
		fm_int MemSize(){return fm->MemSize();};
		fm_int NodeBlockSize(){return fm->NodeBlockSize();};
		fm_int CodeStyleSize(){return fm->CodeStyleSize();};
		double GetRuns(){return fm->GetRuns();};
};
#endif
