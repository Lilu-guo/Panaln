/*============================================
# Filename: Balance_WT.h
# Ver 1.0 2014-06-08
# Copyright (C) 2014 ChenLonggang (chenlonggang.love@163.com)
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: Balabced-shape wavelet tree 
=============================================*/
#ifndef _BALANCE_FM_H
#define _BALANCE_FM_H
#include"ABS_WT.h"

typedef struct balNode_t
{
	unsigned char set[CHAR_SET_SIZE];
	int setSize;
	unsigned char label;
	struct balNode_t * leftChild;
	struct balNode_t * rightChild;
	char code[CODE_MAX_LEN];
}balNode_t;
typedef balNode_t * balanceTree;

class Balance_FM : public ABS_FM
{
	public:
		Balance_FM(const char * filename,int block_size = 1024,int r=16,int D=32):ABS_FM(filename,block_size,r,D){}
		Balance_FM():ABS_FM(filename){} //2025.3.23
		~Balance_FM(){}
	protected:
		int TreeCode();
	private:
		balanceTree CreateBalanceTree();
		int GenerateBalCode(balanceTree balTree);
		int DestroyBalTree(balanceTree balTree);
		int BalanceTreeInit(balNode_t * root);
		int RecurseCreate(balanceTree root);
		int SigCharBalCode(balanceTree tree);
		int MultiCharBalCode(balanceTree root);

};
#endif


