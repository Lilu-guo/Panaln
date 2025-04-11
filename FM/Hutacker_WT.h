/*============================================
# Filename: Hutacker_WT.h
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#ifndef HUTACKER_FM_H
#define HUTACKER_FM_H
#include"ABS_WT.h"
typedef struct hutaNode_t
{
	fm_int freq;
	unsigned char label;
	int level;
	struct hutaNode_t * leftChild;
	struct hutaNode_t * rightChild;
}hutaNode_t;

typedef hutaNode_t * hutackerTree;

class Hutacker_FM : public ABS_FM
{
	public:
		Hutacker_FM(const char * filename,int block_size=1024,int r=16,int D=32):ABS_FM(filename,block_size,r,D){}
		Hutacker_FM():ABS_FM(filename){} //2025.3.23
		~Hutacker_FM(){};
	protected:
			int TreeCode();
	private:
			hutackerTree CreateHutackerTree();
			int HutackerNodesInit(hutaNode_t ** hutNodesPPtr);
			int FindMiniTwoNodes(hutaNode_t ** hutNodesPPtr,int nNodes,int * index1,int * index2);
			int MergeNodes(hutaNode_t ** hutaNodesPPtr,int index1,int index2,int nNodes);
			int GenerateHutackerCode(hutackerTree root);
			int SingleCharHutackerCode(hutackerTree root);
			int MultiCharHutackerCode(hutackerTree root);
			int DestroyHutackerTree(hutackerTree root);
			int GetHutackerTreeHeight(hutackerTree root);
			int ComputerCharDepth(hutackerTree root,int * charDepth);

};
#endif
