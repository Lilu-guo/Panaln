/*============================================
# Filename: Huffman_WT.h
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#ifndef _Huff_FM_H
#define _Huff_FM_H
#include "ABS_WT.h"
typedef struct huffNode_t{
	fm_int freq;
	unsigned char label;
	char code[CODE_MAX_LEN];
	struct huffNode_t * leftChild;
	struct huffNode_t * rightChild;
}huffNode_t;
typedef huffNode_t * huffmanTree;


class Huffman_FM : public ABS_FM
{
	public:
		Huffman_FM(const char * filename,int block_size =256,int r=16,int D=32):ABS_FM(filename,block_size,r,D){}
		Huffman_FM():ABS_FM(filename){} ///2025.3.23
		~Huffman_FM(){}
	protected:
		int TreeCode();
	private:
		huffmanTree CreateHuffTree();
		int HuffNodesInit(huffNode_t ** nodesPPtr);
		int FindMini(int * index,huffNode_t ** nodesPPtr,int nNodes);
		int MergeNode(huffNode_t* node1,huffNode_t *node2,huffNode_t * father);
		int GenerateHuffCode(huffmanTree tree);
		int DestroyHuffTree(huffmanTree tree);
		int SiglCharHuffCode(huffmanTree tree);
		int MultiCharHuffCode(huffmanTree tree);
		int GetHuffmanTreeHeight(huffmanTree root);
};
#endif
