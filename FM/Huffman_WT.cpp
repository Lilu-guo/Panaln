/*============================================
# Filename: Huffman_WT.cpp
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#include"Huffman_WT.h"
int Huffman_FM::TreeCode()                                                                   //总函数（树形编码）
{
	huffmanTree huffTree=NULL;
	huffTree = CreateHuffTree();                                                             //创建Huffman树
	if(!huffTree)
	{
		cout<<"CreateHuffTree failed"<<endl;                                                 //失败
		exit(0);
	}
	int ret = GenerateHuffCode(huffTree);                                                    //遍历树，产生Huffman码表，放到codeTable中
	if(ret < 0)                                                                              //ret为0，代表成功
	{
		cout<<"GenerateHuffTree failed"<<endl;
		DestroyHuffTree(huffTree);
		exit(0);
	}
	DestroyHuffTree(huffTree);                                                               //销毁树
	return 0;
}

huffmanTree Huffman_FM::CreateHuffTree()                                                     //具体构建Huffman树（根据字符频率）
{
	fm_int ret =0;
	huffNode_t *root =NULL;                                                                  //根节点
	
	int nNodes = alphabetsize;                                                               //字符表大小
	int maxnNodes =alphabetsize*2-1;                                                         //最多节点数（Huffman树性质）

	huffNode_t **huffNodesPPtr=(huffNode_t**)malloc(sizeof(huffNode_t*)*maxnNodes);          //节点指针的指针（指针数组）

	ret=HuffNodesInit(huffNodesPPtr);                                                        //节点初始化
	if(ret <0)
	{
		cout<<"huffNodesinit failed"<<endl;
		exit(0);
	}
	if(alphabetsize==1)                                                                      //只有一个字符表
	{
		root = huffNodesPPtr[0];
		free(huffNodesPPtr);
		return root;
	}
	int index1=0;
	int index2=0;
	huffNode_t * mini1 =NULL;                                                               //节点
	huffNode_t * mini2 =NULL;

	while(nNodes < maxnNodes)                                                               //遍历完所有节点，刚好保证构建完Huffman树
	{
		ret = FindMini(&index1,huffNodesPPtr,nNodes);                                       //找频率最小的节点一
		if(index1<0)
		{
			cout<<"findMini failed"<<endl;
			return NULL;
		}
		mini1 = huffNodesPPtr[index1];                                                      //节点一
		huffNodesPPtr[index1]=NULL;                                                         //原数组中，置空

		ret = FindMini(&index2,huffNodesPPtr,nNodes);                                       //找频率最小的节点二
		if(index2<0)
		{
			cout<<"findMini failed"<<endl;
			return NULL;
		}
		mini2 = huffNodesPPtr[index2];
		huffNodesPPtr[index2]=NULL;

		//mere mini1 and mini2 nodes
		ret = MergeNode(mini1,mini2,huffNodesPPtr[nNodes]);                                //合并最小的两个，节点一、节点二
		if(ret < 0)
		{
			cout<<"mergeNode failed"<<endl;
			return NULL;
		}
		nNodes++;                                                                           //节点计数加1
	}
	root = huffNodesPPtr[maxnNodes-1];                                                      //数组中最后一个节点，刚好为root节点
	free(huffNodesPPtr);
	huffNodesPPtr = NULL;                                                                   //数组用完，置空
	return root;                                                                            //返回根节点
}

int Huffman_FM::HuffNodesInit(huffNode_t ** nodesPPtr)                                     //节点初始化
{
	int i=0;
	int curIndex=0;                                                                        //节点下标（CHAR_SET_SIZE为256）
	for(i=0;i<CHAR_SET_SIZE;i++)                                                           //字符集个数
	{
		if(charFreq[i])                                                                    //频率不为0（256中一定有频率为0的）
		{
			nodesPPtr[curIndex]=(huffNode_t*)malloc(sizeof(huffNode_t));                   //初始每个字符为一个节点
			if (!nodesPPtr[curIndex])                                                      //节点指针为空
			{
				cout<<"malloc failed"<<endl;
				exit(0);
			}
			nodesPPtr[curIndex]->freq = charFreq[i];                                       //在Getfile()中已计算
			nodesPPtr[curIndex]->label = i;                                                //字符，赋给label
			nodesPPtr[curIndex]->leftChild = NULL;                                         //左孩子
			nodesPPtr[curIndex]->rightChild = NULL;                                        //右孩子
			memset(nodesPPtr[curIndex]->code,0,sizeof(nodesPPtr[curIndex]->code));         //置0
			curIndex++;                                                                    //下标计数，加1（只统计频率不为0的）
		}
	}
	for(i=curIndex;i<2*curIndex-1;i++)                                                     //又开辟，与不为空字符大小-1的节点（用于非叶子节点，满二叉树性质）
	{
		nodesPPtr[i] = (huffNode_t*)malloc(sizeof(huffNode_t));
		if (!nodesPPtr[i])
		{
			cout<<"malloc failed"<<endl;
			exit(0);
		}
		nodesPPtr[i]->freq = 0;
		nodesPPtr[i]->label = 0;
		nodesPPtr[i]->leftChild =nodesPPtr[i]->rightChild=NULL;
		memset(nodesPPtr[i]->code,0,sizeof(nodesPPtr[i]->code));
	}
	return 0;
}

int Huffman_FM::FindMini(int * index,huffNode_t **nodesPPtr,int nNodes)                  //找到Freq最小的节点下标（nodesPPtr已存在）
{
	if(nNodes <=0)
	{
		cout<<"Wrong paramater"<<endl;
		exit(0);
	}
	int i=0;
	long long int minFreq = ((long long int)1<<63)-1;                                    //初始化为long的最大值
	*index = -1;

	for(i=0;i<nNodes;i++)                                                                //遍历所有节点
	{
		if(nodesPPtr[i]==NULL)                                                           //节点指针为空，则跳过
		{
			continue;
		}
		if(nodesPPtr[i]->freq < minFreq)                                                 //找到比minFreq小的节点
		{
			minFreq = nodesPPtr[i]->freq;                                                //更新minFreq
			*index =i;                                                                   //记录该节点的下标
		}
	}
	return *index;                                                                       //返回Freq最小的节点下标
}

int Huffman_FM::MergeNode(huffNode_t *node1,huffNode_t *node2,huffNode_t * father)      //合并两孩子，产生父节点
{
	father->freq = node1->freq + node2-> freq;                                          //父节点，频率等于左右孩子之和
	father -> leftChild = node1;                                                        //挂为左孩子
	father -> rightChild = node2;                                                       //挂为右孩子
	return 0;
}

int Huffman_FM::GenerateHuffCode(huffmanTree tree)                                      //产生Huffman码表（此时的树已经存在）
{
	memset(codeTable,0,CHAR_SET_SIZE*CODE_MAX_LEN);                                     //置0
	if(tree -> leftChild ==NULL && tree->rightChild==NULL)
	{
		return SiglCharHuffCode(tree);                                                  //单个字符
	}
	else
	{
		int height = GetHuffmanTreeHeight(tree);                                        //获得Huffman树高
		if(height >= CODE_MAX_LEN || height <0)
		{
			cout<<"Codes is too long"<<endl;
			exit(0);
		}
		return MultiCharHuffCode(tree);                                                //具体产生Huffman编码
	}
}


int Huffman_FM::DestroyHuffTree(huffmanTree tree)                       //销毁Huffman树
{
	if(tree==NULL)
	{
		return 0;
	}
	if(tree->leftChild)
	{
		DestroyHuffTree(tree->leftChild);                              //递归
		tree->leftChild = NULL;                                        //指针置空
	}
	if(tree->rightChild)
	{
		DestroyHuffTree(tree->rightChild);
		tree->rightChild = NULL;
	}
	free(tree);
	return 0;
}

int Huffman_FM::SiglCharHuffCode(huffmanTree tree)                     //单字符
{
	codeTable[tree->label][0]='0';
	return 0;
}


int Huffman_FM::MultiCharHuffCode(huffmanTree tree)                    //树已经存在，此时产生codeTable（递归的遍历所有路径）
{
	//for multi character code
	if(tree->leftChild == NULL && tree ->rightChild ==NULL)
	{
		//tree is a single leaf node
		strcpy(codeTable[tree->label],tree->code);                     //找到叶子节点，把当前tree走过的'0'/'1'串，赋给对应字符的codeTable
		return 0;
	}
	int codeLen = strlen(tree->code);
	if(codeLen == CODE_MAX_LEN-1)
	{
		cout<<"the codes is too long"<<endl;
		exit(0);
	}
	strcpy(tree->leftChild->code,tree->code);
	strcpy(tree->rightChild->code,tree->code);

	tree->leftChild->code[codeLen]='0';                                //给codeTable赋值（此时的树已经存在）
	tree->rightChild->code[codeLen]='1';

	MultiCharHuffCode(tree->leftChild);                                //递归更新树根
	MultiCharHuffCode(tree->rightChild);

	return 0;
}

int Huffman_FM::GetHuffmanTreeHeight(huffmanTree root)                 //Huffman树的高度
{
	if(!root)
	{
		cout<<"Tree is empty"<<endl;
		exit(0);
	}
	int leftHeight = 0;
	int rightHeight =0;
	if(root->leftChild)
		leftHeight = GetHuffmanTreeHeight(root->leftChild);            //递归
	if(root -> rightChild)
		rightHeight = GetHuffmanTreeHeight(root->rightChild);
	return (leftHeight > rightHeight? leftHeight : rightHeight)+1;     //在较大树高的孩子基础上，加1
}

