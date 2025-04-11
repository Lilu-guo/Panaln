/*============================================
# Filename: Hutacker_WT.cpp
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#include"Hutacker_WT.h"
int Hutacker_FM::TreeCode()
{
	hutackerTree hutTree=NULL;
	hutTree = CreateHutackerTree();
	if(!hutTree)
	{
		cout<<"Hutacker_FM::TreeCode:CreateHutackerTree  failed"<<endl;
		exit(0);
	}
	int ret = 0;
	ret = GenerateHutackerCode(hutTree);
	if(ret < 0)
	{
		cout<<"Hutacker_FM::TreeCode:GenerateHutackerTree  failed"<<endl;
		DestroyHutackerTree(hutTree);
		exit(0);
	}
	DestroyHutackerTree(hutTree);
	return 0;
}

hutackerTree Hutacker_FM::CreateHutackerTree()
{
	hutackerTree root=NULL;
	int nNodes = this->alphabetsize;
	hutaNode_t **hutackerNodesPPtr=(hutaNode_t**)malloc(sizeof(hutaNode_t*)*nNodes);
	if(!hutackerNodesPPtr)
	{
		cout<<"HutackerTree::CreateHutackerTree:malloc  failed"<<endl;
		exit(0);
	}
	int ret = HutackerNodesInit(hutackerNodesPPtr);
	if(ret < 0)
	{
		cout<<"HutackerTree::CreatehutackerTree:hutackerNodesInit  failed"<<endl;
		exit(0);
	}
	if(alphabetsize==1)
	{
		root = hutackerNodesPPtr[0];
		free(hutackerNodesPPtr);
		return root;
	}

	int index1=0;
	int index2=0;
	//hutaNode_t * node1=NULL;
	//hutaNode_t * node2=NULL;

	int realNodeNum = nNodes;
	while(realNodeNum  > 1)
	{
		ret = FindMiniTwoNodes(hutackerNodesPPtr,nNodes,&index1,&index2);
		if(ret < 0)
		{
			cout<<"HutackerFM::CreateHutackerTree:FindMiniTwoNodes failed"<<endl;
			exit(0);
		}
		ret = MergeNodes(hutackerNodesPPtr,index1,index2,nNodes);
		if(ret < 0)
		{
			cout<<"HutackerFM::CreateHutackerTree:MergeNodes failed"<<endl;
			exit(0);
		}
		realNodeNum--;
	}
	for(int i=0;i<nNodes;i++)
	{
		if(hutackerNodesPPtr[i])
		{
			root = hutackerNodesPPtr[i];
			break;
		}
	}
	return root;
}


int Hutacker_FM::HutackerNodesInit(hutaNode_t ** hutNodesPPtr)
{
	int curIndex = 0;
	for(int i=0;i<CHAR_SET_SIZE;i++)
	{
		if(charFreq[i])
		{
			hutNodesPPtr[curIndex]=(hutaNode_t*)malloc(sizeof(hutaNode_t));
			hutNodesPPtr[curIndex]->label=i;
			hutNodesPPtr[curIndex]->freq=charFreq[i];
			hutNodesPPtr[curIndex]->leftChild=NULL;
			hutNodesPPtr[curIndex]->rightChild=NULL;
			curIndex++;
		}
	}
	return 0;
}


int Hutacker_FM::FindMiniTwoNodes(hutaNode_t ** huNodesPPtr,int nNodes,int * index1,int *index2)
{
	fm_int miniFreq = ((long long int)1<<63)-1;
	*index1 = *index2 = -1;
	for(int i=0;i < nNodes-1;i++)
	{
		if(!huNodesPPtr[i])
		{
			continue;
		}
		for(int j=i+1;j<nNodes;j++)
		{
			if(!huNodesPPtr[j])
			{
				continue;
			}
			if(huNodesPPtr[i]->freq + huNodesPPtr[j]->freq < miniFreq)
			{
				*index1 = i;
				*index2 = j;
				miniFreq = huNodesPPtr[i]->freq + huNodesPPtr[j]->freq;
			}
			if(huNodesPPtr[j]->leftChild==NULL && huNodesPPtr[j]->rightChild==NULL)
				break;

		}
	}
	if(*index1==-1 || *index2 == -1)
	{
		cout<<"Hutacker_FM::FindMiniTwoNode  failed"<<endl;
		exit(0);
	}
	return 0;
}


int Hutacker_FM::MergeNodes(hutaNode_t ** hutaNodesPPtr,int index1,int index2,int nNodes)
{
	hutaNode_t * tmp = (hutaNode_t*)malloc(sizeof(hutaNode_t));
	if(!tmp)
	{
		cout<<"Hutacker_FM::MergeNodes malloc failed"<<endl;
		exit(0);
	}
	
	tmp->freq=hutaNodesPPtr[index1]->freq + hutaNodesPPtr[index2]->freq;
	tmp->leftChild=hutaNodesPPtr[index1];
	tmp->rightChild=hutaNodesPPtr[index2];
	tmp->level=0;//important
	hutaNodesPPtr[index1]=tmp;
	hutaNodesPPtr[index2]=NULL;
	return 0;
}


int Hutacker_FM::GenerateHutackerCode(hutackerTree root)
{
	memset(codeTable,0,CHAR_SET_SIZE*CODE_MAX_LEN);
	int ret=0;
	if(root->leftChild == NULL && root->rightChild == NULL)
		ret = SingleCharHutackerCode(root);
	else
		ret = MultiCharHutackerCode(root);
	return ret;
}


int Hutacker_FM::SingleCharHutackerCode(hutackerTree root)
{
	codeTable[root->label][0]='0';
	return 0;
}


int Hutacker_FM::MultiCharHutackerCode(hutackerTree root)
{
	int charDepth[CHAR_SET_SIZE];
	memset(charDepth,0,sizeof(charDepth));
	int ret = ComputerCharDepth(root,charDepth);
	if(ret < 0)
	{
		cout<<"Hutacker_FM::MutiCharHutackerCode:ComputerCharDepth failed"<<endl;
		exit(0);
	}

	int i=0;
	int j=0;
	char preCode[CODE_MAX_LEN];
	memset(preCode,0,sizeof(preCode));
	for(i=0;i<CHAR_SET_SIZE;i++)
	{
		if(!charDepth[i])
			continue;
		if(strlen(preCode)==0)
		{
			for(j=0;j<charDepth[i];j++)
				preCode[j]='0';
			strcpy(codeTable[i],preCode);
		}
		else
		{
			for(j=strlen(preCode)-1;preCode[j]=='1' && j>=0;j--)
				preCode[j]='\0';
			preCode[strlen(preCode)-1]='1';
			while(strlen(preCode) < (unsigned int)charDepth[i])
				preCode[strlen(preCode)]='0';
			strcpy(codeTable[i],preCode);
		}
	}
	return 0;
}

			
int Hutacker_FM::DestroyHutackerTree(hutackerTree root)
{
	if(!root)
		return 0;
	if(root->leftChild)
	{
		DestroyHutackerTree(root->leftChild);
		root->leftChild = NULL;
	}
	if(root->rightChild)
	{
		DestroyHutackerTree(root->rightChild);
		root->rightChild = NULL;
	}
	free(root);
	return 0;
}


int Hutacker_FM::GetHutackerTreeHeight(hutackerTree root)
{
	if(root -> leftChild == NULL && root -> rightChild == NULL)
		return 0;
	int leftHeight  =0;
	int rightHeight =0;
	
	if(root->leftChild)
		leftHeight = GetHutackerTreeHeight(root->leftChild);
	if(root->rightChild)
		rightHeight =GetHutackerTreeHeight(root -> rightChild);
	return (leftHeight > rightHeight ? leftHeight : rightHeight )+1;
}


int Hutacker_FM::ComputerCharDepth(hutackerTree root,int * charDepth)
{
	int ret = 0;
	if(root->leftChild ==NULL && root->rightChild == NULL)
	{
		//tree leaf 
		charDepth[root->label] = root->level;
		return 0;
	}

	root->leftChild->level = root->level+1;
	root->rightChild->level= root->level+1;
	
	if(root->leftChild)
	{
		ret = ComputerCharDepth(root->leftChild,charDepth);
		if(ret < 0)
			return ret;
	}
	if(root->rightChild)
	{
		ret = ComputerCharDepth(root->rightChild,charDepth);
		if(ret<0)
			return ret;
	}
	return 0;
}
