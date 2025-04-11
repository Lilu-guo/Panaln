/*============================================
# Filename: Balance_WT.cpp
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#include"Balance_WT.h"
int Balance_FM::TreeCode()
{
	balanceTree balTree=NULL;
	balTree = CreateBalanceTree();
	if(!balTree)
	{
		cout<<"CreateBalanceTree failed"<<endl;
		exit(0);
	}
	
	int ret = GenerateBalCode(balTree);
	if(ret < 0)
	{
		cout<<"GenarateBalTree failed"<<endl;
		DestroyBalTree(balTree);
		exit(0);
	}
	DestroyBalTree(balTree);
	return 0;
}


balanceTree Balance_FM::CreateBalanceTree()
{
	balNode_t * root = (balNode_t *)malloc(sizeof(balNode_t));
	if(!root)
	{
		cout<<"malloc failed"<<endl;
		return NULL;
	}
	memset(root,0,sizeof(balNode_t));

	int ret = 0;
	ret = BalanceTreeInit(root);
	if(ret < 0)
	{
		cout<<"balancetreeinit failed"<<endl;
		exit(0);
	}
	if(root->setSize == 1)
	{
		//single char set 
		root->label = root->set[0];
		return root;
	}
	
	ret = RecurseCreate(root);
	if(ret < 0)
	{
		cout<<"recurseCreate"<<endl;
		return NULL;
	}
	return root;
}


int Balance_FM::BalanceTreeInit(balanceTree tree)
{
	for(int i=0;i<CHAR_SET_SIZE;i++)
	{
		if(charMap[i])
		{
			tree->set[tree->setSize] = i;
			tree->setSize++;
		}
	}
	return 0;
}


int Balance_FM::RecurseCreate(balanceTree root)
{
	int ret = 0;
	if(!root || root->setSize == 1)
		return -1;
	balNode_t *leftNode = (balNode_t *)malloc(sizeof(balNode_t));
	if(!leftNode)
	{
		cout<<"Balance_FM::recurseCreate:  Malloc failed"<<endl;
		exit(0);
	}
	
	balNode_t *rightNode =(balNode_t *)malloc(sizeof(balNode_t));
	if(!rightNode)
	{
		cout<<"Balance_FM::recurseCreate:  Malloc failed"<<endl;
		exit(0);
	}

	memset(leftNode,0,sizeof(balNode_t));
	memset(rightNode,0,sizeof(balNode_t));

	leftNode->setSize = root->setSize/2;
	rightNode->setSize = root->setSize - leftNode->setSize;

	int i=0;
	for(i=0; i < leftNode->setSize; i++)
		leftNode->set[i] = root->set[i];
	
	for(i=0; i < rightNode->setSize; i++)
		rightNode->set[i] =root->set[i + root->setSize/2];
	
	root->leftChild = leftNode;
	root->rightChild= rightNode;

	if(leftNode -> setSize > 1)
	{
		ret = RecurseCreate(leftNode);
		if(ret < 0)
		{
			cout<<"Balance_FM::recurseCreate:  resurse Failed"<<endl;
			return ret;
		}
	}
	else
	{
		leftNode->label = leftNode ->set[0];
	}

	if(rightNode -> setSize > 1)
	{
		ret = RecurseCreate(rightNode);
		if(ret < 0)
		{
			cout<<"balance_FM::recurseCreate:  recurse Failed"<<endl;
			return ret;
		}
	}
	else
	{
		rightNode->label = rightNode->set[0];
	}
	return 0;
}


int Balance_FM::GenerateBalCode(balanceTree tree)
{
	memset(codeTable,0,CHAR_SET_SIZE*CODE_MAX_LEN);
	int ret=0;
	
	if(tree -> leftChild == NULL && tree -> rightChild == NULL)
		ret = SigCharBalCode(tree);
	else
		ret = MultiCharBalCode(tree);
	return ret;
}


int Balance_FM::SigCharBalCode(balanceTree tree)
{

	codeTable[tree->label][0]='0';
	return 0;
}


int Balance_FM::MultiCharBalCode(balanceTree root)
{
	int ret = 0;
	int codeLen = strlen(root->code);
	if(root->leftChild == NULL && root->rightChild == NULL)
	{
		strcpy(codeTable[root->label],root->code);
		return 0;
	}

	strcpy(root->leftChild -> code,root -> code);
	strcpy(root->rightChild -> code,root -> code);

	root -> leftChild -> code[codeLen] = '0';
	root -> rightChild -> code[codeLen] = '1';

	ret = MultiCharBalCode(root->leftChild);
	if(ret < 0)
	{
		cout<<"Balance_FM::multiCharBalCode: left branche failed"<<endl;
		exit(0);
	}

	ret = MultiCharBalCode(root -> rightChild);
	if(ret < 0)
	{
		cout<<"Balance_FM::multiCharBalCode: right branche fauiled"<<endl;
		exit(0);
	}

	return 0;
}


int Balance_FM::DestroyBalTree(balanceTree tree)
{
	if(!tree)
	{
		return 0;
	}
	if(tree->leftChild)
	{
		DestroyBalTree(tree->leftChild);
		tree->leftChild=NULL;
	}
	if(tree->rightChild)
	{
		DestroyBalTree(tree->rightChild);
		tree->rightChild=NULL;
	}
	free(tree);
	return 0;
}
