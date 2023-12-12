
#include"Huffman_WT.h"
int Huffman_FM::TreeCode()
{
	huffmanTree huffTree=NULL;
	huffTree = CreateHuffTree();
	if(!huffTree)
	{
		cout<<"CreateHuffTree failed"<<endl;
		exit(0);
	}
	int ret = GenerateHuffCode(huffTree);
	if(ret < 0)
	{
		cout<<"GenerateHuffTree failed"<<endl;
		DestroyHuffTree(huffTree);
		exit(0);
	}
	DestroyHuffTree(huffTree);
	return 0;
}

huffmanTree Huffman_FM::CreateHuffTree()
{
	fm_int ret =0;
	huffNode_t *root =NULL;
	
	int nNodes = alphabetsize;
	int maxnNodes =alphabetsize*2-1;

	huffNode_t **huffNodesPPtr=(huffNode_t**)malloc(sizeof(huffNode_t*)*maxnNodes);

	ret=HuffNodesInit(huffNodesPPtr);
	if(ret <0)
	{
		cout<<"huffNodesinit failed"<<endl;
		exit(0);
	}
	if(alphabetsize==1)
	{
		root = huffNodesPPtr[0];
		free(huffNodesPPtr);
		return root;
	}
	int index1=0;
	int index2=0;
	huffNode_t * mini1 =NULL;
	huffNode_t * mini2 =NULL;

	while(nNodes < maxnNodes)
	{
		ret = FindMini(&index1,huffNodesPPtr,nNodes);
		if(index1<0)
		{
			cout<<"findMini failed"<<endl;
			return NULL;
		}
		mini1 = huffNodesPPtr[index1];
		huffNodesPPtr[index1]=NULL;

		ret = FindMini(&index2,huffNodesPPtr,nNodes);
		if(index2<0)
		{
			cout<<"findMini failed"<<endl;
			return NULL;
		}
		mini2 = huffNodesPPtr[index2];
		huffNodesPPtr[index2]=NULL;

		ret = MergeNode(mini1,mini2,huffNodesPPtr[nNodes]);
		if(ret < 0)
		{
			cout<<"mergeNode failed"<<endl;
			return NULL;
		}
		nNodes++;
	}
	root = huffNodesPPtr[maxnNodes-1];
	free(huffNodesPPtr);
	huffNodesPPtr = NULL;
	return root;
}

int Huffman_FM::HuffNodesInit(huffNode_t ** nodesPPtr)
{
	int i=0;
	int curIndex=0;
	for(i=0;i<CHAR_SET_SIZE;i++)
	{
		if(charFreq[i])
		{
			nodesPPtr[curIndex]=(huffNode_t*)malloc(sizeof(huffNode_t));
			if (!nodesPPtr[curIndex])
			{
				cout<<"malloc failed"<<endl;
				exit(0);
			}
			nodesPPtr[curIndex]->freq = charFreq[i];
			nodesPPtr[curIndex]->label = i;
			nodesPPtr[curIndex]->leftChild = NULL;
			nodesPPtr[curIndex]->rightChild = NULL;
			memset(nodesPPtr[curIndex]->code,0,sizeof(nodesPPtr[curIndex]->code));
			curIndex++;
		}
	}
	for(i=curIndex;i<2*curIndex-1;i++)
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

int Huffman_FM::FindMini(int * index,huffNode_t **nodesPPtr,int nNodes)
{
	if(nNodes <=0)
	{
		cout<<"Wrong paramater"<<endl;
		exit(0);
	}
	int i=0;
	long long int minFreq = ((long long int)1<<63)-1;
	*index = -1;

	for(i=0;i<nNodes;i++)
	{
		if(nodesPPtr[i]==NULL)
		{
			continue;
		}
		if(nodesPPtr[i]->freq < minFreq)
		{
			minFreq = nodesPPtr[i]->freq;
			*index =i;
		}
	}
	return *index;
}

int Huffman_FM::MergeNode(huffNode_t *node1,huffNode_t *node2,huffNode_t * father)
{
	father->freq = node1->freq + node2-> freq;
	father -> leftChild = node1;
	father -> rightChild = node2;
	return 0;
}

int Huffman_FM::GenerateHuffCode(huffmanTree tree)
{
	memset(codeTable,0,CHAR_SET_SIZE*CODE_MAX_LEN);
	if(tree -> leftChild ==NULL && tree->rightChild==NULL)
	{
		return SiglCharHuffCode(tree);
	}
	else
	{
		int height = GetHuffmanTreeHeight(tree);
		if(height >= CODE_MAX_LEN || height <0)
		{
			cout<<"Codes is too long"<<endl;
			exit(0);
		}
		return MultiCharHuffCode(tree);
	}
}


int Huffman_FM::DestroyHuffTree(huffmanTree tree)
{
	if(tree==NULL)
	{
		return 0;
	}
	if(tree->leftChild)
	{
		DestroyHuffTree(tree->leftChild);
		tree->leftChild = NULL;
	}
	if(tree->rightChild)
	{
		DestroyHuffTree(tree->rightChild);
		tree->rightChild = NULL;
	}
	free(tree);
	return 0;
}

int Huffman_FM::SiglCharHuffCode(huffmanTree tree)
{
	codeTable[tree->label][0]='0';
	return 0;
}

int Huffman_FM::MultiCharHuffCode(huffmanTree tree)
{
	if(tree->leftChild == NULL && tree ->rightChild ==NULL)
	{
		strcpy(codeTable[tree->label],tree->code);
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
	tree->leftChild->code[codeLen]='0';
	tree->rightChild->code[codeLen]='1';
	MultiCharHuffCode(tree->leftChild);
	MultiCharHuffCode(tree->rightChild);
	return 0;
}

int Huffman_FM::GetHuffmanTreeHeight(huffmanTree root)
{
	if(!root)
	{
		cout<<"Tree is empty"<<endl;
		exit(0);
	}
	int leftHeight = 0;
	int rightHeight =0;
	if(root->leftChild)
		leftHeight = GetHuffmanTreeHeight(root->leftChild);
	if(root -> rightChild)
		rightHeight = GetHuffmanTreeHeight(root->rightChild);
	return (leftHeight > rightHeight? leftHeight : rightHeight)+1;
}