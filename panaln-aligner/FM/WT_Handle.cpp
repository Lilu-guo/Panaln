
#include"WT_Handle.h"
#include"Huffman_WT.h"
#include"Balance_WT.h"
#include"Hutacker_WT.h"

WT_Handle::WT_Handle():fm(new ABS_FM()),u(){}

WT_Handle::WT_Handle(const char * filename,int block_size,int r,int D,int shape)
{
	if(block_size<=0 || shape<0 || shape >2)
	{
		cout<<"WT_Handle::WT_handle error parmater"<<endl;
		exit(0);
	}

	switch(shape)
	{
		case 0: fm =new Hutacker_FM(filename,block_size,r,D);break;
		case 1: fm =new Huffman_FM(filename,block_size,r,D);break;
		case 2: fm =new Balance_FM(filename,block_size,r,D);break;
		default: fm=new Hutacker_FM(filename,block_size,r,D);break;
	}
	fm->BuildTree();
}

WT_Handle::WT_Handle(const WT_Handle &h):fm(h.fm),u(h.u){}

WT_Handle & WT_Handle::WT_Handle:: operator =(const WT_Handle & h)
{
	if(u.reattach(h.u))
	{
		delete fm;
	}
	fm = h.fm;
	return * this;
}

WT_Handle::~WT_Handle()
{
	if(u.only())
		delete fm;
}
