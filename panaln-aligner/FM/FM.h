
#ifndef FM_H
#define FM_H
#include"loadkit.h"
#include"savekit.h"
#include"ABS_WT.h"
#include"Huffman_WT.h"
#include"Balance_WT.h"
#include"Hutacker_WT.h"
#include"WT_Handle.h"
class FM
{
	public:
		FM(const char * filename,int blocksize,int r,int D,int shape);
		FM();
		~FM(){};
		FM(const FM & h):wt(h.wt){}
		FM& operator =(const FM&h){wt=h.wt;return *this;};
		
		void counting(const char *pattern,fm_int &num);
		fm_int * locating(const char *pattern,fm_int & num);
		unsigned char * extracting(fm_int pos,fm_int len);
		
		fm_int FM_Lookup(fm_int pos){return wt.Fm_Lookup(pos);};
		void Fmocc(fm_int l, fm_int r,unsigned char c ,fm_int *Left,fm_int * Right){wt.Fmocc(l,r,c,Left,Right);};
		void Fm16occ(fm_int l, fm_int r,fm_int *Left,fm_int * Right){wt.Fm16occ(l,r,Left,Right);};

		void Occ1p(unsigned char c, fm_int pos1, fm_int pos2, fm_int &occ1, fm_int &occ2){wt.Occ1p(c,pos1,pos2,occ1,occ2);};
		void Occ8p(unsigned char c, fm_int pos1, fm_int pos2, fm_int *occ81, fm_int *occ82){wt.Occ8p(c,pos1,pos2,occ81,occ82);};
		void Occ16p(fm_int pos1, fm_int pos2, fm_int *occ161, fm_int *occ162){wt.Occ16p(pos1,pos2,occ161,occ162);};

		int load(const char * indexfile);
		int loadIdx();
		int save(const char * indexfile);

		fm_int getN();
		fm_int getAlphabetSize();
		double getRuns();
		fm_int sizeInByte();
		fm_int sizeInByteForCount();
		fm_int sizeInByteForLocate();
		fm_int sizeInByteForExtract();
		fm_int sizeInByteForSAL();
		fm_int sizeInByteForRankL();
		fm_int sizeInByteFornewD();
		double compressRatio();
		double compressRatioForCount();
		double compressRatioForLocate();
		double compressRatioForExtract();
		double compressRatioForSal();
		double compressRatioForRankl();
		double compressRatioFornewD();
		fm_int TreePlainBlock(){return wt.TreePlainBlock();};
		fm_int TreeRLGBlock(){return wt.TreeRLGBlock();};
		fm_int TreeALL01Block(){return wt.TreeALL01Block();};
		fm_int TreeEF01Block(){return wt.TreeEF01Block();};
		fm_int TreePlainSize(){return wt.TreePlainSize();};
		fm_int TreeRLGSize(){return wt.TreeRLGSize();};
		fm_int TreeEF01Size(){return wt.TreeEF01Size();};
		fm_int MemSize(){return wt.MemSize();};
		fm_int NodeBlockSize(){return wt.NodeBlockSize();};
		fm_int CodeStyleSize(){return wt.CodeStyleSize();};

	private:
		WT_Handle wt;
};
#endif