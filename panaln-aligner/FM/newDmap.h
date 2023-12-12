
#include<string.h>
#include"loadkit.h"
#include"savekit.h"
#include"InArray.h"
#include"BaseClass.h"
#include<math.h>
#include<iostream>
#include <bitset>
using namespace std;
class Dmap
{
  public:
	Dmap();
	Dmap(fm_int data_num);
    fm_int GetMemorySize();
	void write(savekit &s);
	void load(loadkit &s);
	~Dmap(void);
	void SetValue(fm_int index, u64 v);
	void constructrank(int Blength, int SBlength);
	fm_int select(fm_int i);
	fm_int BS(InArray *B,fm_int l,fm_int r,fm_int i);
	void test();
	u64 GetValue(fm_int index);
	int getD(fm_int i);
	fm_int rank(fm_int i);
    fm_int rank2(fm_int i);
  private:
	u32 *data;
	fm_int datanum;
	int rankflag;
	InArray *Md;
	InArray *Sd;
	InArray *SBd;
	InArray *Bd;
	int Blength;
	int SBlength;
};