
#include"FM.h"

FM::FM(const char *filename,int blocksize,int r,int D,int shape):wt(filename,blocksize,r,D,shape){}

FM::FM():wt(){}

fm_int FM::getN(){
	return wt.GetN();
}

fm_int FM::getAlphabetSize(){
	return wt.GetAlphabetsize();
}

double FM::getRuns()
{
	return wt.GetRuns();
}

fm_int FM::sizeInByte()
{
	return wt.SizeInByte();
}

fm_int FM::sizeInByteForCount()
{
	return wt.SizeInByte_count();
}

fm_int FM::sizeInByteForLocate()
{
	return wt.SizeInByte_locate();
}

fm_int FM::sizeInByteForExtract()
{
	return wt.SizeInByte_extract();
}

fm_int FM::sizeInByteForSAL()
{
	return wt.SizeinByte_sal();
}

fm_int FM::sizeInByteForRankL()
{
	return wt.SizeinByte_rankl();
}

fm_int FM::sizeInByteFornewD()
{
	return wt.SizeinByte_newD();
}

double FM::compressRatioForCount(){
	return sizeInByteForCount()/(getN()*1.0);
}


double FM::compressRatioForLocate(){
	return sizeInByteForLocate()/(getN()*1.0);
}

double FM::compressRatioForExtract(){
	return sizeInByteForExtract()/(getN()*1.0);
}
double FM::compressRatioForSal(){
	return wt.SizeinByte_sal()/(getN()*1.0);
}
double FM::compressRatioFornewD(){
	return wt.SizeinByte_newD()/(getN()*1.0);
}
double FM::compressRatioForRankl(){
	return wt.SizeinByte_rankl()/(getN()*1.0);
}

double FM::compressRatio(){
	return sizeInByte()/(getN()*1.0);
}

int FM::save(const char * indexfile)
{
	savekit s(indexfile);
	s.writeu64(198809102510);
	wt.Save(s);
	s.close();
	return 0;
}


int FM::load(const char * indexfile)
{
	loadkit s(indexfile);
	unsigned long long int magicnum=0;
	s.loadu64(magicnum);
	if(magicnum!=198809102510)
	{
		cerr<<"Not a FM_Index file"<<endl;
		exit(0);
	}
	wt.Load(s);
	s.close();
	return 0;
}

int FM::loadIdx()
{
	wt.loadIdx();
	return 0;
}

void FM::counting(const char * pattern,fm_int &num)
{
	wt.Counting(pattern,num);
}


fm_int * FM::locating(const char * pattern,fm_int & num)
{
	return wt.Locating(pattern,num);
}


unsigned char * FM::extracting(fm_int pos,fm_int len)
{
	return wt.Extracting(pos,len);
}

