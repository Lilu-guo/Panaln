#include<stdlib.h>
#include<string.h>
#include"FM.h"
#include<ctime>
#include<fstream>
#include<iostream>
#include<vector>
#include<unistd.h>
#include <algorithm>
#define GRCh38
using namespace std;

void AAAloc( FM* fm,string refile,string refile2){
	i64 t=0;
	ifstream infile(refile.c_str(),ios::in);
	ofstream outfile(refile2.c_str(),ios::app);
	string tmp("");
	i64 num = 0;
	i64 cnt = 0;
	while(getline(infile,tmp))
	{
		clock_t s=clock();
    	fm->locating(tmp.c_str(),num);
    	clock_t e=clock();
    	t+=(e-s);
		++cnt;
		tmp = "";
	}
	outfile<<"loc time:		"<<(double)t/cnt/1000<<"ms"<<endl;
	outfile.close();
	infile.close();
}
void getextract( FM* fm,int loop,int patternlength,string refile){
	i64 t=0;
	i64 length=fm->getN();
	ofstream outfile(refile.c_str(),ios::app);
	for(int i=0;i<loop;i++){
		i64 num=0;
		i64 start=rand()%(length-patternlength);
    	auto temp=fm->extracting(start,patternlength);
		#ifdef GRCh38
		while(1)
		{	
			int j = 0;
			for(j;j < patternlength;++j)
			{
				if(temp[j] == 'N')
				{
					start=rand()%(length-patternlength);
					temp=fm->extracting(start,patternlength);
					break;
				}
			}
			if(j == patternlength)break;
		}
		#endif
		outfile<<temp<<endl;
	}
	outfile.close();
}
void testlocate( FM* fm,int loop,int patternlength,string refile){
	i64 t=0;
	i64 length=fm->getN();
	ofstream outfile(refile.c_str(),ios::app);
	for(int i=0;i<loop;i++){
		i64 num=0;
		i64 start=rand()%(length-patternlength);
    	auto temp=fm->extracting(start,patternlength);
		#ifdef GRCh38
		while(1)
		{	
			int j = 0;
			for(j;j < patternlength;++j)
			{
				if(temp[j] == 'N')
				{
					start=rand()%(length-patternlength);
					temp=fm->extracting(start,patternlength);
					break;
				}
			}
			if(j == patternlength)break;
		}
		#endif
		clock_t s=clock();
    	fm->locating((char*)temp,num);
    	clock_t e=clock();
    	t+=(e-s);
	}
	outfile<<"locatetime:"<<(double)t/loop/1000<<"ms"<<endl;
	outfile.close();
}
void testextract( FM* fm,int loop,int patternlength,string refile){
	i64 t=0;
	i64 length=fm->getN();
	ofstream outfile(refile.c_str(),ios::app);
	for(int i=0;i<loop;i++){
		i64 num=0;
		i64 start=rand()%(length-patternlength);
		clock_t s=clock();
    	fm->extracting(start,patternlength);
    	clock_t e=clock();
    	t+=(e-s);
	}
	outfile<<"extracttime:"<<(double)t/loop<<"us"<<endl;
	outfile.close();
}
void testcount( FM* fm,int loop,int patternlength,string refile){
	i64 t=0;
	i64 length=fm->getN();
	ofstream outfile(refile.c_str(),ios::app);
	for(int i=0;i<loop;i++){
		i64 num=0;
		i64 start=rand()%(length-patternlength);
    	unsigned char* tmp = fm->extracting(start,patternlength);
		clock_t s=clock();
		fm->counting((char*)tmp,num);
    	clock_t e=clock();
    	t+=(e-s);
	}
	outfile<<"counttime:"<<(double)t/loop<<"us"<<endl;
	outfile.close();
}

void printQueryUsage(){
	cout << "==========================================================================================" << endl;
	cout << " c pattern: count the number of occurrences of the pattern string in the file" << endl;
	cout << " l pattern: locate the postion of the pattern string in the file" << endl;
	cout << " e pos len: extract the substring of the file from postion pos ,the length of which is len" << endl;
	cout << " s save the index in file filename-blockszie-samplerate-multiple-shape.index" << endl;
	cout << " q quit the query" << endl;
	cout << "==========================================================================================" << endl;
}
void showpos(fm_int *pos,fm_int num){
	sort(pos,pos+num);
	getchar();
	for(fm_int i=0;i<num;i++){
		cout<<pos[i]<<endl;
		if ((i + 1) % 20 == 0)
		{
			char command;
			cout << "===================press enter to continue=======press q to quit==================="<<endl;
			command = getchar();
			if(command=='q'){
				return;
			}
		}
	}
}
void printinformation(FM*fm){
	cout<<"the Number of characters:        "<<fm->getN()<<endl;
	cout<<"the Alphabetsize:                "<<fm->getAlphabetSize()<<endl;
	cout<<"the Runs is:                     "<<fm->getRuns()<<endl;
	cout<<"the size of the Index:           "<<1.0*fm->sizeInByte()/1024/1024<<"MB"<<endl;
	cout<<"the size of the Index(Count)     "<<1.0*fm->sizeInByteForCount()/1024/1024<<"MB"<<endl;
	cout<<"the size of the Index(Locate):   "<<1.0*fm->sizeInByteForLocate()/1024/1024<<"MB"<<endl;
	cout<<"the size of the Index(Extract):  "<<1.0*fm->sizeInByteForExtract()/1024/1024<<"MB"<<endl;
	cout<<"the size of the Wavelet tree:    "<<1.0*fm->sizeInByteForCount()/1024/1024<<"MB"<<endl;
	cout<<"the size of the Sample for SA:   "<<1.0*fm->sizeInByteForSAL()/1024/1024<<"MB"<<endl;
	cout<<"the size of the Sample for ISA:  "<<1.0*fm->sizeInByteForRankL()/1024/1024<<"MB"<<endl;
	cout<<"the size of the Flag array newD: "<<1.0*fm->sizeInByteFornewD()/1024/1024<<"MB"<<endl;
	cout<<"the Ratio of the Index:          "<<fm->compressRatio()*100<<endl;
	cout<<"the Ratio of the Index(Count):   "<<fm->compressRatioForCount()*100<<endl;
	cout<<"the Ratio of the Index(Locate):  "<<fm->compressRatioForLocate()*100<<endl;
	cout<<"the Ratio of the Index(Extract): "<<fm->compressRatioForExtract()*100<<endl;
	cout<<"the Ratio of the Wavelet tree:   "<<fm->compressRatioForCount()*100<<endl;
	cout<<"the Ratio of the Sample for SA:  "<<fm->compressRatioForSal()*100<<endl;
	cout<<"the Ratio of the Sample for ISA: "<<fm->compressRatioForRankl()*100<<endl;
	cout<<"the Ratio of the Flag array newD:"<<fm->compressRatioFornewD()*100<<endl;
}
int main(int argc, char *argv[]){
	string filename="";
	int r=16;
	int blocksize=256;
	int D=128;
	int shape=1;
	int ch;
	FM *fm;
	fm=new FM("/home/lab/gll/panaln/mg-ref/multigenome_bubble_GRCh38.fasta.ref",blocksize,r,32,1);
	cout<<"test................"<<endl;
	printinformation(fm);
	string pattern="CCCC";
	fm_int num=0;
	fm->counting(pattern.c_str(),num);
	fm_int* pos;
	pos=fm->locating(pattern.c_str(),num);
	cout<<pattern<<" occurs "<<num<<" times"<<endl;
	cout<<"occ pos:"<<endl;
	for(int i=0;i<num;i++){
		cout<<pos[i]<<" ";
	} 
	cout<<endl;  
	return 0;
	
	if(argc<2){
		return 0;
	}
	while((ch = getopt(argc, argv,"f:b:d:r:s:")) != -1){
		switch(ch)
		{
			case 'f':
				filename=string(optarg);
				break;
			case 'b':
				blocksize=atoi(optarg);
				break;
			case 'r':
				r=atoi(optarg);
				break;
			case 'd':
				D=atoi(optarg);
				break;
			case 's':
				shape=atoi(optarg);
				break;				
		}
	}
	cerr << "parameters:"<< endl;
	cerr << "	blocksize: "<< blocksize << endl;
	cerr << "	samplerate: "<< D << endl;
	cerr << "	multiple: "<< r << endl;
	cerr << "	shape: " << shape << endl;
	string indexfile=filename+"-"+to_string(blocksize)+"-"+to_string(D)+"-"+to_string(r)+"-"+to_string(shape)+".876790783";
	if(access(indexfile.c_str(),0)==0){
		cerr << "the index already exists "<< endl;
		cerr << "start load the index " << endl;
		fm=new FM();
		fm->load(indexfile.c_str());
		cerr << "load the index finished" << endl;
	}else{
		cerr << "start build the index of "<< filename << endl;
		fm=new FM(filename.c_str(),blocksize,r,D,shape);
		cerr << " index build finished " << endl;		
	}
	string qqqqq("AAAqwe");
	testlocate(fm,1000,20,qqqqq);
	return 0;
#ifdef counttime
	ofstream fout2;
	string outstr("AAAcounttime.txt");
	i64 memsize = fm->MemSize() / 8;
	i64 countsize = fm->sizeInByteForCount();
	fout2.open(outstr.c_str(),ios_base::out | ios_base::app);
	fout2<<"b:		"<<blocksize<<"		r:		"<<r<<endl;
	fout2<<"Radio:		"<<fm->compressRatioForCount()*100<<endl;
	fout2<<"Memsize:	"<<memsize*1.0 / countsize * 100<<endl;
	fout2<<"Newblock:	"<<100 - (memsize*1.0 / countsize * 100)<<endl;
	fout2.close();
	testcount(fm,5000,20,outstr);
	return 0;
#endif
#ifdef spaceradio
	if(blocksize == 256)
	{
		string indexname("AAAindexradio.txt");
		ofstream indexout;
		indexout.open(indexname.c_str(),ios_base::out | ios_base::app);
		indexout<<filename<<endl;
		indexout<<"D:	"<<D<<endl;
		i64 sasize = fm->sizeInByteForSAL();
		i64 dsize = fm->sizeInByteFornewD();
		i64 isasize = fm->sizeInByteForRankL();
		i64 fmsize = fm->sizeInByte();
		indexout<<"IndexRadio:	"<<fm->compressRatio() * 100<<" %"<<endl;
		indexout<<"SA:			"<<sasize * 1.0 / fmsize * 100<<" %"<<endl;
		indexout<<"newD:		"<<dsize * 1.0 / fmsize * 100<<" %"<<endl;
		indexout<<"ISA:			"<<isasize * 1.0 / fmsize * 100<<" %"<<endl;
		indexout<<"TotalRadio:	"<<(sasize + dsize + isasize) * 1.0 / fmsize * 100<<" %"<<endl;
		indexout.close();

		string locname("AAAlcexradio.txt");
		ofstream locout;
		locout.open(locname.c_str(),ios_base::out | ios_base::app);
		locout<<filename<<endl;
		locout<<"D:		"<<D<<endl;
		locout<<"IndexRadio:	"<<fm->compressRatio() * 100<<" %"<<endl;
		locout.close();
		testlocate(fm,5000,20,locname);
		testextract(fm,5000,20,locname);
	}
	return 0;
#endif

#ifdef _auto
	ofstream outfile;
	string statics="statics.txt";
	outfile.open(statics.c_str(),ios::app);

	outfile<<"	filename:  "<<filename<<endl;
	outfile<<"	blocksize: "<<blocksize<<endl;
	outfile<<"	Ratio:	   "<<fm->compressRatio()*100<<endl;
	outfile.close();
return 0;
#endif

#ifdef zhengquexing
	string file3("AAAsources");
	ofstream fout3;
	fout3.open(file3.c_str());
	int i = 0,len3 = fm->getN() - 1;
	for(;i < len3 - 1000;i += 1000)
		fout3<<fm->extracting(i,1000);
	fout3<<fm->extracting(i,len3 - i);
	fout3.close();
	return 0;

	FILE * fp = fopen(filename.c_str(),"r+");
	if(fp==NULL)
	{
		cout<<"Be sure the file is available "<<filename<<endl;
		exit(0);
	}
	fseek(fp,0,SEEK_END);
	i64 lenn = ftell(fp)+1;
	unsigned char * T = new unsigned char[lenn];
	fseek(fp,0,SEEK_SET);

	fm_int e=0;
	fm_int num=0;
	while((e=fread(T+num,sizeof(uchar),lenn-1-num,fp))!=0)
		num = num +e;
	if(num!=lenn-1)
	{
		cout<<num<<"  "<<lenn<<"  Read source file failed"<<endl;
		exit(0);
	}
	T[lenn-1]=0;
	fclose(fp);
	for(i64 qqq = 0;qqq <= lenn - 21;qqq++)
	{
		uchar *substring=fm->extracting(qqq,20);
		for(int www = 0;www < 20;www++)
		{
			if(substring[www] != T[www+qqq])
			{
				cout<<"NO!!!!	:	"<<qqq<<endl;
				exit(1);
			}
		}
		if(qqq % 1000000 == 0)
			cout<<qqq<<endl;
	}
	cout<<"YES!!!!"<<endl;
	return 0;
#endif
	printinformation(fm);
	bool flag=true;
	char command;
	while(flag){
		printQueryUsage();
		cin>>command;
		cout<<command<<endl;
		string pattern;
		string quit("q");
		fm_int num=0;
		fm_int start=0;
		fm_int len=0;
		uchar* substring=NULL;
		fm_int *pos;
		switch(command)
		{
			case 'c':
				cout<<"Counting..."<<endl;
				while(1)
				{
					cin>>pattern;
					if(pattern == quit)break;
					fm->counting(pattern.c_str(),num);
					cout<<pattern<<" occurs "<<num<<" times"<<endl;
				}
				break;
			case 'l':
				cout<<"Locating..."<<endl;
				while(1)
				{
					cin>>pattern;
					if(pattern == quit)break;
					pos=fm->locating(pattern.c_str(),num);
					cout<<pattern<<" occurs "<<num<<" times"<<endl;
					showpos(pos,num);
				}
				break;
			case 'e':
				cout<<"Extracting..."<<endl;
				while(1)
				{
					cin>>start>>len;
					substring=fm->extracting(start,len);
					if(substring!=NULL){
						cout<<"the substring is "<<endl;
						cout<<substring<<endl;
					}
				}
				break;
			case 's':
				cout<<"Saving..."<<endl;
				fm->save(indexfile.c_str());
				break;
			case 'q':
				cout<<"quit!"<<endl;
				flag=false;
				break;
			default:
				cout<<"error input"<<endl;
				printQueryUsage();
				flag=false;
				break;
		}
	}
	return 0;
}
