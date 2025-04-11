#include<stdlib.h>
#include<sys/time.h>
#include"FM.h"
#include<iostream>
using namespace std;
#define  times 10000
int main(int argc,char ** argvs)
{
	
	if(argc!=3)
	{
		cout<<"use it like:   ./my_fm file speedlevel"<<endl;
		exit(0);
	}

	time_t t3 = clock();
	int speedlevel = atoi(argv[2]);
	FM fm(argvs[1], speedlevel);
	time_t t4 = clock();
	cout<<"build time: "<<(t4-t3)/1000000.0<<endl;
		//cout<<"fm is ready"<<endl;
	int n =fm.GetN();
	int size = fm.SizeInByte_count();
	cout<<"szie "<<size/(n*1.0)<<endl;
	unsigned char *p[times];
	for(int i=0;i<times;i++)
	{
		int x = rand()%( n - 50);
		p[i]=fm.Extracting(x,20);
	}
	//cout<<"pattern is ready"<<endl;
	double totnum =0;
	int num = 0;
	struct timeval start;
	struct timeval end;
	gettimeofday(&start,NULL);
	for(int i=0;i<times;i++)
	{
		fm.Counting((const char *)p[i],num);
		//totnum = totnum + num;
	}
	//cout<<totnum<<endl;
	gettimeofday(&end,NULL);
	cout<<"count-time: "<<((end.tv_sec-start.tv_sec)*1000000.0+end.tv_usec - start.tv_usec)/times<<endl;
	for(int i=0;i<times;i++)
		delete [] p[i];
	return 0;
}

