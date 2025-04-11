#include<stdlib.h>
#include<string.h>
#include"FM.h"
#include <time.h>
#include<ctime>
#include<fstream>
#include<iostream>
using namespace std;
void usage();
void helpbuild();
void helpload();
void helpsave();
void helpcount();
void helplocate();
void splitcommand(string command,string result[]);
void showpos( int * pos,int num);
int main(){
FM *fm=new FM("/home/xiaofei/testfiles/test");
	FILE *fp=fopen("/home/xiaofei/testfiles/test","r+");
	fseeko(fp, 0, SEEK_END);
	int n = ftello(fp) + 1;
	unsigned char *T = new unsigned char[n];
	fseeko(fp, 0, SEEK_SET);
	long long e = 0;
	int num = fread(T, sizeof(uchar), n, fp);
	T[n - 1] = 0;
	fclose(fp);
	int k=0;
	int kk=0;
	for(int i=0;i<n-10;i++){
		auto temp=fm->extracting(i,10);
		for(int jj=9;jj>=0;jj--){
			if((int)T[i+jj]!=(int)temp[jj]){
				cout<<jj<<"  "<<temp<<"  "<<T[i+jj]<<endl;
				k++;
			}
		}
		
	}
	cout<<k<<"  "<<n<<endl;
	cout<<"#########################"<<endl;
	auto pat="zxcv";
	int D;
	auto pos=fm->locating("zxcv",D);
	cout<<D<<endl;
	k=0;
	for(int i=0;i<D;i++){
		auto temp=fm->extracting(pos[i],4);
		for(int j=0;j<4;j++){
			if(pat[j]!=temp[j]){
				k++;
				//cout<<pos[i]<<"  "<<temp<<endl;
				break;
			}
		}
	}
	cout<<k<<endl;
}
int mainq(int argc, char* argv[])
{
	usage();
	string command;
	string  result[3];
	char filename[100]={'\0'};
	char indexname[100]={'\0'};
	FM *csa=NULL;
	while(1)
	{
		result[0]="";//the choice of funtion items
		result[1]="";//filename of data file to build function or pattern to count function
		command="";
		cout<<">";
		getline(cin,command);
		command+=" ";
		splitcommand(command,result);
		if(result[0]=="quit")
		{
			break;
		}	
		else if(result[0]=="build")
		{
			if(csa!=NULL)
				delete csa;
			csa=NULL;
			clock_t s = clock();
			int speed = atoi(result[2].c_str());
			csa=new FM(result[1].data(), speed);
			clock_t e = clock();
			cout<<"build time(s):" << (e-s)*1.0/CLOCKS_PER_SEC<<endl;
		}
		else if(result[0]=="count")
		{
			int  num=0;
			if(csa!=NULL)
			{
				csa->counting(result[1].data(),num);
				cout<<"occs: "<<num<<endl;
			}
			else
			{
				cout<<"build a FM first"<<endl;
			}
		}
		else if(result[0]=="locate")
		{
			int * pos;
			int num=0;
			if(csa!=NULL)
			{
				pos=csa->locating(result[1].data(),num);
				showpos(pos,num);
				delete [] pos;
			}
			else
				cout<<"build a FM first"<<endl;
		}
		else if(result[0]=="load")
		{
			if(csa!=NULL)
				delete csa;
			csa=new FM();
			csa->load(result[1].data());
		}
		else if(result[0]=="save")
		{
			csa->save(result[1].data());
		}
		else if(result[0]=="size")
		{
			if(csa!=NULL)
				cout<<csa->compressRatioForCount()<<endl;
			else
				cout<<"build a FM first"<<endl;
		}
		else if(result[0]=="help")
		{
			if(result[1]=="build")
				helpbuild();
			if(result[1]=="count")
				helpcount();
			if(result[1]=="locate")
				helplocate();
			if(result[1]=="load")
				helpload();
			if(result[1]=="save")
				helpsave();
		}
		else
			usage();
	}
	return 0;
}

void showpos(int * pos,int num)
{
	cout<<"occs:"<<num<<endl;
	for(long int i=0;i<num;i++)
	{
		cout<<pos[i]<<endl;
		if((i+1)%20==0)
		{
			char command;
			cout<<"-----------------more---------------------";
			system("stty raw");
			command=getchar();
     		cout<<endl<<'\r';
			system("stty cooked");
			if(command==27)
			{
				cout<<endl;
				return ;
			}
		}
	}
}


void splitcommand(string command,string result[])
{

	int i=0;
	int start=0;
	int len=command.length();
	int j = 0;
	result[0]=command;
	for(i=0;i<len;i++)
	{
        if(command[i] == ' ')
        {
          result[j++] = command.substr(start,i - start);
          start = i + 1;
         }
	}
}

void usage()
{
	cout<<"--------------------------------------------------------------------------------------"<<endl;
	cout<<"The flowing commands are supported "<<endl;
	cout<<"	help XX: show the details for the command XX"<<endl;
	cout<<"	build XX Y:build the index of file XX with speedlevel Y"<<endl;
	cout<<"	load  XX: load the index file XX"<<endl;
	cout<<"	save XX: write the csa to a index file XX"<<endl;
	cout<<"	count XX: count the pattern XX's occs" <<endl;
	cout<<"	locate YY: enum eyery  position of the pattern"<<endl;
	cout<<"	size: the size of the csa"<<endl;
	cout<<"	quit: say goodbye"<<endl;

}
void helpbuild()
{
	cout<<"build XX"<<endl;
	cout<<"	XX:the source file,it's you responsibility to provide a correct path"<<endl;
}
void helpcount()
{
	cout<<"count XX"<<endl;
	cout<<"	XX:the pattern.you have the responsibility to ensure the index csa is nearby,otherwise,nothing you will get"<<endl;
}
void helplocate()
{
	cout<<"locate XX"<<endl;
	cout<<"	XX: the pattern.it's like count, and you can scan the positions of all the occs in a manner like the system command -more-"<<endl;
}
void helpload()
{
	cout<<"load XX"<<endl;
	cout<<"	XX: the FM-index file, the command will read the index file and build a csa secretly"<<endl;
}
void helpsave()
{
	cout<<"save XX"<<endl;
	cout<<"	XX: the FM-index file, the command will save the csa in file XX"<<endl;
}
