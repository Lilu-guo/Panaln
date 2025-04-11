#include<iostream>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
using namespace std;
typedef long long fm_int;
int main()
{
    /*
    FILE *fp;
    fp = fopen("../GRCh38.fna","r");
    fseek(fp,0,SEEK_END);
	long long n = ftell(fp);
	unsigned char * T = new unsigned char[n];
	fseek(fp,0,SEEK_SET);
	fm_int e=0;
	fm_int num=0;
	while((e=fread(T+num,sizeof(unsigned char),n - num,fp))!=0)
		num = num +e;
	if(num!=n)
	{
		cout<<num<<"  "<<n<<"  Read source file failed"<<endl;
		exit(0);
	}
	fclose(fp);
    */
    string fn1("../GRCh38.fna");
    string fn2("./G38");
    ifstream fin;
    ofstream fout;
    fin.open(fn1.c_str(),ios_base::in);
    fout.open(fn2.c_str(),ios_base::out | ios_base::app);
    string tmp("");
    fm_int idx = 0;
    while(getline(fin,tmp))
    {
        if(tmp[0] == '>') continue;
        else
            fout<<tmp;
    }
    fin.close();
    fout.close();
    return 0;
}