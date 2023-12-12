
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <bitset>
#include <vector>

using namespace std;

string *covInfo(const string &data, int col)
{
    int len = data.length();
    string *ss = new string[col];
    int i = 0, num = 0;
    while (i < len)
    {
        char ch = data.at(i);
        if (ch == ' ')
            num++;
        else
            ss[num].push_back(ch);
        i++;
        if (num == col)
            break;
    }
    return ss;
}

int main(int args, char **argv)
{
    cout<<"hello"<<endl;
    string chromes[]={"1.fa","2.fa","3.fa","4.fa","5.fa","6.fa","7.fa","8.fa","9.fa","10.fa","11.fa","12.fa","13.fa","14.fa","15.fa","16.fa","17.fa","18.fa","19.fa","20.fa","21.fa","22.fa","X.fa","Y.fa"};
    string chrnum[]={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
    // string fin = argv[1];
    string fin = "/home/lab/xiaofei/formatsnp/genome.fa";
    string path="/home/lab/xiaofei/formatsnp/chr/";
    ifstream input(fin);
    ofstream output;
    if (input.fail())
    {
        cout << "Fille does not exist!" << endl;
        return 0;
    }
    string data;
    string fout;
    while(getline(input, data)){
        if(data.at(0)=='>'){
            if(output.is_open()){
                output.close();
            }
            string *ss;
            ss = covInfo(data, 2);
            fout=path+"chr_"+ss[0].substr(1)+".fa";
            cout<<fout<<endl;
            output.open(fout);
        }
        output<<data<<endl;
    }

    input.close();
    return 0;
}