
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
#include <unistd.h>

using namespace std;

string *covInfo(const string &data, int col, char delimiter)
{
    int len = data.length();
    string *ss = new string[col];
    int i = 0, num = 0;
    while (i < len)
    {
        char ch = data.at(i);
        if (ch == delimiter)
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
    cout << "hello" << endl;
    string path = "/home/lab/gll/formatsnp/";
    string filename;
    string snpfilename;
    string indelfilename;

    string chromes[] = {"1.fa", "2.fa", "3.fa", "4.fa", "5.fa", "6.fa", "7.fa", "8.fa", "9.fa", "10.fa", "11.fa", "12.fa", "13.fa", "14.fa", "15.fa", "16.fa", "17.fa", "18.fa", "19.fa", "20.fa", "21.fa", "22.fa", "X.fa", "Y.fa", "MT.fa"};
    string chrnum[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"};
    string snpname[] = {"insertion", "single", "deletion"};

    for (int i = 0; i < 25; i++)
    {
        filename = path + chromes[i] + ".vcf";
        snpfilename = path + "/panVcf/SNP.extract.chr" + chrnum[i] + ".data";
        indelfilename = path + "/panVcf/INDEL.extract.chr" + chrnum[i] + ".data";
        ifstream file(filename.c_str());
        ofstream snp(snpfilename.c_str());
        ofstream indel(indelfilename.c_str());
        int num = 0, snpNum = 0, indelNum = 0;
        string temp;
        while (getline(file, temp))
        {
            num++;
            string *ss = covInfo(temp, 5, '\t');
            if (ss[3].length() == ss[4].length())
            {
                snp << atoi(ss[1].c_str()) << "\t" << ss[3] << "\t" << ss[4] << "\t" << 33 << endl;
                snpNum++;
            }
            else
            {
                indel << atoi(ss[1].c_str()) << "\t" << ss[3] << "\t" << ss[4] << "\t" << 33 << endl;
                indelNum++;
            }
        }
        cout << "chr " << chrnum[i] << " num:" << num << " snp:" << snpNum << " indel:" << indelNum;
        if (num == snpNum + indelNum)
            cout << "  Pass!" << endl;
        else
            cout << "  Error!" << endl;
    }
    return 0;
}