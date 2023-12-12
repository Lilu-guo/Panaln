
#include <stdlib.h>
#include <string.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

bool isDNA(string dna)
{
    if (dna[0] == '>')
    {
        return false;
    }
    else
    {
        return true;
    }
}

void splitfile()
{
    ifstream fa("/media/dell198/4fe135dd-6ee8-4bea-8ccd-f450dc5a7100/home/lab/xiaofei/formatsnp/genome.fa");
    string path = "/media/dell198/4fe135dd-6ee8-4bea-8ccd-f450dc5a7100/home/lab/xiaofei/formatsnp/";
    string refile = "";
    string chr = "";
    string temp;
    if (!fa)
    {
        cout << "open error" << endl;
        return;
    }
    getline(fa, temp);
    refile = path + temp.substr(1, 2) + ".fa";
    while (getline(fa, temp))
    {
        if (isDNA(temp))
        {
            chr.append(temp);
        }
        else
        {
            cout << refile << endl;
            ofstream out(refile.c_str());
            out << chr << endl;
            chr = "";
            refile = path + temp.substr(1, 2) + ".fa";
        }
    }
    cout << refile << endl;
    ofstream out(refile.c_str());
    out << chr << endl;
}

void splitstr(string str, vector<string> &res, char ch)
{
    size_t pre = 0;
    res.clear();
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == ch)
        {
            res.push_back(str.substr(pre, i - pre));
            pre = i + 1;
        }
    }
    res.push_back(str.substr(pre));
}

void getfieldnumbers()
{
    ifstream fa("/media/dell198/4fe135dd-6ee8-4bea-8ccd-f450dc5a7100/home/lab/xiaofei/formatsnp/genome.snp");
    vector<string> field;
    vector<string> tempstrs;
    string temp;
    int k = 0;
    while (getline(fa, temp))
    {
        splitstr(temp, tempstrs, '\t');
        if (find(field.begin(), field.end(), tempstrs[1]) == field.end())
        {
            field.push_back(tempstrs[1]);
        }
    }
    for (int i = 0; i < field.size(); i++)
    {
        cout << field[i] << endl;
    }
    cout << field.size() << endl;
}

int main(int argc, char *argv[])
{
    string path = "/home/lab/gll/formatsnp/";
    string filename;
    string refilename;
    string chromes[] = {"1.fa", "2.fa", "3.fa", "4.fa", "5.fa", "6.fa", "7.fa", "8.fa", "9.fa", "10.fa", "11.fa", "12.fa", "13.fa", "14.fa", "15.fa", "16.fa", "17.fa", "18.fa", "19.fa", "20.fa", "21.fa", "22.fa", "X.fa", "Y.fa", "MT.fa"};
    string chrnum[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"};
    string snpname[] = {"insertion", "single", "deletion"};

    for (int i = 0; i < 25; i++)
    {
        filename = path + chromes[i];
        refilename = path + chromes[i] + ".vcf";
        ifstream infile(filename.c_str());
        string chr = "";
        infile >> chr;
        cout << chr.size() << endl;
        ofstream outfile(refilename.c_str());
        ifstream snp("/home/lab/gll/formatsnp/snp144Common.snp");
        string tempsnp;
        while (getline(snp, tempsnp))
        {
            vector<string> tempsnps;
            splitstr(tempsnp, tempsnps, '\t');
            if (tempsnps[2].compare(chrnum[i]) == 0)
            {
                if (tempsnps[1].compare("insertion") == 0)
                {
                    outfile << tempsnps[2] << "\t" << tempsnps[3] << "\t" << tempsnps[0] << "\t" << chr[atoi(tempsnps[3].c_str())] << "\t" << chr[atoi(tempsnps[3].c_str())] + tempsnps[4] << endl;
                }
                else if (tempsnps[1].compare("single") == 0)
                {
                    outfile << tempsnps[2] << "\t" << tempsnps[3] << "\t" << tempsnps[0] << "\t" << chr[atoi(tempsnps[3].c_str())] << "\t" << tempsnps[4] << endl;
                }
                else if (tempsnps[1].compare("deletion") == 0)
                {
                    outfile << tempsnps[2] << "\t" << atoi(tempsnps[3].c_str()) - 1 << "\t" << tempsnps[0] << "\t" << chr.substr(atoi(tempsnps[3].c_str()) - 1, atoi(tempsnps[4].c_str()) + 1) << "\t" << chr[atoi(tempsnps[3].c_str()) - 1] << endl;
                }
            }
        }
    }
    return 0;
}