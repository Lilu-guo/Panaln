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
#include <iomanip>

using namespace std;

string *covInfo(const string &data, int col, char d)
{
    int len = data.length();
    string *ss = new string[col];
    int i = 0, num = 0;
    while (i < len)
    {
        char ch = data.at(i);
        if (ch == d)
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
    string fnvcf = argv[1];
    string fnsnp = argv[2];
    ifstream in(fnvcf);
    ofstream out(fnsnp);
    if (in.fail())
    {
        cout << "Input file does not exist!" << endl;
        return 0;
    }
    if (out.fail())
    {
        cout << "Output file does not exist!" << endl;
        return 0;
    }
    string data, id, type, chr, pos, ref, all, alt;
    while (getline(in, data))
    {
        string *ss;
        ss = covInfo(data, 12, '\t');
        id = ss[4];
        type = ss[11];
        chr = ss[1].substr(3);
        pos = ss[2];
        ref = ss[8];
        all = ss[9];
        string *aa;
        aa = covInfo(all, 2, '/');
        if (type.compare("single") == 0)
        {
            if (aa[0].compare(ref) == 0)
            {
                alt = aa[1];
            }
            else
            {
                alt = aa[0];
            }
        }
        else if (type.compare("insertion") == 0)
        {
            if (aa[0].compare(ref) == 0)
            {
                alt = aa[1];
            }
            else
            {
                alt = aa[0];
            }
        }
        else
        {
            alt = to_string(ref.length());
        }
        out << id << "\t" << type << "\t" << chr << "\t" << pos << "\t" << alt << endl;
    }

    return 0;
}