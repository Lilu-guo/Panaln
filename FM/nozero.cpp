/*
 * @Author: your name
 * @Date: 2021-06-07 21:16:35
 * @LastEditTime: 2021-06-07 22:05:41
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /FM-Adaptive/nozero.cpp
 */
#include<iostream>
#include<fstream>
using namespace std;
int main(void)
{
    ifstream fin;
    fin.open("/home/lab/lab/hzt/bt/coha_txt");
    ofstream fout("/home/lab/lab/hzt/bt/coha_txt.nozero");
    char ch;
    while(fin.get(ch))
    {
        if(ch == 0)continue;
        fout<<ch;
    }
    fin.close();
    fout.close();
    return 0;
}