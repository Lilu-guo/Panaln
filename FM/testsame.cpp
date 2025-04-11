#include<iostream>
#include<fstream>
#include<cstdlib>
#include<unistd.h>
using namespace std;
int main()
{
    FILE *f1,*f2;
    f1 = fopen("AAAqwe","r");
    f2 = fopen("../testfile/dna","r");
    unsigned char *c1,*c2;
    c1 = new unsigned char[403927746 + 10];
    c2 = new unsigned char[403927746 + 10];
    fread(c1,1,403927746,f1);
    fread(c2,1,403927746,f2);
    for(int i = 0;i < 403927746;++i)
    {
        if(c1[i] != c2[i])
        {
            cout<<i<<endl;
            cout<<c1[i]<<"  "<<c2[i]<<endl;
            exit(0);
        }
    }
    fclose(f1);
    fclose(f2);
    cout<<"yes!!!"<<endl;
    /*
    string fn1("AAAqwe"), fn2("../testfile/dna");
    ifstream fin1, fin2;
    fin1.open(fn1.c_str());
    fin2.open(fn2.c_str());
    string tmp1(""), tmp2("");
    long long i = 0;
    while(getline(fin1,tmp1))
    {
        if(tmp1 == "0") 
            cout<<"is 0"<<endl;
    }
    
    while(getline(fin1,tmp1) && getline(fin2,tmp2))
    {
        ++i;
//        if(i >= 1676455)
        if(tmp1 != tmp2)
        {
            cout<<tmp1<<"   "<<tmp2<<endl;
            cout<<i<<endl;
//            exit(0);
            sleep(1);
        }
    }
    
    fin1.close();
    fin2.close();
    */
    return 0;
}