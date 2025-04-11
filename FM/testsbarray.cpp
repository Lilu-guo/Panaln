#include<iostream>
#include<cstdlib>
#include"InArray.h"
int main()
{
    SBArray *q = new SBArray((1 << 21) + 4,33,16,12); 
    long long cnt = 0;
    long long len = (1ull << 33) + 1;
    int step2 = 256, step1 = 256 * 16;
    q->SetValue(0,0,1);
    q->SetValue(0,0,0);
    u64 *sbarray = new u64[(1 << 21) + 1];
    u64 *barray = new u64[(1 << 25 ) + 1];
    sbarray[0] = 0;
    barray[0] = 0;
    for(int i = 1;i < (1 << 21) + 1;++i)
        sbarray[i] = ((1ull << 32) + i);
    for(int i = 1;i < (1 << 25 ) + 1;++i)
        barray[i] = (u32)rand() % 3840;
    i64 sbi = 1, bi = 1, i = 1;
    for(;i < len;++i)
    {
        if(i % step1 == 0)
            q->SetValue(i / step1,sbarray[sbi++],1);
        if(i % step2 == 0)
            q->SetValue(i / step2,barray[bi++],0);
    }
cout<<"construct is ok!!!"<<endl;
    sbi = 0; bi = 0;
    for(i = 0;i < len;++i)
    {
        if(i % step1 == 0)
        {
            u64 tmp = q->GetValue(i / step1,1);
            if(tmp != sbarray[sbi++])
            {
                cout<<"NO!!!"<<endl;
                cout<<tmp<<"    "<<sbarray[sbi - 1];
                exit(0);
            }
        }
        if(i % step2 == 0)
        {
            u64 tmp = q->GetValue(i / step2,0);
            if(tmp != barray[bi++])
            {
                cout<<"NO!!!"<<endl;
                cout<<tmp<<"    "<<barray[bi - 1];
                exit(0);
            }
        }
    }
    delete q;
    delete sbarray;
    delete barray;
    cout<<"is ok!!!"<<endl;
    return 0;
}