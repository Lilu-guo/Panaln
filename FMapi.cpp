#include"FMapi.h"
#include"FM/FM.h"
#ifdef __cplusplus                                   
extern "C" {
#endif
struct FMA                                           
{
    FM *fm = NULL;
};
struct FMA *GetFMapi(char * fn)
{
    struct FMA *tmp = new struct FMA;                
    tmp->fm = new FM(fn);                              
    return tmp;
}
void ReleaseFMapi(struct FMA *ppInstance)
{
    delete ppInstance->fm;
    ppInstance->fm = 0;
}
struct FMA *FmBuild_index(const char * fn)
{
    struct FMA *tmp = new struct FMA;
    tmp->fm = new FM(fn,256,16,32,1);                 
    return tmp;
}
void printinformation(struct FMA *ppInstance){
	cout<<"the Number of characters:        "<<ppInstance->fm->getN()<<endl;
	cout<<"the Alphabetsize:                "<<ppInstance->fm->getAlphabetSize()<<endl;
	cout<<"the Runs is:                     "<<ppInstance->fm->getRuns()<<endl;
}
void Fm_save(const char * fn, struct FMA *ppInstance)
{
    printinformation(ppInstance);
    ppInstance->fm->save(fn);
}
void Fm_load(const char * fn, struct FMA *ppInstance)                   
{
    ppInstance->fm->load(fn);
}
void my_load_idx(struct FMA *ppInstance)
{
    ppInstance->fm->loadIdx();
}
void Fm_occ(long long l, long long r,unsigned char c ,long long *L,long long* R,struct FMA *ppInstance)
{   
    if(l == 1 && r == 0) {*L = 1;*R = 0; return ;}
    ppInstance->fm->Fmocc(l, r ,c, L, R);
}
void Fm_16occ(long long l, long long r,long long *Left,long long *Right,struct FMA *ppInstance)
{
    if(l == 1 && r == 0)                                                                                      
    {
        for(int i = 0;i < 16;++i)
        {
            Left[i] = 1;
            Right[i] = 0;
        }
        return;                                                                                             
    }
    else
    {
        ppInstance->fm->Fm16occ(l, r, Left, Right);                                                          
    }
}
void Occ1p(unsigned char c, fm_int pos1, fm_int pos2, fm_int *occ1, fm_int *occ2, struct FMA *ppInstance)    
{   
    ppInstance->fm->Occ1p(c, pos1, pos2, *occ1, *occ2);
}
void Occ8p(unsigned char c, fm_int pos1, fm_int pos2, fm_int *occ81, fm_int *occ82, struct FMA *ppInstance) 
{
    ppInstance->fm->Occ8p(c, pos1, pos2, occ81, occ82);
}
void Occ16p(fm_int pos1, fm_int pos2, fm_int *occ161, fm_int *occ162, struct FMA *ppInstance)               
{
    ppInstance->fm->Occ16p(pos1, pos2, occ161, occ162);
}
void Fm_GetN(struct FMA *ppInstance, long long *n)
{
    *n = ppInstance->fm->getN();
}
long long Fm_Lookup(struct FMA *ppInstance, long long pos)
{
    return ppInstance->fm->FM_Lookup(pos);
}
unsigned char* Fm_Extract(struct FMA *ppInstance, fm_int pos, fm_int len){
    return ppInstance->fm->extracting(pos,len);
}
#ifdef __cplusplus
};
#endif