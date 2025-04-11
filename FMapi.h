#ifndef _FM_API_H__
#define _FM_API_H__
#ifdef __cplusplus
extern "C" {
#endif
struct FMA; 
extern struct FMA *GetFMapi(char * fn); 
extern void ReleaseFMapi(struct FMA *ppInstance);
extern struct FMA *FmBuild_index(const char * fn);
void printinformation(struct FMA *ppInstance);
extern void Fm_save(const char * fn, struct FMA *ppInstance);
extern void Fm_load(const char * fn, struct FMA *ppInstance);
extern void my_load_idx(struct FMA *ppInstance);                                                                    
extern void Fm_occ(long long l, long long r,unsigned char c,long long *Left,long long *Right,struct FMA *ppInstance);
extern void Fm_16occ(long long l, long long r,long long *Left,long long *Right,struct FMA *ppInstance);
extern void Occ1p(unsigned char c, long long pos1, long long pos2, long long *occ1, long long *occ2, struct FMA *ppInstance);   
extern void Occ8p(unsigned char c, long long pos1, long long pos2, long long *occ81, long long *occ82, struct FMA *ppInstance); 
extern void Occ16p(long long pos1, long long pos2, long long *occ161, long long *occ162, struct FMA *ppInstance);               
extern void Fm_GetN(struct FMA *ppInstance, long long *n);
extern long long Fm_Lookup(struct FMA *ppInstance, long long pos);
extern unsigned char* Fm_Extract(struct FMA *ppInstance, long long pos, long long len);
extern long long* Fm_getC(struct FMA *ppInstance);
extern long long Fm_getprimary(struct FMA *ppInstance);
#ifdef __cplusplus
};
#endif
#endif