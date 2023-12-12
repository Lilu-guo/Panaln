
#ifndef _SAVEKIT
#define _SAVEKIT
#include<stdio.h>
#include<stdlib.h>
#include"BaseClass.h"
#include<fstream>
#include<iostream>
using namespace std;
class savekit
{
	public:
		savekit(const char * file);
		~savekit();
		i32 writei64(i64 value);
		i32 writeu64(u64 value);
		i32 writei32(i32 value);
	    i32 writeu32(u32 value);
		i32 writei16(i16 value);
		i32 writeu16(u16 value);
		i32 writeu8(u8 value);
		i32 writedouble(double value);
		i32 writei64array(i64 * value,i32 len);
		i32 writeu64array(u64 * value,i32 len);
		i32 writei32array(i32 * value,i32 len);
		i32 writeu32array(u32 * value,i32 len);
		i32 writei16array(i16 * value,i32 len);
		i32 writeu16array(u16 * value,i32 len);
		i32 writeu8array(u8 * value,i32 len);
	    void close();
	private:
		// FILE *w;
		ofstream *w;
};
#endif
