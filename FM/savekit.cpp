/*============================================
# Filename: savekit.cpp
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#include"savekit.h"
void savekit::close()
{
	if(w!=NULL)
	{
		// fclose(w);
		w->close();
	}
	w=NULL;
}

savekit::~savekit()
{
	if(w!=NULL){
		// fclose(w);
		w->close();
	}
		
}

savekit::savekit(const char * file)
{
	// this->w=fopen(file,"wb");
	this->w=new ofstream(file,ofstream::binary);
	if(w==NULL)
	{
		cout<<file<<" save Fopen error"<<endl;
		exit(0);
	}
}

i32 savekit::writei64(i64 value)
{
	// fwrite(&value,sizeof(i64),1,w);
	w->write((char*)&value,sizeof(i64));
	return 0;
}

i32 savekit::writeu64(u64 value)
{
	// fwrite(&value,sizeof(u64),1,w);
	w->write((char*)&value,sizeof(u64));
	return 0;
}

i32 savekit::writei32( i32 value)
{
	// fwrite(&value,sizeof( i32),1,w);
	w->write((char*)&value,sizeof(i32));
	return 0;
}
 i32 savekit::writeu32(u32 value)
{
	// fwrite(&value,sizeof(u32),1,w);
	w->write((char*)&value,sizeof(u32));
	return 0;
}
 i32 savekit::writei16(i16 value)
{
	// fwrite(&value,sizeof(i16),1,w);
	w->write((char*)&value,sizeof(i16));
	return 0;
}
 i32 savekit::writeu16(u16 value)
{
	// fwrite(&value,sizeof(u16),1,w);
	w->write((char*)&value,sizeof(u16));
	return 0;
}

i32 savekit::writeu8(u8 value)
{
	// fwrite(&value,sizeof(u8),1,w);
	w->write((char*)&value,sizeof(u8));
	return 0;
}
i32 savekit::writedouble(double value)
{
	// fwrite(&value,sizeof(u8),1,w);
	w->write((char*)&value,sizeof(double));
	return 0;
}
 i32 savekit::writei64array(i64 * value,i32 len)
{
	// fwrite(value,sizeof(i64),len,w);
	w->write((char*)value,len*sizeof(i64));
	return 0;
}
 i32 savekit::writeu64array(u64 * value,i32 len)
{
	// fwrite(value,sizeof(u64),len,w);
	w->write((char*)value,len*sizeof(u64));
	return 0;
}
 i32 savekit::writei32array(i32 * value,i32 len)
{
	// fwrite(value,sizeof(i32),len,w);
	w->write((char*)value,len*sizeof(i32));
	return 0;
}
 i32 savekit::writeu32array(u32* value,i32 len)
{
	// fwrite(value,sizeof(u32),len,w);
	w->write((char*)value,len*sizeof(u32));
	return 0;
}
 i32 savekit::writei16array(i16 * value,i32 len)
{
	// fwrite(value,sizeof(i16),len,w);
	w->write((char*)value,len*sizeof(i16));
	return 0;
}
 i32 savekit::writeu16array(u16 * value,i32 len)
{
	// fwrite(value,sizeof(u16),len,w);
	w->write((char*)value,len*sizeof(u16));
	return 0;
}

i32 savekit::writeu8array(u8 * value,i32 len)
{
	// fwrite(value,sizeof(u8),len,w);
	w->write((char*)value,len*sizeof(u8));
	return 0;
}
















