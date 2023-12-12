
#include"savekit.h"
void savekit::close()
{
	if(w!=NULL)
	{
		w->close();
	}
	w=NULL;
}

savekit::~savekit()
{
	if(w!=NULL){
		w->close();
	}
		
}

savekit::savekit(const char * file)
{
	this->w=new ofstream(file,ofstream::binary);
	if(w==NULL)
	{
		cout<<file<<" save Fopen error"<<endl;
		exit(0);
	}
}

i32 savekit::writei64(i64 value)
{
	w->write((char*)&value,sizeof(i64));
	return 0;
}

i32 savekit::writeu64(u64 value)
{
	w->write((char*)&value,sizeof(u64));
	return 0;
}

i32 savekit::writei32( i32 value)
{
	w->write((char*)&value,sizeof(i32));
	return 0;
}

i32 savekit::writeu32(u32 value)
{
	w->write((char*)&value,sizeof(u32));
	return 0;
}
i32 savekit::writei16(i16 value)
{
	w->write((char*)&value,sizeof(i16));
	return 0;
}
 i32 savekit::writeu16(u16 value)
{
	w->write((char*)&value,sizeof(u16));
	return 0;
}

i32 savekit::writeu8(u8 value)
{
	w->write((char*)&value,sizeof(u8));
	return 0;
}
i32 savekit::writedouble(double value)
{
	w->write((char*)&value,sizeof(double));
	return 0;
}
i32 savekit::writei64array(i64 * value,i32 len)
{
	w->write((char*)value,len*sizeof(i64));
	return 0;
}
i32 savekit::writeu64array(u64 * value,i32 len)
{
	w->write((char*)value,len*sizeof(u64));
	return 0;
}
i32 savekit::writei32array(i32 * value,i32 len)
{
	w->write((char*)value,len*sizeof(i32));
	return 0;
}
i32 savekit::writeu32array(u32* value,i32 len)
{
	w->write((char*)value,len*sizeof(u32));
	return 0;
}
i32 savekit::writei16array(i16 * value,i32 len)
{
	w->write((char*)value,len*sizeof(i16));
	return 0;
}
i32 savekit::writeu16array(u16 * value,i32 len)
{
	w->write((char*)value,len*sizeof(u16));
	return 0;
}

i32 savekit::writeu8array(u8 * value,i32 len)
{
	w->write((char*)value,len*sizeof(u8));
	return 0;
}
