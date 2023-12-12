
#include"loadkit.h"
void loadkit::close()
{
	if(r!=NULL)
		r->close();
	r=NULL;
}

loadkit::~loadkit()
{
	if(r!=NULL)
		r->close();
}

loadkit::loadkit(const char * file)
{
	this->r=new ifstream(file,ifstream::binary);
	if(r==NULL)
	{
		cout<<"Fopen error"<<endl;
		exit(0);
	}
}

i32 loadkit::loadi64(i64 & value)
{
	r->read((char*)&value,sizeof(i64));
	return r->gcount()/sizeof(i64);
}
i32 loadkit::loadu64(u64 & value)
{
	r->read((char*)&value,sizeof(u64));
	return r->gcount()/sizeof(u64);
}
i32 loadkit::loadi32(i32 & value)
{
	r->read((char*)&value,sizeof(i32));
	return r->gcount()/sizeof(i32);
}
i32 loadkit::loadu32(u32 & value)
{
	r->read((char*)&value,sizeof(u32));
	return r->gcount()/sizeof(u32);
}
i32 loadkit::loadi16(i16 & value)
{
	r->read((char*)&value,sizeof(i16));
	return r->gcount()/sizeof(i16);
}
i32 loadkit::loadu16(u16 & value)
{
	r->read((char*)&value,sizeof(u16));
	return r->gcount()/sizeof(u16);
}
i32 loadkit::loadu8(u8 & value)
{
	r->read((char*)&value,sizeof(u8));
	return r->gcount()/sizeof(u8);
}
i32 loadkit::loaddouble(double & value){
	r->read((char*)&value,sizeof(double));
	return r->gcount()/sizeof(double);
}
i32 loadkit::loadi64array(i64 * value,i32 len)
{
	r->read((char*)value,len*sizeof(i64));
	return r->gcount()/sizeof(i64);
}
i32 loadkit::loadu64array(u64 * value,i32 len)
{
	r->read((char*)value,len*sizeof(u64));
	return r->gcount()/sizeof(u64);
}
i32 loadkit::loadi32array(i32 * value,i32 len)
{
	r->read((char*)value,len*sizeof(i32));
	return r->gcount()/sizeof(i32);
}
i32 loadkit::loadu32array(u32 * value,i32 len)
{
	r->read((char*)value,len*sizeof(u32));
	return r->gcount()/sizeof(u32);
}
i32 loadkit::loadi16array(i16 * value,i32 len)
{
	r->read((char*)value,len*sizeof(i16));
	return r->gcount()/sizeof(i16);
}
i32 loadkit::loadu16array(u16 * value,i32 len)
{
	r->read((char*)value,len*sizeof(u16));
	return r->gcount()/sizeof(u16);
}
i32 loadkit::loadu8array(u8 * value,i32 len)
{
	r->read((char*)value,len*sizeof(u8));
	return r->gcount()/sizeof(u8);
}
