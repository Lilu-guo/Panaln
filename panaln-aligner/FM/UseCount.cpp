#include"UseCount.h"
UseCount::UseCount():p(new int(1)){}

UseCount::UseCount(const UseCount &u):p(u.p){++*p;}

UseCount::~UseCount(){if(--*p==0) delete p;}

bool UseCount::only() {return *p==1;}

bool UseCount::reattach(const UseCount & u)
{
	++*u.p;
	if(--*p==0)
	{
		delete p;
		p=u.p;
		return true;
	}
	p=u.p;
	return false;
}
