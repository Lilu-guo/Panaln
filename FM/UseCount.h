/*============================================
# Filename: UseCount.h
# Ver 1.0 2021-09-08
# Copyright (C) Zongtao He, Pengfei Liu, Hongbojiang, Hongwei Huo, and Jeffrey S. Vitter.
#
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 or later of the License.
#
# Description: 
=============================================*/
#ifndef USECOUNT_H
#define USECOUNT_H
class UseCount
{
	public:
		UseCount();
		UseCount(const UseCount &);
		~UseCount();
		bool only();
		bool reattach(const UseCount &);
	private:
		int *p;
		UseCount & operator = (const UseCount &);
};
#endif
