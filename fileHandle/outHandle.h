#include <fstream>
#include <iostream>
#include <string>
#include "filehandle.h"

using namespace fileHandle;

outHandle& outHandle::operator<<(const std::string mesg )
{
    std::cout << mesg << " \t ";
    f << mesg << " \t ";
    return *this;
}

outHandle& outHandle::operator<<(const char mesg )
{
    if(mesg != endl)
	{	
		std::cout << mesg << " \t ";
    	f << mesg << " \t ";
	}else if(mesg == endl)
	{
		std::cout << '\n';
    	f << mesg << '\n';
	}
    return *this;
}

outHandle& outHandle::operator<<(const int mesg )
{
    std::cout << mesg << " \t ";
    f << mesg << " \t ";
    return *this;
}

outHandle& outHandle::operator<<(const double mesg )
{
    std::cout << mesg << " \t ";
    f << mesg << " \t ";
    return *this;
}

void outHandle::open(std::string fname)
{
	f.open(fname.c_str());
}

outHandle::outHandle() {}

outHandle::outHandle(std::string fname) { open(fname); }

outHandle::~outHandle() { f.close(); }
