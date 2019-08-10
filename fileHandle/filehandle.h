#ifndef FILEHANDLE_H
#define FILEHANDLE_H

#include <fstream>
#include <iostream>
#include <string>

namespace fileHandle
{
	const char endl = '\n';
	
	class outHandle
	{
		public:
			void open(std::string fname);
			outHandle();
			outHandle(std::string fname);
			~outHandle();
			
			outHandle& operator<<(const std::string mesg);
			outHandle& operator<<(const char mesg);
			outHandle& operator<<(const int mesg);
			outHandle& operator<<(const double mesg);
			outHandle& operator%(const std::string mesg);		
		private:
			std::ofstream f;		
	};
}


#endif
