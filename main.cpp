#include "hfHandle.h"

int main()
{
	hfHandle hf;
	std::cout << hf.initCalc() << '\n';
	for(int i=0; i<10; i++)
		std::cout << hf.scf() << '\n';
}
