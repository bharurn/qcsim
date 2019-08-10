#ifndef hfHandle_H
#define hfHandle_H

#include <Eigen/Dense>
#include <cmath>

#define SUM(N) N*(N+1)/2

using namespace Eigen;

class hfHandle
{
	public:
	
		void work();
		void init(std::string fname);
		hfHandle();
		hfHandle(std::string fname);
		
	private: 
	
		double E_nuc, E_elec, Eprev, acc;
		MatrixXd overlap, Hcore, two_e, S_inv_sqrt, fock, C, D, Dprev;
		int no, no_e,s_no;
	
		int AtomToIndex(int i, int j, int k, int l)
		{
			k--;
			i--;
		
			int col = SUM(k) + l;
			int row = SUM(i) + j - col + 1;
			int sum = (col-1)*0.5*(2*s_no - col+ 2);
		
			return row+sum-1;
		}
		
		double getEnergy(bool total)
		{
			if(total) return E_elec+E_nuc;
			else return E_elec;
		}
		
		bool converged()
		{
			if(fabs(E_elec - Eprev) <= acc) return true;
			else return false;
		}
		
		void coreSCF();
	
};

#endif
