#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include <sstream>
#include "hfHandle.h"

hfHandle::hfHandle()
{
	
}

hfHandle::hfHandle(std::string fname)
{
	init(fname);
}

void hfHandle::init(std::string fname)
{
	std::ifstream intgrl;
	intgrl.open(fname.c_str());
	
	double val;
	int i, j, k, l;
	
	intgrl >> E_nuc >> no >> no_e >> acc;
	
	if(no_e/2 > no)
	{
		std::cout << "Error: No. of atoms and electons mismatch!";
	}

	intgrl >> i >> j >> val;

	std::string str;
	
	overlap = MatrixXd::Constant(no, no, 0);
		
	Hcore = MatrixXd::Constant(no, no, 0);
	
	while(getline(intgrl, str))
	{
		if (str == "*") break;
	
		std::istringstream ss(str);
			
		ss >> i >> j >> val;
		
		overlap(i-1, j-1) = val;
		
		if(i != j) overlap(j-1, i-1) = val;
	}
		
	while(getline(intgrl, str))
	{
		if (str == "*") break;
		
		std::istringstream ss(str);
			
		ss >> i >> j >> val;
			
		Hcore(i-1, j-1) = val;
			
		if(i != j) Hcore(j-1, i-1) = val;
	}
	
	s_no = SUM(no);
	two_e = MatrixXd::Constant(SUM(s_no), 1, 0);
	
	while(getline(intgrl, str))
	{
		if (str == "*") break;
		
		std::istringstream ss(str);
			
		ss >> i >> j >> k >> l >> val;
		
		two_e(AtomToIndex(i,j,k,l)) = val;
	}	
		
	SelfAdjointEigenSolver<MatrixXd> es_ov(overlap);
		
	S_inv_sqrt = es_ov.operatorInverseSqrt();
	
	fock = S_inv_sqrt.transpose()*Hcore*S_inv_sqrt;
	
	SelfAdjointEigenSolver<MatrixXd> es_fock(fock, ComputeEigenvectors); 
		
	C = S_inv_sqrt*es_fock.eigenvectors();
		
	D = MatrixXd::Constant(no, no, 0);
	
	for(int i=0; i<no; i++)
	
		for(int j=0; j<no; j++)
		{
			val = 0;
			
			for(int m=0; m<no_e/2; m++) val += C(i, m)*C(j, m);
			
			D(i, j) = val;
		}
	
	E_elec = 0;
	
	for(int i=0;i<no; i++)
		for(int j=0; j<no; j++)
			E_elec += D(i,j)*(Hcore(i,j) + fock(i,j));
			
	Eprev = E_elec;
	
	std::cout << E_elec << '\n';
}

void hfHandle::coreSCF()
{
	for(int i=0; i < no; i++)
		for(int j=0; j < no; j++) 
		{
			fock(i, j) = Hcore(i,j);
  	
   			for(int k=0; k < no; k++)
   				for(int l=0; l < no; l++) 
		  			fock(i,j) += D(k,l) * ( 2.0 * two_e( AtomToIndex(i+1,j+1,k+1,l+1) ) - two_e( AtomToIndex(i+1,k+1,j+1,l+1) ) );
			  
  		}
  	
  	double val;
	
	SelfAdjointEigenSolver<MatrixXd> es_fock(fock, ComputeEigenvectors); 
		
	C = S_inv_sqrt*es_fock.eigenvectors();
		
	D = MatrixXd::Constant(no, no, 0);
	
	for(int i=0; i<no; i++)
	
		for(int j=0; j<no; j++)
		{
			val = 0;
			
			for(int m=0; m<no_e/2; m++) val += C(i, m)*C(j, m);
			
			D(i, j) = val;
		}
	
	E_elec = 0;
	
	for(int i=0;i<no; i++)
		for(int j=0; j<no; j++)
			E_elec += D(i,j)*(Hcore(i,j) + fock(i,j));
}

void hfHandle::work()
{
	do
	{
		Eprev = E_elec;
		
		coreSCF();
		
		std::cout << E_elec << " " << E_elec-Eprev << '\n';
	}while(!converged());
}


