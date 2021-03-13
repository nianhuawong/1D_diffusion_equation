#include<iostream>
#include<vector>
#include<fstream>

using namespace std;
#include "diffu_solver.h"

int main()
{
	initialize_parameter();

	flow_initialization();

	for (iter = 0; iter < numberOfTimeSteps; ++iter)
	{
		load_qField();

		//time_marching_FTCS();
		time_marching_full_implicit();
		//time_marching_Crank_Nicolson(); 

		boundary_condition();

		compute_residual();

		output_residual();
	}

	output_results();

	return 0;
}

void output_results()
{
	cout << "dumping results..." << endl;
	fstream file;
	file.open("results.dat", ios_base::out );

	file << "TITLE     = \"results\"" << endl;
	file << "VARIABLES = \"x\", \"qField\"" << endl;

	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		file << xCoordinates[iNode] << "\t" << qField_N1[iNode] << endl;
	}	

	file.close();

	cout << "done!" << endl;
}

void output_residual()
{
	if (iter % 10 == 0)
	{
		cout << "\titer " << "\tresidual" << endl;
	}

	cout << "\t" << iter << "\t" << residual << endl;
}

void load_qField()
{
	qField = qField_N1;
}

void compute_residual()
{
	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		residual += (qField_N1[iNode] - qField[iNode]) * (qField_N1[iNode] - qField[iNode]);
	}

	residual = sqrt(residual / numberOfGridPoints);  //L2
}

void time_marching_FTCS()
{	
	//FTCS
	for (int iNode = 1; iNode < numberOfGridPoints - 1; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] + sigma * (qField[iNode+1]- 2.0 * qField[iNode] + qField[iNode-1]);
	}
}

void time_marching_full_implicit()
{
	a = -sigma;
	b = 1.0 + 2.0 * sigma;
	c = -sigma;

	double u0   = qField[0];
	double rhs1 = qField[1];
	double d1   = rhs1 - a * u0;
	if (u0!= 100)
	{
		int kkk = 1;
	}
	double rhsm = qField[numberOfEquations];
	double ump1 = qField[numberOfEquations + 1];
	double dm   = rhsm - c * ump1;
	
	vector<double> djv(numberOfEquations);
	djv[0] = d1;
	djv[numberOfEquations - 1] = dm;

	for (int iEquation = 1; iEquation < numberOfEquations-1; ++iEquation)
	{
		djv[iEquation] = qField[iEquation];
	}
	
	solve_chase_method(djv);
}

//void solve_chase_method(vector<double>& djv)
//{
//	vector< double > AA(numberOfEquations);
//	vector< double > BB(numberOfEquations);
//
//	double d1 = djv[0];
//	double dm = djv[numberOfEquations-1];
//
//	AA[0] = -c / b;
//	BB[0] = d1 / b;
//
//	for (int iEquation = 1; iEquation < numberOfEquations; ++iEquation)
//	{
//		double dj = djv[iEquation];
//
//		double tmp = b + a * AA[iEquation - 1];
//		AA[iEquation] = -c / tmp;
//		BB[iEquation] = (dj - a * BB[iEquation - 1]) / tmp;
//	}
//
//	qField_N1[numberOfEquations] = (dm - a * BB[numberOfEquations - 1]) / (b + a * AA[numberOfEquations - 1]);
//
//	for (int iEquation = numberOfEquations; iEquation >1; --iEquation)
//	{
//		qField_N1[iEquation-1] = AA[iEquation-1] * qField_N1[iEquation] + BB[iEquation-1];
//	}
//}

void solve_chase_method(vector<double>& djv)
{
	vector<double> qField_tmp(numberOfEquations);

	vector< double > AA(numberOfEquations);
	vector< double > BB(numberOfEquations);
	vector< double > CC(numberOfEquations);
	for (int iEquation = 0; iEquation < numberOfEquations; ++iEquation)
	{
		AA[iEquation] = a;
		BB[iEquation] = b;
		CC[iEquation] = c;
	}

	for (int iEquation = 1; iEquation < numberOfEquations; iEquation++)
	{
		double tmp = AA[iEquation] / BB[iEquation - 1];

		BB[iEquation]  = BB[iEquation]  - tmp * CC[iEquation - 1];
		djv[iEquation] = djv[iEquation] - tmp * djv[iEquation - 1];
	}

	qField_tmp[numberOfEquations - 1] = djv[numberOfEquations - 1] / BB[numberOfEquations - 1];

	for (int iEquation = numberOfEquations - 2; iEquation >= 0; iEquation--)
	{
		qField_tmp[iEquation] = ( djv[iEquation] - CC[iEquation]  * qField_tmp[iEquation+1] )/ BB[iEquation];
	}

	for (int iEquation = 0; iEquation < numberOfEquations; ++iEquation)
	{
		qField_N1[iEquation + 1] = qField_tmp[iEquation];
	}
}

////forward elimination
//
//for (i = 1; i < r; i++)
//{
//	m = A[i][i - 1] / A[i - 1][i - 1];
//	A[i][i] = A[i][i] - (m * A[i - 1][i]);
//	B[i] = B[i] - (m * B[i - 1]);
//}
//
////backward substitution
//
//X[r - 1] = B[r - 1] / A[r - 1][r - 1];
//
//for (i = r - 2; i >= 0; i--)
//{
//	X[i] = (B[i] - A[i][i + 1] * X[i + 1]) / A[i][i];
//}
//printf("\nSolution of X is:\n");
//for (i = 0; i < r; i++)
//printf("%4.2f\n", X[i]);

void time_marching_Crank_Nicolson()
{
	a = 0.5 * sigma;
	b = -1.0 - sigma;
	c = a;

	double u0 = qField[0];
	double u1 = qField[1];
	double u2 = qField[2];

	double rhs1 = - u1 - 0.5 * sigma * (u2 - 2.0 * u1 + u0);
	double d1   = rhs1 - a * u0;

	double um   = qField[numberOfEquations];
	double umm1 = qField[numberOfEquations - 1];
	double ump1 = qField[numberOfEquations + 1];

	double rhsm = - um - 0.5 * sigma * (ump1 - 2.0 * um + umm1);
	double dm   = rhsm - c * ump1;

	vector<double> djv(numberOfEquations);
	djv[0] = d1;
	djv[numberOfEquations - 1] = dm;

	for (int iEquation = 1; iEquation < numberOfEquations-1; ++iEquation)
	{
		double u0 = qField[iEquation - 1];
		double u1 = qField[iEquation    ];
		double u2 = qField[iEquation + 1];
		djv[iEquation] = - u1 - 0.5 * sigma * ( u2 - 2.0 * u1 + u0 );
	}

	solve_chase_method(djv);
}

void boundary_condition()
{
	qField[0] = 100;
	qField[numberOfGridPoints - 1] = 0;

	qField_N1[0] = 100;
	qField_N1[numberOfGridPoints - 1] = 0;
}

void flow_initialization()
{
	//��ֵ��0
	qField.resize(numberOfGridPoints);
	qField_N1.resize(numberOfGridPoints);
	boundary_condition();
}

void initialize_parameter()
{		
	cout << "Enter number of grid points..." << endl;
	cin >> numberOfGridPoints;
	cout << "numberOfGridPoints = " << numberOfGridPoints << endl;

	cout << "Enter totalTime..." << endl;
	cin >> totalTime;
	cout << "totalTime = " << totalTime << endl;

	generate_grid_1D( numberOfGridPoints );

	int iter_min = int( 2.0 * totalTime * beta / ds / ds );

	cout << "Enter number of time steps..." << "iter > " << iter_min << endl;
	cin >> numberOfTimeSteps;
	cout << "numberOfTimeSteps = " << numberOfTimeSteps << endl;

	dt = totalTime / numberOfTimeSteps;
	sigma = beta * dt / ds / ds;
	cout << "sigma = " << sigma << endl;

	//m+2������㣬ֻ��m������
	numberOfEquations = numberOfGridPoints - 2;  //��ʽ��ʽ�������̸���
}

void generate_grid_1D( int numberOfGridPoints )
{
	double startCoord = 0.0; 
	double endCoord = 1.0;
	ds = ( endCoord - startCoord ) / ( numberOfGridPoints - 1 );

	xCoordinates.resize(numberOfGridPoints);
	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		xCoordinates[iNode] = startCoord + ds * iNode;
	}
	//return xCoordinates;
}