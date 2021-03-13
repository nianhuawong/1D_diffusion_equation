#include<iostream>
#include<vector>
#include<fstream>

using namespace std;
#include "diffu_solver.h"

int main()
{
	//double a = 1, b = 2, c = 3;
	//vector<double> djv = { 1, 2, 3 };
	//vector<double> solution(3);
	//solve_chase_method(a, b, c, djv,solution);
	
	initialize_parameter();

	flow_initialization();

	for (iter = 0; iter < numberOfTimeSteps; ++iter)
	{
		load_qField();

		//time_marching_FTCS();
		//time_marching_full_implicit();
		time_marching_Crank_Nicolson(); 

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

	vector<double> solution(numberOfEquations);
	solve_chase_method(a, b, c, djv, solution);

	for (int iEquation = 0; iEquation < numberOfEquations; ++iEquation)
	{
		qField_N1[iEquation + 1] = solution[iEquation];
	}
}

void solve_TDMA_method(double a, double b, double c, vector<double>& djv, vector<double>& solution)
{
	int  numberOfEquations = int( djv.size() );
	vector< double > AA(numberOfEquations);
	vector< double > BB(numberOfEquations);

	double d1 = djv[0];
	double dm = djv[numberOfEquations-1];

	AA[0] = -c / b;
	BB[0] = d1 / b;

	for (int iEquation = 1; iEquation < numberOfEquations; ++iEquation)
	{
		double dj = djv[iEquation];

		double tmp = b + a * AA[iEquation - 1];
		AA[iEquation] = -c / tmp;
		BB[iEquation] = (dj - a * BB[iEquation - 1]) / tmp;
	}

	solution[numberOfEquations-1] = (dm - a * BB[numberOfEquations - 2]) / (b + a * AA[numberOfEquations - 2]);
	for (int iEquation = numberOfEquations-2; iEquation >=0; --iEquation)
	{
		solution[iEquation] = AA[iEquation] * solution[iEquation+1] + BB[iEquation];
	}
}

//void solve_TDMA_method(double a, double b, double c, vector<double>& VD, vector<double>& solution)
//{
//	int  numberOfEquations = VD.size();
//
//	vector< double > VA(numberOfEquations);
//	vector< double > VB(numberOfEquations);
//	vector< double > VC(numberOfEquations);
//	for (int iEquation = 0; iEquation < numberOfEquations; ++iEquation)
//	{
//		VA[iEquation] = a;
//		VB[iEquation] = b;
//		VC[iEquation] = c;
//	}
//	for (int iEquation = 1; iEquation < numberOfEquations; iEquation++)
//	{
//		double tmp = VA[iEquation] / VB[iEquation - 1];
//
//		VB[iEquation] = VB[iEquation] - tmp * VC[iEquation - 1];
//		VD[iEquation] = VD[iEquation] - tmp * VD[iEquation - 1];
//	}
//
//	solution[numberOfEquations - 1] = VD[numberOfEquations - 1] / VB[numberOfEquations - 1];
//
//	for (int iEquation = numberOfEquations - 2; iEquation >= 0; iEquation--)
//	{
//		solution[iEquation] = ( VD[iEquation] - VC[iEquation]  * solution[iEquation+1] )/ VB[iEquation];
//	}
//
//}

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

	vector<double> solution(numberOfEquations);
	solve_chase_method(a, b, c, djv,solution);

	for (int iEquation = 0; iEquation < numberOfEquations; ++iEquation)
	{
		qField_N1[iEquation + 1] = solution[iEquation];
	}
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
	//初值赋0
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

	//m+2个网格点，只有m个方程
	numberOfEquations = numberOfGridPoints - 2;  //隐式格式代数方程个数
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