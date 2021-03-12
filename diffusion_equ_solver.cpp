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

	double u1 = qField[1];
	double d1 = u1 - a * qField[0];

	double um = qField[numberOfGridPoints - 2];
	double dm = um - c * qField[numberOfGridPoints-1];
	
	vector<double> uj(numberOfEquations);
	for (int iEquation = 1; iEquation < numberOfEquations; ++iEquation)
	{
		uj[iEquation] = qField[iEquation];
	}
	
	solve_chase_method(d1,dm,uj);
}

void solve_chase_method(double d1, double dm, vector<double>& uj)
{
	vector< double > AA(numberOfEquations);
	vector< double > BB(numberOfEquations);

	AA[0] = -c / b;
	BB[0] = d1 / b;

	for (int iEquation = 1; iEquation < numberOfEquations; ++iEquation)
	{
		double dj = uj[iEquation];

		double tmp = b + a * AA[iEquation - 1];
		AA[iEquation] = -c / tmp;
		BB[iEquation] = (dj - a * BB[iEquation - 1]) / tmp;
	}

	qField_N1[numberOfEquations] = (dm - a * BB[numberOfEquations - 1]) / (b + a * AA[numberOfEquations - 1]);

	for (int iEquation = numberOfEquations - 1; iEquation >= 1; --iEquation)
	{
		qField_N1[iEquation] = AA[iEquation] * qField_N1[iEquation + 1] + BB[iEquation];
	}
}

void time_marching_Crank_Nicolson()
{
	a = 0.5 * sigma;
	b = -1.0 - sigma;
	c = a;

	double u1 = -qField[1] - 0.5 * sigma * (qField[2] - 2.0 * qField[1] + qField[0]);
	double d1 = u1 - a * qField[0];

	double um = -qField[numberOfGridPoints - 2] - 0.5 * sigma * (qField[numberOfGridPoints - 1] - 2.0 * qField[numberOfGridPoints - 2] + qField[numberOfGridPoints - 3]);
	double dm = um - c * qField[numberOfGridPoints - 1];

	vector<double> uj(numberOfEquations);
	for (int iEquation = 1; iEquation < numberOfEquations; ++iEquation)
	{
		uj[iEquation] = -qField[iEquation] - 0.5 * sigma * (qField[iEquation + 1] - 2.0 * qField[iEquation] + qField[iEquation - 1]);
	}

	solve_chase_method(d1, dm, uj);
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
	double startCoord = 0.0;
	double endCoord = 1.0;
	ds = (endCoord - startCoord) / (numberOfGridPoints - 1);

	int iter_min = int( 200 * beta / ds / ds );

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
	double ds = ( endCoord - startCoord ) / ( numberOfGridPoints - 1 );

	xCoordinates.resize(numberOfGridPoints);
	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		xCoordinates[iNode] = startCoord + ds * iNode;
	}
	//return xCoordinates;
}