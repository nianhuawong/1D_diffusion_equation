//#include "./diffu_solver.h"
#include<iostream>
#include<vector>
#include<fstream>

using namespace std;

const double beta = 1.43e-3;

int iter;
int numberOfGridPoints;
int numberOfTimeSteps;

double totalTime;
double dt;
double ds;
double sigma;

vector<double> qField;
vector<double> qField_N1;
double residual;

vector<double> xCoordinates;

void initialize_parameter();
void flow_initialization();
void load_qField();
void time_marching_FTCS();
void time_marching_full_implicit();
void time_marching_Crank_Nicolson();
void boundary_condition();
void compute_residual();
void output_residual();
void output_results();

void generate_grid_1D(int numberOfGridPoints);

int main()
{
	initialize_parameter();

	flow_initialization();

	for (iter = 0; iter < numberOfTimeSteps; ++iter)
	{
		load_qField();

		//time_marching_FTCS();
		time_marching_full_implicit();
		//time_marching_Crank_Nicolson

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
	double a = - sigma; 
	double b = 1.0 + 2.0 * sigma;
	double c = - sigma;
	//m+2个网格点，只有m个方程
	int numberOfEquations = numberOfGridPoints - 2;  //代数方程个数

	double d1 = qField[1] - a * qField[0];
	double dm = qField[numberOfGridPoints - 2] - c * qField[numberOfGridPoints-1];

	vector< double > AA(numberOfEquations);
	vector< double > BB(numberOfEquations);

	AA[0]= -c / b;
	BB[0] = d1 / b;
	
	for (int iEquation = 1; iEquation < numberOfEquations; ++iEquation)
	{
		double dj = qField[iEquation];

		double tmp = b + a * AA[iEquation - 1];
		AA[iEquation] = - c / tmp;
		BB[iEquation] = (dj - a * BB[iEquation - 1]) / tmp;
	}
	
	qField_N1[numberOfEquations] = (dm - a * BB[numberOfEquations - 1]) / (b + a * AA[numberOfEquations - 1]);

	for (int iEquation = numberOfEquations-1; iEquation >= 1; --iEquation)
	{
		qField_N1[iEquation] = AA[iEquation] * qField_N1[iEquation+1] + BB[iEquation];
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