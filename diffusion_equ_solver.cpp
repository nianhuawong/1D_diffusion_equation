#include<iostream>
#include<vector>
#include<fstream>
#include "diffu_solver.h"
using namespace std;

int main()
{
	initialize_parameter();

	flow_initialization();

	for (int iter = 0; iter < numberOfTimeSteps; ++iter)
	{
		load_qField();

		time_marching();

		boundary_condition();

		compute_residual();

		output_residual();
	}

	output_results();

}

void output_results()
{
	fstream file;
	file.open("results.dat", ios_base::out | ios_base::app);

	file << "TITLE     = \"results\"" << endl;
	file << "VARIABLES = \"x\", \"qField\"" << endl;

	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		file << xCoordinates[iNode] << "\t" << qField_N1[iNode] << endl;
	}	

	file.close();
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

void time_marching()
{	
	//FTCS
	for (int iNode = 1; iNode < numberOfGridPoints - 1; ++iNode)
	{
		qField_N1[iNode] = qField[iNode] + sigma * (qField[iNode+1]- 2.0 * qField[iNode] + qField[iNode-1]);
	}
}

void boundary_condition()
{
	qField[0] = 100;
	qField[numberOfGridPoints - 1] = 0;
}

void flow_initialization()
{
	//³õÖµ¸³0
	qField.resize(numberOfGridPoints);
	boundary_condition();

	qField_N1 = qField;
}

void initialize_parameter()
{		
	cout << "Enter number of grid points..." << endl;
	cin >> numberOfGridPoints;
	cout << "numberOfGridPoints = " << numberOfGridPoints << endl;
	
	cout << "Enter number of time steps..." << endl;
	cin >> numberOfTimeSteps;
	cout << "numberOfTimeSteps = " << numberOfTimeSteps << endl;

	cout << "Enter totalTime..." << endl;
	cin >> totalTime;
	cout << "totalTime = " << totalTime << endl;

	xCoordinates = generate_grid_1D( numberOfGridPoints );
	double startCoord = 0.0;
	double endCoord = 1.0;
	ds = (endCoord - startCoord) / (numberOfGridPoints - 1);

	dt = totalTime / numberOfTimeSteps;
	sigma = beta * dt / ds / ds;
}

vector<double> generate_grid_1D( int numberOfGridPoints )
{
	double startCoord = 0.0; 
	double endCoord = 1.0;
	double ds = ( endCoord - startCoord ) / ( numberOfGridPoints - 1 );

	vector<double> xCoordinates(numberOfGridPoints);
	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		xCoordinates[iNode] = startCoord + ds * iNode;
	}
	return xCoordinates;
}