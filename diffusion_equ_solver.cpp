#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>

using namespace std;
#include "diffu_solver.h"

int main()
{	
	//验证三对角矩阵追赶法求解是否正确
	//double a = 1, b = 2, c = 3;			//A={1,2,3}, {3,2,5}//b={1,2,3}
	//vector<double> d = { 1, 2, 3 };		//x={-4,3,0}, {-0.8462,0.5385,0.6923}
	//vector<double> solution(3);
	//solve_TDMA_method(a, b, c, d, solution);

	initialize_parameter();

	flow_initialization();

	for (iter = 0; iter < numberOfTimeSteps; ++iter)
	{
		load_qField();

		time_marching();
		//time_marching_FTCS();
		//time_marching_full_implicit();
		//time_marching_Crank_Nicolson(); 

		boundary_condition();

		compute_residual();

		output_residual();

		physicalTime += dt;
	}

	output_results(outFile);

	return 0;
}

void output_results(string outFile)
{
	cout << "dumping results..." << endl;
	fstream file;
	file.open(outFile, ios_base::out );

	file << "TITLE     = \"results\"" << endl;
	file << "VARIABLES = \"x\", \"qField\"" << endl;

	file << setiosflags(ios::right);
//	file << setiosflags(ios::scientific);
	file << setprecision(8);
	for (int iNode = 0; iNode < numberOfGridPoints; ++iNode)
	{
		file << xCoordinates[iNode] << "\t" << qField_N1[iNode] << endl;
	}	

	file.close();

	cout << "done!" << endl;
}

void output_residual()
{	
	if (iter % 1000 == 0)
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
		djv[iEquation] = qField[iEquation+1];//这个地方是iEquation+1，找了好久才找到这个Bug:(，仔细看公式下标
	}

	vector<double> solution(numberOfEquations);
	solve_TDMA_method(a, b, c, djv, solution);

	for (int iEquation = 0; iEquation < numberOfEquations; ++iEquation)
	{
		qField_N1[iEquation + 1] = solution[iEquation];
	}
}

//求解三对角矩阵的追赶法，Tri-Diagonal Matrix Algorithm, 参考中科院教材做法
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

//求解三对角矩阵的追赶法，Tri-Diagonal Matrix Algorithm, 参考网络开源代码
void solve_TDMA_method_version2(double a, double b, double c, vector<double>& VD, vector<double>& solution)
{
	int  numberOfEquations = int( VD.size() );

	vector< double > VA(numberOfEquations);
	vector< double > VB(numberOfEquations);
	vector< double > VC(numberOfEquations);
	for (int iEquation = 0; iEquation < numberOfEquations; ++iEquation)
	{
		VA[iEquation] = a;
		VB[iEquation] = b;
		VC[iEquation] = c;
	}
	for (int iEquation = 1; iEquation < numberOfEquations; iEquation++)
	{
		double tmp = VA[iEquation] / VB[iEquation - 1];

		VB[iEquation] = VB[iEquation] - tmp * VC[iEquation - 1];
		VD[iEquation] = VD[iEquation] - tmp * VD[iEquation - 1];
	}

	solution[numberOfEquations - 1] = VD[numberOfEquations - 1] / VB[numberOfEquations - 1];

	for (int iEquation = numberOfEquations - 2; iEquation >= 0; iEquation--)
	{
		solution[iEquation] = ( VD[iEquation] - VC[iEquation]  * solution[iEquation+1] )/ VB[iEquation];
	}
}

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
		double u0 = qField[iEquation ];		
		double u1 = qField[iEquation + 1];		//注意此处下标，仔细对照公式
		double u2 = qField[iEquation + 2];
		djv[iEquation] = - u1 - 0.5 * sigma * ( u2 - 2.0 * u1 + u0 );
	}

	vector<double> solution(numberOfEquations);
	solve_TDMA_method(a, b, c, djv,solution);

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

	set_time_march_method();

	generate_grid_1D( numberOfGridPoints );

	int iter_min = int( 2.0 * totalTime * beta / ds / ds );

	cout << "Enter number of time steps..." << "for FTCS, iter > " << iter_min << endl;
	cin >> numberOfTimeSteps;
	cout << "numberOfTimeSteps = " << numberOfTimeSteps << endl;

	dt = totalTime / numberOfTimeSteps;
	sigma = beta * dt / ds / ds;
	cout << "sigma = " << sigma << endl;

	//m+2个网格点，只有m个方程
	numberOfEquations = numberOfGridPoints - 2;  //隐式格式代数方程个数
}

void set_time_march_method()
{
	cout << "1--FTCS;\t2--fully implict;\t3--Crank_Nicolson, please choose!" << endl;
	int time_march_method;
	cin >> time_march_method;

	if (time_march_method == 1)
	{
		time_marching = &time_marching_FTCS;
		cout << "time marching method is FTCS!" << endl;
		outFile = "results-explicit.dat";
	}
	else if (time_march_method == 2)
	{
		time_marching = &time_marching_full_implicit;
		cout << "time marching method is fully implict!" << endl;
		outFile = "results-full.dat";
	}
	else if (time_march_method == 3)
	{
		time_marching = &time_marching_Crank_Nicolson;
		cout << "time marching method is Crank_Nicolson!" << endl;
		outFile = "results-CN.dat";
	}
	else
	{
		cout << "invalid time marching method, program ends!" << endl;
		exit(1);
	}
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