#pragma once
const double beta = 1.43e-3;

int numberOfGridPoints;
int numberOfTimeSteps;
double totalTime;
double dt;
double ds;
double sigma;
vector<double> qField;
vector<double> qField_N1;
double residual;
int iter;
vector<double> xCoordinates;

void initialize_parameter();
void compute_residual();
vector<double> generate_grid_1D( int numberOfGridPoints);

