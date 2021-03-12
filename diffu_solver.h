#pragma once
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
void time_marching();
void boundary_condition();
void compute_residual();
void output_residual();
void output_results();

void generate_grid_1D( int numberOfGridPoints);

