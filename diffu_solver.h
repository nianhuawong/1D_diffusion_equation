#pragma once

const double beta = 1.43e-3;

int iter;
int numberOfGridPoints;
int numberOfTimeSteps;
int numberOfEquations;

double a, b, c;
double totalTime;
double dt;
double ds;
double sigma;

double residual;

vector< double > qField;
vector< double > qField_N1;
vector< double > xCoordinates;

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
void solve_TDMA_method(double a, double b, double c, vector<double>& VD, vector<double>& solution);
