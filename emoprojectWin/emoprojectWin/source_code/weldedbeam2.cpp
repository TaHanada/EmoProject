#include "weldedbeam2.h"
#include <math.h>

void WELDEDBEAM::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	tmp_num_obj = 2;
	tmp_num_cons = 4;
	tmp_dim = 4;
	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	dim = tmp_dim;

}

void WELDEDBEAM::range_setting(int swit, vector<vector<double>>& range) {
	range.resize(dim);
	for (int i = 0; i < dim; ++i) {
		range[i].resize(2);
	}
	range[0][0] = 0.123;
	range[0][1] = 5.0;
	range[1][0] = 0.1;
	range[1][1] = 10.0;
	range[2][0] = 0.1;
	range[2][1] = 10.0;
	range[3][0] = 0.123;
	range[3][1] = 5.0;
}

individual WELDEDBEAM::problems(int swit, individual ind) {
	ind = weldedbeamproblem(ind);
	return ind;

}

//制約値は，負のときに制約違反しているとみなす．(>=0の形にせねば)
individual WELDEDBEAM::weldedbeamproblem(individual ind) {
	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;
	vector<long double> constraint = ind.constraint;
	
	long double r2 = 6000 * (14.0 + 0.5 * x[1]) * sqrt(0.25 * (pow(x[1], 2) + pow(x[0] + x[2], 2))) / (2 * (sqrt(2) * x[0] * x[2] * (pow(x[1], 2) / 12 + 0.25 * pow(x[0] + x[2], 2))));
	long double r1 = 6000 / (sqrt(2) * x[0] * x[1]);
	long double r = sqrt(pow(r1, 2) + pow(r2, 2) + x[1] * r1 * r2 / sqrt(0.25 * (pow(x[1], 2) + pow(x[0] + x[2], 2))));

	f[0] = 1.10471 * pow(x[0], 2) * x[1] + 0.04811 * x[2] * x[3] * (14.0 + x[1]);
	f[1] = 2.1952 / (x[3] * pow(x[2], 3));
	
	/*
	long double r2 = 6000 * (14.0 + 0.5 * x[1]) * sqrt(0.25 * (pow(x[1], 2) + pow(x[0] + x[2], 2))) / (2 * (0.707 * x[0] * x[1] * (pow(x[1], 2) / 12 + 0.25 * (pow(x[1], 2) + pow(x[0] + x[2], 2)))));
	long double r1 = 6000 / (sqrt(2) * x[0] * x[1]);
	long double r = sqrt(pow(r1, 2) + pow(r2, 2) + x[1] * r1 * r2 / sqrt(0.25 * (pow(x[1], 2) + pow(x[0] + x[2], 2))));

	f[0] = 1.10471 * pow(x[0], 2) * x[1] + 0.04811 * x[2] * x[3] * (14.0 + x[2]);
	f[1] = 2.1952 / (x[3] * pow(x[2], 3));
	*/

	constraint[0] = 13600 - r;
	constraint[1] = 30000 - 504000 / (pow(x[2], 2) * x[3]);
	constraint[2] = x[3] - x[0];
	// constraint[3] = 64764.022 * (1 - 0.0282346 * x[2]) * x[2] * pow(x[3], 3) - 60000;
	constraint[3] = 64764.022 * (1 - 0.0282346 * x[2]) * x[2] * pow(x[3], 3) - 6000;

	ind.cost = f;
	ind.constraint = constraint;
	return ind;

}
