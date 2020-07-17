#include "speedreducerdesign.h"
#include <math.h>

void SPEEDREDUCER::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	tmp_num_obj = 2;
	tmp_num_cons = 11;
	tmp_dim = 7;
	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	dim = tmp_dim;

}

void SPEEDREDUCER::range_setting(int swit, vector<vector<double>>& range) {
	range.resize(dim);
	for (int i = 0; i < dim; ++i) {
		range[i].resize(2);
	}
	range[0][0] = 2.6;
	range[0][1] = 3.6;
	range[1][0] = 0.7;
	range[1][1] = 0.8;
	range[2][0] = 17;
	range[2][1] = 28;
	range[3][0] = 7.3;
	range[3][1] = 8.3;
	range[4][0] = 7.3;
	range[4][1] = 8.3;
	range[5][0] = 2.9;
	range[5][1] = 3.9;
	range[6][0] = 5.0;
	range[6][1] = 5.5;
}

individual SPEEDREDUCER::problems(int swit, individual ind) {
	ind = speedreducerproblem(ind);
	return ind;

}

//制約値は，負のときに制約違反しているとみなす．(>=0の形にせねば)
individual SPEEDREDUCER::speedreducerproblem(individual ind) {
	vector<long double> x = ind.var;
	x[2] = floor(x[2]);
	vector<long double> f = ind.cost;
	vector<long double> constraint = ind.constraint;

	f[0] = 0.7854 * x[0] * pow(x[1], 2) * ((10 * pow(x[2], 2)) / 3 + 14.933 * x[2] - 43.0934) - 1.508 * x[0] * (pow(x[5], 2) + pow(x[6], 2)) + 7.477 * (pow(x[5], 3) + pow(x[6], 3)) + 0.7854 * (x[3] * pow(x[5], 2) + x[4] * pow(x[6], 2));
	f[1] = pow(pow((745 * x[3]) / (x[1] * x[2]), 2) + 16900000, 0.5) / (0.1 * pow(x[5], 3));

	constraint[0] = 1.0 / 27 - 1.0 / (x[0] * pow(x[1], 2) * x[2]);
	constraint[1] = 1.0 / 397.5 - 1.0 / (x[0] * pow(x[1], 2) * pow(x[2], 2));
	constraint[2] = 1.0 / 1.93 - pow(x[3], 3) / (x[1] * x[2] * pow(x[5], 4));
	constraint[3] = 1.0 / 1.93 - pow(x[4], 3) / (x[1] * x[2] * pow(x[6], 4));
	constraint[4] = 40 - x[1] * x[2];
	constraint[5] = 12 - x[0] / x[1];
	constraint[6] = x[0] / x[1] - 5.0;
	constraint[7] = x[3] - 1.5 * x[5] - 1.9;
	constraint[8] = x[4] - 1.1 * x[6] - 1.9;
	constraint[9] = 1300 - f[1];
	constraint[10] = 1100 - pow(pow((745 * x[4]) / (x[1] * x[2]), 2) + 127500000, 0.5) / (0.1 * pow(x[6], 3));

	ind.cost = f;
	ind.constraint = constraint;
	return ind;

}
