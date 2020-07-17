#include "twobartrussdesign.h"
#include <math.h>

void TWOBARTRUSS::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	tmp_num_obj = 2;
	tmp_num_cons = 3;
	tmp_dim = 3;
	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	dim = tmp_dim;

}

void TWOBARTRUSS::range_setting(int swit, vector<vector<double>>& range) {
	range.resize(dim);
	for (int i = 0; i < dim; ++i) {
		range[i].resize(2);
	}
	range[0][0] = 0.000825;
	range[0][1] = 0.023865;
	range[1][0] = 0.0011317;
	range[1][1] = 0.0683065;
	range[2][0] = 1.0;
	range[2][1] = 3.0;
}

individual TWOBARTRUSS::problems(int swit, individual ind) {
	ind = twobartrussproblem(ind);
	return ind;

}

//制約値は，負のときに制約違反しているとみなす．(>=0の形にせねば)
individual TWOBARTRUSS::twobartrussproblem(individual ind) {
	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;
	vector<long double> constraint = ind.constraint;

	f[0] = x[0] * pow(16 + pow(x[2], 2), 0.5) + x[1] * pow(1 + pow(x[2], 2), 0.5);
	f[1] = (20 * pow(16 + pow(x[2], 2), 0.5)) / (x[2] * x[0]);

	constraint[0] = 0.1 - f[0];
	constraint[1] = 100000 - f[1];
	constraint[2] = 100000 - (80 * pow(1 + pow(x[2], 2), 0.5)) / (x[2] * x[1]);

	ind.cost = f;
	ind.constraint = constraint;
	return ind;

}
