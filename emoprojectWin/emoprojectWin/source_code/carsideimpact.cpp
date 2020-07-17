#include "carsideimpact.h"
#include <math.h>

void CARSIDE::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	tmp_num_obj = 3;
	tmp_num_cons = 10;
	tmp_dim = 7;
	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	dim = tmp_dim;

}

void CARSIDE::range_setting(int swit, vector<vector<double>>& range) {
	range.resize(dim);
	for (int i = 0; i < dim; ++i) {
		range[i].resize(2);
	}
	range[0][0] = 0.5;
	range[0][1] = 1.5;
	range[1][0] = 0.45;
	range[1][1] = 1.35;
	range[2][0] = 0.5;
	range[2][1] = 1.5;
	range[3][0] = 0.5;
	range[3][1] = 1.5;
	range[4][0] = 0.875;
	range[4][1] = 2.625;
	range[5][0] = 0.4;
	range[5][1] = 1.2;
	range[6][0] = 0.4;
	range[6][1] = 1.2;
}

individual CARSIDE::problems(int swit, individual ind) {
	ind = carsideproblem(ind);
	return ind;

}

// object function ans constrains
individual CARSIDE::carsideproblem(individual ind) {
	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;
	vector<long double> constraint = ind.constraint;

	long double v_mbp = 10.58 - 0.674 * x[0] * x[1] - 0.67275 * x[1];
	long double v_fd = 16.45 - 0.489 * x[2] * x[6] - 0.843 * x[4] * x[5];

	f[0] = 1.98 + 4.9 * x[0] + 6.67 * x[1] + 6.98 * x[2] + 4.01 * x[3] + 1.78 * x[4] + 0.00001 * x[5] + 2.73 * x[6];
	f[1] = 4.72 - 0.5 * x[3] - 0.19 * x[1] * x[2];
	f[2] = 0.5 * (v_mbp + v_fd);

	constraint[0] = 1.0 - 1.16 + 0.3717 * x[1] * x[3] + 0.0092928 * x[2];
	constraint[1] = 0.32 - 0.261 + 0.0159 * x[0] * x[1] + 0.06486 * x[0] + 0.019 * x[1] * x[6] - 0.0144 * x[2] * x[4] - 0.0154464 * x[5];
	constraint[2] = 0.32 - 0.214 - 0.00817 * x[4] + 0.045195 * x[0] + 0.0135168 * x[0] - 0.03099 * x[1] * x[5] + 0.018 * x[1] * x[6] - 0.007176 * x[2] - 0.023232 * x[2] + 0.00364 * x[4] * x[5] + 0.018 * pow(x[1], 2);
	constraint[3] = 0.32 - 0.74 + 0.61 * x[1] + 0.031296 * x[2] + 0.031872 * x[6] - 0.277 * pow(x[1], 2);
	constraint[4] = 32 - 28.98 - 3.818 * x[2] + 4.2 * x[0] * x[1] - 1.27296 * x[5] + 2.68065 * x[6];
	constraint[5] = 32 - 33.86 - 2.95 * x[2] + 5.057 * x[0] * x[1] + 3.795 * x[1] + 3.4431 * x[6] - 1.45728;
	constraint[6] = 32 - 46.36 + 9.9 * x[1] + 4.4505 * x[0];
	constraint[7] = 4.0 - f[1];
	constraint[8] = 9.9 - v_mbp;
	constraint[9] = 15.7 - v_fd;

	ind.cost = f;
	ind.constraint = constraint;
	return ind;

}
