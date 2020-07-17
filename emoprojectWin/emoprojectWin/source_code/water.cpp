#include "water.h"
#include <math.h>

void WATER::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	tmp_num_obj = 5;
	tmp_num_cons = 7;
	tmp_dim = 3;
	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	dim = tmp_dim;

}

void WATER::range_setting(int swit, vector<vector<double>>& range) {
	range.resize(dim);
	for (int i = 0; i < dim; ++i) {
		range[i].resize(2);
		if (i < 1) {
			range[i][0] = 0.01;
			range[i][1] = 1.45;
		}
		else {
			range[i][0] = 0.01;
			range[i][1] = 0.10;
		}

	}

}

individual WATER::problems(int swit, individual ind) {
	ind = waterproblem(ind);
	return ind;

}

//§–ñ’l‚ÍC•‰‚Ì‚Æ‚«‚É§–ñˆá”½‚µ‚Ä‚¢‚é‚Æ‚Ý‚È‚·D
individual WATER::waterproblem(individual ind) {
	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;
	vector<long double> constraint = ind.constraint;

	f[0] = 106780.37 * (x[1] + x[2]) + 61704.67;
	f[1] = 3000.0 * x[0];
	f[2] = 305700.0 * 2289.0 * x[1] / pow(0.06 * 2289.0, 0.65);
	f[3] = 250.0 * 2289.0 * exp(-39.75 * x[1] + 9.9 * x[2] + 2.74);
	f[4] = 25.0 * (1.39 / (x[0] * x[1]) + 4940.0 * x[2] - 80);

	constraint[0] = 1 - (0.00139 / (x[0] * x[1]) + 4.94 * x[2] - 0.08);
	constraint[1] = 1 - (0.000306 / (x[0] * x[1]) + 1.082 * x[2] - 0.0986);
	constraint[2] = 50000.0 - (12.307 / (x[0] * x[1]) + 49408.24 * x[2] + 4051.02);
	constraint[3] = 16000.0 - (2.098 / (x[0] * x[1]) + 8046.33 * x[2] - 696.71);
	constraint[4] = 10000.0 - (2.138 / (x[0] * x[1]) + 7883.39 * x[2] - 705.04);
	constraint[5] = 2000.0 - (0.417 * x[0] * x[1] + 1721.26 * x[2] - 136.54);
	constraint[6] = 550.0 - (0.164 / (x[0] * x[1]) + 631.13 * x[2] - 54.48);

	ind.cost = f;
	ind.constraint = constraint;
	return ind;
}