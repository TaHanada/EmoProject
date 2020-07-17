#include "discbrakedesign.h"
#include <math.h>

void DISCBRAKE::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	tmp_num_obj = 2;
	tmp_num_cons = 5;
	tmp_dim = 4;
	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	dim = tmp_dim;

}

void DISCBRAKE::range_setting(int swit, vector<vector<double>>& range) {
	range.resize(dim);
	for (int i = 0; i < dim; ++i) {
		range[i].resize(2);
	}
	switch (swit) {
	case 1:
		range[0][0] = 55;
		range[0][1] = 80;
		range[1][0] = 75;
		range[1][1] = 110;
		range[2][0] = 1000;
		range[2][1] = 3000;
		range[3][0] = 2;
		range[3][1] = 20.999999;
		break;

	case 2: //x_3 = 3000 ‚ÉŒÅ’è
		range[0][0] = 55;
		range[0][1] = 80;
		range[1][0] = 75;
		range[1][1] = 110;
		range[2][0] = 3000;
		range[2][1] = 3000;
		range[3][0] = 2;
		range[3][1] = 20.999999;
		break;
	}

}

individual DISCBRAKE::problems(int swit, individual ind) {
	ind = discproblem(ind);
	return ind;

}

//§–ñ’l‚ÍC•‰‚Ì‚Æ‚«‚É§–ñˆá”½‚µ‚Ä‚¢‚é‚Æ‚Ý‚È‚·D
individual DISCBRAKE::discproblem(individual ind) {
	vector<long double> x = ind.var;
	x[3] = floor(x[3]);
	vector<long double> f = ind.cost;
	vector<long double> constraint = ind.constraint;

	f[0] = 4.9 * pow(10, -5) * (pow(x[1], 2) - pow(x[0], 2)) * (x[3] - 1);
	f[1] = (9.82 * pow(10, 6) * (pow(x[1], 2) - pow(x[0], 2))) / (x[2] * x[3] * (pow(x[1], 3) - pow(x[0], 3)));

	constraint[0] = (x[1] - x[0]) - 20;
	constraint[1] = 30 - 2.5 * (x[3] + 1);
	constraint[2] = 0.4 - x[2] / (3.14 * (pow(x[1], 2) - pow(x[0], 2)));
	constraint[3] = 1 - (2.22 * pow(10, -3) * x[2] * (pow(x[1], 3) - pow(x[0], 3))) / pow((pow(x[1], 2) - pow(x[0], 2)), 2);
	constraint[4] = (2.66 * pow(10, -2) * x[2] * x[3] * (pow(x[1], 3) - pow(x[0], 3))) / (pow(x[1], 2) - pow(x[0], 2)) - 900;

	ind.cost = f;
	ind.constraint = constraint;
	return ind;

}