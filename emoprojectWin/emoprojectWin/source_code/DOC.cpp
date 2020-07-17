#include"DOC.h"
#define _USE_MATH_DEFINES
#include<math.h>


void DOC::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	switch (swit) {
	case 1:
		tmp_num_obj = 2;
		tmp_num_cons = 7;
		tmp_dim = 6;
		break;

	case 2:
		tmp_num_obj = 2;
		tmp_num_cons = 7;
		tmp_dim = 16;
		break;

	case 3:
		tmp_num_obj = 2;
		tmp_num_cons = 10;
		tmp_dim = 10;
		break;

	case 4:
		tmp_num_obj = 2;
		tmp_num_cons = 6;
		tmp_dim = 8;
		break;

	case 5:
		tmp_num_obj = 2;
		tmp_num_cons = 9;
		tmp_dim = 8;
		break;

	case 6:
		tmp_num_obj = 2;
		tmp_num_cons = 10;
		tmp_dim = 11;
		break;

	case 7:
		tmp_num_obj = 2;
		tmp_num_cons = 6;
		tmp_dim = 11;
		break;

	case 8:
		tmp_num_obj = 3;
		tmp_num_cons = 7;
		tmp_dim = 10;
		break;

	case 9:
		tmp_num_obj = 3;
		tmp_num_cons = 14;
		tmp_dim = 11;
		break;

	default:
		break;
	}


	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	dim = tmp_dim;
}


void DOC::range_setting(int swit, vector<vector<double>>& range) {
	range.resize(dim);

	switch (swit) {
	case 1: {
		double lower[6] = { 0,78,33,27,27,27 };
		double upper[6] = { 1,102,45,45,45,45 };
		for (int i = 0; i < dim; ++i) {
			range[i].push_back(lower[i]), range[i].push_back(upper[i]);
		}
		break;
	}
	case 2: {
		double lower[16] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
		double upper[16] = { 1,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10 };
		for (int i = 0; i < dim; ++i) {
			range[i].push_back(lower[i]), range[i].push_back(upper[i]);
		}
		break;
	}

	case 3: {
		double lower[10] = { 0,0,0,0,0,0,0,0,0,0.01 };
		double upper[10] = { 1,1,300,100,200,100,1,100,200,0.03 };
		for (int i = 0; i < dim; ++i) {
			range[i].push_back(lower[i]), range[i].push_back(upper[i]);
		}
		break;
	}

	case 4: {
		double lower[8] = { 0,-10,-10,-10,-10,-10,-10,-10 };
		double upper[8] = { 1,10,10,10,10,10,10,10 };
		for (int i = 0; i < dim; ++i) {
			range[i].push_back(lower[i]), range[i].push_back(upper[i]);
		}
		break;
	}

	case 5: {
		double lower[8] = { 0,0,0,0,100,6.3,5.9,4.5 };
		double upper[8] = { 1,1000,40,40,300,6.7,6.4,6.25 };
		for (int i = 0; i < dim; ++i) {
			range[i].push_back(lower[i]), range[i].push_back(upper[i]);
		}
		break;
	}

	case 6: {
		double lower[11] = { 0,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10 };
		double upper[11] = { 1,10,10,10,10,10,10,10,10,10,10 };
		for (int i = 0; i < dim; ++i) {
			range[i].push_back(lower[i]), range[i].push_back(upper[i]);
		}
		break;
	}

	case 7: {
		double lower[11] = { 0,0,0,0,0,0,0,0,0,0,0 };
		double upper[11] = { 1,10,10,10,10,10,10,10,10,10,10 };
		for (int i = 0; i < dim; ++i) {
			range[i].push_back(lower[i]), range[i].push_back(upper[i]);
		}
		break;
	}

	case 8: {
		double lower[10] = { 0,0,500,1000,5000,100,100,100,100,100 };
		double upper[10] = { 1,1,1000,2000,6000,500,500,500,500,500 };
		for (int i = 0; i < dim; ++i) {
			range[i].push_back(lower[i]), range[i].push_back(upper[i]);
		}
		break;
	}

	case 9: {
		double lower[11] = { 0,-1,-1,-1,-1,-1,-1,-1 - 1,-1,-1 };
		double upper[11] = { 1,1,10,10,10,10,10,10,10,10,10 };
		for (int i = 0; i < dim; ++i) {
			range[i].push_back(lower[i]), range[i].push_back(upper[i]);
		}
		break;
	}

	default:
		break;
	}


}


individual DOC::problems(int swit, individual ind) {
	switch (swit) {
	case 1:
		return DOC1(ind);
		break;

	case 2:
		return DOC2(ind);
		break;

	case 3:
		return DOC3(ind);
		break;

	case 4:
		return DOC4(ind);
		break;

	case 5:
		return DOC5(ind);
		break;

	case 6:
		return DOC6(ind);
		break;

	case 7:
		return DOC7(ind);
		break;

	case 8:
		return DOC8(ind);
		break;

	case 9:
		return DOC9(ind);
		break;

	default:
		break;
	}
}


individual DOC::DOC1(individual ind) {
	double g, g_temp;
	vector<long double> x = ind.var, y = ind.cost, c = ind.constraint;
	// basic multi - objective problem
	g_temp = 5.3578547 * pow(x[3],2) + 0.8356891 * x[1]*x[5] + 37.293239 * x[1] - 40792.141;
	g = g_temp + 30665.5386717834 + 1;
	y[0] = x[0];
	y[1] = g*(1 - sqrt(y[0])/ g);

	//constraints in objective space
	c[0] = max(-1*(pow(y[0],2) + pow(y[1],2) - 1), (long double)0.0);

	//constraints in decision space
	c[1] = 85.334407 + 0.0056858 * x[2]*x[5] + 0.0006262 * x[1]*x[4] - 0.0022053 * x[3]*x[5] - 92;
	c[2] = -85.334407 - 0.0056858 * x[2]*x[5] - 0.0006262 * x[1]*x[4] + 0.0022053 * x[3]*x[5];
	c[3] = 80.51249 + 0.0071317 * x[2]*x[5] + 0.0029955 * x[1]*x[2] + 0.0021813 *pow( x[3], 2) - 110;
	c[4] = -80.51249 - 0.0071317 * x[2]*x[5] - 0.0029955 * x[1]*x[2] - 0.0021813 * pow(x[3],2) + 90;
	c[5] = 9.300961 + 0.0047026 * x[3]*x[5] + 0.0012547 * x[1]*x[3] + 0.0019085 * x[3]*x[4] - 25;
	c[6] = -9.300961 - 0.0047026 * x[3]*x[5] - 0.0012547 * x[1]*x[3] - 0.0019085 * x[3]*x[4] + 20;

	ind.cost = y;
	for (int i = 0; i < num_cons; ++i) {
		ind.constraint[i] = -1*c[i];
	}

	return ind;
}


individual DOC::DOC2(individual ind) {
	double g, g_temp;
	vector<long double> x = ind.var, y = ind.cost, c = ind.constraint;

	double a[10][5] = { { -16,2,0,1,0 },
	{ 0,-2,0,0.4,2 },
	{ -3.5,0,2,0,0 },
	{ 0,-2,0,-4,-1 },
	{ 0,-9,-2,1,-2.8 },
	{ 2,0,-4,0,0 },
	{ -1,-1,-1,-1,-1 },
	{ -1,-2,-3,-2,-1 },
	{ 1,2,3,4,5 },
	{ 1,1,1,1,1 } };
	double b[10] = { -40,-2,-0.25,-4,-4,-1,-40,-60,5,1 };
	double c1[5][5] = { { 30,-20,-10,32,-10 },
	{ -20,39,-6,-31,32 },
	{ -10,-6,10,-6,-10 },
	{ 32,-31,-6,39,-20 },
	{ -10,32,-10,-20,30 } };
	double d[5] = { 4,8,10,6,2 };
	double e[5] = { -15,-27,-36,-18,-12 };

	//basic multi - objective problem
	g_temp = 0;
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			g_temp += c1[i][j] * x[10 + i] * x[10 + j];
		}
	}
	for (int i = 0; i < 5; ++i) {
		g_temp += 2 * d[i] * pow(x[10 + i], 3);
	}
	for (int i = 0; i < 10; ++i) {
		g_temp -= b[i] * x[1+i];
	}

	g = (g_temp - 32.6555929502) + 1;
	y[0] = x[0];
	y[1] = g*(1 - pow(y[0], 1.0 / 3.0) / g);

	//constraints in objective space
	double d1[3];
	c[0] = max(-1 * (sqrt(y[0]) + y[1] - 1), (long double)0.0);
	d1[0] = max((pow(y[0] - 1.0 / 8.0, 2) + pow(y[1] - 1 + sqrt(1.0 / 8.0), 2) - 0.15*0.15), (long double)0.0);
	d1[1] = max((pow(y[0] - 1.0 / 2.0, 2) + pow(y[1] - 1 + sqrt(1.0 / 2.0), 2) - 0.15*0.15), (long double)0.0);
	d1[2] = max((pow(y[0] - 7.0 / 8.0, 2) + pow(y[1] - 1 + sqrt(7.0 / 8.0), 2) - 0.15*0.15), (long double)0.0);
	c[1] = min(min(d1[0], d1[1]), d1[2]);


	//constraints in decision space
	for (int i = 0; i < 5; ++i) {
		c[2+i] = 0;
		for (int j = 0; j < 5; ++j) {
			c[2+i] += -2 * c1[j][i] * x[10 + j];
		}
		c[2 + i] += -3 * d[i] * pow(x[10 + i], 2) - e[i];
		for (int j = 0; j < 10; ++j) {
			c[2 + i] += a[j][i] * x[1+j];
		}
	}


	ind.cost = y;
	for (int i = 0; i < num_cons; ++i) {
		ind.constraint[i] = -1 * c[i];
	}

	return ind;
}


individual DOC::DOC3(individual ind) {
	double g, g_temp;
	vector<long double> x = ind.var, y = ind.cost, c = ind.constraint;
	// basic multi - objective problem
	g_temp = -9 * x[5] - 15 * x[8] + 6 * x[1] + 16 * x[2] + 10.* (x[6] + x[7]);
	g = (g_temp + 400.0551) + 1;
	y[0] = x[0];
	y[1] = g*(1 - (y[0]) / g);

	//constraints in objective space
	c[0] = max(-1 * (pow(y[0], 2) + pow(y[1], 2) - 1), (long double)0.0);
	c[1] = max(-1 * (abs((-1 * y[0] + y[1] - 0.5) / sqrt(2)) - 0.1 / sqrt(2)), (long double)0.0);
	c[2] = max(-1 * (abs((-1 * y[0] + y[1] - 0) / sqrt(2)) - 0.1 / sqrt(2)), (long double)0.0);
	c[3] = max(-1 * (abs((-1 * y[0] + y[1] + 0.5) / sqrt(2)) - 0.1 / sqrt(2)), (long double)0.0);

	//constraints in decision space
	c[4] = x[9] * x[3] + 0.02*x[6] - 0.025*x[5];
	c[5] = x[9] * x[4] + 0.02*x[7] - 0.015*x[8];
	c[6] = abs(x[1] + x[2] - x[3] - x[4]) - 0.0001;
	c[7] = abs(0.03*x[1] + 0.01* x[2] - x[9] * (x[3] + x[4])) - 0.0001;
	c[8] = abs(x[3] + x[6] - x[5]) - 0.0001;
	c[9] = abs(x[4] + x[7] - x[8]) - 0.0001;

	ind.cost = y;
	for (int i = 0; i < num_cons; ++i) {
		ind.constraint[i] = -1 * c[i];
	}

	return ind;
}


individual DOC::DOC4(individual ind) {
	double g, g_temp;
	vector<long double> x = ind.var, y = ind.cost, c = ind.constraint;
	// basic multi - objective problem
	g_temp = pow(x[1] - 10, 2) + 5 * pow(x[2] - 12, 2) + pow(x[3], 4) + 3 * pow(x[4] - 11, 2) + 10 * pow(x[5], 6) +
		7 * pow(x[6], 2) + pow(x[7], 4) - 4 * x[6] * x[7] - 10 * x[6] - 8 * x[7];

	g = g_temp - 680.6300573745 + 1;
	y[0] = x[0];
	y[1] = g*(1 - sqrt(y[0]) / g);

	//constraints in objective space
	c[0] = max(-1 * (y[0] + y[1] - 1), (long double)0.0);
	c[1] = max(-1 * (y[0] + y[1] - 1 - abs(sin(10 * M_PI*(y[0] - y[1] + 1)))), (long double)0.0);

	//constraints in decision space
	c[2] = -127 + 2 * pow(x[1], 2) + 3 * pow(x[2], 4) + x[3] + 4 * pow(x[4], 2) + 5 * x[5];
	c[3] = -282 + 7 * x[1] + 3 * x[2] + 10 * pow(x[3], 2) + x[4] - x[5];
	c[4] = -196 + 23 * x[1] + pow(x[2], 2) + 6 * pow(x[6], 2) - 8 * x[7];
	c[5] = 4 * pow(x[1], 2) + pow(x[2], 2) - 3 * x[1] * x[2] + 2 * pow(x[3], 2) + 5 * x[6] - 11 * x[7];

	ind.cost = y;
	for (int i = 0; i < num_cons; ++i) {
		ind.constraint[i] = -1 * c[i];
	}

	return ind;
}


individual DOC::DOC5(individual ind) {
	double g, g_temp;
	vector<long double> x = ind.var, y = ind.cost, c = ind.constraint;
	// basic multi - objective problem
	g_temp = x[1];
	g = g_temp - 193.724510070035 + 1;
	y[0] = x[0];
	y[1] = g*(1 - sqrt(y[0])/ g);

	//constraints in objective space
	c[0] = max(-1 * (y[0] + y[1] - 1), (long double)0.0);
	c[1] = max(-1 * (y[0] + y[1] - 1 - abs(sin(10 * M_PI*(y[0] - y[1] + 1)))), (long double)0.0);
	c[2] = max((y[0] - 0.8)*(y[1] - 0.6), (long double)0.0);

	//constraints in decision space
	c[3] = -1 * x[1] + 35 * pow(x[2], 0.6) + 35 * pow(x[3], 0.6);
	c[4] = abs(-300 * x[3] + 7500 * x[5] - 7500 * x[6] - 25 * x[4] * x[5] + 25 * x[4] * x[6] + x[3] * x[4]) - 0.0001;
	c[5] = abs(100 * x[2] + 155.365 * x[4] + 2500 * x[7] - x[2] * x[4] - 25 * x[4] * x[7] - 15536.5) - 0.0001;
	c[6] = abs(-1 * x[5] + log(-1 * x[4] + 900)) - 0.0001;
	c[7] = abs(-1 * x[6] + log(x[4] + 300)) - 0.0001;
	c[8] = abs(-1 * x[7] + log(-2 * x[4] + 700)) - 0.0001;

	ind.cost = y;
	for (int i = 0; i < num_cons; ++i) {
		ind.constraint[i] = -1 * c[i];
	}

	return ind;
}


individual DOC::DOC6(individual ind) {
	double g, g_temp;
	vector<long double> x = ind.var, y = ind.cost, c = ind.constraint;
	// basic multi - objective problem
	g_temp = pow(x[1], 2) + pow(x[2], 2) + x[1] * x[2] - 14 * x[1] - 16 * x[2] + pow(x[3] - 10, 2) + 4 * pow(x[4] - 5, 2) +
		pow(x[5] - 3, 2) + 2 * pow(x[6] - 1, 2) + 5 * pow(x[7], 2) + 7 * pow(x[8] - 11, 2) + 2 * pow(x[9] - 10, 2) + pow(x[10] - 7, 2) + 45;
	g = g_temp - 24.3062090681 + 1;
	y[0] = x[0];
	y[1] = g*(1 - sqrt(y[0]) / g);

	//constraints in objective space
	c[0] = max(-1 * (y[0] + y[1] - 1), (long double)0.0);
	c[1] = max(-1 * (y[0] - 0.5)*(y[0] + y[1] - 1 - abs(sin(10 * M_PI*(y[0] - y[1] + 1)))), (long double)0.0);

	//constraints in decision space
	c[2] = -105 + 4 * x[1] + 5 * x[2] - 3 * x[7] + 9 * x[8];
	c[3] = 10 * x[1] - 8 * x[2] - 17 * x[7] + 2 * x[8];
	c[4] = -8 * x[1] + 2 * x[2] + 5 * x[9] - 2 * x[10] - 12;
	c[5] = 3 * pow(x[1] - 2, 2) + 4 * pow(x[2] - 3, 2) + 2 * pow(x[3], 2) - 7 * x[4] - 120;
	c[6] = 5 * pow(x[1], 2) + 8 * x[2] + pow(x[3] - 6, 2) - 2 * x[4] - 40;
	c[7] = pow(x[1], 2) + 2 * pow(x[2] - 2,2) - 2 * x[1]*x[2] + 14 * x[5] - 6 * x[6];
	c[8] = 0.5 * pow(x[1] - 8, 2) + 2 * pow(x[2] - 4, 2) + 3 * pow(x[5], 2) - x[6] - 30;
	c[9] = -3 * x[1] + 6 * x[2] + 12 * pow(x[9] - 8, 2) - 7 * x[10];

	ind.cost = y;
	for (int i = 0; i < num_cons; ++i) {
		ind.constraint[i] = -1 * c[i];
	}

	return ind;
}


individual DOC::DOC7(individual ind) {
	double g, g_temp;
	vector<long double> x = ind.var, y = ind.cost, c = ind.constraint;
	double c1[10] = { -6.089, -17.164, -34.054, -5.914, -24.721, -14.986, -24.1, -10.708, -26.662, -22.179 };
	// basic multi - objective problem
	g_temp = 0;
	double sum = pow(10,-20);
	for (int i = 0; i < 10; ++i) {
		sum += x[1 + i];
	}
	for (int i = 0; i < 10; ++i) {
		g_temp += x[1 + i] * (c[i] + log(x[1 + i] / sum));
	}
	g = g_temp + 47.7648884595 + 1;
	y[0] = x[0];
	y[1] = g*(1 - sqrt(y[0]) / g);

	//constraints in objective space
	c[0] = max(-1 * (y[0] + y[1] - 1), (long double)0.0);
	c[1] = max(-1 * (y[0] - 0.5)*(y[0] + y[1] - 1 - abs(sin(10 * M_PI*(y[0] - y[1] + 1)))), (long double)0.0);
	c[2] = max(-1 * (abs(-1 * y[0] + y[1]) / sqrt(2) - 0.1 / sqrt(2)), (long double)0.0);

	//constraints in decision space
	c[3] = abs(x[1] + 2 * x[2] + 2 * x[3] + x[6] + x[10] - 2) - 0.0001;
	c[4] = abs(x[4] + 2 * x[5] + x[6] + x[7] - 1) - 0.0001;
	c[5] = abs(x[3] + x[7] + x[8] + 2 * x[9] + x[10] - 1) - 0.0001;

	ind.cost = y;
	for (int i = 0; i < num_cons; ++i) {
		ind.constraint[i] = -1 * c[i];
	}

	return ind;
}


individual DOC::DOC8(individual ind) {
	double g, g_temp;
	vector<long double> x = ind.var, y = ind.cost, c = ind.constraint;
	// basic multi - objective problem
	g_temp = x[2] + x[3] + x[4];
	g = g_temp - 7049.2480205286 + 1;
	y[0] = (x[0]*x[1])*g;
	y[1] = (x[0]*(1 - x[1]))*g;
	y[2] = (1 - x[0])*g;

	//constraints in objective space
	c[0] = max(-1*(y[2] - 0.4)*(y[2] - 0.6), (long double)0.0);

	//constraints in decision space
	c[1] = -1 + 0.0025 * (x[5] + x[7]);
	c[2] = -1 + 0.0025 * (x[6] + x[8] - x[5]);
	c[3] = -1 + 0.01 * (x[9] - x[6]);
	c[4] = -1 * x[2] * x[7] + 833.33252 * x[5] + 100 * x[2] - 83333.333;
	c[5] = -1 * x[3] * x[8] + 1250 * x[6] + x[3] * x[5] - 1250 * x[5];
	c[6] = -1 * x[4] * x[9] + 1250000 + x[4] * x[6] - 2500 * x[6];

	ind.cost = y;
	for (int i = 0; i < num_cons; ++i) {
		ind.constraint[i] = -1 * c[i];
	}

	return ind;
}


individual DOC::DOC9(individual ind) {
	double g, g_temp;
	vector<long double> x = ind.var, y = ind.cost, c = ind.constraint;
	// basic multi - objective problem
	g_temp = -0.5 * (x[2]*x[5] - x[3]*x[4] + x[4]*x[10] - x[6]*x[10] + x[6]*x[9] - x[7]*x[8]);
	g = g_temp + 0.8660254038 + 1;
	y[0] = cos(0.5*M_PI*x[0])*cos(0.5*M_PI*x[1])*g;
	y[1] = cos(0.5*M_PI*x[0])*sin(0.5*M_PI*x[1])*g;
	y[2] = sin(0.5*M_PI*x[0])*g;

	//constraints in objective space
	c[0] = max(-1 * (pow(y[0], 2) + pow(y[1], 2) - 1), (long double)0.0);

	//constraints in decision space
	c[1] = pow(x[4], 2) + pow(x[5], 2) - 1;
	c[2] = pow(x[10], 2) - 1;
	c[3] = pow(x[6], 2) + pow(x[7], 2) - 1;
	c[4] = pow(x[2], 2) + pow(x[3] - x[10], 2) - 1;
	c[5] = pow(x[2] - x[6], 2) + pow(x[3] - x[7], 2) - 1;
	c[6] = pow(x[2] - x[8], 2) + pow(x[3] - x[9], 2) - 1;
	c[7] = pow(x[4] - x[6], 2) + pow(x[5] - x[7], 2) - 1;
	c[8] = pow(x[4] - x[8], 2) + pow(x[5] - x[9], 2) - 1;
	c[9] = pow(x[8], 2) + pow(x[9] - x[10], 2) - 1;
	c[10] = x[3] * x[4] - x[2] * x[5];
	c[11] = -1*x[4]*x[10];
	c[12] = x[6]*x[10];
	c[13] = x[7]*x[8] - x[6]*x[9];

	ind.cost = y;
	for (int i = 0; i < num_cons; ++i) {
		ind.constraint[i] = -1 * c[i];
	}

	return ind;
}