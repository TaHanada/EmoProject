//	cec09.cpp

#include "cec09.h"
#include <math.h>
#define PI  3.1415926535897932384626433832795
#define MYSIGN(x) ((x)>0?1.0:-1.0)

/****************************************************************************/
// constraint test instances
/****************************************************************************/
void CEC09::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	switch (swit) {
	case 1:
		tmp_num_obj = 2;
		tmp_num_cons = 2;
		tmp_dim = 4;
		break;

	case 2:
		tmp_num_obj = 2;
		tmp_num_cons = 1;
		tmp_dim = 4;
		break;

	case 3:
		tmp_num_obj = 2;
		tmp_num_cons = 1;
		tmp_dim = 4;
		break;

	case 4:
		tmp_num_obj = 2;
		tmp_num_cons = 1;
		tmp_dim = 4;
		break;

	case 5:
		tmp_num_obj = 2;
		tmp_num_cons = 1;
		tmp_dim = 4;
		break;

	case 6:
		tmp_num_obj = 2;
		tmp_num_cons = 1;
		tmp_dim = 4;
		break;

	case 7:
		tmp_num_obj = 2;
		tmp_num_cons = 1;
		tmp_dim = 4;
		break;

	case 8:
		tmp_num_obj = 2;
		tmp_num_cons = 2;
		tmp_dim = 4;
		break;

	case 9:
		tmp_num_obj = 2;
		tmp_num_cons = 2;
		tmp_dim = 4;
		break;

	default:
		break;
	}


	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	dim = tmp_dim;
}


void CEC09::range_setting(int swit, vector<vector<double>>& range) {
	range.resize(dim);

	switch (swit) {
	case 1: 
		for (int i = 0; i < dim; ++i) {
			range[i].resize(2);
			if (i < 1) {
				range[i][0] = 0;
				range[i][1] = 1;
			}
			else {
				range[i][0] = -5;
				range[i][1] = 5;
			}

		}
		break;

	case 2: 
		for (int i = 0; i < dim; ++i) {
			range[i].resize(2);
			if (i < 1) {
				range[i][0] = 0;
				range[i][1] = 1;
			}
			else {
				range[i][0] = -5;
				range[i][1] = 5;
			}
		}
		break;
	

	case 3: 
		for (int i = 0; i < dim; ++i) {
			range[i].resize(2);
			if (i < 1) {
				range[i][0] = 0;
				range[i][1] = 1;
			}
			else {
				range[i][0] = -5;
				range[i][1] = 5;
			}
		}
		break;
	

	case 4: 
		for (int i = 0; i < dim; ++i) {
			range[i].resize(2);
			if (i < 1) {
				range[i][0] = 0;
				range[i][1] = 1;
			}
			else {
				range[i][0] = -5;
				range[i][1] = 5;
			}
		}
		break;
	

	case 5: 
		for (int i = 0; i < dim; ++i) {
			range[i].resize(2);
			if (i < 1) {
				range[i][0] = 0;
				range[i][1] = 1;
			}
			else {
				range[i][0] = -5;
				range[i][1] = 5;
			}
		}
		break;
	

	case 6: 
		for (int i = 0; i < dim; ++i) {
			range[i].resize(2);
			if (i < 1) {
				range[i][0] = 0;
				range[i][1] = 1;
			}
			else {
				range[i][0] = -5;
				range[i][1] = 5;
			}
		}
		break;
	

	case 7: 
		for (int i = 0; i < dim; ++i) {
			range[i].resize(2);
			if (i < 1) {
				range[i][0] = 0;
				range[i][1] = 1;
			}
			else {
				range[i][0] = -5;
				range[i][1] = 5;
			}
		}
		break;
	

	case 8: 
		for (int i = 0; i < dim; ++i) {
			range[i].resize(2);
			if (i < 1) {
				range[i][0] = 0;
				range[i][1] = 1;
			}
			else {
				range[i][0] = -5;
				range[i][1] = 5;
			}
		}
		break;
	

	case 9: 
		for (int i = 0; i < dim; ++i) {
			range[i].resize(2);
			if (i < 1) {
				range[i][0] = 0;
				range[i][1] = 1;
			}
			else {
				range[i][0] = -5;
				range[i][1] = 5;
			}
		}
		break;

	default:
		break;
	}


}



individual CEC09::problems(int swit, individual ind) {

	switch (swit) {
	case 1:
		ind = CTP1(ind);
		break;

	case 2:
		ind = CTP2(ind);
		break;

	case 3:
		ind = CTP3(ind);
		break;

	case 4:
		ind = CTP4(ind);
		break;

	case 5:
		ind = CTP5(ind);
		break;

	case 6:
		ind = CTP6(ind);
		break;

	case 7:
		ind = CTP7(ind);
		break;

	case 8:
		ind = CTP8(ind);
		break;

	case 9:
		ind = CTP9(ind);
		break;

	default:
		break;
	}

	return ind;
}



long double CEC09::rastrigin(vector<long double> var) {
	double value = 0;

	//Rastrigin's function
	for (int i = 1; i < var.size(); ++i) {
		value += pow(var[i], 2) - 10.0*cos(4 * PI*var[i]) + 10.0;
	}

	value += 1;
	return value;
}

individual CEC09::CTP1(individual ind)
{
	double a1 = 0.858, a2 = 0.728;
	double b1 = 0.541, b2 = 0.295;

	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;
	int nx = ind.var.size();
	long double ras = rastrigin(x);

	f[0] = x[0];
	f[1] = ras*exp(-1.0*f[0] / ras);

	ind.constraint[0] = f[1] - a1*exp(-1.0*b1*f[0]);
	ind.constraint[1] = f[1] - a2*exp(-1.0*b2*f[0]);

	ind.cost = f;
	return ind;
}

individual CEC09::CTP2(individual ind)
{
	double theta = -0.2*PI, a = 0.2;
	double b = 10, e = 1;
	double c = 1, d = 6;

	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;
	int nx = ind.var.size();
	double ras = rastrigin(x);

	f[0] = x[0];
	f[1] = ras*(1.0 - sqrt(f[0] / ras));
	ind.constraint[0] = cos(theta)*(f[1] - e) - sin(theta)*f[0] - a*pow((abs(sin(b*PI*pow(sin(theta)*(f[1] - e) + cos(theta)*f[0], c)))), d);
	ind.cost = f;

	return ind;
}

individual CEC09::CTP3(individual ind)
{
	double theta = -0.2*PI, a = 0.1;
	double b = 10, e = 1;
	double c = 1, d = 0.5;

	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;
	int nx = ind.var.size();
	long double ras = rastrigin(x);

	f[0] = x[0];
	f[1] = ras*(1.0 - sqrt(f[0] / ras));
	ind.constraint[0] = cos(theta)*(f[1] - e) - sin(theta)*f[0] - a*pow(abs(sin(b*PI*pow(sin(theta)*(f[1] - e) + cos(theta)*f[0], c))), d);

	ind.cost = f;

	return ind;
}

individual CEC09::CTP4(individual ind)
{
	double theta = -0.2*PI, a = 0.75;
	double b = 10, e = 1;
	double c = 1, d = 0.5;

	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;
	int nx = ind.var.size();
	long double ras = rastrigin(x);

	f[0] = x[0];
	f[1] = ras*(1.0 - sqrt(f[0] / ras));
	ind.constraint[0] = cos(theta)*(f[1] - e) - sin(theta)*f[0] - a*pow(abs(sin(b*PI*pow(sin(theta)*(f[1] - e) + cos(theta)*f[0], c))), d);

	ind.cost = f;

	return ind;
}

individual CEC09::CTP5(individual ind)
{
	double theta = -0.2*PI, a = 0.1;
	double b = 10, e = 1;
	double c = 2, d = 0.5;

	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;

	int nx = ind.var.size();
	long double ras = rastrigin(x);

	f[0] = x[0];
	f[1] = ras*(1.0 - sqrt(f[0] / ras));
	ind.constraint[0] = cos(theta)*(f[1] - e) - sin(theta)*f[0] - a*pow(abs(sin(b*PI*pow(sin(theta)*(f[1] - e) + cos(theta)*f[0], c))), d);

	ind.cost = f;

	return ind;
}

individual CEC09::CTP6(individual ind)
{
	double theta = 0.1*PI, a = 40;
	double b = 0.5, e = -2;
	double c = 1, d = 2;

	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;
	int nx = ind.var.size();
	long double ras = rastrigin(x);

	f[0] = x[0];
	f[1] = ras*(1.0 - sqrt(f[0] / ras));
	ind.constraint[0] = cos(theta)*(f[1] - e) - sin(theta)*f[0] - a*pow(abs(sin(b*PI*pow(sin(theta)*(f[1] - e) + cos(theta)*f[0], c))), d);

	ind.cost = f;

	return ind;
}

individual CEC09::CTP7(individual ind)
{
	double theta = -0.05*PI, a = 40;
	double b = 5, e = 0;
	double c = 1, d = 6;

	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;
	int nx = ind.var.size();
	long double ras = rastrigin(x);

	f[0] = x[0];
	f[1] = ras*(1.0 - sqrt(f[0] / ras));
	ind.constraint[0] = cos(theta)*(f[1] - e) - sin(theta)*f[0] - a*pow(abs(sin(b*PI*pow(sin(theta)*(f[1] - e) + cos(theta)*f[0], c))), d);

	ind.cost = f;

	return ind;
}

individual CEC09::CTP8(individual ind)
{
	double theta1 = 0.1*PI, a1 = 40;
	double b1 = 0.5, e1 = -2;
	double c1 = 1, d1 = 2;

	double theta2 = -0.05*PI, a2 = 40;
	double b2 = 2, e2 = 0;
	double c2 = 1, d2 = 6;

	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;
	double con[2];
	int nx = ind.var.size();
	double ras = rastrigin(x);

	f[0] = x[0];
	f[1] = ras*(1.0 - sqrt(f[0] / ras));
	ind.constraint[0] = cos(theta1)*(f[1] - e1) - sin(theta1)*f[0] - a1*pow(abs(sin(b1*PI*pow(sin(theta1)*(f[1] - e1) + cos(theta1)*f[0], c1))), d1);
	ind.constraint[1] = cos(theta2)*(f[1] - e2) - sin(theta2)*f[0] - a2*pow(abs(sin(b2*PI*pow(sin(theta2)*(f[1] - e2) + cos(theta2)*f[0], c2))), d2);

	ind.cost = f;

	return ind;
}

individual CEC09::CTP9(individual ind)
{
	double theta1 = 0.4*PI, a1 = 4000;
	double b1 = 0.5, e1 = -2;
	double c1 = 1, d1 = 2;

	double theta2 = -0.05*PI, a2 = 4000;
	double b2 = 2, e2 = 0;
	double c2 = 1, d2 = 6;

	vector<long double> x = ind.var;
	vector<long double> f = ind.cost;
	double con[2];
	int nx = ind.var.size();
	long double ras = rastrigin(x);

	f[0] = x[0];
	f[1] = ras*(1.0 - sqrt(f[0] / ras));
	ind.constraint[0] = cos(theta1)*(f[1] - e1) - sin(theta1)*f[0] - a1*pow(abs(sin(b1*PI*pow(sin(theta1)*(f[1] - e1) + cos(theta1)*f[0], c1))), d1);
	ind.constraint[1] = cos(theta2)*(f[1] - e2) - sin(theta2)*f[0] - a2*pow(abs(sin(b2*PI*pow(sin(theta2)*(f[1] - e2) + cos(theta2)*f[0], c2))), d2);

	ind.cost = f;

	return ind;
}