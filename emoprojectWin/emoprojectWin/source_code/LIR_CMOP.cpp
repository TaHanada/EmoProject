//	LIR_CMOP.cpp

#include "LIR_CMOP.h"
#include <math.h>
#include<fstream>
#include<string>
#include<sstream>
#define M_PI  3.1415926535897932384626433832795
#define MYSIGN(x) ((x)>0?1.0:-1.0)

/****************************************************************************/
// constraint test instances
/****************************************************************************/
void LIR_CMOP::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	
	tmp_num_obj = 2;
	tmp_dim = 30;
	if (swit==3||swit==4||swit == 7 || swit == 8) {
		tmp_num_cons = 3;
	}
	else {
		tmp_num_cons = 2;
	}


	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	dim = tmp_dim;
	problemID = swit;
}



void LIR_CMOP::range_setting(int swit, vector<vector<double>>& range) {
	range.resize(dim);
	for (int i = 0; i < dim; ++i) {
		range[i].resize(2);
		range[i][0] = 0;
		range[i][1] = 1;
	}
}



individual LIR_CMOP::problems(int swit, individual ind) {

	switch (swit) {
	case 1:
		ind = LIR_CMOP1(ind);
		break;

	case 2:
		ind = LIR_CMOP2(ind);
		break;

	case 3:
		ind = LIR_CMOP3(ind);
		break;

	case 4:
		ind = LIR_CMOP4(ind);
		break;

	case 5:
		ind = LIR_CMOP5(ind);
		break;

	case 6:
		ind = LIR_CMOP6(ind);
		break;

	case 7:
		ind = LIR_CMOP7(ind);
		break;

	case 8:
		ind = LIR_CMOP8(ind);
		break;

	case 9:
		ind = LIR_CMOP9(ind);
		break;

	case 10:
		ind = LIR_CMOP10(ind);
		break;

	case 11:
		ind = LIR_CMOP11(ind);
		break;

	case 12:
		ind = LIR_CMOP12(ind);
		break;

	case 13:
		ind = LIR_CMOP13(ind);
		break;

	case 14:
		ind = LIR_CMOP14(ind);
		break;

	default:
		break;
	}

	return ind;
}



individual LIR_CMOP::LIR_CMOP1(individual ind) {

	double f1 = 0;
	double f2 = 0;
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double a = 0.51;
	double b = 0.5;


	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*M_PI*x[0]), 2);
	}

	// f1
	f1 = x[0] + g1;

	//f2
	f2 = 1.0 - pow(x[0], 2) + g2;

	ind.constraint[0] = (a - g1)*(g1 - b);
	ind.constraint[1] = (a - g2)*(g2 - b);

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP2(individual ind) {

	double f1 = 0;
	double f2 = 0;
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double a = 0.51;
	double b = 0.5;


	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*M_PI*x[0]), 2);
	}

	// f1
	f1 = x[0] + g1;

	//f2
	f2 = 1.0 - sqrt(x[0]) + g2;

	ind.constraint[0] = (a - g1)*(g1 - b);
	ind.constraint[1] = (a - g2)*(g2 - b);

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP3(individual ind) {

	double f1 = 0;
	double f2 = 0;
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double a = 0.51;
	double b = 0.5;
	double c = 20;

	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*M_PI*x[0]), 2);
	}

	// f1
	f1 = x[0] + g1;

	//f2
	f2 = 1.0 - pow(x[0], 2) + g2;

	ind.constraint[0] = (a - g1)*(g1 - b);
	ind.constraint[1] = (a - g2)*(g2 - b);
	ind.constraint[2] = sin(c*M_PI*x[0])-0.5;

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP4(individual ind) {

	double f1 = 0;
	double f2 = 0;
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double a = 0.51;
	double b = 0.5;
	double c = 20;

	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*M_PI*x[0]), 2);
	}

	// f1
	f1 = x[0] + g1;

	//f2
	f2 = 1.0 - sqrt(x[0]) + g2;

	ind.constraint[0] = (a - g1)*(g1 - b);
	ind.constraint[1] = (a - g2)*(g2 - b);
	ind.constraint[2] = sin(c*M_PI*x[0])-0.5;

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP5(individual ind) {

	double f1 = ind.cost[0];
	double f2 = ind.cost[1];
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double p[2] = { 1.6, 2.5 };
	double q[2] = { 1.6,2.5 };
	double a[2] = { 2,2 };
	double b[2] = { 4,8 };
	double r = 0.1;
	double theta = -0.25*M_PI;
	int k = 2;

	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}

	// f1
	f1 = x[0] + 10 * g1 + 0.7057;

	//f2
	f2 = 1.0 - sqrt(x[0]) + 10 * g2 + 0.7057;

	for (int i = 0; i < k; ++i) {
		ind.constraint[i] = pow((f1 - p[i])*cos(theta) - (f2 - q[i])*sin(theta), 2) / pow(a[i], 2)
			+ pow((f1 - p[i])*sin(theta) + (f2 - q[i])*cos(theta), 2) / pow(b[i], 2)
			- r;
	}

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP6(individual ind) {

	double f1 = ind.cost[0];
	double f2 = ind.cost[1];
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double p[2] = { 1.8, 2.8 };
	double q[2] = { 1.8,2.8 };
	double a[2] = { 2,2 };
	double b[2] = { 8,8 };
	double r = 0.1;
	double theta = -0.25*M_PI;
	int k = 2;

	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}

	// f1
	f1 = x[0] + 10 * g1 + 0.7057;

	//f2
	f2 = 1.0 - pow(x[0], 2) + 10 * g2 + 0.7057;

	for (int i = 0; i < k; ++i) {
		ind.constraint[i] = pow((f1 - p[i])*cos(theta) - (f2 - q[i])*sin(theta), 2) / pow(a[i], 2)
			+ pow((f1 - p[i])*sin(theta) + (f2 - q[i])*cos(theta), 2) / pow(b[i], 2)
			- r;
	}

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP7(individual ind) {

	double f1 = ind.cost[0];
	double f2 = ind.cost[1];
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double p[3] = { 1.2, 2.25, 3.5 };
	double q[3] = { 1.2,2.25,3.5 };
	double a[3] = { 2,2.5,2.5 };
	double b[3] = { 6,12,10 };
	double r = 0.1;
	double theta = -0.25*M_PI;
	int k = 3;

	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}

	// f1
	f1 = x[0] + 10 * g1 + 0.7057;

	//f2
	f2 = 1.0 - sqrt(x[0]) + 10 * g2 + 0.7057;


	for (int i = 0; i < k; ++i) {
		ind.constraint[i] = pow((f1 - p[i])*cos(theta) - (f2 - q[i])*sin(theta), 2) / pow(a[i], 2)
			+ pow((f1 - p[i])*sin(theta) + (f2 - q[i])*cos(theta), 2) / pow(b[i], 2)
			- r;
	}

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP8(individual ind) {
	
	double f1 = ind.cost[0];
	double f2 = ind.cost[1];
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double p[3] = { 1.2, 2.25, 3.5 };
	double q[3] = { 1.2,2.25,3.5 };
	double a[3] = { 2,2.5,2.5 };
	double b[3] = { 6,12,10 };
	double r = 0.1;
	double theta = -0.25*M_PI;
	int k = 3;

	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}

	// f1
//	f1 = x[0] + 10 * g1 + 0.7057;

	//f2
//	f2 = 1.0 - pow(x[0], 2) + 10 * g2 + 0.7057;

	for (int i = 0; i < k; ++i) {
		ind.constraint[i] = pow((f1 - p[i])*cos(theta) - (f2 - q[i])*sin(theta), 2) / pow(a[i], 2)
			+ pow((f1 - p[i])*sin(theta) + (f2 - q[i])*cos(theta), 2) / pow(b[i], 2)
			- r;
	}

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP9(individual ind) {

	double f1 = ind.cost[0];
	double f2 = ind.cost[1];
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double p = 1.4;
	double q = 1.4;
	double a = 1.5;
	double b = 6.0;
	double r = 0.1;
	double alpha = 0.25*M_PI;
	double theta = -0.25*M_PI;
	int k = 3;

	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}

	// f1
//	f1 = 1.7057*x[0] * (10 * g1 + 1.0);

	//f2
//	f2 = 1.7057*(1.0 - pow(x[0],2)) * (10 * g2 + 1.0);


	//c1
	ind.constraint[0] = pow((f1 - p)*cos(theta) - (f2 - q)*sin(theta), 2) / pow(a, 2)
		+ pow((f1 - p)*sin(theta) + (f2 - q)*cos(theta), 2) / pow(b, 2)
		- r;

	//c2
	ind.constraint[1] = f1*sin(alpha) + f2*cos(alpha)
		- sin(4.0*M_PI*(f1*cos(alpha) - f2*sin(alpha))) - 2.0;

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP10(individual ind) {

	double f1 = ind.cost[0];
	double f2 = ind.cost[1];
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double p = 1.1;
	double q = 1.2;
	double a = 2.0;
	double b = 4.0;
	double r = 0.1;
	double alpha = 0.25*M_PI;
	double theta = -0.25*M_PI;
	int k = 3;

	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}

	// f1
//	f1 = 1.7057*x[0] * (10 * g1 + 1.0);

	//f2
//	f2 = 1.7057*(1.0 - sqrt(x[0])) * (10 * g2 + 1.0);


	//c1
	ind.constraint[0] = pow((f1 - p)*cos(theta) - (f2 - q)*sin(theta), 2) / pow(a, 2)
		+ pow((f1 - p)*sin(theta) + (f2 - q)*cos(theta), 2) / pow(b, 2)
		- r;

	//c2
	ind.constraint[1] = f1*sin(alpha) + f2*cos(alpha)
		- sin(4.0*M_PI*(f1*cos(alpha) - f2*sin(alpha))) - 1.0;

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP11(individual ind) {

	double f1 = ind.cost[0];
	double f2 = ind.cost[1];
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double p = 1.2;
	double q = 1.2;
	double a = 1.5;
	double b = 5.0;
	double r = 0.1;
	double alpha = 0.25*M_PI;
	double theta = -0.25*M_PI;
	int k = 3;

	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}

	// f1
//	f1 = 1.7057*x[0] * (10 * g1 + 1.0);

	//f2
//	f2 = 1.7057*(1.0 - sqrt(x[0])) * (10 * g2 + 1.0);


	//c1
	ind.constraint[0] = pow((f1 - p)*cos(theta) - (f2 - q)*sin(theta), 2) / pow(a, 2)
		+ pow((f1 - p)*sin(theta) + (f2 - q)*cos(theta), 2) / pow(b, 2)
		- r;

	//c2
	ind.constraint[1] = f1*sin(alpha) + f2*cos(alpha)
		- sin(4.0*M_PI*(f1*cos(alpha) - f2*sin(alpha))) - 2.1;

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP12(individual ind) {

	double f1 = ind.cost[0];
	double f2 = ind.cost[1];
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double p = 1.6;
	double q = 1.6;
	double a = 1.5;
	double b = 6.0;
	double r = 0.1;
	double alpha = 0.25*M_PI;
	double theta = -0.25*M_PI;
	int k = 3;

	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}

	// f1
//	f1 = 1.7057*x[0] * (10 * g1 + 1.0);

	//f2
//	f2 = 1.7057*(1.0 - pow(x[0],2)) * (10 * g2 + 1.0);


	//c1
	ind.constraint[0] = pow((f1 - p)*cos(theta) - (f2 - q)*sin(theta), 2) / pow(a, 2)
		+ pow((f1 - p)*sin(theta) + (f2 - q)*cos(theta), 2) / pow(b, 2)
		- r;

	//c2
	ind.constraint[1] = f1*sin(alpha) + f2*cos(alpha)
		- sin(4.0*M_PI*(f1*cos(alpha) - f2*sin(alpha))) - 2.5;

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP13(individual ind) {

	double f1 = 0;
	double f2 = 0;
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double p = 1.6;
	double q = 1.6;
	double a = 1.5;
	double b = 6.0;
	double r = 0.1;
	double alpha = 0.25*M_PI;
	double theta = -0.25*M_PI;
	int k = 3;

	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}

	// f1
//	f1 = 1.7057*x[0] * (10 * g1 + 1.0);

	//f2
//	f2 = 1.7057*(1.0 - pow(x[0], 2)) * (10 * g2 + 1.0);


	//c1
	ind.constraint[0] = pow((f1 - p)*cos(theta) - (f2 - q)*sin(theta), 2) / pow(a, 2)
		+ pow((f1 - p)*sin(theta) - (f2 - q)*cos(theta), 2) / pow(b, 2)
		- r;

	//c2
	ind.constraint[1] = f1*sin(alpha) + f2*cos(alpha)
		- sin(4.0*M_PI*(f1*cos(alpha) - f2*sin(alpha))) - 2.5;

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


individual LIR_CMOP::LIR_CMOP14(individual ind) {

	double f1 = 0;
	double f2 = 0;
	double g1 = 0;
	double g2 = 0;
	double c1 = 0;
	double c2 = 0;

	double p = 1.6;
	double q = 1.6;
	double a = 1.5;
	double b = 6.0;
	double r = 0.1;
	double alpha = 0.25*M_PI;
	double theta = -0.25*M_PI;
	int k = 3;

	vector<double> x = ind.var;

	//g
	vector<int> j1;
	vector<int> j2;
	for (int i = 1; i < ind.var.size(); ++i) {
		if ((i + 1) % 2 == 0) {
			j2.push_back(i);
		}
		else {
			j1.push_back(i);
		}
	}

	//g1
	for (int i = 0; i < j1.size(); ++i) {
		g1 += pow(x[j1[i]] - sin(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}
	//g2
	for (int i = 0; i < j2.size(); ++i) {
		g2 += pow(x[j2[i]] - cos(0.5*(j1[i] + 1) / 30.0*M_PI*x[0]), 2);
	}

	// f1
	f1 = 1.7057*x[0] * (10 * g1 + 1.0);

	//f2
	f2 = 1.7057*(1.0 - pow(x[0], 2)) * (10 * g2 + 1.0);


	//c1
	ind.constraint[0] = pow((f1 - p)*cos(theta) - (f2 - q)*sin(theta), 2) / pow(a, 2)
		+ pow((f1 - p)*sin(theta) - (f2 - q)*cos(theta), 2) / pow(b, 2)
		- r;

	//c2
	ind.constraint[1] = f1*sin(alpha) + f2*cos(alpha)
		- sin(4.0*M_PI*(f1*cos(alpha) - f2*sin(alpha))) - 2.5;

	ind.cost[0] = f1;
	ind.cost[1] = f2;

	return ind;
}


void LIR_CMOP::fileout(int swit, int gen, int reps, population pop) {

	switch (swit) {
	case 1:
	{
		string name1 = to_string(reps) + "_" + to_string(gen) + "gen.dat";
		ofstream fout("./" + to_string(problemID) + "/tch/offsprings/" + name1);

		for (int k = 0; k < pop.size; ++k) {
			//gen 
			fout << gen << "\t";

			//feasible
			if (pop.ind[k].total_c == 0) {
				fout << "1" << "\t";
			}
			else {
				fout << "0" << "\t";
			}

			//function value
			for (int i = 0; i < num_obj - 1; ++i) {
				fout << max(pop.ind[k].cost[i], -1 * pop.ind[k].cost[i]) << "\t";
			}
			fout << max(pop.ind[k].cost[num_obj - 1], -1 * pop.ind[k].cost[num_obj - 1]) << endl;
		}
	}
	break;

	case 2:
	{
		string name1 = to_string(reps) + "_" + to_string(gen) + "gen.dat";
		ofstream fout1("./" + to_string(problemID) + "/tch/island1/curents/" + name1);
		ofstream fout2("./" + to_string(problemID) + "/tch/island2/curents/" + name1);

		for (int k = 0; k < pop.size; ++k) {

			//feasible
			if (pop.ind[k].total_c == 0) {
				for (int i = 0; i < num_obj-1; ++i) {
					fout1 <<max(pop.ind[k].cost[i], -1 * pop.ind[k].cost[i]) << "\t";
				}
				fout1 << max(pop.ind[k].cost[num_obj - 1], -1 * pop.ind[k].cost[num_obj - 1]) << endl;
			}
			else {
				for (int i = 0; i < num_obj - 1; ++i) {
					fout2 << max(pop.ind[k].cost[i], -1 * pop.ind[k].cost[i]) << "\t";
				}
				fout2 << max(pop.ind[k].cost[num_obj - 1], -1 * pop.ind[k].cost[num_obj - 1]) << endl;
			}
		}
	}

	default:
		break;
	}

}


void LIR_CMOP::constrained_landscape() {
	double f1[12] = { 2,2,1.5,1.5,5,5,6,6,3,3,3,3.5 };
	double f2[12] = { 1.6,1.6,1.6,1.5,4.5,5,6,6,3,3,3,3.5 };

	for (int prob = 12; prob < 13; ++prob) {
//		ofstream fout("ConstrainedLandscape_" + to_string(prob)+".csv");
		ofstream fout("InfeasibleRegion_" + to_string(prob) + ".csv");

		int tmp_num_obj, tmp_num_cons, tmp_dim;
		initialize(prob, tmp_num_obj, tmp_num_cons, tmp_dim);
		individual ind;
		ind.var.resize(dim);
		ind.cost.resize(num_obj);
		ind.constraint.resize(num_cons);

		for (int i = 100; i > 0; --i) {
			for (int j = 0; j < 100; ++j) {
				ind.cost[0] = f1[prob-1] * (double)j / 100.0;
				ind.cost[1] = f2[prob-1] * (double)i / 100.0;
				ind = problems(prob, ind);
				ind.total_c = 100 * num_cons;
/*				for (int c = 0; c < num_cons; ++c) {
					ind.total_c += pow(max(0.0,-1.0*ind.constraint[c]),2)-10*cos(100*M_PI*max(0.0, -1.0*ind.constraint[c]));
					ind.total_c += pow(max(0.0, -1.0*ind.constraint[c]), 2)/100 - 100 * cos(100 * M_PI*max(0.0, -1.0*ind.constraint[c]));
				}*/
				ind.total_c = 0;
				for (int i = 0; i < num_cons; ++i) {
					if (ind.constraint[i] < 0) {
						ind.total_c += -1.0*ind.constraint[i];
					}
				}

			//	cout << ind.cost[0] << " " << ind.cost[1] << " " << ind.constraint[0] << " " << ind.constraint[1] <<" "<< ind.total_c << endl;
//				fout << ind.total_c << ",";
				if (ind.total_c > 0) fout << f1[prob-1] * (double)j / 100.0 << "," << f2[prob-1] * (double)i / 100.0 << endl;
			}
//			fout << endl;
		}

	}


	cout << "fin" << endl;
	cin >> dim;
}


individual LIR_CMOP::CLIR_constraint(individual ind) {
	double totalC = 10 * num_cons;

	for (int c = 0; c < num_cons; ++c) {
		totalC += pow(min(0.0, ind.constraint[c]), 2) / 100 - 10 * cos(100 * M_PI*min(0.0, ind.constraint[c]));
		ind.constraint[c] = 0;
	}
	ind.constraint[0] = -1.0*totalC;

	return ind;
}