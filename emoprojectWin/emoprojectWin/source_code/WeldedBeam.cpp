#include"WeldedBeam.h"
#define _USE_MATH_DEFINES
#include<math.h>


individual WB::WeldedBeam(individual ind) {
	//h=0, b=1, l=2, t=3
	double h = ind.var[0], b = ind.var[1], l = ind.var[2], t = ind.var[3];
	double tau_, delta_, pc_;
	double tau1_ = 6000.0 / sqrt(2) / h / l;
	double tau2_ = (6000.0 * (14.0 + 0.5*l)*sqrt(0.25*(pow(l, 2) + pow(h + t, 2)))) / (2 * (0.707*h*l*(pow(l, 2) / 12.0 + 0.25*pow(h + t, 2))));

	tau_ = sqrt(pow(tau1_, 2) + pow(tau2_, 2) + (l*tau1_*tau2_) / sqrt(0.25*(pow(l, 2) + pow(h + t, 2))));
	delta_ = 504000.0 / pow(t, 2) / b;
	pc_ = 64746.022*(1.0 - 0.0282346*t)*t*pow(b, 3);


	//objectives
	ind.cost[0] = 1.10471*h * h * l + 0.04811*t * b * (14.0 + l);
	ind.cost[1] = 2.1952 / pow(t, 3) / b;

	//constraints
	ind.constraint[0] = 13600.0 - tau_;
	ind.constraint[1] = 30000.0 - delta_;
	ind.constraint[2] = b - h;
	ind.constraint[3] = pc_ - 60000.0;

	return ind;
}


void WB::check() {
	individual x;
	x.var.resize(12);//size of phenotype
	x.cost.resize(num_obj);//number of tasks

	for (int i = 0; i < 12; ++i) {
		x.var[i] = 0.5;
	}
	x.var[0] = 0;
	x.var[1] = 1;



	x = WeldedBeam(x);


	cout << x.cost[0] << " " << x.cost[1] << " " << x.cost[2] << endl;
}


individual WB::problems(int swit, individual ind) {

	switch (swit) {
	case 1:
		ind = WeldedBeam(ind);
		break;


	default:
		break;
	}

	return ind;
}