//	cec09.cpp

#include "CDTLZ.h"
#include <math.h>
#include<fstream>
#include<string>
#include<sstream>
#define M_PI  3.1415926535897932384626433832795
#define MYSIGN(x) ((x)>0?1.0:-1.0)

/****************************************************************************/
// constraint test instances
/****************************************************************************/
void CDTLZ::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	switch (swit) {
	case 1://2 C1-DTLZ1
		tmp_num_obj = 2;
		tmp_num_cons = tmp_num_obj;
		tmp_dim = tmp_num_obj + 5 - 1;
		break;

	case 2://2 C1-DTLZ3
		tmp_num_obj = 2;
		tmp_num_cons = 1;
		tmp_dim = tmp_num_obj + 10 - 1;
		break;

	case 3://2 C2-DTLZ2
		tmp_num_obj = 2;
		tmp_num_cons = 1;
		tmp_dim = tmp_num_obj + 10 - 1;
		break;

	case 4://2 C3-DTLZ1
		tmp_num_obj = 2;
		tmp_num_cons = 1;
		tmp_dim = tmp_num_obj + 5 - 1;
		break;

	case 5://2 C3-DTLZ4
		tmp_num_obj = 2;
		tmp_num_cons = 1;
		tmp_dim = tmp_num_obj + 10 - 1;
		break;

	case 6://2 C3-DTLZ4
		tmp_num_obj = 2;
		tmp_num_cons = 2;
		tmp_dim = tmp_num_obj + 10 - 1;
		break;

	default:
		break;
	}


	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	dim = tmp_dim;
	problemID = swit;
}


void CDTLZ::range_setting(int swit, vector<vector<double>>& range) {
	range.resize(dim);
	for (int i = 0; i < dim; ++i) {
		range[i].resize(2);
		range[i][0] = 0;
		range[i][1] = 1;
	}
}



individual CDTLZ::problems(int swit, individual ind) {

	switch (swit) {
	case 1:
		ind = C1_DTLZ1(ind);
		break;

	case 2:
		ind = C1_DTLZ3(ind);
		break;

	case 3:
		ind = C2_DTLZ2(ind);
		break;

	case 4:
		ind = C3_DTLZ1(ind);
		break;

	case 5:
		ind = C3_DTLZ4(ind);
		break;

	case 6:
		ind = C1C3_DTLZ3(ind);
		break;

	default:
		break;
	}

	return ind;
}


individual CDTLZ::DTLZ1(individual ind) {
	int k = 5;
	double g = 0;
	/*
	for (int i = 0; i < ind.var.size(); ++i) {
	ind.var[i] = 0.5;
	}*/

	//g function
	for (int i = 0; i < k; ++i) {
		g += (ind.var[i + num_obj - 1] - 0.5)*(ind.var[i + num_obj - 1] - 0.5) - cos(20.0*M_PI*(ind.var[i + num_obj - 1] - 0.5));
	}


	g = 100.0 * (k + g);


	//f1
	ind.cost[0] = 0.5*(1.0 + g);
	for (int i = 0; i < num_obj - 1; ++i) {
		ind.cost[0] = ind.cost[0] * ind.var[i];
	}


	//f
	for (int i = 1; i < num_obj; ++i) {
		ind.cost[i] = 0.5*(1.0 + g);
		for (int j = 0; j < num_obj - 1 - i; ++j) {
			ind.cost[i] = ind.cost[i] * ind.var[j];
		}
		ind.cost[i] = ind.cost[i] * (1.0 - ind.var[num_obj - 1 - i]);
	}

	return ind;
}


individual CDTLZ::DTLZ2(individual ind) {
	int k = 10;
	double g = 0;

	//g function
	for (int i = 0; i < k; ++i) {
		g += pow(ind.var[i + num_obj - 1] - 0.5, 2);
	}


	//f1
	ind.cost[0] = 1.0 + g;
	for (int i = 0; i < num_obj - 1; ++i) {
		ind.cost[0] = ind.cost[0] * cos(M_PI / 2.0*ind.var[i]);
	}


	//f
	for (int i = 1; i < num_obj; ++i) {
		ind.cost[i] = 1.0 + g;
		for (int j = 0; j < num_obj - 1 - i; ++j) {
			ind.cost[i] = ind.cost[i] * cos(M_PI / 2.0*ind.var[j]);

		}
		ind.cost[i] = ind.cost[i] * sin(M_PI / 2.0*ind.var[num_obj - 1 - i]);
	}


	return ind;
}


individual CDTLZ::DTLZ3(individual ind) {
	int k = 10;
	double g = 0;


	//g function
	for (int i = 0; i < k; ++i) {
		g += (ind.var[i + num_obj - 1] - 0.5)*(ind.var[i + num_obj - 1] - 0.5) - cos(20.0*M_PI*(ind.var[i + num_obj - 1] - 0.5));
	}
	g = 100 * (k + g);


	//f1
	ind.cost[0] = 1.0 + g;
	for (int i = 0; i < num_obj - 1; ++i) {
		ind.cost[0] = ind.cost[0] * cos(M_PI / 2.0*ind.var[i]);
	}


	//f
	for (int i = 1; i < num_obj; ++i) {
		ind.cost[i] = 1.0 + g;
		for (int j = 0; j < num_obj - 1 - i; ++j) {
			ind.cost[i] = ind.cost[i] * cos(M_PI / 2.0*ind.var[j]);
		}
		ind.cost[i] = ind.cost[i] * sin(M_PI / 2.0*ind.var[num_obj - 1 - i]);
	}


	return ind;
}


individual CDTLZ::DTLZ4(individual ind) {
	int k = 10;
	double g = 0;
	double alpha = 100;

	//g function
	for (int i = 0; i < k; ++i) {
		g += pow(ind.var[i + num_obj - 1] - 0.5, 2);
	}


	//f1
	ind.cost[0] = 1.0 + g;
	for (int i = 0; i < num_obj - 1; ++i) {
		ind.cost[0] = ind.cost[0] * cos(M_PI / 2.0*pow(ind.var[i], alpha));
	}


	//f
	for (int i = 1; i < num_obj; ++i) {
		ind.cost[i] = 1.0 + g;
		for (int j = 0; j < num_obj - 1 - i; ++j) {
			ind.cost[i] = ind.cost[i] * cos(M_PI / 2.0*pow(ind.var[j], alpha));
		}
		ind.cost[i] = ind.cost[i] * sin(M_PI / 2.0*pow(ind.var[num_obj - 1 - i], 100));
	}


	return ind;
}


individual CDTLZ::C1_DTLZ1(individual ind) {

	ind = DTLZ1(ind);
	double con = 1;

	for (int i = 0; i < num_obj; ++i) {
		if (i == num_obj - 1) {
			ind.constraint[i] -= ind.cost[i] / 0.6;
		}
		else {
			ind.constraint[i] -= ind.cost[i] / 0.5;
		}
	}


	return ind;
}


individual CDTLZ::C1_DTLZ3(individual ind) {
	double r[] = { 0, 0, 9, 9, 0, 12.5, 0, 0, 12.5, 0, 15, 15 };

	ind = DTLZ3(ind);

	double pow = 0;
	double left = 0;
	double right = 0;

	for (int i = 0; i < num_obj; ++i) {
		pow += ind.cost[i] * ind.cost[i];
	}
	left = pow - 16;
	right = pow - r[num_obj] * r[num_obj];


	ind.constraint[0] = left*right;

	return ind;
}


individual CDTLZ::C2_DTLZ2(individual ind) {
	double r[] = { 0, 0, 0.15, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };

	ind = DTLZ2(ind);
	ind.constraint[0] = 0;

	double tmp = INT_MAX;
	double right = 0;
	double left = 0;


	for (int i = 0; i < num_obj; ++i) {
		right = pow(ind.cost[i] - 1.0, 2);

		for (int j = 0; j < num_obj; ++j) {
			if (i != j) {
				right += pow(ind.cost[j], 2);
			}
		}
		right -= r[num_obj] * r[num_obj];

		if (right < tmp) {
			tmp = right;
		}
	}
	right = tmp;

	for (int i = 0; i < num_obj; ++i) {
		left += pow(ind.cost[i] - 1.0 / sqrt(num_obj), 2);
	}
	left -= r[num_obj] * r[num_obj];


	if (right <= left) {
		ind.constraint[0] = right;
	}
	else {
		ind.constraint[0] = left;
	}


	return ind;
}


individual CDTLZ::C3_DTLZ1(individual ind) {
	double tmp = 0;
	ind = DTLZ1(ind);
	ind.constraint[0] = 0;


	for (int i = 0; i < num_obj; ++i) {
		tmp = ind.cost[i] - 1;

		for (int j = 0; j < num_obj; ++j) {
			if (i != j) {
				tmp += ind.cost[j] / 0.5;
			}
		}

		if (-1.0*tmp >= 0) {
			ind.constraint[0] += tmp;
		}
	}

	return ind;
}


individual CDTLZ::C3_DTLZ4(individual ind) {
	double tmp = 0;
	ind = DTLZ4(ind);
	ind.constraint[0] = 0;


	for (int i = 0; i < num_obj; ++i) {
		tmp = pow(ind.cost[i], 2) / 4.0 - 1;

		for (int j = 0; j < num_obj; ++j) {
			if (i != j) {
				tmp += pow(ind.cost[j], 2);
			}
		}

		if (-1.0*tmp >= 0) {
			ind.constraint[0] += tmp;
		}
	}



	return ind;
}


individual CDTLZ::C1C3_DTLZ3(individual ind) {
	double r[] = { 0, 0, 12.5, 9, 0, 12.5, 0, 0, 12.5, 0, 15, 15 };

	ind = DTLZ3(ind);

	double pow = 0;
	double left = 0;
	double right = 0;

	for (int i = 0; i < num_obj; ++i) {
		pow += ind.cost[i] * ind.cost[i];
	}
	left = pow - 36;
	right = pow - r[num_obj] * r[num_obj];


	ind.constraint[0] = left*right;
	ind.constraint[1] = pow - 4;

//	ind.constraint[0] = 0;
//	ind.constraint[1] = 0;

	return ind;
}


void CDTLZ::fileout(int swit, int gen, int reps, population pop) {

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
//			if(pow(pop.ind[k].cost[0],2)+ pow(pop.ind[k].cost[1], 2)>4){
				for (int i = 0; i < num_obj - 1; ++i) {
					fout1 << max(pop.ind[k].cost[i], -1 * pop.ind[k].cost[i]) << "\t";
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