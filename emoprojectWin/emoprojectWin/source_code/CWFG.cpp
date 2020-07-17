//	CWFG.cpp

#include "CWFG.h"
#include <math.h>
#include<fstream>
#include<string>
#include<sstream>
#define M_PI  3.1415926535897932384626433832795
#define MYSIGN(x) ((x)>0?1.0:-1.0)

/****************************************************************************/
// constraint test instances
/****************************************************************************/
void CWFG::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	switch (swit) {
	case 1://2 CWFG1_correlation
		tmp_num_obj = 2;
		tmp_num_cons = 2;
		tmp_dim = 24;
		break;

	default:
		break;
	}


	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	c = 5;
	dim = tmp_dim;
	problemID = swit;
}


void CWFG::range_setting(int swit, vector<vector<double>>& range) {
	range.resize(dim);
	for (int i = 0; i < dim; ++i) {
		range[i].resize(2);
		range[i][0] = 0;
		range[i][1] = 1;
	}
}



individual CWFG::problems(int swit, individual ind) {

	switch (swit) {
	case 1:
		ind = CWFG1_correlation(ind);
		break;

	default:
		break;
	}

	return ind;
}


individual CWFG::WFG1(individual ind) {
	int k = num_obj - 1;
	vector<double> z;
	double g = 0;

	for (int i = 0; i < ind.var.size()-c; ++i) {
		z.push_back(ind.var[i] * 2.0*(i + 1));
	}
	ind.cost = WF1(z, k, num_obj);

	return ind;
}



individual CWFG::CWFG1_correlation(individual ind) {
	int k = num_obj - 1;
	double del_ = 10.0;
	int c = 5;
	vector<double> z;


	for (int i = 0; i < ind.var.size(); ++i) {
		z.push_back(ind.var[i] * 2.0*(i + 1));
	}

	//obj
	ind = WFG1(ind);

	//constrained var
	vector< double > cvar = get_constrained_var_WFG1(z, k, num_obj, c, num_cons);

	//correlation constrained function
	double ym = cvar[num_cons] + 0.55 / del_;
	cout << ym << endl;
	for (int i = 0; i < num_cons; ++i) {

		if (i+1 == 1) ind.constraint[i]=cos(del_*ym * M_PI);
		else {
			ind.constraint[i] = cos(del_*ym*M_PI - 1.0 / (double)i*sin(del_*ym*M_PI));
		}
	}
	
	
	return ind;
}


void CWFG::fileout(int swit, int gen, int reps, population pop) {

	switch (swit) {
	case 1:
	{
		string name1 = to_string(reps) + "_" + to_string(gen) + "gen.dat";
		ofstream fout("./" + to_string(problemID) + "/tch/offsprings/" + name1, ios::out | ios::app);

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
			for (int i = 0; i < num_obj; ++i) {
				fout << max(pop.ind[k].cost[i], -1 * pop.ind[k].cost[i]) << "\t";
			}
			for (int i = 0; i < dim - 1; ++i) {
				fout << pop.ind[k].var[i] << "\t";
			}
			fout << pop.ind[k].var[dim - 1] << endl;
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