
#include "divide.h"
#include "competition.h"
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <direct.h>
#include <iomanip>
#define PI  3.1415926535897932384626433832795
#define MYSIGN(x) ((x)>0?1.0:-1.0)

void competition::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	switch (swit) {
	case 1:
		tmp_num_obj = 1;
		tmp_num_cons = 22;
		tmp_dim = 32;
		break;
	case 2:
		tmp_num_obj = 5;
		tmp_num_cons = 22;
		tmp_dim = 32;
		break;
	}

	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	dim = tmp_dim;
}

void competition::range_setting(int swit, vector<vector<double>>& range) {
	range.resize(dim);
	for (int i = 0; i < range.size(); i++) {
		range[i].resize(2);
	}

	for (int i = 0; i < 4; ++i) {
		range[i][0] = 1.0;
		range[i][1] = 5.3;
	}
	range[4][0] = 0.1;
	range[4][1] = 0.3;
	for (int i = 5; i < 9; ++i) {
		range[i][0] = -5.0;
		range[i][1] = 30.0;
	}
	for (int i = 9; i < 19; ++i) {
		range[i][0] = 0.005;
		range[i][1] = 0.2;
	}
	for (int i = 19; i < 22; ++i) {
		range[i][0] = -6.3;
		range[i][1] = 0;
	}
	range[22][0] = 6.0;
	range[22][1] = 14.0;
	range[23][0] = 6.0;
	range[23][1] = 20.0;
	range[24][0] = 50;
	range[24][1] = 80;
	range[25][0] = 20;
	range[25][1] = 70;
	for (int i = 26; i < 29; i++) {
		range[i][0] = 3.87;
		range[i][1] = 6.3;
	}
	for (int i = 29; i < 32; i++) {
		range[i][0] = 0.005;
		range[i][1] = 0.1;
	}
}

individual competition::problems(int swit, individual ind) {
	switch (swit) {
	case 1:
		ind = singleobjective(ind);
		break;
	case 2:
		ind = multiobjective(ind);
		break;
	}
	return ind;
	
}

individual competition::singleobjective(individual ind) {
	return ind;
}

individual competition::multiobjective(individual ind) {
	vector<long double> normalize_x = ind.var;
	vector<long double> f = ind.cost;
	vector<long double> z = ind.constraint;
	Divide* obj = Divide::getInstance();

	//normalization
	for (int i = 0; i < normalize_x.size(); i++) {
		if (0 <= i && i <= 3) {
			normalize_x[i] = (normalize_x[i] - 1.0) / (5.3 - 1.0);
		}
		else if (i == 4) {
			normalize_x[i] = (normalize_x[i] - 0.1) / (0.3 - 0.1);
		}
		else if (5 <= i && i <= 8) {
			normalize_x[i] = (normalize_x[i] + 5.0) / (30.0 + 5.0);
		}
		else if (9 <= i && i <= 18) {
			normalize_x[i] = (normalize_x[i] - 0.005) / (0.2 - 0.005);
		}
		else if (19 <= i && i <= 21) {
			normalize_x[i] = (normalize_x[i] + 6.3) / 6.3;
		}
		else if (i == 22) {
			normalize_x[i] = (normalize_x[i] - 6.0) / (14.0 - 6.0);
		}
		else if (i == 23) {
			normalize_x[i] = (normalize_x[i] - 6.0) / (20.0 - 6.0);
		}
		else if (i == 24) {
			normalize_x[i] = (normalize_x[i] - 50.0) / (80.0 - 50.0);
		}
		else if (i == 25) {
			normalize_x[i] = (normalize_x[i] - 20.0) / (70.0 - 20.0);
		}
		else if (26 <= i && i <= 28) {
			normalize_x[i] = (normalize_x[i] - 3.87) / (6.3 - 3.87);
		}
		else if (29 <= i && i <= 31) {
			normalize_x[i] = (normalize_x[i] - 0.005) / (0.1 - 0.005);
		}
	}
	//string extermination = "python windturbine_mop.py ind_data_" + to_string(obj->get_file_divide());
	ofstream outputfile;
	outputfile.open("./ind_data_" + to_string(obj->get_file_divide()) + "/pop_vars_eval.txt");
	outputfile << setprecision(17) << normalize_x[0];
	for (int i = 1; i < normalize_x.size(); i++) {
		outputfile << " \t" << setprecision(17) << normalize_x[i];
	}
	outputfile.close();

	//確認用
	/*ifstream vars_file;
	vars_file.open("./ind_data_" + to_string(obj->get_file_divide()) + "/pop_vars_eval.txt");
	string temp2;
	int index2 = 0;
	while (getline(vars_file, temp2, '\t')) {
		normalize_x[index2] = stod(temp2);
		index2++;
	}
	vars_file.close();
	char dir[255];
	_getcwd(dir, 255);
	cout << "Current Directory : " << dir << endl;*/

	//int ret;
	//ret = system(extermination.c_str());

	//python実行後に出るtxtファイルから，目的関数値と制約関数値を読み込む．
	ifstream objs_file;
	objs_file.open("./ind_data_" + to_string(obj->get_file_divide()) + "/pop_objs_eval.txt");
	string temp;
	int index = 0;
	while (getline(objs_file, temp, '\t')) {
		f[index] = stold(temp);
		index++;
	}
	objs_file.close();
	ifstream cons_file;
	cons_file.open("./ind_data_" + to_string(obj->get_file_divide()) + "/pop_cons_eval.txt");
	index = 0;
	while (getline(cons_file, temp, '\t')) {
		z[index] = stold(temp);
		index++;
	}
	cons_file.close();

	ind.cost = f;
	ind.constraint = z;

	return ind;
}