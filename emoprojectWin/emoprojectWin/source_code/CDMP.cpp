#include"CDMP.h"
#define _USE_MATH_DEFINES
#include<math.h>

float sign(vector<long double> p1, vector<double> p2,vector<double> p3) {
	
	return (p1[0] - p3[0])*(p2[1] - p3[1]) - (p2[0] - p3[0])*(p1[1] - p3[1]);
}

void CDMP::initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim) {
	switch (swit) {
	case 1:
		tmp_num_obj = 3;
		tmp_num_cons = 100;
		tmp_dim = 2;
		break;

	case 2:
		tmp_num_obj = 3;
		tmp_num_cons = 100;
		tmp_dim = 2;
		break;

	case 3:
		tmp_num_obj = 3;
		tmp_num_cons = 100;
		tmp_dim = 2;
		break;

	case 4:
		tmp_num_obj = 3;
		tmp_num_cons = 200;
		tmp_dim = 2;
		break;

	case 5:
		tmp_num_obj = 3;
		tmp_num_cons = 200;
		tmp_dim = 2;
		break;

	case 6:
		tmp_num_obj = 3;
		tmp_num_cons = 200;
		tmp_dim = 2;
		break;

	case 7:
		tmp_num_obj = 3;
		tmp_num_cons = 300;
		tmp_dim = 2;
		break;

	case 8:
		tmp_num_obj = 3;
		tmp_num_cons = 300;
		tmp_dim = 2;
		break;

	case 9:
		tmp_num_obj = 3;
		tmp_num_cons = 300;
		tmp_dim = 2;
		break;


	default:
		break;
	}


	num_obj = tmp_num_obj;
	num_cons = tmp_num_cons;
	dim = tmp_dim;

	//define a, v
	r = 5;
	center.push_back(-35);
	center.push_back(-35);
	a.resize(num_obj);
	for (int i = 0; i < num_obj; ++i) {
		a[i].resize(2);
		a[i][0] = 0 * cos(2 * M_PI*(i - 1) / (double)num_obj) + r * (-1 * sin(2 * M_PI*(i - 1) / (double)num_obj)) + center[0];
		a[i][1] = 0 * sin(2 * M_PI*(i - 1) / (double)num_obj) + r * cos(2 * M_PI*(i - 1) / (double)num_obj) + center[1];
	}
	v.resize(2);
	v[0].resize(dim);
	v[1].resize(dim);
	for (int i = 0; i < dim; ++i) {
		v[0][i] = (i + 1) % 2;
		v[1][i] = i % 2;
	}

	//fout
	ofstream fout("Points.dat");
	for (int i = 0; i < num_obj; ++i) {
		fout << a[i][0] << "\t" << a[i][1] << endl;
	}

	//Mixed constraints
	if (swit <= 9) {
		mixed_points.resize(num_cons);

		//input lower and upper
		ifstream fin("./2D/2D_random"+to_string(num_cons)+"_ver"+to_string((swit-1)%3+1)+".csv");
		cout << "./2D/2D_random" + to_string(num_cons) + "_ver" + to_string((swit-1) % 3 + 1) + ".csv" << endl;
		if (fin.fail()) {
			cout << "入力ファイルをオープンできません!!!!!" << endl;
		}


		int input_counter = 0;
		for (string line_in; getline(fin, line_in);) {
			if (line_in.size() == 0) continue;
			string token;
			istringstream ss(line_in);
			while (getline(ss, token, ',')) {
				mixed_points[input_counter].push_back(atof(token.c_str()));
			}
			input_counter++;
		}
	}

//	check(swit);
	
}


individual CDMP::Objectives(individual ind) {
	for (int i = 0; i < num_obj; ++i) {
		ind.cost[i] = 0;
		for (int j = 0; j < dim; ++j) {
			ind.cost[i] += pow(ind.var[j] - (a[i][0] * v[0][j] + a[i][1] * v[1][j]), 2);
		}
		ind.cost[i] = sqrt(ind.cost[i]);
	}

	return ind;
}


void CDMP::check(int swit) {
	ofstream fout("Infeasible"+to_string(swit)+".dat");
	ofstream fout1("Feasible"+to_string(swit)+".dat");

	individual x;
	x.var.resize(12);//size of phenotype
	x.cost.resize(num_obj);//number of tasks
	x.constraint.resize(num_cons);//number of tasks

	for (int i = 0; i <= 100; ++i) {
		for (int j = 0; j <= 100; ++j) {
			x.var[0] = 100 * (double)i / 100.0 - 50;
			x.var[1] = 100 * (double)j / 100.0 - 50;
			x = Mixed(x);

			int infeasible = 0;
			for (int i = 0; i < num_cons; ++i) {
				if (x.constraint[i] < 0) {
					fout << x.var[0] << "\t" << x.var[1] << endl;
					infeasible = 1;
					break;
				}
			}
			if (infeasible == 0) {
				fout1 << x.var[0] << "\t" << x.var[1] << endl;
			}
		}
	}
}


individual CDMP::Constraints(int swit, individual ind) {

	switch (swit) {
	case 1:
		ind = Mixed(ind);
		break;

	case 2:
		ind = Mixed(ind);
		break;

	case 3:
		ind = Mixed(ind);
		break;

	case 4:
		ind = Mixed(ind);
		break;

	case 5:
		ind = Mixed(ind);
		break;

	case 6:
		ind = Mixed(ind);
		break;

	case 7:
		ind = Mixed(ind);
		break;

	case 8:
		ind = Mixed(ind);
		break;

	case 9:
		ind = Mixed(ind);
		break;


	default:
		break;
	}

	return ind;
}


individual CDMP::Mixed(individual ind) {
	//
	if (is_ParetoSolution(ind) == 1) {
		for (int i = 0; i < num_cons; ++i) {
			ind.constraint[i] = 0;
		}
	}
	else {
		for (int i = 0; i < num_cons; ++i) {
			ind.constraint[i] = sqrt(pow(ind.var[0] - mixed_points[i][0], 2) + pow(ind.var[1] - mixed_points[i][1], 2)) - mixed_points[i][2];
		}
	}
	
	return ind;
}


individual CDMP::Checker(individual ind) {

	return ind;
}

individual CDMP::Vertices(individual ind) {

	return ind;
}

individual CDMP::Moat(individual ind) {

	return ind;
}

individual CDMP::Center(individual ind) {

	return ind;
}


int CDMP::is_ParetoSolution(individual ind) {
	double x[2] = { ind.var[0] - center[0], ind.var[1] - center[1] };

	bool b1, b2, b3;



	for (int i = 0; i < num_obj; ++i) {
		b1 = sign(ind.var, center, a[i%num_obj]) < 0.0f;
		b2 = sign(ind.var, a[i%num_obj], a[(i+1)%num_obj]) < 0.0f;
		b3 = sign(ind.var, a[(i+1)%num_obj], center) < 0.0f;
		if ((b1 == b2) && (b2 == b3) == 1) return 1;
	}

	return 0;
}


