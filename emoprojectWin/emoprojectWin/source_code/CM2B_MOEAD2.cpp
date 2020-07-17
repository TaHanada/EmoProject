#include"CM2B_MOEAD2.h"
#include"weight_vec.h"
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<stdio.h>
#include <iomanip>

struct rank_no {
	int no;
	double distance;

	bool operator<(const rank_no& another) const {
		if (distance != another.distance) {
			return distance < another.distance;
		}
		return no < another.no;
	}
};

/*
void CM2B_MOEAD2::range_setting(int ProblemID) {
	switch (ProblemID)
	{
	case 1:
	{
		//input lower and upper
		ifstream fin("decisionspace.txt");
		if (fin.fail()) {
			cout << "入力ファイルをオープンできません" << endl;
		}

		range.resize(dim);
		int input_counter = 0;
		for (string line_in; getline(fin, line_in);) {
			if (line_in.size() == 0) continue;

			string token;
			istringstream ss(line_in);
			while (getline(ss, token, '\t')) {
				range[input_counter].push_back(atof(token.c_str()));

			}
			input_counter++;
		}
	}
	break;

	case 2:
		range.resize(dim);
		for (int i = 0; i < dim; ++i) {
			range[i].resize(2);
		}
		range[0][0] = 0.125, range[0][1] = 5.0;
		range[1][0] = 0.125, range[1][1] = 5.0;
		range[2][0] = 0.1, range[2][1] = 10.0;
		range[3][0] = 0.1, range[3][1] = 10.0;
		break;

	default:
		break;
	}


}*/


void CM2B_MOEAD2::Parameter_setting(string algorithm, int ProblemID, int ProblemSubID) {
	//Problem specific parameters
	problems.set_ProblemID(ProblemID,ProblemSubID);
	problems.parameter_setting(num_obj, num_const, dim, range);


	//num_vec
	switch (num_obj) {
	case 2:
		num_vec = 100;
		H1 = num_vec - 1, H2 = 0;//分割数 12(M=3) 6(M=5) 3,2(M=8) 3,2(M=10)
		break;

	case 3:
		num_vec = 300;
		H1 = 23, H2 = 0;//分割数 12(M=3) 6(M=5) 3,2(M=8) 3,2(M=10)
		break;
	default:
		//only competition problem
		if (ProblemID == 6 || ProblemID == 7) {
			num_vec = 210;
			H1 = 6, H2 = 0; //このとき，分割数は6．
		}
		break;
	}


	//others
	pop_m = 10;
	pop_c = 1;
	//only competition 個体群サイズ210だと，世代数は48で端数が80．
	if (ProblemID == 6 || ProblemID == 7) {
		gen = 47;
	}
	else {
		gen = 100;
	}
	crossover_probability = 1.0;
	mutation_probability = 1.0 / (double)dim;
	weight.resize(num_vec, num_obj, H1, H2, 1, num_vec);//コンストラクタの代わり
	weight.initialize();//重みベクトルの生成と距離の計算
	scalar= "tch";
	cout << mutation_probability << endl;
}


//consol output population
void CM2B_MOEAD2::cout_population(int swit) {
	switch (swit) {
	case 1:
		for (int i = 0; i < num_vec; ++i) {
			curent[i].cout_population();
		}
		break;

	case 2:
		for (int i = 0; i < num_vec; ++i) {
			offspring[i].cout_population();
		}
		break;

	default:
		break;
	}
}


//consol output property
void CM2B_MOEAD2::cout_property(int swit) {
	switch (swit) {
	case 1:
		for (int i = 0; i <num_vec; ++i) {
			curent[i].cout_property();
		}
		break;

	case 2:
		for (int i = 0; i < num_vec; ++i) {
			offspring[i].cout_property();
		}
		break;

	default:
		break;
	}
}


//weighted sum
double CM2B_MOEAD2::weighted_sum(individual ind) {
	double value = 0;

	for (int i = 0; i < num_obj; ++i) {
		value += ind.cost[i] * weight.vec[ind.weight_no][i];
	}

	return value;
}


//Tchebycheff no absolute
double CM2B_MOEAD2::tchebycheff_notabsolute(individual ind) {
	double max = 0;
	double g = 0;

	for (int i = 0; i < num_obj; ++i) {
		g = weight.vec[ind.weight_no][i] * ((reference_point[i] - ind.cost[i]));
		if (max < g) {
			max = g;
		}
	}

	return max;//最小化されている
}


//Tchebycheff
double CM2B_MOEAD2::tchebycheff(individual ind) {
	double max = 0;
	double g = 0;

	for (int i = 0; i < num_obj; ++i) {
		if (weight.vec[ind.weight_no][i] == 0) {
			g = pow(10, -6) * abs((0.0 - (ind.cost[i] - reference_point[i]) / (nadir_point[i] - reference_point[i])));//normalization
		}
		else {
			g = weight.vec[ind.weight_no][i] * abs((0.0 - (ind.cost[i] - reference_point[i]) / (nadir_point[i] - reference_point[i])));//normalization
		}

		if (max < g) {
			max = g;
		}
	}
/*	for (int i = 0; i < num_obj; ++i) {
		if (weight.vec[ind.weight_no][i] == 0) {
			g = pow(10, -6) * abs(reference_point[i] - ind.cost[i]);//normalization
		}
		else {
			g = weight.vec[ind.weight_no][i] * abs(reference_point[i] - ind.cost[i]);//normalization
		}

		if (max < g) {
			max = g;
		}
	}*/
	return max;
}


//PBI
double CM2B_MOEAD2::PBI(individual ind) {
	double theta = 5.0;//パラメータ
	double value = 0;
	double dt = 0;
	double dn = 0;
	double ramda = 0;

/*
	//dtの計算
	for (int i = 0; i < num_obj; ++i) {
		//分子
		dt += (reference_point[i] - ind.cost[i])*weight.vec[ind.weight_no][i];

		//分母
		ramda += pow(weight.vec[ind.weight_no][i], 2);
	}
	dt = dt / sqrt(ramda);



	//dnの計算
	for (int i = 0; i < num_obj; ++i) {
		dn += pow((reference_point[i] - ind.cost[i]) - dt*weight.vec[ind.weight_no][i] / sqrt(ramda), 2);
	}
	dn = sqrt(dn);

	value = abs(dt) + theta*dn;
	*/
	return value;
}


//IPBI
double CM2B_MOEAD2::IPBI(individual ind) {
	double theta = 0.1;//パラメータ
	double value = 0;
	double dt = 0;
	double dn = 0;
	double ramda = 0;


	//dtの計算
	for (int i = 0; i < num_obj; ++i) {
		//分子
		dt += (nadir_point[i] - ind.cost[i])*weight.vec[ind.weight_no][i];

		//分母
		ramda += pow(weight.vec[ind.weight_no][i], 2);
	}
	dt = dt / sqrt(ramda);



	//dnの計算
	for (int i = 0; i < num_obj; ++i) {
		dn += pow((nadir_point[i] - ind.cost[i]) - dt*weight.vec[ind.weight_no][i] / sqrt(ramda), 2);
	}
	dn = sqrt(dn);

	value = abs(dt) - theta*dn;

	return value;
}


//AOF
double CM2B_MOEAD2::AOF(individual ind) {
	int k = ind.weight_no;

	double p = -1.0*ind.cost[1] + k/4;
	double value=0;
	
	/*	double theta = 0.1;//パラメータ
	double value = 0;

	double ramda = 0;


	//dtの計算
	for (int i = 0; i < num_obj; ++i) {
		//分子
		dt += (nadir_point[i] - ind.cost[i])*weight.vec[ind.weight_no][i];

		//分母
		ramda += pow(weight.vec[ind.weight_no][i], 2);
	}
	dt = dt / sqrt(ramda);



	//dnの計算
	for (int i = 0; i < num_obj; ++i) {
		dn += pow((nadir_point[i] - ind.cost[i]) - dt*weight.vec[ind.weight_no][i] / sqrt(ramda), 2);
	}
	dn = sqrt(dn);

	value = abs(dt) - theta*dn;

	//penalty
	value -= ind.total_c*penalty_para[island];

	*/
	return value;
}


//scalarizing function
double CM2B_MOEAD2::scalarizing_function(individual ind) {
	if (scalar == "WS") {
		return weighted_sum(ind);
	}
	else if (scalar == "tch") {
		return tchebycheff(ind);
	}
	else if (scalar == "PBI") {
		return  PBI(ind);
	}
	else if (scalar == "IPBI") {
		return  IPBI(ind);
	}
	else if (scalar == "IPBI") {
		return  AOF(ind);
	}
	else {
		cout << "ERROR: " << scalar << " is not defined!" << endl;
		return 0;
	}
}


//calcurate total constraints
void CM2B_MOEAD2::total_constraints(individual& ind) {
	ind.total_c = 0;

	for (int i = 0; i < num_const; ++i) {
		if (ind.constraint[i] < 0) {
			ind.total_c += -1.0*ind.constraint[i];
		}
	}
}


//initialize reference point and nadir point
void CM2B_MOEAD2::initialize_reference_nadir_point() {
	reference_point.resize(num_obj);
	nadir_point.resize(num_obj);
	for (int j = 0; j < num_obj; ++j) {
		reference_point[j] = 1000000000;
		nadir_point[j] = -1000000000;
	}

	for (int i = 0; i < num_obj; ++i) {
		//generate solution
		individual tmp_ind;
		weight_vec weight;
		tmp_ind.var.resize(dim);
		tmp_ind.cost.resize(num_obj);
		tmp_ind.constraint.resize(num_const);

		for (int j = 0; j < dim; ++j) {
			tmp_ind.var[j] = range[j][0] + (range[j][1] - range[j][0])*nextDouble();
		}

		//evaluation
		tmp_ind = problems.evaluation(tmp_ind);


		for (int j = 0; j < num_obj; ++j) {
			if (reference_point[j] > tmp_ind.cost[j]) {
				reference_point[j] = tmp_ind.cost[j];//set reference point
			}
			if (nadir_point[j] < tmp_ind.cost[j]) {
				nadir_point[j] = tmp_ind.cost[j];//set nadir point
			}
		}
	}

}


//initialize population
void CM2B_MOEAD2::initialize_population() {
	vector<individual> update_population;


	for (int i = 0; i < num_vec; ++i) {
		for (int pp = 0; pp < pop_c; ++pp) {
			for (int j = 0; j < dim; ++j) {
				offspring[i].ind[pp].var[j] = range[j][0] + (range[j][1] - range[j][0])*nextDouble();
			}

			//calculate costs
			offspring[i].ind[pp] = problems.evaluation(offspring[i].ind[pp]);
			total_constraints(offspring[i].ind[pp]);//cal total constraints value


			//update reference point and nadir point
			for (int j = 0; j < num_obj; ++j) {
				if (reference_point[j] > offspring[i].ind[pp].cost[j]) {
					reference_point[j] = offspring[i].ind[pp].cost[j];
				}
				if (nadir_point[j] < offspring[i].ind[pp].cost[j]) {
					nadir_point[j] = offspring[i].ind[pp].cost[j];
				}
			}

			update_population.push_back(offspring[i].ind[pp]);
		}
	}

	//file out solutions
	vector<double> rf;
	for (int i = 0; i < num_vec; ++i) {
		rf.push_back(0);
	}
	file_allsolutions(rep, 0, scalar, offspring, rf);

	for (int i = 0; i < num_vec; ++i) {
		for (int j = 0; j < update_population.size(); ++j) {
			//adjustment of weight vector
			update_population[j].weight_no = i;
			update_population[j].fitness = scalarizing_function(update_population[j]);
		}
	
		update_NSGAII(update_population, i);
	}	
}


//curent solution's objective values
void CM2B_MOEAD2::file_objectives(int k, int island) {
	string name = "data_" + to_string(k) + ".dat";
	ofstream fout("./CM2B_MOEAD/" + to_string(ProblemSubID) + "/" + scalar + "/island" + to_string(island + 1) + "/" + name);

	if (island == 0) {
		for (int i = 0; i < num_vec; ++i) {
			for (int j = 0; j < pop_m; ++j) {
				if (curent[i].ind[j].total_c == 0) {
					for (int k = 0; k < num_obj; ++k) {
						//fout << setprecision(7) << max(curent[i].ind[j].cost[k], -1 * curent[i].ind[j].cost[k]) << "\t";
						fout << setprecision(10) << setprecision(7) << curent[i].ind[j].cost[k] << "\t";
					}
					for (int k = 0; k < num_const; ++k) {
						fout << setprecision(10) << curent[i].ind[j].constraint[k] << "\t";
					}
					for (int k = 0; k < dim - 1; ++k) {
						fout << setprecision(10) << curent[i].ind[j].var[k] << "\t";
					}
					fout << setprecision(10) << curent[i].ind[j].var[dim - 1] << endl;
					//fout << setprecision(7) << max(curent[i].ind[j].cost[num_obj - 1], -1* curent[i].ind[j].cost[num_obj - 1]) << endl;
				}
			}
		}
	}
	else if (island > 0) {
		for (int i = 0; i < num_vec; ++i) {
			for (int j = 0; j < pop_m; ++j) {
				if (curent[i].ind[j].total_c != 0) {
					for (int k = 0; k < num_obj-1; ++k) {
						//fout << setprecision(7) << max(curent[i].ind[j].cost[k], -1 * curent[i].ind[j].cost[k]) << "\t";
						fout << setprecision(10) << curent[i].ind[j].cost[k] << "\t";
					}for (int k = 0; k < num_const; ++k) {
						fout << setprecision(10) << curent[i].ind[j].constraint[k] << "\t";
					}
					for (int k = 0; k < dim - 1; ++k) {
						fout << setprecision(10) << curent[i].ind[j].var[k] << "\t";
					}
					fout << setprecision(10) << curent[i].ind[j].var[dim - 1] << endl;
					//fout << setprecision(7) << max(curent[i].ind[j].cost[num_obj - 1], -1 * curent[i].ind[j].cost[num_obj - 1]) << endl;
				}
			}
		}
	}
}


//evolution (mu, mu) 
void CM2B_MOEAD2::evolution_mu(int curent_gen) {
	int n1 = 0, n2 = 0, n3=0, n4=0;//parent ID
	int i1, i2, i3, i4;//parent ID

	int equal_var = 1;
	vector<double> ras;
	ras.resize(num_vec);

	//generate permutation of weight vector
	vector<int> permutation;
	for (int i = 0; i < num_vec; ++i) {
		permutation.push_back(i);
	}
	random_shuffle(permutation.begin(), permutation.end());


	for (int ii = 0; ii < num_vec; ++ii) {
		int i = permutation[ii];
		//check the number of solutions that have better scalarizing function value than minID
		int minI = get_minID(curent[i]);
		double rasv = 0;
		for (int pp = 0; pp < pop_m; ++pp) {
			if (curent[i].ind[pp].fitness < curent[i].ind[minI].fitness) rasv = rasv + 1;
		}
		ras[i] = rasv / (double)pop_m;


		for (int pp = 0; pp < pop_c; ++pp) {
			equal_var = 1;
			while (equal_var == 1) {
				parentselection(i, i1, i2, n1, n2);


				//SBX
				if (nextDouble() <= crossover_probability) {
					offspring[i].ind[pp] = SBX(i1, i2, n1, n2);
				}
				else {
					offspring[i].ind[pp] = curent[i1].ind[n1];
				}


				//polynomial mutation
				PM(i, pp);


				//repair
				for (int j = 0; j < dim; ++j) {
					if (offspring[i].ind[pp].var[j] < range[j][0]) {
						offspring[i].ind[pp].var[j] = range[j][0];
					}
					else if (offspring[i].ind[pp].var[j] > range[j][1]) {
						offspring[i].ind[pp].var[j] = range[j][1];
					}
				}

				//judge equal
				equal_var = 0;
				for (int j = 0; j < num_vec; ++j) {
					for (int r = 0; r < pop_m; ++r) {
						if (curent[j].ind[r].var == offspring[i].ind[pp].var) {
							equal_var = 1;
						}
					}
				}
				for (int j = 0; j <= i; ++j) {
					for (int r = 0; r <= pp; ++r) {
						if (j == i&&r == pp) {

						}
						else if (offspring[j].ind[r].var == offspring[i].ind[pp].var) {
							equal_var = 1;
						}
					}
				}
			}//while (equal_var == 1)
		}//for (int pp = 0; pp < pop_c; ++pp)
	}//for (int ii = 0; ii < num_vec; ++ii)



	//evaluation and update
	for (int i = 0; i < num_vec; ++i) {
		for (int pp = 0; pp < pop_c; ++pp) {
			//evaluation
			offspring[i].ind[pp] = problems.evaluation(offspring[i].ind[pp]);
			total_constraints(offspring[i].ind[pp]);


			//update ref and nadir
			for (int j = 0; j < num_obj; ++j) {
				if (reference_point[j] > offspring[i].ind[pp].cost[j]) {
					reference_point[j] = offspring[i].ind[pp].cost[j];
				}
				if (nadir_point[j] < offspring[i].ind[pp].cost[j]) {
					nadir_point[j] = offspring[i].ind[pp].cost[j];
				}
			}
		}
	}


	//file output
	file_allsolutions(rep, curent_gen, scalar, offspring, ras);

	//update each sub-population
	vector<individual> update_population;
	update_population.resize(pop_c*num_vec + pop_m);
	for (int i = 0; i < num_vec; ++i) {
		for (int pp = 0; pp < pop_c; ++pp) {
			update_population[pop_m + i] = offspring[i].ind[pp];
		}
	}

	for (int i = 0; i < num_vec; ++i) {
		for (int j = 0; j < pop_m; ++j) {
			update_population[j] = curent[i].ind[j];
		}

		for (int j = 0; j < pop_c*num_vec + pop_m; ++j) {
			update_population[j].weight_no = i;
			update_population[j].fitness=scalarizing_function(update_population[j]);
		}

		update_NSGAII(update_population, i);
	}
}


//run
void CM2B_MOEAD2::run_algorithm(int reps) {
	rep = reps;
	srand(reps);
	init_genrand(reps);
	//generate files
	cout << rep << endl;
	file_initialization(rep, gen);
	//file_initialization(rep + 1, gen);

	//weight vector
	weight.reselect_nighbor();
	weight.shuffle();


	//initialize ref. and nadir points
	initialize_reference_nadir_point();


	//inititalize population
	initialize_population();

	update_constraint_max(offspring);

	//CDMP
	for (int island = 0; island < 2; ++island) {
		file_curents(reps, 0, scalar, island);
	}

//	cout << "initialize population" << endl;

	//evolution
	for (int j = 0; j < gen; ++j) {
		if (j%10==0) cout << "generation: " << j << endl;

		//curent　population
		evolution_mu(j + 1);

		//CDMP
		for (int island = 0; island < 2; ++island) {
			file_curents(reps, j+1, scalar, island);
		}
	}


	//file output of final population
	for (int island = 0; island < 2; ++island) {
		file_objectives(reps, island);
		//file_objectives(reps + 1, island);
	}
}


//initialization files
void CM2B_MOEAD2::file_initialization(int  k, int gen) {
	int count = 0;
	int tmp_gen = 0;

	tmp_gen = gen;
	//offspirngs
	{
		string			name1 = to_string(k) + "_" + to_string(0) + "gen.dat";
		ofstream fout("./CM2B_MOEAD/" + to_string(ProblemSubID) + "/" + scalar + "/offsprings/" + name1);
		fout << "#gen	Feasible	";
		for (int i = 0; i < num_obj; ++i) {
			fout << "f" << to_string(i+1) <<"	";
		}
		for (int i = 0; i < num_const; ++i) {
			fout << "c" << to_string(i + 1) << "	";
		}
		for (int i = 0; i < dim - 1; ++i) {
			fout << "x" << to_string(i + 1) << "	";
		}
		fout << "x" << to_string(dim) << "	";
		fout << "ras" << endl;


		for (int i = 1; i <= tmp_gen; ++i) {
			string name1 = to_string(k) + "_" + to_string(i) + "gen.dat";
			ofstream fout("./CM2B_MOEAD/" + to_string(ProblemSubID) + "/" + scalar + "/offsprings/" + name1);
			fout << "#gen	Feasible	";
			for (int i = 0; i < num_obj; ++i) {
				fout << "f" << to_string(i + 1) << "	";
			}
			for (int i = 0; i < num_const; ++i) {
				fout << "c" << to_string(i + 1) << "	";
			}
			for (int i = 0; i < dim - 1; ++i) {
				fout << "x" << to_string(i + 1) << "	";
			}
			fout << "x" << to_string(dim) << "	";
			fout << "ras" << endl;

		}
	}


	//現個体
	for (int island = 0; island < 2; ++island) {
		string name1 = to_string(k) + "_" + to_string(0) + "gen.dat";
		ofstream fout("./CM2B_MOEAD/" + to_string(ProblemSubID) + "/" + scalar + "/island" + to_string(island + 1) + "/curents/" + name1);
		fout << "#gen	Feasible	";
		for (int i = 0; i < num_obj; ++i) {
			fout << "f" << to_string(i + 1) << "	";
		}
		for (int i = 0; i < num_const; ++i) {
			fout << "c" << to_string(i + 1) << "	";
		}
		for (int i = 0; i < dim - 1; ++i) {
			fout << "x" << to_string(i + 1) << "	";
		}
		fout << "x" << to_string(dim) << endl;
		//fout << "f1	f2	c1	c2	c3	c4	c5	c6	c7	c8	c9	c10	c11	c12	c13	c14	c15	c16	c17	c18	c19	c20	c21	c22	c23	c24	c25	c26	c27	c28	c29	c30	c31	c32	c33	c34	c35	c36	c37	c38	c39	c40	c41	c42	c43	c44	c45	c46	c47	c48	c49	c50	c51	c52	c53	c54	x1	x2	x3	x4	x5	x6	x7	x8	x9	x10	x11	x12	x13	x14	x15	x16	x17	x18	x19	x20	x21	x22	x23	x24	x25	x26	x27	x28	x29	x30	x31	x32	x33	x34	x35	x36	x37	x38	x39	x40	x41	x42	x43	x44	x45	x46	x47	x48	x49	x50	x51	x52	x53	x54	x55	x56	x57	x58	x59	x60	x61	x62	x63	x64	x65	x66	x67	x68	x69	x70	x71	x72	x73	x74	x75	x76	x77	x78	x79	x80	x81	x82	x83	x84	x85	x86	x87	x88	x89	x90	x91	x92	x93	x94	x95	x96	x97	x98	x99	x100	x101	x102	x103	x104	x105	x106	x107	x108	x109	x110	x111	x112	x113	x114	x115	x116	x117	x118	x119	x120	x121	x122	x123	x124	x125	x126	x127	x128	x129	x130	x131	x132	x133	x134	x135	x136	x137	x138	x139	x140	x141	x142	x143	x144	x145	x146	x147	x148	x149	x150	x151	x152	x153	x154	x155	x156	x157	x158	x159	x160	x161	x162	x163	x164	x165	x166	x167	x168	x169	x170	x171	x172	x173	x174	x175	x176	x177	x178	x179	x180	x181	x182	x183	x184	x185	x186	x187	x188	x189	x190	x191	x192	x193	x194	x195	x196	x197	x198	x199	x200	x201	x202	x203	x204	x205	x206	x207	x208	x209	x210	x211	x212	x213	x214	x215	x216	x217	x218	x219	x220	x221	x222" << endl;
		for (int i = 1; i <= tmp_gen; ++i) {
			string name1 = to_string(k) + "_" + to_string(i) + "gen.dat";
			ofstream fout("./CM2B_MOEAD/" + to_string(ProblemSubID) + "/" + scalar + "/island" + to_string(island + 1) + "/curents/" + name1);
			fout << "#gen	Feasible	";
			for (int i = 0; i < num_obj; ++i) {
				fout << "f" << to_string(i + 1) << "	";
			}
			for (int i = 0; i < num_const; ++i) {
				fout << "c" << to_string(i + 1) << "	";
			}
			for (int i = 0; i < dim - 1; ++i) {
				fout << "x" << to_string(i + 1) << "	";
			}
			fout << "x" << to_string(dim) << endl;
			//fout << "#gen	Feasible	f1	f2	c1	c2	c3	c4	c5	c6	c7	c8	c9	c10	c11	c12	c13	c14	c15	c16	c17	c18	c19	c20	c21	c22	c23	c24	c25	c26	c27	c28	c29	c30	c31	c32	c33	c34	c35	c36	c37	c38	c39	c40	c41	c42	c43	c44	c45	c46	c47	c48	c49	c50	c51	c52	c53	c54	x1	x2	x3	x4	x5	x6	x7	x8	x9	x10	x11	x12	x13	x14	x15	x16	x17	x18	x19	x20	x21	x22	x23	x24	x25	x26	x27	x28	x29	x30	x31	x32	x33	x34	x35	x36	x37	x38	x39	x40	x41	x42	x43	x44	x45	x46	x47	x48	x49	x50	x51	x52	x53	x54	x55	x56	x57	x58	x59	x60	x61	x62	x63	x64	x65	x66	x67	x68	x69	x70	x71	x72	x73	x74	x75	x76	x77	x78	x79	x80	x81	x82	x83	x84	x85	x86	x87	x88	x89	x90	x91	x92	x93	x94	x95	x96	x97	x98	x99	x100	x101	x102	x103	x104	x105	x106	x107	x108	x109	x110	x111	x112	x113	x114	x115	x116	x117	x118	x119	x120	x121	x122	x123	x124	x125	x126	x127	x128	x129	x130	x131	x132	x133	x134	x135	x136	x137	x138	x139	x140	x141	x142	x143	x144	x145	x146	x147	x148	x149	x150	x151	x152	x153	x154	x155	x156	x157	x158	x159	x160	x161	x162	x163	x164	x165	x166	x167	x168	x169	x170	x171	x172	x173	x174	x175	x176	x177	x178	x179	x180	x181	x182	x183	x184	x185	x186	x187	x188	x189	x190	x191	x192	x193	x194	x195	x196	x197	x198	x199	x200	x201	x202	x203	x204	x205	x206	x207	x208	x209	x210	x211	x212	x213	x214	x215	x216	x217	x218	x219	x220	x221	x222" << endl;
		}
	}
/*	for (int island = 0; island < 2; ++island) {

		string name1 = to_string(k) + "_" + to_string(0) + "gen.dat";
		ofstream fout("./" + scalar + "/island" + to_string(island + 1) + "/curents/" + name1);
		fout << "#gen	Feasible	f1	f2	c1	c2	c3	c4	c5	c6	c7	c8	c9	c10	c11	c12	c13	c14	c15	c16	c17	c18	c19	c20	c21	c22	c23	c24	c25	c26	c27	c28	c29	c30	c31	c32	c33	c34	c35	c36	c37	c38	c39	c40	c41	c42	c43	c44	c45	c46	c47	c48	c49	c50	c51	c52	c53	c54	x1	x2	x3	x4	x5	x6	x7	x8	x9	x10	x11	x12	x13	x14	x15	x16	x17	x18	x19	x20	x21	x22	x23	x24	x25	x26	x27	x28	x29	x30	x31	x32	x33	x34	x35	x36	x37	x38	x39	x40	x41	x42	x43	x44	x45	x46	x47	x48	x49	x50	x51	x52	x53	x54	x55	x56	x57	x58	x59	x60	x61	x62	x63	x64	x65	x66	x67	x68	x69	x70	x71	x72	x73	x74	x75	x76	x77	x78	x79	x80	x81	x82	x83	x84	x85	x86	x87	x88	x89	x90	x91	x92	x93	x94	x95	x96	x97	x98	x99	x100	x101	x102	x103	x104	x105	x106	x107	x108	x109	x110	x111	x112	x113	x114	x115	x116	x117	x118	x119	x120	x121	x122	x123	x124	x125	x126	x127	x128	x129	x130	x131	x132	x133	x134	x135	x136	x137	x138	x139	x140	x141	x142	x143	x144	x145	x146	x147	x148	x149	x150	x151	x152	x153	x154	x155	x156	x157	x158	x159	x160	x161	x162	x163	x164	x165	x166	x167	x168	x169	x170	x171	x172	x173	x174	x175	x176	x177	x178	x179	x180	x181	x182	x183	x184	x185	x186	x187	x188	x189	x190	x191	x192	x193	x194	x195	x196	x197	x198	x199	x200	x201	x202	x203	x204	x205	x206	x207	x208	x209	x210	x211	x212	x213	x214	x215	x216	x217	x218	x219	x220	x221	x222" << endl;


		for (int i = 1; i <= tmp_gen; ++i) {
			string name1 = to_string(k) + "_" + to_string(i) + "gen.dat";
			ofstream fout("./" + scalar + "/island" + to_string(island + 1) + "/curents/" + name1);
			fout << "#gen	Feasible	f1	f2	c1	c2	c3	c4	c5	c6	c7	c8	c9	c10	c11	c12	c13	c14	c15	c16	c17	c18	c19	c20	c21	c22	c23	c24	c25	c26	c27	c28	c29	c30	c31	c32	c33	c34	c35	c36	c37	c38	c39	c40	c41	c42	c43	c44	c45	c46	c47	c48	c49	c50	c51	c52	c53	c54	x1	x2	x3	x4	x5	x6	x7	x8	x9	x10	x11	x12	x13	x14	x15	x16	x17	x18	x19	x20	x21	x22	x23	x24	x25	x26	x27	x28	x29	x30	x31	x32	x33	x34	x35	x36	x37	x38	x39	x40	x41	x42	x43	x44	x45	x46	x47	x48	x49	x50	x51	x52	x53	x54	x55	x56	x57	x58	x59	x60	x61	x62	x63	x64	x65	x66	x67	x68	x69	x70	x71	x72	x73	x74	x75	x76	x77	x78	x79	x80	x81	x82	x83	x84	x85	x86	x87	x88	x89	x90	x91	x92	x93	x94	x95	x96	x97	x98	x99	x100	x101	x102	x103	x104	x105	x106	x107	x108	x109	x110	x111	x112	x113	x114	x115	x116	x117	x118	x119	x120	x121	x122	x123	x124	x125	x126	x127	x128	x129	x130	x131	x132	x133	x134	x135	x136	x137	x138	x139	x140	x141	x142	x143	x144	x145	x146	x147	x148	x149	x150	x151	x152	x153	x154	x155	x156	x157	x158	x159	x160	x161	x162	x163	x164	x165	x166	x167	x168	x169	x170	x171	x172	x173	x174	x175	x176	x177	x178	x179	x180	x181	x182	x183	x184	x185	x186	x187	x188	x189	x190	x191	x192	x193	x194	x195	x196	x197	x198	x199	x200	x201	x202	x203	x204	x205	x206	x207	x208	x209	x210	x211	x212	x213	x214	x215	x216	x217	x218	x219	x220	x221	x222" << endl;

		}
	}*/

}


//file output offsprings
void CM2B_MOEAD2::file_allsolutions(int k, int gen, string scalar, vector<population> pop, vector<double> ras) {
	{
		string name1 = to_string(k) + "_" + to_string(gen) + "gen.dat";
		ofstream fout("./CM2B_MOEAD/"+to_string(ProblemSubID)+"/" + scalar + "/offsprings/" + name1, ios::out | ios::app);
		for (int k = 0; k < num_vec; ++k) {
			for (int j = 0; j < pop_c; j++) {
				//gen 
				fout << gen << "\t";

				//feasible
				if (pop[k].ind[j].total_c == 0) {
					fout << "1" << "\t";
				}
				else {
					fout << "0" << "\t";
				}

				//function value
				for (int i = 0; i < num_obj; ++i) {
					fout << setprecision(10) << pop[k].ind[j].cost[i] << "\t";
				}
				for (int i = 0; i < num_const; i++) {
					fout << setprecision(10) << pop[k].ind[j].constraint[i] << "\t";
				}
				for (int i = 0; i < dim - 1; ++i) {
					fout << setprecision(10) << pop[k].ind[j].var[i] << "\t";
				}
				fout << setprecision(10) << pop[k].ind[j].var[dim - 1] << "\t";
				//constraints
		/*		for (int i = 0; i < 54; ++i) {
					fout << ind.constraint[i] << "\t";
				}*/

				//var
		/*		for (int i = 0; i <221; i++)
				{
					fout << ind.var[i] << "\t";
				}
				fout << ind.var[221] << "\t";*/
				fout << ras[k] << endl;
			}
		}
	}
}


//目的関数値のプロット(現個体)
void CM2B_MOEAD2::file_curents(int k, int gen, string scalar, int island) {

	string name1 = to_string(k) + "_" + to_string(gen) + "gen.dat";
	ofstream fout("./CM2B_MOEAD/" + to_string(ProblemSubID) + "/" + scalar + "/island" + to_string(island + 1) + "/curents/" + name1, ios::out | ios::app);

	if (island == 0) {
		for (int j = 0; j < num_vec; ++j) {
			double min_totalC = 1000000000;
			double min_fitness = 1000000000;
			int minID = 0;
			for (int i = 0; i < pop_m; ++i) {
				if (min_totalC > curent[j].ind[i].total_c) {
					min_totalC = curent[j].ind[i].total_c;
					minID = i;
				}
				else if (min_totalC == curent[j].ind[i].total_c&&min_fitness > curent[j].ind[i].fitness) {
					min_fitness = curent[j].ind[i].fitness;
					minID = i;
				}
			}

			//gen 
			fout << gen << "\t";

			//feasible
			if (curent[j].ind[minID].total_c == 0) {
				fout << "1" << "\t";
			}
			else {
				fout << "0" << "\t";
			}


			
			//function value
			//fout << curent[j].ind[minID].cost[0] << "\t" << max(curent[j].ind[minID].cost[1] , -1*curent[j].ind[minID].cost[1]) << "\t";
			for (int i = 0; i < num_obj; ++i) {
				fout << setprecision(10) << curent[j].ind[minID].cost[i] << "\t";
			}
			//constraints
			for (int i = 0; i < num_const; ++i) {
				fout << setprecision(10) << curent[j].ind[minID].constraint[i] << "\t";
			}

			//var
			for (int i = 0; i < dim-1; i++)
			{
				fout << setprecision(10) << curent[j].ind[minID].var[i] << "\t";
			}
			fout << setprecision(10) << curent[j].ind[minID].var[dim-1] << endl;
		}
	}
	else if (island > 0) {
		for (int k = 0; k < num_vec; ++k) {
			for (int j = 1; j < pop_m; ++j) {
				//gen 
				fout << gen << "\t";

				//feasible
				if (curent[k].ind[j].total_c == 0) {
					fout << "1" << "\t";
				}
				else {
					fout << "0" << "\t";
				}

				//function value
				//fout << curent[k].ind[j].cost[0] << "\t" << -1.0*curent[k].ind[j].cost[1] << "\t";
				for (int i = 0; i < num_obj; ++i) {
					fout << setprecision(10) << curent[k].ind[j].cost[i] << "\t";
				}

				//constraints
				for (int i = 0; i < num_const; ++i) {
					fout << setprecision(10) << curent[k].ind[j].constraint[i] << "\t";
				}

				//var
				for (int i = 0; i < dim-1; i++)
				{
					fout << setprecision(10) << curent[k].ind[j].var[i] << "\t";
				}
				fout << setprecision(10) << curent[k].ind[j].var[dim-1] << endl;
			}
		}
	}
}


int CM2B_MOEAD2::calc_cossim(individual offspring1) {
	double a, b, ab;
	double min_cos = 0;
	int min_vec = 0;


	for (int i = 0; i < pop_m; ++i) {
		a = 0;
		b = 0;
		ab = 0;
		for (int j = 0; j < num_obj; ++j) {
			a += pow(offspring1.cost[j], 2);
			b += pow(weight.vec[i][j], 2);
			ab += offspring1.cost[j] * weight.vec[i][j];
		}

		if (ab / a / b > min_cos) {
			min_cos = ab / a / b;
			min_vec = i;
		}

	}

	return min_vec;
};


int CM2B_MOEAD2::tournament_PRandCV(int tournamentSize, population pop) {
	vector<int> candidates;
	int candidate=0;
	double min_totalC = 1000000;
	int min_rank = 10000000;
	int min_solID = 0;

	for (int i = 0; i < tournamentSize; ++i) {
		candidate = genrand_int32() % pop.size;
		if (min_rank > pop.ind[candidate].rank) {
			min_totalC = pop.ind[candidate].total_c;
			min_rank = pop.ind[candidate].rank;
			min_solID = candidate;
		}
		else if (min_rank == pop.ind[candidate].rank&&min_totalC > pop.ind[candidate].total_c) {
			min_totalC = pop.ind[candidate].total_c;
			min_rank = pop.ind[candidate].rank;
			min_solID = candidate;
		}
	}

	return min_solID;
}


int CM2B_MOEAD2::tournament_totalC(int tournamentSize, population pop) {
	vector<int> candidates;
	int candidate = 0;
	double min_totalC = 1000000;
	int min_rank = 10000000;
	int min_solID = 0;

	for (int i = 0; i < tournamentSize; ++i) {
		candidate = genrand_int32() % pop.size;
		if (min_totalC > pop.ind[candidate].total_c) {
			min_totalC = pop.ind[candidate].total_c;
			min_solID = candidate;
		}
	}

	return min_solID;
}


int CM2B_MOEAD2::tournament(int tournamentSize, population pop) {
	vector<int> candidates;
	int candidate = 0;
	double max_dist = 1000000;
	int min_rank = 10000000;
	int min_solID = 0;

	for (int i = 0; i < tournamentSize; ++i) {
		candidate = genrand_int32() % pop.size;
		if (min_rank > pop.ind[candidate].rank) {
			max_dist = pop.ind[candidate].dist;
			min_rank = pop.ind[candidate].rank;
			min_solID = candidate;
		}
		else if (min_rank == pop.ind[candidate].rank&&max_dist < pop.ind[candidate].dist) {
			max_dist = pop.ind[candidate].dist;
			min_rank = pop.ind[candidate].rank;
			min_solID = candidate;
		}
	}

	return min_solID;
}


int CM2B_MOEAD2::tournament_fitness(int tournamentSize, population pop) {
	vector<int> candidates;
	int candidate = 0;
	double min_fitness = 1000000;
	int min_solID = 0;

	for (int i = 0; i < tournamentSize; ++i) {
		candidate = genrand_int32() % pop.size;
		if (min_fitness > pop.ind[candidate].fitness) {
			min_fitness = pop.ind[candidate].fitness;
			min_solID = candidate;
		}
	}

	return min_solID;
}


int CM2B_MOEAD2::tournament_inverse(int tournamentSize, population pop) {
	vector<int> candidates;
	int candidate = 0;
	double min_totalC = -10000000;
	int min_solID = 0;

	for (int i = 0; i < tournamentSize; ++i) {
		candidate = genrand_int32() % pop.size;
		if (pop.ind[candidate].total_c > min_totalC) {
			min_totalC = pop.ind[candidate].total_c;
			min_solID = candidate;
		}
	}

	return min_solID;
}


individual CM2B_MOEAD2::SBX(int vec1, int vec2, int n1, int n2) {
	individual offspring1, offspring2;
	double mu = 20.0;//distribution index
	double SBX_u = 0;
	double beta1, beta2;
	double alpha;
	double beta;

	offspring1 = curent[vec1].ind[n1];
	offspring2 = curent[vec2].ind[n2];

	for (int j = 0; j < dim; ++j) {
//		n1 = tournament_totalC(2, curent[vec1]);
//		n2 = tournament_totalC(2, curent[vec2]);
		offspring1.var[j] = curent[vec1].ind[n1].var[j];
		offspring2.var[j] = curent[vec2].ind[n2].var[j];

		if (nextDouble() <= 0.5) {
			if (fabs(curent[vec1].ind[n1].var[j] - curent[vec2].ind[n2].var[j]) > pow(10, -14)) {

				if (curent[vec1].ind[n1].var[j] < curent[vec2].ind[n2].var[j]) {
					beta1 = 1.0 + 2.0*(curent[vec1].ind[n1].var[j] - range[j][0]) / (curent[vec2].ind[n2].var[j] - curent[vec1].ind[n1].var[j]);
					beta2 = 1.0 + 2.0*(range[j][1] - curent[vec2].ind[n2].var[j]) / (curent[vec2].ind[n2].var[j] - curent[vec1].ind[n1].var[j]);

					if (nextDouble() <= 0.5) {
						//offspring1
						SBX_u = nextDouble();
						alpha = 2.0 - pow(beta1, -1.0*(mu + 1.0));
						if (SBX_u <= 1.0 / alpha) {
							beta = pow(SBX_u*alpha, 1.0 / (mu + 1.0));
						}
						else {
							beta = pow(2.0 - SBX_u*alpha, -1.0 / (mu + 1.0));
						}
						offspring1.var[j] = 0.5*((1.0 + beta)*curent[vec1].ind[n1].var[j] + (1.0 - beta)*curent[vec2].ind[n2].var[j]);
					}
					else {
						//offspring2
						SBX_u = nextDouble();
						alpha = 2.0 - pow(beta2, -1.0*(mu + 1.0));
						if (SBX_u <= 1.0 / alpha) {
							beta = pow(SBX_u*alpha, 1.0 / (mu + 1.0));
						}
						else {
							beta = pow(2.0 - SBX_u*alpha, -1.0 / (mu + 1.0));
						}
						offspring2.var[j] = 0.5*((1.0 + beta)*curent[vec2].ind[n2].var[j] + (1.0 - beta)*curent[vec1].ind[n1].var[j]);

						offspring1.var[j] = offspring2.var[j];
					}


				}
				else if (curent[vec1].ind[n1].var[j] > curent[vec2].ind[n2].var[j]) {
					beta1 = 1.0 + 2.0*(curent[vec2].ind[n2].var[j] - range[j][0]) / (curent[vec1].ind[n1].var[j] - curent[vec2].ind[n2].var[j]);
					beta2 = 1.0 + 2.0*(range[j][1] - curent[vec1].ind[n1].var[j]) / (curent[vec1].ind[n1].var[j] - curent[vec2].ind[n2].var[j]);

					if (nextDouble() <= 0.5) {
						//offspring1
						SBX_u = nextDouble();
						alpha = 2.0 - pow(beta1, -1.0*(mu + 1.0));
						if (SBX_u <= 1.0 / alpha) {
							beta = pow(SBX_u*alpha, 1.0 / (mu + 1.0));
						}
						else {
							beta = pow(2.0 - SBX_u*alpha, -1.0 / (mu + 1.0));
						}


						offspring1.var[j] = 0.5*((1.0 + beta)*curent[vec2].ind[n2].var[j] + (1.0 - beta)*curent[vec1].ind[n1].var[j]);

					}
					else {
						//offspring2
						SBX_u = nextDouble();
						alpha = 2.0 - pow(beta2, -1.0*(mu + 1.0));
						if (SBX_u <= 1.0 / alpha) {
							beta = pow(SBX_u*alpha, 1.0 / (mu + 1.0));
						}
						else {
							beta = pow(2.0 - SBX_u*alpha, -1.0 / (mu + 1.0));
						}
						offspring2.var[j] = 0.5*((1.0 + beta)*curent[vec1].ind[n1].var[j] + (1.0 - beta)*curent[vec2].ind[n2].var[j]);
						offspring1.var[j] = offspring2.var[j];
					}

				}
			}//if (fabs(curent.ind[n1].var[j] - curent.ind[n2].var[j]) > pow(10, -14))
		}//if (nextDouble() <= 0.5)
	}

	return offspring1;
}


individual CM2B_MOEAD2::DE(int vec1, int vec2, int vec3, int vec4, int n1, int n2, int n3, int n4) {
	individual offspring1, offspring2;
	double mu = 20.0;//distribution index
	double SBX_u = 0;
	double beta1, beta2;
	double alpha;
	double beta;
	double F[3] = { 0.6,0.8,1.0 };

	offspring1 = curent[vec1].ind[n1];
	offspring2 = curent[vec2].ind[n2];
	


	for (int i = 0; i < dim; ++i) {
		offspring1.var[i] = curent[vec1].ind[n1].var[i] + F[genrand_int32() % 3] * (curent[vec2].ind[n2].var[i] - curent[vec1].ind[n1].var[i]) + F[genrand_int32() % 3] * (curent[vec3].ind[n3].var[i] - curent[vec4].ind[n4].var[i]);
	}

	return offspring1;
}


void CM2B_MOEAD2::PM(int vec, int pp) {
	double mum = 20.0;//index of polynomial mutation
	double del = 0;
	double p = 0;

	p = nextDouble();
	if (p < 0.5) {
		del = pow(2 * p, 1.0 / (double)(1 + mum)) - 1.0;
	}
	else {
		del = 1 - pow(2 - 2 * p, 1.0 / (double)(1 + mum));
	}
	for (int j = 0; j < dim; ++j) {
		if (nextDouble() < mutation_probability) {
			offspring[vec].ind[pp].var[j] += del*(range[j][1] - range[j][0]);
		}
	}
}


void CM2B_MOEAD2::update_NSGAII_const(vector<individual> update_population, int vec) {
	vector<int> frontID;
	vector<vector<int>> flag2;
	int flag_dominated = 0;
	frontID.resize(update_population.size());
	int num_pushback=0;

	//solutions whose constraints violation is higher than epsilon
	vector<int> vio_pop;
/*	for (int i = 0; i < update_population.size(); ++i) {

		if (update_population[i].total_c > epsilon_) {
			vio_pop.push_back(i);
			num_pushback++;
			frontID[i] = 1;
		}
	}*/

	while  (num_pushback<update_population.size()){
		vector<int> flag1;
		for (int i = 0; i < update_population.size(); ++i) {
			if (frontID[i] == 0) {
				flag_dominated = 0;
				for (int j = 0; j < update_population.size(); ++j) {
					if (frontID[j] == 0) {
						if (update_population[i].fitness == update_population[j].fitness&&update_population[i].total_c == update_population[j].total_c) {
							frontID[j] = 1;
							num_pushback++;
						}
						else if (update_population[i].fitness >= update_population[j].fitness&&update_population[i].total_c >= update_population[j].total_c) {
							flag_dominated = 1;
							break;
						}
					}
				}

				if (flag_dominated == 0) {
					flag1.push_back(i);
					num_pushback++;
				}
			}
		}

		for (int i = 0; i < flag1.size(); ++i) {
			frontID[flag1[i]] = 1;
		}
		flag2.push_back(flag1);
	}
	flag2.push_back(vio_pop);//solutions whose constraints violation is higher than epsilon

	int pp = 0;
	int sum = 0;
	for (int i = 0; i < flag2.size(); ++i) {
		sum += flag2[i].size();
		if (sum <= pop_m) {
			for (int j = 0; j < flag2[i].size(); ++j) {
				curent[vec].ind[pp] = update_population[flag2[i][j]];
				curent[vec].ind[pp].rank = i;
				pp++;
			}
		}
		else {
			double min_c = 1000000;
			int minID = 0;
			while (pp < pop_m) {
				min_c = 1000000;
				minID = 0;
				for (int j = 0; j < flag2[i].size(); ++j) {
					if (min_c > update_population[flag2[i][j]].total_c) {
						min_c = update_population[flag2[i][j]].total_c;
						minID = j;
					}
				}
				curent[vec].ind[pp] = update_population[flag2[i][minID]];
				curent[vec].ind[pp].rank = i;
				update_population[flag2[i][minID]].total_c = 10000000;
				pp++;
			}
			break;
		}
	}

	/*
	for (int i = 0; i < flag2.size(); ++i) {
		ofstream fout("test"+to_string(i)+".csv");
		for (int j = 0; j < flag2[i].size(); ++j) {
			fout << update_population[flag2[i][j]].fitness << "\t" << update_population[flag2[i][j]].total_c << endl;
		}
	}
	cout << "fin" << endl;
	cin >> flag_dominated;
	*/
}


void CM2B_MOEAD2::update_NSGAII(vector<individual> update_population, int vec) {
	vector<int> frontID;
	vector<vector<int>> flag2;
	int flag_dominated = 0;
	frontID.resize(update_population.size());
	int num_pushback = 0;

	//solutions whose constraints violation is higher than epsilon
/*	vector<int> vio_pop;
	for (int i = 0; i < update_population.size(); ++i) {
		if (update_population[i].total_c > epsilon_) {
			vio_pop.push_back(i);
			num_pushback++;
			frontID[i] = 1;
		}
	}
	int coun = 0;
	*/

	while (num_pushback<pop_m) {
		vector<int> flag1;
		for (int i = 0; i < update_population.size(); ++i) {
			if (frontID[i] == 0) {
				flag_dominated = 0;
				for (int j = 0; j < update_population.size(); ++j) {
					if (frontID[j] == 0) {
						if (update_population[i].fitness == update_population[j].fitness&&update_population[i].total_c == update_population[j].total_c) {
/*							frontID[j] = 1;
							num_pushback++;
							coun++;*/
						}
						else if (update_population[i].fitness >= update_population[j].fitness&&update_population[i].total_c >= update_population[j].total_c) {
							flag_dominated = 1;
							break;
						}
					}
				}

				if (flag_dominated == 0) {
					flag1.push_back(i);
					num_pushback++;
				}
			}
		}

		for (int i = 0; i < flag1.size(); ++i) {
			frontID[flag1[i]] = 1;
		}
		flag2.push_back(flag1);
	}
//	flag2.push_back(vio_pop);//solutions whose constraints violation is higher than epsilon


	//crowding distance
	double max = 0, min = 1000000000;
	//initialize
	for (int i = 0; i < update_population.size(); ++i) {
		update_population[i].dist = 0;
	}

	//calc crowding distance for fitness 
	for (int i = 0; i < flag2.size(); ++i) {
		if (flag2[i].size() != 0) {
			//		cout << i <<" "<< flag2[i].size() <<" "<<flag2.size()<<" "<<update_population.size()<<" "<<coun<< endl;
			max = update_population[flag2[i][0]].fitness, min = update_population[flag2[i][0]].fitness;
			vector<rank_no> cd;
			cd.resize(flag2[i].size());
			for (int j = 0; j < flag2[i].size(); ++j) {
				if (max < update_population[flag2[i][j]].fitness) max = update_population[flag2[i][j]].fitness;
				if (min > update_population[flag2[i][j]].fitness) min = update_population[flag2[i][j]].fitness;
				cd[j].distance = update_population[flag2[i][j]].fitness;
				cd[j].no = j;
			}
			sort(cd.begin(), cd.end());
			update_population[flag2[i][cd[0].no]].dist += 1;
			update_population[flag2[i][cd[flag2[i].size() - 1].no]].dist += 1;
			for (int j = 1; j < flag2[i].size() - 1; ++j) {
				update_population[flag2[i][cd[j].no]].dist += abs(cd[j - 1].distance - cd[j + 1].distance) / (max - min);
			}
		}
	}

	//calc crowding distance for total_c 
	for (int i = 0; i < flag2.size(); ++i) {
		if (flag2[i].size() != 0) {
			max = update_population[flag2[i][0]].total_c, min = update_population[flag2[i][0]].total_c;
			vector<rank_no> cd;
			cd.resize(flag2[i].size());
			for (int j = 0; j < flag2[i].size(); ++j) {
				if (max < update_population[flag2[i][j]].total_c) max = update_population[flag2[i][j]].total_c;
				if (min > update_population[flag2[i][j]].total_c) min = update_population[flag2[i][j]].total_c;
				cd[j].distance = update_population[flag2[i][j]].total_c;
				cd[j].no = j;
			}
			sort(cd.begin(), cd.end());
			update_population[flag2[i][cd[0].no]].dist += 2;
			update_population[flag2[i][cd[flag2[i].size() - 1].no]].dist += 1;
			for (int j = 1; j < flag2[i].size() - 1; ++j) {
				update_population[flag2[i][cd[j].no]].dist += abs(cd[j - 1].distance - cd[j + 1].distance) / (max - min);
			}
		}
	}




	int pp = 0;
	int sum = 0;
	for (int i = 0; i < flag2.size(); ++i) {
//		ofstream fout("test" + to_string(i) + ".csv");
		sum += flag2[i].size();
		if (sum <= pop_m) {
			for (int j = 0; j < flag2[i].size(); ++j) {
				curent[vec].ind[pp] = update_population[flag2[i][j]];
				curent[vec].ind[pp].rank = i;
//				fout << curent[vec].ind[pp].fitness << "\t" << curent[vec].ind[pp].total_c << endl;
				pp++;
			}
		}
		else {
			double min_c = 0;
			int minID = 0;

			while (pp < pop_m) {
				min_c = 0;
				minID = 0;
				for (int j = 0; j < flag2[i].size(); ++j) {
					if (min_c < update_population[flag2[i][j]].dist) {
						min_c = update_population[flag2[i][j]].dist;
						minID = j;
					}
				}
				curent[vec].ind[pp] = update_population[flag2[i][minID]];
				curent[vec].ind[pp].rank = i;
				update_population[flag2[i][minID]].dist = -1;
//				fout << curent[vec].ind[pp].fitness << "\t" << curent[vec].ind[pp].total_c << endl;
				pp++;
			}
			
			break;
		}
	}
}


int CM2B_MOEAD2::get_minID(population pop) {
	double min_totalC = 1000000000;
	double min_fitness = 1000000000;
	int minID = 0;
	for (int i = 0; i < pop_m; ++i) {
		if (min_totalC > pop.ind[i].total_c) {
			min_totalC = pop.ind[i].total_c;
			minID = i;
		}
		else if (min_totalC == pop.ind[i].total_c&&min_fitness > pop.ind[i].fitness) {
			min_fitness = pop.ind[i].fitness;
			minID = i;
		}
	}

	return minID;
}


void CM2B_MOEAD2::parentselection(int vec,int& i1, int& i2, int& n1, int& n2) {
	double delta = 0.5;


/*	if (nextDouble() < delta) {
		i1 = vec;
		i2 = vec;
		n1 = tournament_totalC(2, curent[i1]);
		n2 = tournament_totalC(2, curent[i2]);
	}
	else {
		if (nextDouble() < 0.5) {
			i1 = vec;
			i2 = genrand_int32() % num_vec;
			n1 = get_minID(curent[i1]);
			n2 = tournament_totalC(2, curent[i2]);
		}
		else {
			i1 = genrand_int32() % num_vec;
			i2 = vec;
			n1 = tournament_totalC(2, curent[i2]);
			n2 = get_minID(curent[i1]);
		}
	}*/

	if (nextDouble() < delta) {
		i1 = vec;
		i2 = vec;
		n1 = tournament_totalC(2, curent[i1]);
		n2 = tournament_totalC(2, curent[i2]);
	}
	else {
		if (nextDouble() < 0.5) {
			i1 = vec;
			i2 = genrand_int32() % num_vec;
		}
		else {
			i1 = genrand_int32() % num_vec;
			i2 = vec;
		}
	}
	n1 = tournament_totalC(2, curent[i1]);
	n2 = tournament_totalC(2, curent[i2]);
}

void CM2B_MOEAD2::update_constraint_max(vector<population> pop) {
	constraint_max.resize(num_const);
}