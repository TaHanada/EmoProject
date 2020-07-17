#include"MOEAD_ACDP.h"
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
void CM2B_MOEAD::range_setting(int ProblemID) {
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


void MOEAD_ACDP::Parameter_setting(string algorithm, int ProblemID,int ProblemSubID) {
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
		if (ProblemID == 6) {
			num_vec = 210;
			H1 = 6, H2 = 0; //このとき，分割数は6．
		}
		break;
	}


	//others
	pop_m = num_vec;
	if (ProblemID == 6) {
		gen = 47;
	}
	else if (ProblemID != 6) {
		gen = 100;
	}
	crossover_probability = 1.0;
	mutation_probability = 1.0 / (double)dim;
	weight.resize(num_vec, num_obj, H1, H2, 10, 10);//コンストラクタの代わり
	weight.initialize();//重みベクトルの生成と距離の計算
	scalar= "tch";
	theta0_ = M_PI / 2.0 / (double)pop_m;
}


//consol output population
void MOEAD_ACDP::cout_population(int swit) {
	switch (swit) {
	case 1:
		curent.cout_population();
		break;

	case 2:
		offspring.cout_population();
		break;

	default:
		break;
	}
}


//consol output property
void MOEAD_ACDP::cout_property(int swit) {
	switch (swit) {
	case 1:
		curent.cout_property();
		break;

	case 2:
		offspring.cout_property();
		break;

	default:
		break;
	}
}


//weighted sum
double MOEAD_ACDP::weighted_sum(individual ind) {
	double value = 0;

	for (int i = 0; i < num_obj; ++i) {
		value += ind.cost[i] * weight.vec[ind.weight_no][i];
	}

	return value;
}


//Tchebycheff no absolute
double MOEAD_ACDP::tchebycheff_notabsolute(individual ind) {
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
double MOEAD_ACDP::tchebycheff(individual ind) {
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
double MOEAD_ACDP::PBI(individual ind) {
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
double MOEAD_ACDP::IPBI(individual ind) {
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
double MOEAD_ACDP::AOF(individual ind) {
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
double MOEAD_ACDP::scalarizing_function(individual ind) {
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
void MOEAD_ACDP::total_constraints(individual& ind) {
	ind.total_c = 0;

	for (int i = 0; i < num_const; ++i) {
		if (ind.constraint[i] < 0) {
			ind.total_c += -1.0*ind.constraint[i];
		}
	}
}


//initialize reference point and nadir point
void MOEAD_ACDP::initialize_reference_nadir_point() {
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
void MOEAD_ACDP::initialize_population() {

	for (int i = 0; i < num_vec; ++i) {
		for (int j = 0; j < dim; ++j) {
			curent.ind[i].var[j] = range[j][0] + (range[j][1] - range[j][0])*nextDouble();
		}

		//calculate costs
		curent.ind[i] = problems.evaluation(curent.ind[i]);
		curent.ind[i].weight_no = i;
		curent.ind[i].fitness = scalarizing_function(curent.ind[i]);//cal scalarizing function value
		total_constraints(curent.ind[i]);//cal total constraints value


		//update reference point and nadir point
		for (int j = 0; j < num_obj; ++j) {
			if (reference_point[j] > curent.ind[i].cost[j]) {
				reference_point[j] = curent.ind[i].cost[j];
			}
			if (nadir_point[j] < curent.ind[i].cost[j]) {
				nadir_point[j] = curent.ind[i].cost[j];
			}
		}
	}

	//file out solutions
	file_allsolutions(rep, 0, scalar, curent);
}


//curent solution's objective values
void MOEAD_ACDP::file_objectives(int k, int island) {
	string name = "data_" + to_string(k) + ".dat";
	ofstream fout("./MOEAD_ACDP/" + to_string(ProblemSubID) + "/" + scalar + "/island" + to_string(island + 1) + "/" + name);

	if (island == 0) {
		for (int i = 0; i < num_vec; ++i) {
			if (curent.ind[i].total_c == 0) {
				for (int k = 0; k < num_obj; ++k) {
					fout << setprecision(10) << curent.ind[i].cost[k] << "\t";
				}
				for (int k = 0; k < num_const; ++k) {
					fout << setprecision(10) << curent.ind[i].constraint[k] << "\t";
				}
				for (int k = 0; k < dim - 1; ++k) {
					fout << setprecision(10) << curent.ind[i].var[k] << "\t";
				}
				fout << setprecision(10) << curent.ind[i].var[dim - 1] << endl;
			}
		}
	}
	else if (island > 0) {
		for (int i = 0; i < num_vec; ++i) {
			if (curent.ind[i].total_c != 0) {
				for (int k = 0; k < num_obj; ++k) {
					fout << setprecision(10) << curent.ind[i].cost[k] << "\t";
				}
				for (int k = 0; k < num_const; ++k) {
					fout << setprecision(10) << curent.ind[i].constraint[k] << "\t";
				}
				for (int k = 0; k < dim - 1; ++k) {
					fout << setprecision(10) << curent.ind[i].var[k] << "\t";
				}
				fout << setprecision(10) << curent.ind[i].var[dim - 1] << endl;
			}
		}
	}
}


//evolution (mu, 1) 
void MOEAD_ACDP::evolution_one(int gen) {
	int n1 = 0, n2 = 0;//parent ID
	int i1, i2;//parent ID
		int count = 0;
	//generate permutation of weight vector
	vector<int> permutation;
	for (int i = 0; i < num_vec; ++i) {
		permutation.push_back(i);
	}
	random_shuffle(permutation.begin(), permutation.end());

	for (int ii = 0; ii < num_vec; ++ii) {
		int i = permutation[ii];
		//check the number of solutions that have better scalarizing function value than minID
		if (nextDouble() < 1.0) {
			i1 = weight.neighbor[i][genrand_int32() % weight.select_size].no;
			i2 = weight.neighbor[i][genrand_int32() % weight.select_size].no;
		}
		else {
			i1 = genrand_int32() % num_vec;
			i2 = genrand_int32() % num_vec;
		}

		//SBX
		if (nextDouble() <= 1) {
			offspring.ind[i] = SBX(curent.ind[i1], curent.ind[i2]);
		}
		else {
			offspring.ind[i] = curent.ind[i1];
		}

		//repair
		for (int j = 0; j < dim; ++j) {
			if (offspring.ind[i].var[j] < range[j][0]) {
				offspring.ind[i].var[j] = range[j][0];
			}
			else if (offspring.ind[i].var[j] > range[j][1]) {
				offspring.ind[i].var[j] = range[j][1];
			}
		}


		//polynomial mutation
		offspring.ind[i] = PM(offspring.ind[i]);


		//repair
		for (int j = 0; j < dim; ++j) {
			if (offspring.ind[i].var[j] < range[j][0]) {
				offspring.ind[i].var[j] = range[j][0];
			}
			else if (offspring.ind[i].var[j] > range[j][1]) {
				offspring.ind[i].var[j] = range[j][1];
			}
		}


		//evaluation and update
		offspring.ind[i] = problems.evaluation(offspring.ind[i]);
		total_constraints(offspring.ind[i]);

		//update ref and nadir
		for (int j = 0; j < num_obj; ++j) {
			if (reference_point[j] > offspring.ind[i].cost[j]) {
				reference_point[j] = offspring.ind[i].cost[j];
			}
			if (nadir_point[j] < offspring.ind[i].cost[j]) {
				nadir_point[j] = offspring.ind[i].cost[j];
			}
		}

		int count = 0;
		double pf=0;
		for (int i = 0; i < pop_m; ++i) {
			if (curent.ind[i].total_c == 0) count++;
		}
		pf = (double)count / (double)pop_m;
		for (int j = 0; j < weight.update_size; ++j) {
			int update_vec = weight.neighbor[i][j].no;
			if (update_vec == 51) count++;
			//assign weight vector
			curent.ind[update_vec].weight_no = update_vec;
			offspring.ind[i].weight_no = update_vec;

			//calc scalarizing function value
			curent.ind[update_vec].fitness = scalarizing_function(curent.ind[update_vec]);
			offspring.ind[i].fitness = scalarizing_function(offspring.ind[i]);

			//compare by CDP manner
			if (curent.ind[update_vec].total_c == 0 && offspring.ind[i].total_c == 0) {
				if (curent.ind[update_vec].fitness > offspring.ind[i].fitness) curent.ind[update_vec] = offspring.ind[i];
			}
			else if (get_angle(curent.ind[update_vec],offspring.ind[i])<=theta_) {
				if (curent.ind[update_vec].total_c > offspring.ind[i].total_c) {
					curent.ind[update_vec] = offspring.ind[i];
				}
			}
			else if (get_angle(curent.ind[update_vec], offspring.ind[i]) > theta_&&nextDouble() < pf) {
				if (curent.ind[update_vec].fitness > offspring.ind[i].fitness) {
					curent.ind[update_vec] = offspring.ind[i];
				}
			}
		}
	}//for (int ii = 0; ii < num_vec; ++ii)

	//file output
	file_allsolutions(rep, gen, scalar, offspring);
}


//run
void MOEAD_ACDP::run_algorithm(int reps) {
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

	//CDMP
	for (int island = 0; island < 2; ++island) {
		curent_fout(reps, 0, island);
	}

	//evolution
	for (int j = 0; j < gen; ++j) {
		if (j%10==0) cout << "generation: " << j << endl;

		//update theta
		update_theta(j);

		//curent　population
		evolution_one(j + 1);

		//CDMP
		for (int island = 0; island < 2; ++island) {
			curent_fout(reps, j + 1, island);
		}
	}


	//file output of final population
	for (int island = 0; island < 2; ++island) {
		file_objectives(reps, island);
		//file_objectives(reps + 1, island);
	}
}



//initialization files
void MOEAD_ACDP::file_initialization(int  k, int gen) {
	int count = 0;
	int tmp_gen = 0;

	tmp_gen = gen;
	//offspirngs
	{
		string			name1 = to_string(k) + "_" + to_string(0) + "gen.dat";
		ofstream fout("./MOEAD_ACDP/" + to_string(ProblemSubID) + "/" + scalar + "/offsprings/" + name1);
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
		fout << "x" << to_string(dim) << endl;
		//fout << "f" << to_string(num_obj) << endl;


		for (int i = 1; i <= tmp_gen; ++i) {
			string name1 = to_string(k) + "_" + to_string(i) + "gen.dat";
			ofstream fout("./MOEAD_ACDP/" + to_string(ProblemSubID) + "/" + scalar + "/offsprings/" + name1);
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
			//fout << "f" << to_string(num_obj) << endl;

		}
	}


	//現個体
	for (int island = 0; island < 2; ++island) {
		string name1 = to_string(k) + "_" + to_string(0) + "gen.dat";
		ofstream fout("./MOEAD_ACDP/" + to_string(ProblemSubID) + "/" + scalar + "/island" + to_string(island + 1) + "/curents/" + name1);
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
			ofstream fout("./MOEAD_ACDP/" + to_string(ProblemSubID) + "/" + scalar + "/island" + to_string(island + 1) + "/curents/" + name1);
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
void MOEAD_ACDP::file_allsolutions(int k, int gen, string scalar, population pop) {
	{
		string name1 = to_string(k) + "_" + to_string(gen) + "gen.dat";
		ofstream fout("./MOEAD_ACDP/" + to_string(ProblemSubID) + "/" + scalar + "/offsprings/" + name1, ios::out | ios::app);
		for (int k = 0; k < num_vec; ++k) {
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
				fout << setprecision(10) << pop.ind[k].cost[i] << "\t";
			}
			for (int i = 0; i < num_const; i++) {
				fout << setprecision(10) << pop.ind[k].constraint[i] << "\t";
			}
			for (int i = 0; i < dim - 1; ++i) {
				fout << setprecision(10) << pop.ind[k].var[i] << "\t";
			}
			fout << setprecision(10) << pop.ind[k].var[dim - 1] << endl;
		}
		
	}
}


int MOEAD_ACDP::calc_cossim(individual offspring1) {
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


individual MOEAD_ACDP::SBX(individual ind1, individual ind2) {
	individual offspring1, offspring2;
	double mu = 20.0;//distribution index
	double SBX_u = 0;
	double beta1, beta2;
	double alpha;
	double beta;

	offspring1 = ind1;
	offspring2 = ind2;

	for (int j = 0; j < dim; ++j) {


		if (nextDouble() <= 0.5) {
			if (fabs(ind1.var[j] - ind2.var[j]) > pow(10, -14)) {

				if (ind1.var[j] < ind2.var[j]) {
					beta1 = 1.0 + 2.0*(ind1.var[j] - range[j][0]) / (ind2.var[j] - ind1.var[j]);
					beta2 = 1.0 + 2.0*(range[j][1] - ind2.var[j]) / (ind2.var[j] - ind1.var[j]);

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
						offspring1.var[j] = 0.5*((1.0 + beta)*ind1.var[j] + (1.0 - beta)*ind2.var[j]);
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
						offspring2.var[j] = 0.5*((1.0 + beta)*ind2.var[j] + (1.0 - beta)*ind1.var[j]);

						offspring1.var[j] = offspring2.var[j];
					}


				}
				else if (ind1.var[j] > ind2.var[j]) {
					beta1 = 1.0 + 2.0*(ind2.var[j] - range[j][0]) / (ind1.var[j] - ind2.var[j]);
					beta2 = 1.0 + 2.0*(range[j][1] - ind1.var[j]) / (ind1.var[j] - ind2.var[j]);

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


						offspring1.var[j] = 0.5*((1.0 + beta)*ind2.var[j] + (1.0 - beta)*ind1.var[j]);

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
						offspring2.var[j] = 0.5*((1.0 + beta)*ind1.var[j] + (1.0 - beta)*ind2.var[j]);
						offspring1.var[j] = offspring2.var[j];
					}

				}
			}//if (fabs(curent.ind[n1].var[j] - curent.ind[n2].var[j]) > pow(10, -14))
		}//if (nextDouble() <= 0.5)
	}

	return offspring1;
}


individual MOEAD_ACDP::PM(individual ind) {
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
			ind.var[j] += del*(range[j][1] - range[j][0]);
		}
	}

	return ind;
}


double MOEAD_ACDP::get_angle(individual ind1, individual ind2) {
	double angle = 0;
	double upper = 0, bottom1 = 0, bottom2 = 0;

	for (int i = 0; i < num_obj; ++i) {
		upper += (ind1.cost[i] - reference_point[i])*(ind2.cost[i] - reference_point[i]);
		bottom1 += pow(ind1.cost[i] - reference_point[i], 2);
		bottom2 += pow(ind2.cost[i] - reference_point[i], 2);
	}
	angle = acos(upper / (sqrt(bottom1)*sqrt(bottom2) + pow(10, -5)));

	return angle;
}

void MOEAD_ACDP::update_theta(int curent_gen) {
	double cp_ = log(M_PI / 2.0 / theta0_) / log(1 + alpha_);

	if (curent_gen < alpha_*gen) {
		theta_ = theta0_*pow(1 + (double)curent_gen / (double)gen,cp_);
	}
	else {
		theta_ = M_PI / 2.0;
	}
}

void MOEAD_ACDP::curent_fout(int k, int j, int island) {
	string name1 = to_string(k) + "_" + to_string(j) + "gen.dat";
	ofstream fout("./MOEAD_ACDP/" + to_string(ProblemSubID) + "/" + scalar + "/island" + to_string(island + 1) + "/curents/" + name1, ios::out | ios::app);

	if (island == 0) {
		for (int i = 0; i < pop_m; ++i) {
			if (curent.ind[i].total_c == 0) {
				fout << j << "\t" << "1" << "\t";
				for (int k = 0; k < num_obj; ++k) {
					fout << setprecision(10) << curent.ind[i].cost[k] << "\t";
				}
				for (int k = 0; k < num_const; ++k) {
					fout << setprecision(10) << curent.ind[i].constraint[k] << "\t";
				}
				for (int k = 0; k < dim - 1; ++k) {
					fout << setprecision(10) << curent.ind[i].var[k] << "\t";
				}
				fout << setprecision(10) << curent.ind[i].var[dim - 1] << endl;
			}
		}
	}
	else if (island > 0) {
		for (int i = 0; i < pop_m; ++i) {
			if (curent.ind[i].total_c != 0) {
				fout << j << "\t" << "0" << "\t";
				for (int k = 0; k < num_obj; ++k) {
					fout << setprecision(10) << curent.ind[i].cost[k] << "\t";
				}
				for (int k = 0; k < num_const; ++k) {
					fout << setprecision(10) << curent.ind[i].constraint[k] << "\t";
				}
				for (int k = 0; k < dim - 1; ++k) {
					fout << setprecision(10) << curent.ind[i].var[k] << "\t";
				}
				fout << setprecision(10) << curent.ind[i].var[dim - 1] << endl;
			}
		}
	}
}