#include"DG_MOEAD.h"
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
			cout << "error?" << endl;
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


void DG_MOEAD::Parameter_setting(string algorithm, int ProblemID, int ProblemSubID) {
	//Problem specific parameters
	problems.set_ProblemID(ProblemID,ProblemSubID);
	problems.parameter_setting(num_obj, num_const, dim, range);


	//num_vec
	switch (num_obj) {
	case 2:
		num_vec = 100;
		H1 = num_vec - 1, H2 = 0;//������ 12(M=3) 6(M=5) 3,2(M=8) 3,2(M=10)
		break;

	case 3:
		num_vec = 300;
		H1 = 23, H2 = 0;//������ 12(M=3) 6(M=5) 3,2(M=8) 3,2(M=10)
		break;
	default:
		break;
	}


	//others
	pop_m = 2;
	pop_c = 1;
	gen = 500;
	crossover_probability = 1.0;
	mutation_probability = 1.0 / (double)dim;
	weight.resize(num_vec, num_obj, H1, H2, num_vec/10, num_vec/10);//�R���X�g���N�^�̑���
	weight.initialize();//�d�݃x�N�g���̐����Ƌ����̌v�Z
	scalar= "tch";
	cout << mutation_probability << endl;
}


//consol output population
void DG_MOEAD::cout_population(int swit) {
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
void DG_MOEAD::cout_property(int swit) {
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
double DG_MOEAD::weighted_sum(individual ind) {
	double value = 0;

	for (int i = 0; i < num_obj; ++i) {
		value += ind.cost[i] * weight.vec[ind.weight_no][i];
	}

	return value;
}


//Tchebycheff no absolute
double DG_MOEAD::tchebycheff_notabsolute(individual ind) {
	double max = 0;
	double g = 0;

	for (int i = 0; i < num_obj; ++i) {
		g = weight.vec[ind.weight_no][i] * ((reference_point[i] - ind.cost[i]));
		if (max < g) {
			max = g;
		}
	}

	return max;//�ŏ�������Ă���
}


//Tchebycheff
double DG_MOEAD::tchebycheff(individual ind) {
	double max = 0;
	double g = 0;

	for (int i = 0; i < num_obj; ++i) {
		if (weight.vec[ind.weight_no][i] == 0) {
			g = pow(10, -6) * abs((0.0 - (ind.cost[i] - reference_point[i]) / (nadir_point[i] - reference_point[i])));//normalization
		}
		else {
	//		cout << ind.weight_no << " " << ind.cost[i] << " " << reference_point[i] << " " << nadir_point[i] << endl;
			g = weight.vec[ind.weight_no][i] * abs((0.0 - (ind.cost[i] - reference_point[i]) / (nadir_point[i] - reference_point[i])));//normalization
		}

		if (max < g) {
			max = g;
		}
	}
	
/*	for (int i = 0; i < num_obj; ++i) {
		if (weight.vec[ind.weight_no][i] == 0) {
			g = pow(10, -6) * abs(reference_point[i] - ind.cost[i]);
		}
		else {
			g = weight.vec[ind.weight_no][i] * abs(reference_point[i] - ind.cost[i]);
		}

		if (max < g) {
			max = g;
		}
	}*/
	return max;
}


//PBI
double DG_MOEAD::PBI(individual ind) {
	double theta = 5.0;//�p�����[�^
	double value = 0;
	double dt = 0;
	double dn = 0;
	double ramda = 0;


	//dt�̌v�Z
	for (int i = 0; i < num_obj; ++i) {
		//���q
		dt += (reference_point[i] - ind.cost[i])*weight.vec[ind.weight_no][i];

		//����
		ramda += pow(weight.vec[ind.weight_no][i], 2);
	}
	dt = dt / sqrt(ramda);



	//dn�̌v�Z
	for (int i = 0; i < num_obj; ++i) {
		dn += pow((reference_point[i] - ind.cost[i]) - dt*weight.vec[ind.weight_no][i] / sqrt(ramda), 2);
	}
	dn = sqrt(dn);

	value = abs(dt) + theta*dn;
	
	return value;
}


//IPBI
double DG_MOEAD::IPBI(individual ind) {
	double theta = 0.1;//�p�����[�^
	double value = 0;
	double dt = 0;
	double dn = 0;
	double ramda = 0;


	//dt�̌v�Z
	for (int i = 0; i < num_obj; ++i) {
		//���q
		dt += (nadir_point[i] - ind.cost[i])*weight.vec[ind.weight_no][i];

		//����
		ramda += pow(weight.vec[ind.weight_no][i], 2);
	}
	dt = dt / sqrt(ramda);



	//dn�̌v�Z
	for (int i = 0; i < num_obj; ++i) {
		dn += pow((nadir_point[i] - ind.cost[i]) - dt*weight.vec[ind.weight_no][i] / sqrt(ramda), 2);
	}
	dn = sqrt(dn);

	value = abs(dt) - theta*dn;

	return value;
}


//AOF
double DG_MOEAD::AOF(individual ind) {
	int k = ind.weight_no;

	double p = -1.0*ind.cost[1] + k/4;
	double value=0;
	
	/*	double theta = 0.1;//�p�����[�^
	double value = 0;

	double ramda = 0;


	//dt�̌v�Z
	for (int i = 0; i < num_obj; ++i) {
		//���q
		dt += (nadir_point[i] - ind.cost[i])*weight.vec[ind.weight_no][i];

		//����
		ramda += pow(weight.vec[ind.weight_no][i], 2);
	}
	dt = dt / sqrt(ramda);



	//dn�̌v�Z
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
double DG_MOEAD::scalarizing_function(individual ind) {
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
void DG_MOEAD::total_constraints(individual& ind) {
	ind.total_c = 0;

	for (int i = 0; i < num_const; ++i) {
		if (ind.constraint[i] < 0) {
			ind.total_c += -1.0*ind.constraint[i];
		}
	}
}


//initialize reference point and nadir point
void DG_MOEAD::initialize_reference_nadir_point() {
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
void DG_MOEAD::initialize_population() {
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


	//file output
	for (int i = 0; i < update_population.size(); ++i) {
		temp_offspring.ind[i] = update_population[i];
	}

	problems.fileout(1, 0, rep, temp_offspring);


	for (int j = 0; j < update_population.size(); ++j) {
		//adjustment of weight vector
		curent[j].ind[0] = update_population[j];
		curent[j].ind[1] = update_population[j];
	}
}


//curent solution's objective values
void DG_MOEAD::file_objectives(int k, int island) {
	string name = "data_" + to_string(k) + ".dat";
	ofstream fout("./" + to_string(ProblemSubID) + "/tch/island" + to_string(island + 1) + "/" + name);

	if (island == 0) {
		for (int i = 0; i < num_vec; ++i) {
			for (int j = 0; j < pop_m; ++j) {
				if (curent[i].ind[j].total_c == 0) {
					for (int k = 0; k < num_obj-1; ++k) {
						fout << setprecision(7) << max(curent[i].ind[j].cost[k],-1* curent[i].ind[j].cost[k]) << "\t";
					}
					fout << setprecision(7) << max(curent[i].ind[j].cost[num_obj - 1], -1* curent[i].ind[j].cost[num_obj - 1]) << endl;
				}
			}
		}
	}
	else if (island > 0) {
		for (int i = 0; i < num_vec; ++i) {
			for (int j = 0; j < pop_m; ++j) {
				if (curent[i].ind[j].total_c != 0) {
					for (int k = 0; k < num_obj-1; ++k) {
						fout << setprecision(7) << max(curent[i].ind[j].cost[k], -1 * curent[i].ind[j].cost[k]) << "\t";
					}
					fout << setprecision(7) << max(curent[i].ind[j].cost[num_obj - 1], -1 * curent[i].ind[j].cost[num_obj - 1]) << endl;
				}
			}
		}
	}
}


//evolution (mu, mu) 
void DG_MOEAD::evolution_one(int curent_gen) {
	int n1 = 0, n2 = 0, n3=0, n4=0;//parent ID
	int i1, i2, i3, i4;//parent ID

	//generate permutation of weight vector
	vector<int> permutation;
	for (int i = 0; i < num_vec; ++i) {
		permutation.push_back(i);
	}
	random_shuffle(permutation.begin(), permutation.end());


	for (int ii = 0; ii < num_vec; ++ii) {
		int i = permutation[ii];


		for (int pp = 0; pp < pop_c; ++pp) {
				//crossover
				if (nextDouble() <= 1) {
//					offspring[i].ind[pp] = SBX(i1, i2, n1, n2);
					offspring[i].ind[pp] = DE(i, curent_gen);
				}
				else {
					offspring[i].ind[pp] = curent[i].ind[n1];
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


				//update each sub-population
				vector<int> perm_update;
				for (int j = 0; j < weight.update_size; ++j) {
					perm_update.push_back(weight.neighbor[i][j].no);
				}
				random_shuffle(perm_update.begin(), perm_update.end());

				for (int j = 0; j < weight.update_size; ++j) {
					int x = perm_update[j];

					vector<individual> update_population;
					update_population.push_back(offspring[i].ind[pp]);

					for (int k = 0; k < pop_m; ++k) {
						update_population.push_back(curent[x].ind[k]);
					}

					for (int k = 0; k < update_population.size(); ++k) {
						update_population[k].weight_no = x;
						update_population[k].fitness = scalarizing_function(update_population[k]);
					}

					select_next_ind(update_population, x);
				}
		}//for (int pp = 0; pp < pop_c; ++pp)
	}//for (int ii = 0; ii < num_vec; ++ii)


	for (int i = 0; i < num_vec; ++i) {
		for (int j = 0; j < pop_c; ++j) {
			temp_offspring.ind[i*pop_c + j] = offspring[i].ind[j];
		}
	}


	//file output
	problems.fileout(1, curent_gen, rep, temp_offspring);
}


//run
void DG_MOEAD::run_algorithm(int reps) {
	rep = reps;
	srand(reps);
	init_genrand(reps);
	//generate files
	cout << rep << endl;
	file_initialization(rep + 1, gen);

	//weight vector
	weight.reselect_nighbor();
	weight.shuffle();


	//initialize ref. and nadir points
	initialize_reference_nadir_point();


	//inititalize population
	initialize_population();

//	cout << "initialize population" << endl;

	//evolution
	for (int j = 0; j < gen; ++j) {
		if (j%100==0) cout << "generation: " << j << endl;

		//curent�@population
		evolution_one(j + 1);

		for (int i = 0; i < num_vec; ++i) {
			for (int j = 0; j < pop_m; ++j) {
				temp_curent.ind[i*pop_m + j] = curent[i].ind[j];
			}
		}
		problems.fileout(2, j + 1, rep, temp_curent);
	}


	//file output of final population
	for (int island = 0; island < 2; ++island) {
		file_objectives(reps + 1, island);
	}
}


//initialization files
void DG_MOEAD::file_initialization(int  k, int gen) {
	int count = 0;
	int tmp_gen = 0;

	tmp_gen = gen;
	//offspirngs
	{
		string			name1 = to_string(k) + "_" + to_string(0) + "gen.dat";
		ofstream fout("./" + to_string(ProblemSubID) + "/tch/offsprings/" + name1);
		fout << "#gen	Feasible	";
		for (int i = 0; i < num_obj; ++i) {
			fout << "f" << to_string(i+1) <<"	";
		}
		fout << "ras" << endl;


		for (int i = 1; i <= tmp_gen; ++i) {
			string name1 = to_string(k) + "_" + to_string(i) + "gen.dat";
			ofstream fout("./" + to_string(ProblemSubID) + "/" + scalar + "/offsprings/" + name1);
			fout << "#gen	Feasible	";
			for (int i = 0; i < num_obj; ++i) {
				fout << "f" << to_string(i + 1) << "	";
			}
			fout << "ras" << endl;

		}
	}


	//����
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
void DG_MOEAD::file_allsolutions(int k, int gen, string scalar, vector<population> pop, vector<double> ras) {
	{
		string name1 = to_string(k) + "_" + to_string(gen) + "gen.dat";
		ofstream fout("./"+to_string(ProblemSubID)+"/tch/offsprings/" + name1, ios::out | ios::app);
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
					fout << max(pop[k].ind[j].cost[i], -1*pop[k].ind[j].cost[i]) << "\t";
				}

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


//�ړI�֐��l�̃v���b�g(����)
void DG_MOEAD::file_curents(int k, int gen, string scalar, int island) {

	string name1 = to_string(k) + "_" + to_string(gen) + "gen.dat";
	ofstream fout("./" + scalar + "/island" + to_string(island + 1) + "/curents/" + name1, ios::out | ios::app);

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
			fout << curent[j].ind[minID].cost[0] << "\t" << max(curent[j].ind[minID].cost[1] , -1*curent[j].ind[minID].cost[1]) << "\t";


			//constraints
			for (int i = 0; i < 54; ++i) {
				fout << curent[j].ind[minID].constraint[i] << "\t";
			}

			//var
			for (int i = 0; i < 221; i++)
			{
				fout << curent[j].ind[minID].var[i] << "\t";
			}
			fout << curent[j].ind[minID].var[221] << endl;
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
				fout << curent[k].ind[j].cost[0] << "\t" << -1.0*curent[k].ind[j].cost[1] << "\t";


				//constraints
				for (int i = 0; i < 54; ++i) {
					fout << curent[k].ind[j].constraint[i] << "\t";
				}

				//var
				for (int i = 0; i < 221; i++)
				{
					fout << curent[k].ind[j].var[i] << "\t";
				}
				fout << curent[k].ind[j].var[221] << endl;
			}
		}
	}
}


int DG_MOEAD::calc_cossim(individual offspring1) {
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


individual DG_MOEAD::SBX(int vec1, int vec2, int n1, int n2) {
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


individual DG_MOEAD::DE(int vec, int curent_gen) {
	individual offspring1, ind1, ind2, ind3;
	double F = 0.5;//distribution index
	double CR = 1.0;
	int jrand = genrand_int32() % dim;

	//parent selection
	ind1 = curent[vec].ind[genrand_int32() % pop_m];
	ind2= curent[weight.neighbor[vec][genrand_int32() % weight.select_size].no].ind[genrand_int32() % pop_m];
	ind3 = curent[weight.neighbor[vec][genrand_int32() % weight.select_size].no].ind[genrand_int32() % pop_m];

	//DE operator
	offspring1 = ind1;
	for (int j = 0; j < dim; ++j) {
		if (nextDouble() < CR || j == jrand) offspring1.var[j] = ind1.var[j] + F*(ind2.var[j] - ind3.var[j]);
	}

	return offspring1;
}



void DG_MOEAD::PM(int vec, int pp) {
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


void DG_MOEAD::select_next_ind(vector<individual> update_population, int vec) {
	int patternA = 0, patternB = 0, patternC = 0;
	

	//patternA
	for (int i = 1; i < update_population.size(); ++i) {
		if (update_population[i].total_c < update_population[patternA].total_c) {
			patternA = i;
		}
		else if (update_population[i].total_c == update_population[patternA].total_c) {
			if (update_population[i].fitness < update_population[patternA].fitness) {
				patternA = i;
			}
		}
	}

	//update!!!!
	curent[vec].ind[0] = update_population[patternA];
	patternB = (patternA + 1) % 3;
	patternC = (patternA + 2) % 3;

	if (DR2(update_population[patternA], update_population[patternB])==1 && DR2(update_population[patternA], update_population[patternC])==1) {
		curent[vec].ind[1]= update_population[patternA];
	}
	else if (DR2(update_population[patternA], update_population[patternB])==1 || DR2(update_population[patternC], update_population[patternB])==1) {
		curent[vec].ind[1] = update_population[patternC];
	}
	else if (DR2(update_population[patternA], update_population[patternC]) == 1 || DR2(update_population[patternB], update_population[patternC]) == 1) {
		curent[vec].ind[1] = update_population[patternB];
	}
	else {
		if (nextDouble() < 0.5) curent[vec].ind[1] = update_population[patternB];
		else curent[vec].ind[1] = update_population[patternC];
	}
}

int DG_MOEAD::DR2(individual a,individual b) {

	if (a.fitness == b.fitness&&a.total_c == b.total_c) {
		return 0;
	}
	else if (a.fitness <= b.fitness&&a.total_c <= b.total_c) {
		return 1;
	}
	else if (a.fitness >= b.fitness&&a.total_c >= b.total_c) {
		return -1;
	}
	else {
		return 0;
	}

}

int DG_MOEAD::get_minID(population pop) {
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