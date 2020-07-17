#include"DNEA_CDP.h"
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
cout << "���̓t�@�C�����I�[�v���ł��܂���" << endl;
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


void DNEA_CDP::Parameter_setting(string algorithm, int ProblemID, int ProblemSubID) {
	//Problem specific parameters
	problems.set_ProblemID(ProblemID, ProblemSubID);
	problems.parameter_setting(num_obj, num_const, dim, range);


	//num_vec
	switch (num_obj) {
	case 2:
		pop_m = 100;
		break;

	case 3:
		pop_m = 200;
		break;
	default:
		pop_m = 200;
		break;
	}

	//others
	gen = 49;
	crossover_probability = 0.9;
	mutation_probability = 1.0 / (double)dim;
	niche_obj = 1.0;
	niche_var = 1.0;

	switch (ProblemID) {
	case 8://disc_brake_design_problem
		pop_m = 100;
		gen = 49;
		crossover_probability = 0.9;
		mutation_probability = 1.0 / (double)dim;
		break;
	}
}


//consol output population
void DNEA_CDP::cout_population(int swit) {
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
void DNEA_CDP::cout_property(int swit) {
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


//calcurate total constraints
void DNEA_CDP::total_constraints(individual& ind) {
	ind.total_c = 0;

	for (int i = 0; i < num_const; ++i) {
		if (ind.constraint[i] < 0) {
			ind.total_c += -1.0 * ind.constraint[i];
		}
	}
}


//initialize population
void DNEA_CDP::initialize_population() {

	for (int i = 0; i < pop_m; ++i) {

		if (ProblemID == 3) {
			for (int j = 0; j < dim; ++j) {
				curent.ind[i].var[j] = 30 + 10 * nextDouble();
			}
		}
		else {
			for (int j = 0; j < dim; ++j) {
				curent.ind[i].var[j] = range[j][0] + (range[j][1] - range[j][0]) * nextDouble();
			}
		}


		//calculate costs
		curent.ind[i] = problems.evaluation(curent.ind[i]);
		total_constraints(curent.ind[i]);//cal total constraints value
	}

	//calc. fitness
	evaluation(curent);


	//file out solutions
	file_allsolutions(rep, 0, curent);
}


//current solution's objective values
void DNEA_CDP::file_objectives(int k, int island) {
	string name = "data_" + to_string(k) + ".dat";
	ofstream fout("./DNEA_CDP/" + to_string(ProblemSubID) + "/" + scalar + "/island" + to_string(island + 1) + "/" + name);
	/*fout << "#gen	Feasible	";
	for (int i = 0; i < num_obj; ++i) {
		fout << "f" << to_string(i + 1) << "	";
	}
	for (int i = 0; i < num_const; ++i) {
		fout << "c" << to_string(i + 1) << "	";
	}
	for (int i = 0; i < dim - 1; ++i) {
		fout << "x" << to_string(i + 1) << "	";
	}
	fout << "x" << to_string(dim) << endl;*/

	vector<double> maxcost(num_obj, 0);
	vector<double> mincost(num_obj, 1000000);
	if (island == 0) {
		for (int i = 0; i < pop_m; ++i) {
			if (curent.ind[i].total_c == 0) {
				for (int k = 0; k < num_obj; ++k) {
					//fout << setprecision(7) << max(curent.ind[i].cost[k], -1 * curent.ind[i].cost[k]) << "\t";
					fout << setprecision(10) << curent.ind[i].cost[k] << "\t";
					if (curent.ind[i].cost[k] > maxcost[k]) {
						maxcost[k] = curent.ind[i].cost[k];
					}
					if (curent.ind[i].cost[k] < mincost[k]) {
						mincost[k] = curent.ind[i].cost[k];
					}
				}
				for (int k = 0; k < dim; ++k) {
					fout << setprecision(10) << curent.ind[i].var[k] << "\t";
				}
				for (int k = 0; k < num_const - 1; ++k) {
					fout << setprecision(10) << curent.ind[i].constraint[k] << "\t";
				}
				fout << setprecision(10) << curent.ind[i].constraint[num_const - 1] << endl;
			}
		}
	}
	else if (island > 0) {
		for (int i = 0; i < pop_m; ++i) {
			if (curent.ind[i].total_c != 0) {
				for (int k = 0; k < num_obj; ++k) {
					//fout << setprecision(7) << max(curent.ind[i].cost[k], -1 * curent.ind[i].cost[k]) << "\t";
					fout << setprecision(10) << curent.ind[i].cost[k] << "\t";
					if (curent.ind[i].cost[k] > maxcost[k]) {
						maxcost[k] = curent.ind[i].cost[k];
					}
					if (curent.ind[i].cost[k] < mincost[k]) {
						mincost[k] = curent.ind[i].cost[k];
					}
				}
				for (int k = 0; k < dim; ++k) {
					fout << setprecision(10) << curent.ind[i].var[k] << "\t";
				}
				for (int k = 0; k < num_const - 1; ++k) {
					fout << setprecision(10) << curent.ind[i].constraint[k] << "\t";
				}
				fout << setprecision(10) << curent.ind[i].constraint[num_const - 1] << endl;
			}
		}
	}
	for (int k = 0; k < num_obj; ++k) {
		cout << k + 1 << "max:" << maxcost[k] << endl;
		cout << k + 1 << "min:" << mincost[k] << endl;
	}
}


//evolution (mu, 1) 
void DNEA_CDP::evolution(int gen) {
	int i1 = 0, i2 = 0;//parent ID
	int count = 0;

	for (int i = 0; i < pop_m; ++i) {

		//check the number of solutions that have better scalarizing function value than minID
		//select
		i1 = tournament(2, curent);
		i2 = tournament(2, curent);

		//SBX
		if (nextDouble() <= crossover_probability) {
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
	}

	//merge
	for (int i = 0; i < pop_m; ++i) {
		merge.ind[i] = curent.ind[i];
	}
	for (int i = 0; i < pop_m; ++i) {
		merge.ind[i + pop_m] = offspring.ind[i];
	}
	evaluation(merge);

	//file output
	file_allsolutions(rep, gen, offspring);
}


//run
void DNEA_CDP::run_algorithm(int reps) {
	rep = reps;
	srand(reps);
	init_genrand(reps);
	//generate files
	cout << rep << endl;
	file_initialization(rep, gen);
	//file_initialization(rep + 1, gen);


	//inititalize population
	initialize_population();


	//CDMP
	for (int island = 0; island < 2; ++island) {
		curent_fout(reps, 0, island);
	}


	//evolution
	for (int j = 0; j < gen; ++j) {
		if (j % 10 == 0) cout << "generation: " << j << endl;

		//curent�@population
		evolution(j + 1);

		//CDMP
		for (int island = 0; island < 2; ++island) {
			curent_fout(reps, j + 1, island);
		}
	}


	//file output of final population
	for (int island = 0; island < 2; ++island) {
		file_objectives(reps, island);
		//file_objectives(reps+1, island);
	}
}



//initialization files
void DNEA_CDP::file_initialization(int  k, int gen) {
	int count = 0;
	int tmp_gen = 0;

	tmp_gen = gen;
	//offspirngs
	{
		string			name1 = to_string(k) + "_" + to_string(0) + "gen.dat";

		ofstream fout("./DNEA_CDP/" + to_string(ProblemSubID) + "/" + scalar + "/offsprings/" + name1);
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


		for (int i = 1; i <= tmp_gen; ++i) {
			string name1 = to_string(k) + "_" + to_string(i) + "gen.dat";
			ofstream fout("./DNEA_CDP/" + to_string(ProblemSubID) + "/" + scalar + "/offsprings/" + name1);
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
		}
	}


	//����
	for (int island = 0; island < 2; ++island) {
		string name1 = to_string(k) + "_" + to_string(0) + "gen.dat";
		ofstream fout("./DNEA_CDP/" + to_string(ProblemSubID) + "/" + scalar + "/island" + to_string(island + 1) + "/curents/" + name1);
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
			ofstream fout("./DNEA_CDP/" + to_string(ProblemSubID) + "/" + scalar + "/island" + to_string(island + 1) + "/curents/" + name1);
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

}


//file output offsprings
void DNEA_CDP::file_allsolutions(int k, int gen, population pop) {
	{
		string name1 = to_string(k) + "_" + to_string(gen) + "gen.dat";
		ofstream fout("./DNEA_CDP/" + to_string(ProblemSubID) + "/" + scalar + "/offsprings/" + name1, ios::out | ios::app);
		for (int k = 0; k < pop_m; ++k) {
			//gen 
			fout << gen << "\t";

			//feasible
			if (pop.ind[k].total_c == 0) {
				fout << "1" << "\t";
			}
			else {
				fout << "0" << "\t";
			}

			//function value , constraint value , variable value
			for (int i = 0; i < num_obj; ++i) {
				//fout << max(pop.ind[k].cost[i], -1 * pop.ind[k].cost[i]) << "\t";
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


individual DNEA_CDP::SBX(individual ind1, individual ind2) {
	individual offspring1, offspring2;
	double mu = 20.0;//distribution index 30
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
					beta1 = 1.0 + 2.0 * (ind1.var[j] - range[j][0]) / (ind2.var[j] - ind1.var[j]);
					beta2 = 1.0 + 2.0 * (range[j][1] - ind2.var[j]) / (ind2.var[j] - ind1.var[j]);

					if (nextDouble() <= 0.5) {
						//offspring1
						SBX_u = nextDouble();
						alpha = 2.0 - pow(beta1, -1.0 * (mu + 1.0));
						if (SBX_u <= 1.0 / alpha) {
							beta = pow(SBX_u * alpha, 1.0 / (mu + 1.0));
						}
						else {
							beta = pow(2.0 - SBX_u * alpha, -1.0 / (mu + 1.0));
						}
						offspring1.var[j] = 0.5 * ((1.0 + beta) * ind1.var[j] + (1.0 - beta) * ind2.var[j]);
					}
					else {
						//offspring2
						SBX_u = nextDouble();
						alpha = 2.0 - pow(beta2, -1.0 * (mu + 1.0));
						if (SBX_u <= 1.0 / alpha) {
							beta = pow(SBX_u * alpha, 1.0 / (mu + 1.0));
						}
						else {
							beta = pow(2.0 - SBX_u * alpha, -1.0 / (mu + 1.0));
						}
						offspring2.var[j] = 0.5 * ((1.0 + beta) * ind2.var[j] + (1.0 - beta) * ind1.var[j]);

						offspring1.var[j] = offspring2.var[j];
					}


				}
				else if (ind1.var[j] > ind2.var[j]) {
					beta1 = 1.0 + 2.0 * (ind2.var[j] - range[j][0]) / (ind1.var[j] - ind2.var[j]);
					beta2 = 1.0 + 2.0 * (range[j][1] - ind1.var[j]) / (ind1.var[j] - ind2.var[j]);

					if (nextDouble() <= 0.5) {
						//offspring1
						SBX_u = nextDouble();
						alpha = 2.0 - pow(beta1, -1.0 * (mu + 1.0));
						if (SBX_u <= 1.0 / alpha) {
							beta = pow(SBX_u * alpha, 1.0 / (mu + 1.0));
						}
						else {
							beta = pow(2.0 - SBX_u * alpha, -1.0 / (mu + 1.0));
						}


						offspring1.var[j] = 0.5 * ((1.0 + beta) * ind2.var[j] + (1.0 - beta) * ind1.var[j]);

					}
					else {
						//offspring2
						SBX_u = nextDouble();
						alpha = 2.0 - pow(beta2, -1.0 * (mu + 1.0));
						if (SBX_u <= 1.0 / alpha) {
							beta = pow(SBX_u * alpha, 1.0 / (mu + 1.0));
						}
						else {
							beta = pow(2.0 - SBX_u * alpha, -1.0 / (mu + 1.0));
						}
						offspring2.var[j] = 0.5 * ((1.0 + beta) * ind1.var[j] + (1.0 - beta) * ind2.var[j]);
						offspring1.var[j] = offspring2.var[j];
					}

				}
			}//if (fabs(curent.ind[n1].var[j] - curent.ind[n2].var[j]) > pow(10, -14))
		}//if (nextDouble() <= 0.5)
	}

	return offspring1;
}


individual DNEA_CDP::PM(individual ind) {
	double mum = 30.0;//index of polynomial mutation
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
			ind.var[j] += del * (range[j][1] - range[j][0]);
		}
	}

	return ind;
}


void DNEA_CDP::evaluation(population pop) {
	//calc. ParetoRank and CrowdingDistance
	//assign rank with ParetoRank and CrowindDistance
	//select top pop_m solutions for next curent

	//calc. ParetoRank
	vector<int> frontID;
	vector<vector<int>> ParetoRankID;
	int flag_dominated = 0;
	frontID.resize(pop.size);
	int num_pushback = 0;


	while (num_pushback < pop.size) {
		vector<int> flag1;
		for (int i = 0; i < pop.size; ++i) {
			if (frontID[i] == 0) {
				flag_dominated = 0;
				for (int j = 0; j < pop.size; ++j) {
					if (frontID[j] == 0) {
						if (check_dominance(pop.ind[i], pop.ind[j]) == -1) {

							flag_dominated = 1;
							break;
						}
					}
				}
				if (flag_dominated == 0) {//if solution i is not dominated by solution in the archive
					flag1.push_back(i);
					num_pushback++;
				}
			}
		}
		for (int i = 0; i < flag1.size(); ++i) {
			frontID[flag1[i]] = 1;
		}
		if (flag1.size() != 0) {
			ParetoRankID.push_back(flag1);
		}
	}



	//assign Pareto Rank
	for (int i = 0; i < ParetoRankID.size(); ++i) {
		for (int j = 0; j < ParetoRankID[i].size(); ++j) {
			pop.ind[ParetoRankID[i][j]].fitness = 0;
			for (int k = 0; k < i; ++k) {
				pop.ind[ParetoRankID[i][j]].fitness += ParetoRankID[k].size();
			}
		}
	}

	//calc. Crowding Distance
	for (int i = 0; i < ParetoRankID.size(); ++i) {
		if (ParetoRankID[i].size() != 1) pop = assign_CrowdingDistance(ParetoRankID[i], pop);
	}

	//select top pop_m solutions
	for (int i = 0; i < pop.size; ++i) {
		if (pop.ind[i].fitness < pop_m) {
			curent.ind[pop.ind[i].fitness] = pop.ind[i];
		}
	}
}


int DNEA_CDP::tournament(int tournamentSize, population pop) {
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


int DNEA_CDP::check_dominance(individual a, individual b) {
	//if a dominate b, then return 1
	//if a is dominated by b, then return -1
	//otherwise return 0

	int i;
	int flag1;
	int flag2;
	flag1 = 0;
	flag2 = 0;
	if (a.total_c > 0 && b.total_c > 0)
	{
		if (a.total_c < b.total_c)
		{
			return (1);
		}
		else
		{
			if (a.total_c > b.total_c)
			{
				return (-1);
			}
			else
			{
				return (0);
			}
		}
	}
	else
	{
		if (a.total_c > 0 && b.total_c == 0)
		{
			return (-1);
		}
		else
		{
			if (a.total_c == 0 && b.total_c > 0)
			{
				return (1);
			}
			else
			{
				for (i = 0; i < num_obj; i++)
				{
					if (a.cost[i] < b.cost[i])
					{
						flag1 = 1;

					}
					else
					{
						if (a.cost[i] > b.cost[i])
						{
							flag2 = 1;
						}
					}
				}
				if (flag1 == 1 && flag2 == 0)
				{
					return (1);
				}
				else
				{
					if (flag1 == 0 && flag2 == 1)
					{
						return (-1);
					}
					else
					{
						return (0);
					}
				}
			}
		}
	}
}


population DNEA_CDP::assign_CrowdingDistance(vector<int> IDs, population pop) {
	vector<vector<rank_no>> cost, cost_ds;
	vector<rank_no> dist, arc, arc_ds;
	vector<vector<double>> arc2, arc2_ds;
	cost.resize(IDs.size(), vector<rank_no>(num_obj));
	cost_ds.resize(IDs.size(), vector<rank_no>(dim));
	dist.resize(IDs.size());
	arc.resize(IDs.size());
	arc_ds.resize(IDs.size());
	arc2.resize(IDs.size(), vector<double>(IDs.size()));
	arc2_ds.resize(IDs.size(), vector<double>(IDs.size()));
	double rnd = nextDouble();
	int fk = IDs.size();
	int nb = min(3, max(fk - 1, 0));

	//initialize
	for (int i = 0; i < IDs.size(); ++i) {
		dist[i].distance = 0;
	}

	// �ړI�֐���Ԃƌ���ϐ���Ԃł��ׂẲ��Ƃ̋����̑��a�I�Ȃ��̂����
	// no�ƖړI�֐��l,����ϐ��l�̊���
	for (int i = 0; i < IDs.size(); ++i) {
		for (int j = 0; j < num_obj; ++j) {
			cost[i][j].no = IDs[i];
			cost[i][j].distance = pop.ind[IDs[i]].cost[j];
		}
		for (int j = 0; j < dim; ++j) {
			cost_ds[i][j].no = IDs[i];
			cost_ds[i][j].distance = pop.ind[IDs[i]].var[j];
		}
	}
	/*
	// �ÓIver
	// �j�b�`���a�͐ݒ肳��Ă���p�����[�^
	for (int i = 0; i < IDs.size(); ++i) {
		double f_ds = 0;
		for (int j = 0; j < IDs.size(); ++j) {
			double d_obj = 0;
			double d_var = 0;
			for (int s = 0; s < num_obj; ++s) {
				d_obj += pow(cost[i][s].distance - cost[j][s].distance, 2);
			}
			for (int s = 0; s < dim; ++s) {
				d_var += pow(cost_ds[i][s].distance - cost_ds[j][s].distance, 2);
			}
			double sh = max(0.0, 1 - sqrt(d_obj) / niche_obj);
			double sh_ds = max(0.0, 1 - sqrt(d_var) / niche_var);
			f_ds += sh + sh_ds;
		}
		// dist�ɕۑ�
		dist[i].no = IDs[i];
		dist[i].distance = f_ds;
	}
	*/
	// ���Iver
	// ���������ׂČv�Z���Ă���j�b�`���a�o���Ȃ��Ⴂ���Ȃ�
	double nichea_obj = 0;
	double nichea_var = 0;
	for (int i = 0; i < IDs.size(); ++i) {
		double f_ds = 0;
		for (int j = 0; j < IDs.size(); ++j) {
			double d_obj = 0;
			double d_var = 0;
			for (int s = 0; s < num_obj; ++s) {
				d_obj += pow(cost[i][s].distance - cost[j][s].distance, 2);
			}
			for (int s = 0; s < dim; ++s) {
				d_var += pow(cost_ds[i][s].distance - cost_ds[j][s].distance, 2);
			}
			// i��j�ɑ΂��鋗�����i�[�����̂�nb�̕��ς����߂�
			arc[j].no = IDs[j];
			arc[j].distance = sqrt(d_obj);
			arc_ds[j].no = IDs[j];
			arc_ds[j].distance = sqrt(d_var);
			arc2[i][j] = sqrt(d_obj);
			arc2_ds[i][j] = sqrt(d_var);
		}
		sort(arc.begin(), arc.end());
		sort(arc_ds.begin(), arc_ds.end());
		double da_obj = 0;
		double da_var = 0;
		for (int j = 0; j < nb; ++j) {
			da_obj += arc[j].distance;
			da_var += arc_ds[j].distance;
		}
		nichea_obj += da_obj / nb;
		nichea_var += da_var / nb;
	}
	nichea_obj = nichea_obj / fk;
	nichea_var = nichea_var / fk;

	for (int i = 0; i < IDs.size(); ++i) {
		double f_ds = 0;
		for (int j = 0; j < IDs.size(); ++j) {
			double sh = max(0.0, 1 - arc2[i][j] / nichea_obj);
			double sh_ds = max(0.0, 1 - arc2_ds[i][j] / nichea_var);
			f_ds += sh + sh_ds;
		}
		// cout << f_ds << endl;
		// dist�ɕۑ�
		dist[i].no = IDs[i];
		dist[i].distance = f_ds;
	}
	
	sort(dist.begin(), dist.end());
	for (int i = 0; i < IDs.size(); ++i) {
		pop.ind[dist[i].no].fitness += i;
	}

	return pop;
}

void DNEA_CDP::curent_fout(int k, int j, int island) {
	string name1 = to_string(k) + "_" + to_string(j) + "gen.dat";
	ofstream fout("./DNEA_CDP/" + to_string(ProblemSubID) + "/" + scalar + "/island" + to_string(island + 1) + "/curents/" + name1, ios::out | ios::app);

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