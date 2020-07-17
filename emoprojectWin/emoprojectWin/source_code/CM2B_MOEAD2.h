#ifndef _MOEAD2_
#define _MOEAD2_
#include<fstream>
#include"population.h"
#include"MT.h"
#include"weight_vec.h"
#include"Problems.h"


class CM2B_MOEAD2 {
private:
	int pop_m;//population size
	int pop_c;//number of offsprings
	int gen;//generation count
	int num_obj;//number of objects
	int num_const;//number of constraints
	int dim;//dimention
	double crossover_probability;
	double mutation_probability;
	string scalar;
	int rep;
	vector<double> reference_point;//éQè∆ì_
	vector<double> nadir_point;//ç≈à´ì_
	vector<double> constraint_max;

	vector<population> curent;
	vector<population> offspring;

	weight_vec weight;
	int H1, H2;
	int num_vec;
	vector<vector<double>> range;
	Problems problems;
	int ProblemSubID;

public:

	//constractor
	CM2B_MOEAD2(int tmp_rep, int ProblemID, int tmp_ProblemSubID) {
		Parameter_setting("CM2B_MOEAD2", ProblemID, tmp_ProblemSubID);
		rep = tmp_rep;
		ProblemSubID = tmp_ProblemSubID;

		//set populations
		curent.resize(num_vec);
		offspring.resize(num_vec);
		for (int i = 0; i < num_vec; ++i) {
			curent[i].initialize(pop_m, num_obj, dim, num_const);
			offspring[i].initialize(pop_c, num_obj, dim, num_const);
		}
	};


//	void range_setting(int);


	void Parameter_setting(string, int, int);


	void cout_population(int);//console output population


	void cout_property(int);//console output property of solutions


	void initialize_reference_nadir_point();//initialize ref. and nadir points


	void initialize_population();//initialize population


	double weighted_sum(individual ind);//weighed sum


	double tchebycheff_notabsolute(individual ind);//Tchebycheff no absolution


	double tchebycheff(individual ind);//Tchebycheff


	double PBI(individual ind);//PBI


	double IPBI(individual ind);//IPBI


	double AOF(individual ind);//IPBI


	double scalarizing_function(individual ind);

	void total_constraints(individual&);


	void file_objectives(int, int );//output objectives


	void evolution_mu(int gen);//evolution


	void run_algorithm(int reps);//run algorithm


	void file_initialization(int k, int gen);//file initialization


	void file_allsolutions(int k, int gen, string scalar, vector<population> pop, vector<double> rf);//file allsolutions


	void file_curents(int k, int gen, string scalar, int island);//file curents


	int calc_cossim(individual offspring1);


	int tournament_totalC(int A, population pop);


	int tournament_PRandCV(int A, population pop);


	int tournament(int A, population pop);


	int tournament_fitness(int A, population pop);


	int tournament_inverse(int A, population pop);


	individual SBX(int vec1, int vec2, int n1, int n2);


	individual DE(int vec1, int vec2, int vec3, int vec4, int n1, int n2, int n3, int n4);


	void PM(int vec, int pp);


	void update_NSGAII_const(vector<individual> update_population, int vec);


	void update_NSGAII(vector<individual> update_population, int vec);


	int get_minID(population pop);


	void parentselection(int vec, int& i1, int& i2, int& n1, int& n2);

	void update_constraint_max(vector<population> pop);
};




#endif _MOEAD2_// !_MOEAD2_
