#ifndef _NSGAII_IEpsilon_
#define _NSGAII_IEpsilon_
#include<fstream>
#include"population.h"
#include"MT.h"
#include"weight_vec.h"
#include"Problems.h"


class NSGAII_IEpsilon {
private:
	int pop_m;//population size
	int gen;//generation count
	int num_obj;//number of objects
	int num_const;//number of constraints
	int dim;//dimention
	double crossover_probability;
	double mutation_probability;
	int rep;
	string scalar;

	population curent;
	population offspring;
	population merge;

	//IEpsilon parameter
	double epsilon_;
	double epsilon0_;
	double epsilonMax_;

	int G_;
	int Tc_;
	double cp_;
	double theta_;
	double alpha_;
	double tau_;

	vector<vector<double>> range;
	Problems problems;
	int ProblemID;
	int ProblemSubID;

public:

	//constractor
	NSGAII_IEpsilon(int tmp_rep, int tmp_ProblemID, int tmp_ProblemSubID) {
		Parameter_setting("MOEAD_CDP", tmp_ProblemID, tmp_ProblemSubID);
		rep = tmp_rep;
		ProblemID = tmp_ProblemID;
		ProblemSubID = tmp_ProblemSubID;
		scalar = "tch";//good luck charm

		//set populations
		curent.initialize(pop_m, num_obj, dim, num_const);
		offspring.initialize(pop_m, num_obj, dim, num_const);
		merge.initialize(2*pop_m, num_obj, dim, num_const);
	};


	//	void range_setting(int);


	void Parameter_setting(string, int, int);


	void cout_population(int);//console output population


	void cout_property(int);//console output property of solutions


	void initialize_population();//initialize population


	void total_constraints(individual&);


	void file_objectives(int, int);//output objectives


	void curent_fout(int, int,int);//output objectives


	void evolution(int gen);//evolution


	void run_algorithm(int reps);//run algorithm


	void file_initialization(int k, int gen);//file initialization


	void file_allsolutions(int k, int gen, population pop);//file allsolutions


	individual SBX(individual ind1, individual ind2);


	individual PM(individual ind1);


	void evaluation(population pop);


	int tournament(int tournamentSize, population pop);


	int check_dominance(individual a, individual b);


	void initial_epsilon();


	void update_epsilon();


	population assign_CrowdingDistance(vector<int> ID, population pop);
};




#endif _NSGAII_IEpsilon_// !_NSGAII_IEpsilon_
