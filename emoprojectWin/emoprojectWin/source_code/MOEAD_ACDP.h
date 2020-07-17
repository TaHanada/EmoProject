#ifndef _MOEAD_ACDP_
#define _MOEAD_ACDP_
#include<fstream>
#include"population.h"
#include"MT.h"
#include"weight_vec.h"
#include"Problems.h"
#include<math.h>


class MOEAD_ACDP {
private:
	int pop_m;//population size
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

	double theta_ = 0;
	double theta0_ = 0;
	double alpha_ = 0.8;

	population curent;
	population offspring;

	weight_vec weight;
	int H1, H2;
	int num_vec;
	vector<vector<double>> range;
	Problems problems;
	int ProblemSubID;

public:

	//constractor
	MOEAD_ACDP(int tmp_rep, int ProblemID, int tmp_ProblemSubID) {
		Parameter_setting("MOEAD_CDP", ProblemID, tmp_ProblemSubID);
		rep = tmp_rep;
		ProblemSubID = tmp_ProblemSubID;


		//set populations
		curent.initialize(pop_m, num_obj, dim, num_const);
		offspring.initialize(pop_m, num_obj, dim, num_const);
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


	void evolution_one(int gen);//evolution


	void run_algorithm(int reps);//run algorithm


	void file_initialization(int k, int gen);//file initialization


	void file_allsolutions(int k, int gen, string scalar, population pop);//file allsolutions


	int calc_cossim(individual offspring1);


	individual SBX(individual ind1, individual ind2);


	individual PM(individual ind1);


	double get_angle(individual ind1, individual ind2);


	void update_theta(int curent_gen);


	void curent_fout(int, int, int);//output objectives
};




#endif _MOEAD_ACDP_// !_MOEAD_
