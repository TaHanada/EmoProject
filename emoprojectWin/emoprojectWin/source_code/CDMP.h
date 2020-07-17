#define _USE_MATH_DEFINES
#include "individual.h"
#include"population.h"
#include"MT.h"
#include <cmath>
#include<fstream>
#include<string>
#include<sstream>


class CDMP {
private:
	int problemID;
	int num_obj;
	int num_cons;
	int dim;
	double r;
	vector<double> center;
	vector<vector<double>> a, v; 
	vector<vector<double>> mixed_points;
	
public:
	void check(int);
	void initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim);

	individual Mixed(individual ind);
	individual Mixed_noexception(individual ind);
	individual Checker(individual ind);
	individual Vertices(individual ind);
	individual Center(individual ind);
	individual Moat(individual ind);



	individual Objectives(individual ind);
	individual Constraints(int swit, individual ind);

	void fileout(int swit, int gen, int reps, population pop);


	int is_ParetoSolution(individual ind);
};