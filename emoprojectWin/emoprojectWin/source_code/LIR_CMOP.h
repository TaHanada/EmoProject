#pragma once
#include"population.h"
#include "individual.h"

class LIR_CMOP
{
private:
	int num_obj;
	int num_cons;
	int dim;
	int problemID;


public:


	void initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim);
	void range_setting(int swit, vector<vector<double>>& range);

	individual problems(int swit, individual ind);


	individual LIR_CMOP1(individual ind);


	individual LIR_CMOP2(individual ind);
	

	individual LIR_CMOP3(individual ind);


	individual LIR_CMOP4(individual ind);


	individual LIR_CMOP5(individual ind);


	individual LIR_CMOP6(individual ind);


	individual LIR_CMOP7(individual ind);


	individual LIR_CMOP8(individual ind);


	individual LIR_CMOP9(individual ind);


	individual LIR_CMOP10(individual ind);


	individual LIR_CMOP11(individual ind);


	individual LIR_CMOP12(individual ind);


	individual LIR_CMOP13(individual ind);


	individual LIR_CMOP14(individual ind);

	void fileout(int swit, int gen, int reps, population pop);

	void constrained_landscape();

	individual CLIR_constraint(individual ind);
};