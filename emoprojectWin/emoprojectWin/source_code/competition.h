#pragma once
#define _USE_MATH_DEFINES
#include "individual.h"

class competition
{
private:
	int num_obj;
	int num_cons;
	int dim;

public:
	void initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim);
	void range_setting(int swit, vector<vector<double>>& range);
	individual problems(int swit, individual ind);

	individual singleobjective(individual ind);
	individual multiobjective(individual ind);

};