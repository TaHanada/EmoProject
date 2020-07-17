#pragma once
#include "individual.h"
class CARSIDE
{
private:
	int num_obj;
	int num_cons;
	int dim;

public:
	void initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim);
	void range_setting(int swit, vector<vector<double>>& range);
	individual problems(int swit, individual ind);

	individual carsideproblem(individual ind);
};
