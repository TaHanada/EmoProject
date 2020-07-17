#pragma once
#define _USE_MATH_DEFINES
#include "individual.h"
#include <cmath>
#include<fstream>
#include<string>
#include<sstream>


class DOC {
private:
	int num_obj;
	int num_cons;
	int dim;


public:
	//constractor
	void initialize(int swit, int& tmp_num_obj, int& tmp_num_cons, int& tmp_dim);
	void range_setting(int swit, vector<vector<double>>& range);

	individual DOC1(individual ind);
	individual DOC2(individual ind);
	individual DOC3(individual ind);
	individual DOC4(individual ind);
	individual DOC5(individual ind);
	individual DOC6(individual ind);
	individual DOC7(individual ind);
	individual DOC8(individual ind);
	individual DOC9(individual ind);


	individual problems(int swit, individual ind);

};