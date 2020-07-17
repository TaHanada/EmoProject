#pragma once
#define _USE_MATH_DEFINES
#include "individual.h"
#include <cmath>


class WB {
private:
	int num_obj;


public:
	//constractor
	void initialize(int tmp_num_obj) {
		num_obj = tmp_num_obj;
	}


	individual WeldedBeam(individual ind);


	void check();


	individual problems(int swit, individual ind);

};