#ifndef _individual_
#define _individual_

#include<iostream>
#include<stdio.h>
#include<vector>
#include<algorithm>

using namespace std;

class individual {
public:
	vector<long double> var;//phenotype
	vector<long double> cost;//�e�ړI�֐���cost
	int weight_no;//MOEAD�̏d�݃x�N�^
	double fitness;//scalar fitness
	vector<long double> constraint;
	double total_c;
	int feasible;
	int rank;
	double dist;

	individual &operator=(const individual &a) {
		var = a.var;
		cost = a.cost;
		weight_no = a.weight_no;
		fitness = a.fitness;
		constraint = a.constraint;
		total_c = a.total_c;
		feasible = a.feasible;
		rank = a.rank;
		dist = a.dist;

		return *this;
	};//�R�s�[���Z�q
};

#endif _individual_