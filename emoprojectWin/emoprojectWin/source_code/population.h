#ifndef _population_
#define _population_
#include<vector>
#include<iostream>
#include"individual.h"
#include"MT.h"

using namespace std;

class population {
public:
	int size;//size of population
	vector<individual> ind;//population

	void initialize(int pop_m,int num_obj, int num_items, int num_const) {
		size = pop_m;
		ind.resize(size);
		for (int i = 0; i < size; ++i) {
			ind[i].var.resize(num_items);//size of phenotype
			ind[i].cost.resize(num_obj);//number of tasks
			ind[i].constraint.resize(num_const);//number of constraint
			ind[i].feasible = 0;
			ind[i].dist = 0;
		}
	}//������

	void cout_property();//�̌Q�̏o��

	void cout_population();//�̂̏o��

	void erase();//�S�̂̍폜

	void operator=(const population &a);//�R�s�[���Z�q

	void operator+=(const population &a);//�ǉ����Z�q
};

#endif _population_