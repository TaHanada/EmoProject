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
	}//初期化

	void cout_property();//個体群の出力

	void cout_population();//個体の出力

	void erase();//全個体の削除

	void operator=(const population &a);//コピー演算子

	void operator+=(const population &a);//追加演算子
};

#endif _population_