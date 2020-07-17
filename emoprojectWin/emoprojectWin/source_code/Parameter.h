#ifndef _Parameter_
#define _Parameter_

int pop_m = 100;
int dim = 4;
int num_const = 4;
int num_obj = 2;
int num_vec = 300;
int gen = 49;
int reps = 21;
int seed = 1;
double crossover_probability = 0.9;
double mutation_probability = 1.0 / (double)dim;
int pop_c = 300 / num_vec;
int pop_a = 210;
double niche_obj = 1.0;
double niche_var = 1.0;

void Parameter_setting() {



}

#endif _Parameter_