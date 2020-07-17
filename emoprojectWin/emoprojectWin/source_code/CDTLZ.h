// cec09.h
//
// C++ source codes of the test instances for CEC 2009 MOO competition
//
// If the source codes are not consistant with the report, please use the version in the report
//
// History
// 	  v1   Sept.8  2008
// 	  v1.1 Sept.22 2008: add R2_DTLZ2_M5 R3_DTLZ3_M5 WFG1_M5
//    v1.2 Oct.2  2008: fix the bugs in CF1-CF4, CF6-CF10, thank Qu Boyang for finding the bugs
//    v1.3 Oct.8  2008: fix a bug in R2_DTLZ2_M5, thank Santosh Tiwari
//    v1.4 Nov.26 2008: fix the bugs in CF4-CF5, thank Xueqiang Li for finding the bugs
//    v1.5 Dec.2  2008: fix the bugs in CF5 and CF7, thank Santosh Tiwari for finding the bugs
#pragma once
#include "individual.h"
#include"population.h"


class CDTLZ
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


	individual DTLZ1(individual ind);
	individual DTLZ2(individual ind);
	individual DTLZ3(individual ind);
	individual DTLZ4(individual ind);

	individual C1_DTLZ1(individual);
	individual C1_DTLZ3(individual);
	individual C2_DTLZ2(individual);
	individual C3_DTLZ1(individual);
	individual C3_DTLZ4(individual);
	individual C1C3_DTLZ3(individual);

	void fileout(int swit, int gen, int reps, population pop);
};