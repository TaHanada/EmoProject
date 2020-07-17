#pragma once
#define _USE_MATH_DEFINES
#include "individual.h"
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<stdio.h>
#include <iomanip>
#include"WeldedBeam.h"
#include "mazda_mop.h"
#include"CDMP.h"
#include"DOC.h"
#include"cec09.h"
#include"competition.h"
#include"water.h"
#include"discbrakedesign.h"
#include"twobartrussdesign.h"
#include"speedreducerdesign.h"
#include"carsideimpact.h"
#include"weldedbeam2.h"

using namespace std;
namespace bp = benchmark_problem;


class Problems {
public:
	int ProblemID, ProblemSubID;
	WB WBprob;
	CDMP CDMPprob;
	DOC DOCprob;
	CEC09 CEC09prob;
	bp::MazdaMop mop;
	competition compeprob;
	WATER waterprob;
	DISCBRAKE discprob;
	TWOBARTRUSS twobartrussprob;
	SPEEDREDUCER speedprob;
	CARSIDE carsideprob;
	WELDEDBEAM weldedbeamprob;



	void set_ProblemID(int tmp_ProblemID, int tmp_ProblemSubID) {
		ProblemID = tmp_ProblemID;
		ProblemSubID = tmp_ProblemSubID;
	}


	void parameter_setting(int& num_obj, int& num_const, int& dim, vector<vector<double>>& range) {
		switch (ProblemID) {
		case 1://MAZDA
			num_obj = 2;
			num_const = 55;
			dim = 222;
			range_setting(ProblemID, range, dim);
			break;

		case 2://Welded Beam
			num_obj = 2;
			num_const = 4;
			dim = 4;
			range_setting(ProblemID, range, dim);
			break;

		case 3://CDMP
			CDMPprob.initialize(ProblemSubID, num_obj, num_const, dim);
			range_setting(ProblemID, range, dim);
			break;

		case 4:
			DOCprob.initialize(ProblemSubID, num_obj, num_const, dim);
			range_setting(ProblemID, range, dim);
			break;

		case 5:
			CEC09prob.initialize(ProblemSubID, num_obj, num_const, dim);
			range_setting(ProblemID, range, dim);
			break;

		case 6://2019_competition
			compeprob.initialize(ProblemSubID, num_obj, num_const, dim);
			range_setting(ProblemID, range, dim);
			break;

		case 7://water_problem
			waterprob.initialize(ProblemSubID, num_obj, num_const, dim);
			range_setting(ProblemID, range, dim);
			break;

		case 8://disc_brake_design_problem
			discprob.initialize(ProblemSubID, num_obj, num_const, dim);
			range_setting(ProblemID, range, dim);
			break;

		case 9://two_bar_truss_design_problem
			twobartrussprob.initialize(ProblemSubID, num_obj, num_const, dim);
			range_setting(ProblemID, range, dim);
			break;

		case 10://speedreducer_design_problem
			speedprob.initialize(ProblemSubID, num_obj, num_const, dim);
			range_setting(ProblemID, range, dim);
			break;

		case 11://car_side_impact_design_problem
			carsideprob.initialize(ProblemSubID, num_obj, num_const, dim);
			range_setting(ProblemID, range, dim);
			break;

		case 12://welded_beam_problem(RWLPs)
			weldedbeamprob.initialize(ProblemSubID, num_obj, num_const, dim);
			range_setting(ProblemID, range, dim);
			break;

		default:
			break;
		}


	}


	void range_setting(int ProblemID, vector<vector<double>>& range, int dim) {
		switch (ProblemID)
		{
		case 1:
		{
			//input lower and upper
			ifstream fin("decisionspace.txt");
			if (fin.fail()) {
				cout << "入力ファイルをオープンできません" << endl;
			}

			range.resize(dim);
			int input_counter = 0;
			for (string line_in; getline(fin, line_in);) {
				if (line_in.size() == 0) continue;

				string token;
				istringstream ss(line_in);
				while (getline(ss, token, '\t')) {
					range[input_counter].push_back(atof(token.c_str()));

				}
				input_counter++;
			}
		}
		break;

		case 2:
			range.resize(dim);
			for (int i = 0; i < dim; ++i) {
				range[i].resize(2);
			}
			range[0][0] = 0.125, range[0][1] = 5.0;
			range[1][0] = 0.125, range[1][1] = 5.0;
			range[2][0] = 0.1, range[2][1] = 10.0;
			range[3][0] = 0.1, range[3][1] = 10.0;
			break;

		case 3:
			range.resize(dim);
			for (int i = 0; i < dim; ++i) {
				range[i].push_back(-50);
				range[i].push_back(50);
			}
			break;

		case 4:
			DOCprob.range_setting(ProblemSubID, range);
			break;

		case 5:
			CEC09prob.range_setting(ProblemSubID, range);
			break;

		case 6:
			compeprob.range_setting(ProblemSubID, range);
			break;

		case 7:
			waterprob.range_setting(ProblemSubID, range);
			break;

		case 8:
			discprob.range_setting(ProblemSubID, range);
			break;

		case 9:
			twobartrussprob.range_setting(ProblemSubID, range);
			break;

		case 10:
			speedprob.range_setting(ProblemSubID, range);
			break;

		case 11:
			carsideprob.range_setting(ProblemSubID, range);
			break;
		
		case 12:
			weldedbeamprob.range_setting(ProblemSubID, range);
			break;

		default:
			break;
		}


	}



	individual evaluation(individual ind) {
		switch (ProblemID) {
		case 1:
			mop.evaluate(ind.var, ind.cost, ind.constraint);
			break;

		case 2:
			ind = WBprob.WeldedBeam(ind);
			break;

		case 3:
			ind = CDMPprob.Objectives(ind);
			ind = CDMPprob.Constraints(ProblemSubID, ind);
			break;

		case 4:
			ind = DOCprob.problems(ProblemSubID, ind);
			break;

		case 5:
			ind = CEC09prob.problems(ProblemSubID, ind);
			break;

		case 6:
			ind = compeprob.problems(ProblemSubID, ind);
			break;

		case 7:
			ind = waterprob.problems(ProblemSubID, ind);
			break;

		case 8:
			ind = discprob.problems(ProblemSubID, ind);
			break;

		case 9:
			ind = twobartrussprob.problems(ProblemSubID, ind);
			break;

		case 10:
			ind = speedprob.problems(ProblemSubID, ind);
			break;

		case 11:
			ind = carsideprob.problems(ProblemSubID, ind);
			break;

		case 12:
			ind = weldedbeamprob.problems(ProblemSubID, ind);
			break;

		default:
			break;
		}

		return ind;
	}
};