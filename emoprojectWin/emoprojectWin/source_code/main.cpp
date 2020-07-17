#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include "divide.h"
#include"Parameter.h"
#include"weight_vec.h"
#include"CM2B_MOEAD.h"
#include"MOEAD_CDP.h"
#include"NSGAII_CDP.h"
#include"MOEAD_IEpsilon.h"
#include"NSGAII_SP.h"
#include"MOEAD_ACDP.h"
#include"CM2B_MOEAD2.h"
#include"CM2B_MOEAD3.h"
#include"DNEA_CDP.h"
#include"NSGAII_IEpsilon.h"
#include <string.h>


using namespace std;

int main(int argc, char *argv[]) {
	//argv[1] : problem_name
	//argv[2] : algorithm_name
	//argv[3] : seed and file_divide

	string str_seed = argv[3];
	int change_seed = stoi(str_seed);
	seed = change_seed;
	Divide* obj = Divide::getInstance();
	obj->change_file_divide(seed);

	cout << "start!" << endl;

	//parameter///////////////////////////////////////////////////////////////////////////////
	//problem_setting
	vector<int> ProblemIDs;
	vector<int> ProblemSubIDs;
	if (strcmp(argv[1],"MAZDA") == 0) {
		ProblemIDs.push_back(1);
	}
	else if (strcmp(argv[1],"Welded Beam") == 0) {
		ProblemIDs.push_back(2);
	}
	else if (strcmp(argv[1],"CDMP") == 0) {
		ProblemIDs.push_back(3);
	}
	else if (strcmp(argv[1],"DOC") == 0) {
		ProblemIDs.push_back(4);
	}
	else if (strcmp(argv[1],"CTP") == 0) {
		ProblemIDs.push_back(5);
		ProblemSubIDs.push_back(1);
		ProblemSubIDs.push_back(2);
		ProblemSubIDs.push_back(3);
		ProblemSubIDs.push_back(4);
		ProblemSubIDs.push_back(5);
		ProblemSubIDs.push_back(6);
		ProblemSubIDs.push_back(7);
		ProblemSubIDs.push_back(8);
		ProblemSubIDs.push_back(9);
	}
	else if (strcmp(argv[1], "competition") == 0) {
		ProblemIDs.push_back(6);
		ProblemSubIDs.push_back(2);
	}
	else if (strcmp(argv[1], "water") == 0) {
		ProblemIDs.push_back(7);
		ProblemSubIDs.push_back(1);
	}
	else if (strcmp(argv[1], "disc") == 0) {
		ProblemIDs.push_back(8);
		ProblemSubIDs.push_back(1);
		ProblemSubIDs.push_back(2);
	}
	else if (strcmp(argv[1], "twobartruss") == 0) {
		ProblemIDs.push_back(9);
		ProblemSubIDs.push_back(1);
	}
	else if (strcmp(argv[1], "speedreducer") == 0) {
		ProblemIDs.push_back(10);
		ProblemSubIDs.push_back(1);
	}
	else if (strcmp(argv[1], "carside") == 0) {
		ProblemIDs.push_back(11);
		ProblemSubIDs.push_back(1);
	}
	else if (strcmp(argv[1], "weldedbeam") == 0) {
		ProblemIDs.push_back(12);
		ProblemSubIDs.push_back(1);
	}

	//algorithm_setting
	if (strcmp(argv[2],"NSGAII_CDP") == 0) {
		for (int j = 0; j < ProblemIDs.size(); ++j) {
			for (int k = 0; k < ProblemSubIDs.size(); ++k) {
				for (int i = 0; i < reps; ++i) {
					NSGAII_CDP data(seed+i, ProblemIDs[j], ProblemSubIDs[k]);//ƒAƒ‹ƒSƒŠƒYƒ€‚ÌéŒ¾
					data.run_algorithm(seed+i);
				}
			}
		}
	}
	else if (strcmp(argv[2],"NSGAII_SP") == 0) {
		for (int j = 0; j < ProblemIDs.size(); ++j) {
			for (int k = 0; k < ProblemSubIDs.size(); ++k) {
				for (int i = 0; i < reps; ++i) {
					NSGAII_SP data(seed+i, ProblemIDs[j], ProblemSubIDs[k]);//ƒAƒ‹ƒSƒŠƒYƒ€‚ÌéŒ¾
					data.run_algorithm(seed+i);
				}
			}
		}
	}
	else if (strcmp(argv[2],"MOEAD_CDP") == 0) {
		for (int j = 0; j < ProblemIDs.size(); ++j) {
			for (int k = 0; k < ProblemSubIDs.size(); ++k) {
				for (int i = 0; i < reps; ++i) {
					MOEAD_CDP data(seed+i, ProblemIDs[j], ProblemSubIDs[k]);//ƒAƒ‹ƒSƒŠƒYƒ€‚ÌéŒ¾
					data.run_algorithm(seed+i);
				}
			}
		}
	}
	else if (strcmp(argv[2],"MOEAD_ACDP") == 0) {
		for (int j = 0; j < ProblemIDs.size(); ++j) {
			for (int k = 0; k < ProblemSubIDs.size(); ++k) {
				for (int i = 0; i < reps; ++i) {
					MOEAD_ACDP data(seed+i, ProblemIDs[j], ProblemSubIDs[k]);//ƒAƒ‹ƒSƒŠƒYƒ€‚ÌéŒ¾
					data.run_algorithm(seed+i);
				}
			}
		}
	}
	else if (strcmp(argv[2],"MOEAD_IEpsilon") == 0) {
		for (int j = 0; j < ProblemIDs.size(); ++j) {
			for (int k = 0; k < ProblemSubIDs.size(); ++k) {
				for (int i = 0; i < reps; ++i) {
					MOEAD_IEpsilon data(seed+i, ProblemIDs[j], ProblemSubIDs[k]);//ƒAƒ‹ƒSƒŠƒYƒ€‚ÌéŒ¾
					data.run_algorithm(seed+i);
				}
			}
		}
	}
	else if (strcmp(argv[2],"CM2B_MOEAD") == 0) {
		for (int j = 0; j < ProblemIDs.size(); ++j) {
			for (int k = 0; k < ProblemSubIDs.size(); ++k) {
				for (int i = 0; i < reps; ++i) {
					CM2B_MOEAD data(seed+i, ProblemIDs[j], ProblemSubIDs[k]);//ƒAƒ‹ƒSƒŠƒYƒ€‚ÌéŒ¾
					data.run_algorithm(seed+i);
				}
			}
		}
	}
	else if (strcmp(argv[2], "CM2B_MOEAD2") == 0) {
		for (int j = 0; j < ProblemIDs.size(); ++j) {
			for (int k = 0; k < ProblemSubIDs.size(); ++k) {
				for (int i = 0; i < reps; ++i) {
					CM2B_MOEAD2 data(seed + i, ProblemIDs[j], ProblemSubIDs[k]);//ƒAƒ‹ƒSƒŠƒYƒ€‚ÌéŒ¾
					data.run_algorithm(seed + i);
				}
			}
		}
	}
	else if (strcmp(argv[2], "CM2B_MOEAD3") == 0) {
		for (int j = 0; j < ProblemIDs.size(); ++j) {
			for (int k = 0; k < ProblemSubIDs.size(); ++k) {
				for (int i = 0; i < reps; ++i) {
					CM2B_MOEAD3 data(seed + i, ProblemIDs[j], ProblemSubIDs[k]);//ƒAƒ‹ƒSƒŠƒYƒ€‚ÌéŒ¾
					data.run_algorithm(seed + i);
				}
			}
		}
	}
	else if (strcmp(argv[2], "DNEA_CDP") == 0) {
		for (int j = 0; j < ProblemIDs.size(); ++j) {
			for (int k = 0; k < ProblemSubIDs.size(); ++k) {
				for (int i = 0; i < reps; ++i) {
					DNEA_CDP data(seed + i, ProblemIDs[j], ProblemSubIDs[k]);//ƒAƒ‹ƒSƒŠƒYƒ€‚ÌéŒ¾
					data.run_algorithm(seed + i);
				}
			}
		}
	}
	else if (strcmp(argv[2], "NSGAII_IEpsilon") == 0) {
		for (int j = 0; j < ProblemIDs.size(); ++j) {
			for (int k = 0; k < ProblemSubIDs.size(); ++k) {
				for (int i = 0; i < reps; ++i) {
					NSGAII_IEpsilon data(seed + i, ProblemIDs[j], ProblemSubIDs[k]);//ƒAƒ‹ƒSƒŠƒYƒ€‚ÌéŒ¾
					data.run_algorithm(seed + i);
				}
			}
		}
	}


	int finish;
	cout << "finish!" << endl;
	cin >> finish;
	return 0;
}
