#include "mazda_mop.h"
namespace bp = benchmark_problem;

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

#include "g01_SUV_FrFL.h"
#include "g02_SUV_ODB.h"
#include "g03_SUV_ODB.h"
#include "g04_SUV_ODB.h"
#include "g05_SUV_SIDE.h"
#include "g06_SUV_SIDE.h"
#include "g07_SUV_SIDE.h"
#include "g08_SUV_REAR.h"
#include "g09_SUV_REAR.h"
#include "g10_SUV_LEV.h"
#include "g11_SUV_LEV.h"
#include "g12_SUV_LEV.h"
#include "g13_SUV_BS.h"
#include "g14_SUV_BS.h"
#include "g19_CDW_FrFL.h"
#include "g20_CDW_ODB.h"
#include "g21_CDW_ODB.h"
#include "g22_CDW_ODB.h"
#include "g23_CDW_SIDE.h"
#include "g24_CDW_SIDE.h"
#include "g25_CDW_SIDE.h"
#include "g26_CDW_REAR.h"
#include "g27_CDW_REAR.h"
#include "g28_CDW_LEV.h"
#include "g29_CDW_LEV.h"
#include "g30_CDW_LEV.h"
#include "g31_CDW_BS.h"
#include "g32_CDW_BS.h"
#include "g37_C5H_FrFL.h"
#include "g38_C5H_ODB.h"
#include "g39_C5H_ODB.h"
#include "g40_C5H_ODB.h"
#include "g41_C5H_SIDE.h"
#include "g42_C5H_SIDE.h"
#include "g43_C5H_SIDE.h"
#include "g44_C5H_REAR.h"
#include "g45_C5H_REAR.h"
#include "g46_C5H_LEV.h"
#include "g47_C5H_LEV.h"
#include "g48_C5H_LEV.h"
#include "g49_C5H_BS.h"
#include "g50_C5H_BS.h"
#include "m01_SUV_Mass.h"
#include "m02_CDW_Mass.h"
#include "m03_C5H_Mass.h"

bp::MazdaMop::MazdaMop()
{
  ncar = 3;

  nobj = 5;
  ncon = 54; // (14 + 4) * ncar = 54;
  nvar = 222;

  g01_SUV_FrFL::initialization();
  g02_SUV_ODB::initialization();
  g03_SUV_ODB::initialization();
  g04_SUV_ODB::initialization();
  g05_SUV_SIDE::initialization();
  g06_SUV_SIDE::initialization();
  g07_SUV_SIDE::initialization();
  g08_SUV_REAR::initialization();
  g09_SUV_REAR::initialization();
  g10_SUV_LEV::initialization();
  g11_SUV_LEV::initialization();
  g12_SUV_LEV::initialization();
  g13_SUV_BS::initialization();
  g14_SUV_BS::initialization();
  g19_CDW_FrFL::initialization();
  g20_CDW_ODB::initialization();
  g21_CDW_ODB::initialization();
  g22_CDW_ODB::initialization();
  g23_CDW_SIDE::initialization();
  g24_CDW_SIDE::initialization();
  g25_CDW_SIDE::initialization();
  g26_CDW_REAR::initialization();
  g27_CDW_REAR::initialization();
  g28_CDW_LEV::initialization();
  g29_CDW_LEV::initialization();
  g30_CDW_LEV::initialization();
  g31_CDW_BS::initialization();
  g32_CDW_BS::initialization();
  g37_C5H_FrFL::initialization();
  g38_C5H_ODB::initialization();
  g39_C5H_ODB::initialization();
  g40_C5H_ODB::initialization();
  g41_C5H_SIDE::initialization();
  g42_C5H_SIDE::initialization();
  g43_C5H_SIDE::initialization();
  g44_C5H_REAR::initialization();
  g45_C5H_REAR::initialization();
  g46_C5H_LEV::initialization();
  g47_C5H_LEV::initialization();
  g48_C5H_LEV::initialization();
  g49_C5H_BS::initialization();
  g50_C5H_BS::initialization();
  m01_SUV_Mass::initialization();
  m02_CDW_Mass::initialization();
  m03_C5H_Mass::initialization();
}


void bp::MazdaMop::evaluate(const vector<long double> &var, vector<long double> &obj, vector<long double> &con)
{
	const int offset = nvar / ncar;

	vector<double> varSUV(&var[0], &var[0] + offset);
	vector<double> varCDW(&var[0] + offset, &var[0] + offset * 2);
	vector<double> varC5H(&var[0] + offset * 2, &var[0] + offset * 3);

	obj[0] = m01_SUV_Mass::evaluate(varSUV) + m02_CDW_Mass::evaluate(varCDW) + m03_C5H_Mass::evaluate(varC5H); // minimize

	obj[1] = -1.0*evaluate_common_parts_num(varSUV, varCDW, varC5H); // maximize

//  obj[2] = m01_SUV_Mass::evaluate(varSUV); // SUV weight
//  obj[3] = m02_CDW_Mass::evaluate(varCDW); // CDW weight
//  obj[4] = m03_C5H_Mass::evaluate(varC5H); // C5H weight

#ifdef MINIMIZE_ALL
  obj[1] = obj[1] != 0 ? 1.0/obj[1] : 100.0; // minimize 1/common_parts_num
#endif
  con[0] = g01_SUV_FrFL::evaluate(varSUV);
  con[1] = g02_SUV_ODB::evaluate(varSUV);
  con[2] = g03_SUV_ODB::evaluate(varSUV);
  con[3] = g04_SUV_ODB::evaluate(varSUV);
  con[4] = g05_SUV_SIDE::evaluate(varSUV);
  con[5] = g06_SUV_SIDE::evaluate(varSUV);
  con[6] = g07_SUV_SIDE::evaluate(varSUV);
  con[7] = g08_SUV_REAR::evaluate(varSUV);
  con[8] = g09_SUV_REAR::evaluate(varSUV);
  con[9] = g10_SUV_LEV::evaluate(varSUV);
  con[10] = g11_SUV_LEV::evaluate(varSUV);
  con[11] = g12_SUV_LEV::evaluate(varSUV);
  con[12] = g13_SUV_BS::evaluate(varSUV);
  con[13] = g14_SUV_BS::evaluate(varSUV);
  con[14] = varSUV[13] - varSUV[12]; // g15
  con[15] = varSUV[15] - varSUV[14]; // g16
  con[16] = varSUV[12] - varSUV[63]; // g17
  con[17] = varSUV[14] - varSUV[63]; // g18
  
  con[18] = g19_CDW_FrFL::evaluate(varCDW);
  con[19] = g20_CDW_ODB::evaluate(varCDW);
  con[20] = g21_CDW_ODB::evaluate(varCDW);
  con[21] = g22_CDW_ODB::evaluate(varCDW);
  con[22] = g23_CDW_SIDE::evaluate(varCDW);
  con[23] = g24_CDW_SIDE::evaluate(varCDW);
  con[24] = g25_CDW_SIDE::evaluate(varCDW);
  con[25] = g26_CDW_REAR::evaluate(varCDW);
  con[26] = g27_CDW_REAR::evaluate(varCDW);
  con[27] = g28_CDW_LEV::evaluate(varCDW);
  con[28] = g29_CDW_LEV::evaluate(varCDW);
  con[29] = g30_CDW_LEV::evaluate(varCDW);
  con[30] = g31_CDW_BS::evaluate(varCDW);
  con[31] = g32_CDW_BS::evaluate(varCDW);
  con[32] = varCDW[13] - varCDW[12]; // g33
  con[33] = varCDW[15] - varCDW[14]; // g34
  con[34] = varCDW[12] - varCDW[63]; // g35
  con[35] = varCDW[14] - varCDW[63]; // g36
  
  con[36] = g37_C5H_FrFL::evaluate(varC5H);
  con[37] = g38_C5H_ODB::evaluate(varC5H);
  con[38] = g39_C5H_ODB::evaluate(varC5H);
  con[39] = g40_C5H_ODB::evaluate(varC5H);
  con[40] = g41_C5H_SIDE::evaluate(varC5H);
  con[41] = g42_C5H_SIDE::evaluate(varC5H);
  con[42] = g43_C5H_SIDE::evaluate(varC5H);
  con[43] = g44_C5H_REAR::evaluate(varC5H);
  con[44] = g45_C5H_REAR::evaluate(varC5H);
  con[45] = g46_C5H_LEV::evaluate(varC5H);
  con[46] = g47_C5H_LEV::evaluate(varC5H);
  con[47] = g48_C5H_LEV::evaluate(varC5H);
  con[48] = g49_C5H_BS::evaluate(varC5H);
  con[49] = g50_C5H_BS::evaluate(varC5H);
  con[50] = varC5H[13] - varC5H[12]; // g51
  con[51] = varC5H[15] - varC5H[14]; // g52
  con[52] = varC5H[12] - varC5H[63]; // g53
  con[53] = varC5H[14] - varC5H[63]; // g54
}


double bp::MazdaMop::evaluate_common_parts_num(vector<double> &v1, vector<double> &v2, vector<double> &v3)
{
  const double eps = 0.05;

  double commonPartsNum = 0;

  const size_t N = v1.size();
  for (size_t i = 0; i<N; ++i) {
    double ub = max(v1[i], max(v2[i], v3[i]));
    double lb = min(v1[i], min(v2[i], v3[i]));
    if (ub - lb < eps) ++commonPartsNum; // ith parts is the common part.
  }

  return commonPartsNum;
}



void bp::MazdaMop::evaluate_con(std::vector<long double> &var, std::vector<long double> &con, std::vector<int> &whitch_con) {
	const int offset = nvar / ncar;

	vector<double> varSUV(&var[0], &var[0] + offset);
	vector<double> varCDW(&var[0] + offset, &var[0] + offset * 2);
	vector<double> varC5H(&var[0] + offset * 2, &var[0] + offset * 3);

	for (int i = 0; i < whitch_con.size(); ++i) {
		if (whitch_con[i] == 1) {
			switch (i) {

			case 0:
				con[0] = g01_SUV_FrFL::evaluate(varSUV);
				break;

			case 1:
				con[1] = g02_SUV_ODB::evaluate(varSUV);
				break;

			case 2:
				con[2] = g03_SUV_ODB::evaluate(varSUV);
				break;

			case 3:
				con[3] = g04_SUV_ODB::evaluate(varSUV);
				break;

			case 4:
				con[4] = g05_SUV_SIDE::evaluate(varSUV);
				break;

			case 5:
				con[5] = g06_SUV_SIDE::evaluate(varSUV);
				break;

			case 6:
				con[6] = g07_SUV_SIDE::evaluate(varSUV);
				break;

			case 7:
				con[7] = g08_SUV_REAR::evaluate(varSUV);
				break;

			case 8:
				con[8] = g09_SUV_REAR::evaluate(varSUV);
				break;

			case 9:
				con[9] = g10_SUV_LEV::evaluate(varSUV);
				break;

			case 10:
				con[10] = g11_SUV_LEV::evaluate(varSUV);
				break;

			case 11:
				con[11] = g12_SUV_LEV::evaluate(varSUV);
				break;

			case 12:
				con[12] = g13_SUV_BS::evaluate(varSUV);
				break;

			case 13:
				con[13] = g14_SUV_BS::evaluate(varSUV);
				break;

			case 14:
				con[14] = varSUV[13] - varSUV[12]; // g15
				break;

			case 15:
				con[15] = varSUV[15] - varSUV[14]; // g16
				break;

			case 16:
				con[16] = varSUV[12] - varSUV[63]; // g17
				break;

			case 17:
				con[17] = varSUV[14] - varSUV[63]; // g18
				break;

			case 18:
				con[18] = g19_CDW_FrFL::evaluate(varCDW);
				break;

			case 19:
				con[19] = g20_CDW_ODB::evaluate(varCDW);
				break;

			case 20:
				con[20] = g21_CDW_ODB::evaluate(varCDW);
				break;

			case 21:
				con[21] = g22_CDW_ODB::evaluate(varCDW);
				break;

			case 22:
				con[22] = g23_CDW_SIDE::evaluate(varCDW);
				break;

			case 23:
				con[23] = g24_CDW_SIDE::evaluate(varCDW);
				break;

			case 24:
				con[24] = g25_CDW_SIDE::evaluate(varCDW);
				break;

			case 25:
				con[25] = g26_CDW_REAR::evaluate(varCDW);
				break;

			case 26:
				con[26] = g27_CDW_REAR::evaluate(varCDW);
				break;

			case 27:
				con[27] = g28_CDW_LEV::evaluate(varCDW);
				break;

			case 28:
				con[28] = g29_CDW_LEV::evaluate(varCDW);
				break;

			case 29:
				con[29] = g30_CDW_LEV::evaluate(varCDW);
				break;

			case 30:
				con[30] = g31_CDW_BS::evaluate(varCDW);
				break;

			case 31:
				con[31] = g32_CDW_BS::evaluate(varCDW);
				break;

			case 32:
				con[32] = varCDW[13] - varCDW[12]; // g33
				break;

			case 33:
				con[33] = varCDW[15] - varCDW[14]; // g34
				break;

			case 34:
				con[34] = varCDW[12] - varCDW[63]; // g35
				break;

			case 35:
				con[35] = varCDW[14] - varCDW[63]; // g36
				break;

			case 36:
				con[36] = g37_C5H_FrFL::evaluate(varC5H);
				break;

			case 37:
				con[37] = g38_C5H_ODB::evaluate(varC5H);
				break;

			case 38:
				con[38] = g39_C5H_ODB::evaluate(varC5H);
				break;

			case 39:
				con[39] = g40_C5H_ODB::evaluate(varC5H);
				break;

			case 40:
				con[40] = g41_C5H_SIDE::evaluate(varC5H);
				break;

			case 41:
				con[41] = g42_C5H_SIDE::evaluate(varC5H);
				break;

			case 42:
				con[42] = g43_C5H_SIDE::evaluate(varC5H);
				break;

			case 43:
				con[43] = g44_C5H_REAR::evaluate(varC5H);
				break;

			case 44:
				con[44] = g45_C5H_REAR::evaluate(varC5H);
				break;

			case 45:
				con[45] = g46_C5H_LEV::evaluate(varC5H);
				break;

			case 46:
				con[46] = g47_C5H_LEV::evaluate(varC5H);
				break;

			case 47:
				con[47] = g48_C5H_LEV::evaluate(varC5H);
				break;

			case 48:
				con[48] = g49_C5H_BS::evaluate(varC5H);
				break;

			case 49:
				con[49] = g50_C5H_BS::evaluate(varC5H);
				break;

			case 50:
				con[50] = varC5H[13] - varC5H[12]; // g51
				break;

			case 51:
				con[51] = varC5H[15] - varC5H[14]; // g52
				break;

			case 52:
				con[52] = varC5H[12] - varC5H[63]; // g53
				break;

			case 53:
				con[53] = varC5H[14] - varC5H[63]; // g54
				break;

			default:
				break;
			}

		}

	}

}