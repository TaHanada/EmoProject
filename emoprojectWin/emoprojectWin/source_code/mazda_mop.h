#ifndef MAZDA_MOP_H
#define MAZDA_MOP_H

#include <string>
#include <vector>
#include <utility>

namespace benchmark_problem {

class MazdaMop {
private:
  int ncar;

public:
  int nobj;
  int ncon;
  int nvar;

  std::vector< std::pair<double, double> > varRanges;

public:
  MazdaMop();

  void evaluate(const std::vector<long double> &var,
                      std::vector<long double> &obj,
                      std::vector<long double> &con);

  //original
  void evaluate_con(std::vector<long double> &var, std::vector<long double> &con, std::vector<int> &whitch_con);


protected:
  double evaluate_common_parts_num(std::vector<double> &v1,
                                   std::vector<double> &v2,
                                   std::vector<double> &v3);

};

}

#endif // MAZDA_MOP_H
