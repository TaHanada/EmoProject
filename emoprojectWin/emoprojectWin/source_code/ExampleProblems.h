/*
 * Copyright © 2005 The Walking Fish Group (WFG).
 *
 * This material is provided "as is", with no warranty expressed or implied.
 * Any use is at your own risk. Permission to use or copy this software for
 * any purpose is hereby granted without fee, provided this notice is
 * retained on all copies. Permission to modify the code and to distribute
 * modified code is granted, provided a notice that the code was modified is
 * included with the above copyright notice.
 *
 * http://www.wfg.csse.uwa.edu.au/
 */


/*
 * ExampleProblems.h
 *
 * Defines the problems described in the EMO 2005 paper, namely WFG1--WFG9
 * and I1--I5. For the specifics of each problem, refer to the EMO 2005 paper
 * (available from the WFG web site).
 */


#ifndef EXAMPLE_PROBLEMS_H
#define EXAMPLE_PROBLEMS_H


//// Standard includes. /////////////////////////////////////////////////////

#include <vector>


//// Definitions/namespaces. ////////////////////////////////////////////////

namespace WFG
{

namespace Toolkit
{

namespace Examples
{

namespace Problems
{

/*
 * For all problems, the first "k" elements of "z" are the position-related
 * parameters, and "M" is the number of objectives.
 */

//** The WFG1 problem. ******************************************************
std::vector< double > WF1
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The WFG2 problem. ******************************************************
std::vector< double > WF2
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The WFG3 problem. ******************************************************
std::vector< double > WF3
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The WFG4 problem. ******************************************************
std::vector< double > WF4
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The WFG5 problem. ******************************************************
std::vector< double > WF5
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The WFG6 problem. ******************************************************
std::vector< double > WF6
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The WFG7 problem. ******************************************************
std::vector< double > WF7
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The WFG8 problem. ******************************************************
std::vector< double > WF8
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The WFG9 problem. ******************************************************
std::vector< double > WF9
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The I1 problem. ********************************************************
std::vector< double > I1
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The I2 problem. ********************************************************
std::vector< double > I2
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The I3 problem. ********************************************************
std::vector< double > I3
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The I4 problem. ********************************************************
std::vector< double > I4
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** The I5 problem. ********************************************************
std::vector< double > I5
(
  const std::vector< double >& z,
  const int k,
  const int M
);

//** constrained variobles of CWFG *****************************************
std::vector<double> get_constrained_var_WFG1(const std::vector< double >& z, const int k, const int M,const int c,const int num_cons);

}  // Problems namespace

}  // Examples namespace

}  // Toolkit namespace

}  // WFG namespace

#endif