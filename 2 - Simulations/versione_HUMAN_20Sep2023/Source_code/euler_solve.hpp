#ifndef EULER_SOLVE_H
#define EULER_SOLVE_H

#include <vector>

void euler_solve(std::vector<std::vector<double>>& X,
                 std::fstream& fV,
                 std::fstream& fT,
                 std::fstream& fStates,
                 int WIDTH,
                 int LENGTH,
                 int nCells,
                 double sim_time,
                 const std::vector<int>& cell_type,
                 const std::vector<int>& V_idx,
                 const std::vector<double>& Cm,
                 const std::vector<std::vector<double>>& gJ,
                 const std::vector<std::vector<double>>& rand_g,  // Function prototype, its declaration
	 	 const std::string &output_states_string);

#endif
