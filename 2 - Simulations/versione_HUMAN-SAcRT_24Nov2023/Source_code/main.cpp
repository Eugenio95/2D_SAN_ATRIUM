// -*- mode: c++; c-basic-offset: 8; -*-

// Libraries
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string.h>
#include <vector>
#include <array>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <functional>

// Function declarations
#include "euler_solve.hpp" // Esplicit Euler (fixed step)
#include "const_def.hpp"   // Definition of global constants
#include "mat2D.hpp"       // 2D matrix
#include "configuration.hpp"

double now(void)
{
    return ((double)clock()) / CLOCKS_PER_SEC;
}

const char *usage_message = "Usage:\n\t\
./2D_parallel_SAN_atrium [simulation-params-location/] [params-file-suffix]\n";

int main(int argc, char *argv[])
{
    std::string sim_params_dir;

    // Suffix used to find parameter files
    // Viable options include:
    //     - "_1e-03s_GRAD_SEP_F_101x101_16-May-2022.txt"
    //     - "_1e-03s_GRAD_SEP_F_201x201_12-Sep-2022.txt"
    //     - "_1e-03s_GRAD_SEP_F_71x71_16-May-2022.txt"
    const char* param_file_suffix;

    if (argc < 3) {
        std::cerr << "You must provide the location of your simulation parameter files and a valid params file suffix...\n\n"
                << usage_message;
        return EXIT_FAILURE;
    }

    sim_params_dir = argv[1];
    param_file_suffix = argv[2];

//    Configuration config(sim_params_dir, param_file_suffix);

    std::cout << sim_params_dir << "\n";

    if (sim_params_dir[sim_params_dir.size() - 1] != '/') {
        sim_params_dir += '/';
    }

    Configuration config(sim_params_dir, param_file_suffix); 
    Configuration output_config("Results/", param_file_suffix);

    // Sim parameters 
    std::vector<double> sim_param(sim_param_lgth);
    std::fstream input_parameters = config.get_sim_param_file();

    for (int i = 0; i < sim_param_lgth; ++i)
    { // sim_param.size()
        input_parameters >> sim_param[i];
    }
    input_parameters.close();

    const int WIDTH = (int)sim_param[4];
    const int LENGTH = (int)sim_param[5];
    const double sim_time = sim_param[6];
    const int nCells = WIDTH * LENGTH;
    std::cout << "WIDTH = " << WIDTH << "\n"
         << "LENGTH = " << LENGTH << "\n"
         << "nCells = " << nCells << "\n"
         << "nStates = " << nStates << "\n"
         << "sim_time = " << sim_time << "\n";

	// Initial conditions
	const std::vector<double>
                Y0_SAN{ 
9.40955466945073e-10,
5.12151293481067e-09,
0.844789199516772,
0.155209801130264,
0.273249985872226,
0.165356005672784,
0.113657488564386,
0.0241385808823832,
0.330736929410558,
0.591220749142250,
0.326178030517752,
0.405442066853697,
6.53624830996597e-05,
0.000123909404628413,
-54.5792330377543,
5,
0.000148722415066278,
0.752066181975311,
0.909166753920778,
0.0490707813590274,
0.295933272357133,
0.00263517623841142,
0.00275294121481809,
0.552413899643282,
0.821187134047542,
0.103503074914848,
0.00329822864316888,
0.837948026756321,
0.0421788200039755,
0.181096786989486,
0.00854397970001008,
0.626550885781769,
0.00710587690817621,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0
};

	const std::vector<double>
                Y0_atrium{ 
0.620719406044090,
0.609406980537085,
0.592291226256379,
0.575591322152835,
0.000135710466668922,
0.000138348843595388,
0.000144403571940399,
0.000156364223917474,
0.000162116214501502,
1.05984109023633e-05,
0.998857617215077,
0.998863334145518,
0.974382567759986,
0.0562858539581593,
4.18521440100245e-05,
0.00410213465228339,
0.000310955671231170,
0.975097757230022,
0.904008627935403,
0.904064283232700,
0.00277498208142155,
0.000958915466218931,
0.954360587503485,
-75.3650171682924,
134.743859979303,
0.192830560807385,
0.201265472120290,
0.216457519533416,
0.245627260166889,
0.999378177115573,
0.999511232500315,
0.999561036106774,
0.999971439786227,
9.45128738919411e-05,
7.75793114825318e-05,
5.68354785299318e-05,
3.98571607319244e-05,
0.00466085221641181,
0.00453366886835807,
0.00434703388018010,
0.00427041242214423,
9.26631346151888,
8.67693179582787
};

	std::vector<double> Y0_fat(nStates);

	const std::vector<double>
                Y0_fibro{ 
			-46.6791233290693,
			0.185201718950183,
			0.261543602018550,
			0.198824499759182,
			0.987485898712744,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
                        0,
                        0,
			0
 };

	// Define a unique initial condition matrix and index vectors for voltage, membrane capacitance and cell type
	std::vector< std::vector<double> > Y0(nCells, std::vector<double>(nStates, 0));

//	std::vector<std::vector<double>> Y_n(nCells, std::vector<double>(nStates, 0));

	// Cell phenotypes
	std::vector<int> V_idx(nCells);
	std::vector<double> Cm(nCells);
	// Load geometry
	std::vector<int> cell_type(nCells);

    	std::fstream input_geometry  = config.get_atrial_geometry_file();

	for (int i = 0; i < nCells; ++i) {
		input_geometry >> cell_type[i];
	}
	input_geometry.close();

	int numSAN = 0;
	int numATR = 0;
	int numFAT = 0;
	int numFIBRO = 0;

	for (int i = 0; i < nCells; i++) {
		if ( cell_type[i] == idx_san) {
			V_idx[i] = 14;
			Cm[i] = Cm_SAN;
			numSAN++;

			for (int j = 0; j < nStates; j++) {
				Y0[i][j] = Y0_SAN[j];
			}
		}
		else if ( cell_type[i] == idx_atr ) {
			V_idx[i] = 23;
			Cm[i] = Cm_atrium;
			numATR++;

			for (int j = 0; j < nStates; j++) {
				Y0[i][j] = Y0_atrium[j];
			}
		}
		else if (cell_type[i] == idx_fat) {
			V_idx[i] = 0;
			Cm[i] = Cm_atrium; // I can not divide by 0!!!
			numFAT++;

			for (int j = 0; j < nStates; j++) {
				Y0[i][j] = Y0_fat[j];
			}
		}
		else if (cell_type[i] == idx_fibro) {
			V_idx[i] = 0;
			Cm[i] = Cm_fibro;
			numFIBRO++;

			for (int j = 0; j < nStates; j++) {
				Y0[i][j] = Y0_fibro[j];
			}
		}
	};

	//// Load random conductances
	std::vector<std::vector<double>> rand_g(nCells, std::vector<double>(num_g_rand, 0));
	std::fstream input_conductances = config.get_rand_cond_file();

	for (int i = 0; i < nCells; ++i) {
		for (int j = 0; j < num_g_rand; ++j) {
			input_conductances >> rand_g[i][j];
		}
	}
	input_conductances.close();

	///  Load gap junctions
	std::vector< std::vector < double > > gJ(nCells, std::vector<double>(4, 0));
	std::fstream input_gap_junction = config.get_gap_junc_file();

	for (int i = 0; i < nCells; i++) {
		for (int j = 0; j < 4; j++) {
			input_gap_junction >> gJ[i][j];
		}
	}
	input_gap_junction.close();

	std::fstream fV = output_config.get_Vm_file();
	std::fstream ft = output_config.get_t_file();
	std::fstream fStates = output_config.get_saved_states_file();

	const char* output_states_pre = "saved_states";
	
	char* output_states_string;
	output_states_string = (char*) malloc( strlen(output_states_pre) + strlen(param_file_suffix) );
	strcpy(output_states_string, output_states_pre);
	strcat(output_states_string, param_file_suffix);

	// Start simulation
    	const double start = now();
    	euler_solve(Y0, fV, ft, fStates, WIDTH, LENGTH, nCells, sim_time, cell_type, V_idx, Cm, gJ, rand_g, output_states_string);
    	const double end = now();

    	// Measure time
    	const double ode_time = end - start;
    	std::cout << "Elapsed simulation time: " << ode_time << " s\n";

    	fV.close();
    	ft.close();
    	fStates.close();
	free(output_states_string);

    	return 0;
}
