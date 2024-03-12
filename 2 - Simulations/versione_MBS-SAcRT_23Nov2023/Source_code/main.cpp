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
                        0.00710587690817621};

	const std::vector<double>
                Y0_atrium{ 
                        -75.8020781109938,
                        0.00627822142185011,
                        0.828100420624812,
                        0.828329449762932,
                        9.95212610778580e-06,
                        0.998912579235043,
                        0.926915254483311,
                        0.000927468829338745,
                        0.955681288633480,
                        0.000298196222285850,
                        0.962523547624839,
                        0.00482701758423975,
                        3.95381366832499e-05,
                        0.0584381048968852,
                        0.000488490292989204,
                        0.758559635332726,
                        0.177010526088852,
                        0.00186483114202541,
                        0.179386049546111,
                        0.140668156867733,
                        0.00776483646752613,
                        0.00768976069147392,
                        7.82770477140038,
                        7.86122335997478,
                        0.000168477350034983,
                        0.000172512414083110,
                        0.609141524109881,
                        0.608937740971786,
                        0,
                        0,
                        0,
                        0,
                        0};

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
                        0};

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
			V_idx[i] = 0;
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
