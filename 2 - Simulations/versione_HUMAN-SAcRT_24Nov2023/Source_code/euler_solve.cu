// -*- mode: c++; c-basic-offset: 8; -*-

/*
    EULER SOLVE function
    This functions performs numerical integration of the ODEs using the esplicit euler method, cycling in time and on the cells.

*/


#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <iomanip>

#include "const_def.hpp"
#include "euler_solve.hpp"
#include "f_sys.hpp"
#include "hpc.h"

using namespace std;

// n. of threads per block (must be an integer multiple of 32)
#define BLKDIM 640

void euler_solve(vector<vector<double>>& X,            // state vector matrix(size: nCells x nubmer of ODEs per cell)
                 fstream& fV,                          // file to save membrane voltage data of the cells
                 fstream& fT,                          // file to save time vector
                 fstream& fStates,                     // file to save all of the cells' state variables every N steps
                 int WIDTH,                            // width of the 2D tissue matrix (size: scalar)
                 int LENGTH,                           // length of the 2D tissue matrix (size: scalar)
                 int nCells,                           // total number of cells in the tissue matrix (size: scalar)
                 double sim_time,                      // time of the simulation in seconds (size: scalar)
                 const vector<int>& cell_type,         // array containing logicals expressing the type of the cell (size: nCells)
                 const vector<int>& V_idx,             // array containing index of membrane voltage state variable in X, different between cell types (size:nCells)
                 const vector<double>& Cm,             // array containing membrane capacitances, different between cell types (size:nCells)
                 const vector<vector<double>>& gJ,     // matrix containing gap junctional resistance values among cells (size: nCells x 4, every cells has Rgap with its 4 neighbours)
                 const vector<vector<double>>& rand_g, // matrix containing randomized ionic maximal conductances values (size: nCells x 12)
		 const string &output_states_string
		)
{
        vector< vector<double> > X_n(nCells, vector<double>(nStates, 0)); // FIXME: l'ho tolto dai parametri della funzione e dichiarato qui

	    // PARAMETERS
        double t = 0;
        const int sim_steps = (int) (1.2 / integration_step); // sim_time
	int last_2s;
const vector<int>
	        front_cells{23836,
		23837,
		23838,
		23839,
		23840,
		23841,
		23842,
		23843,
		23844,
		23845,
		23956,
		23957,
		23958,
		23959,
		23960,
		23961,
		23962,
		23963,
		23964,
		23965,
		24466,
		24467,
		24468,
		24469,
		24470,
		24471,
		24472,
		24473,
		24474,
		24475,
		24496,
		24497,
		24498,
		24499,
		24500,
		24501,
		24502,
		24503,
		24504,
		24505,
		24526,
		24527,
		24528,
		24529,
		24530,
		24531,
		24532,
		24533,
		24534,
		24535};
        
	/* We need to copy X, X_n, gJ and rand_g into properly-sized
           C-style arrays */
        const size_t X_ROWS = X.size();
        const size_t X_COLS = X[0].size();
        const size_t X_n_ROWS = X_n.size();
        const size_t X_n_COLS = X_n[0].size();
        const size_t gJ_ROWS = gJ.size();
        const size_t gJ_COLS = gJ[0].size();
        const size_t rand_g_ROWS = rand_g.size();
        const size_t rand_g_COLS = rand_g[0].size();

        assert(X_ROWS == (size_t)nCells);
        assert(X_COLS == (size_t)nStates);
        assert(X_n_ROWS == (size_t)nCells);
        assert(X_n_COLS == (size_t)nStates);
        assert(gJ_ROWS == (size_t)nCells);
        assert(gJ_COLS == (size_t)4);
        assert(rand_g_ROWS == (size_t)nCells);
        assert(rand_g_COLS == (size_t)num_g_rand);

        double *X_arr = new double[X_ROWS * X_COLS];
        double *X_n_arr = new double[X_n_ROWS * X_n_COLS];
        double *gJ_arr = new double[gJ_ROWS * gJ_COLS];
        double *rand_g_arr = new double[rand_g_ROWS * rand_g_COLS];
        // Each cell requires an array w1[] of nStates elements.  To
        // simplify porting the code to CUDA, we allocate a matrix
        // with nCells * nStates elements; cell j will operate on row
        // j of that matrix.
        double *w1_arr = new double[nCells * nStates]; assert(w1_arr);
        const int *V_idx_arr = V_idx.data();

        /* Copy data in */
        for (size_t i=0; i<X_ROWS; i++) {
                for (size_t j=0; j<X_COLS; j++) {
                        X_arr[i*X_COLS + j] = X[i][j];
                }
                for (size_t j=0; j<gJ_COLS; j++) {
                        gJ_arr[i*gJ_COLS + j] = gJ[i][j];
                }
                for (size_t j=0; j<rand_g_COLS; j++) {
                        rand_g_arr[i*rand_g_COLS + j] = rand_g[i][j];
                }
        }

        // Definition of device copies of the input vectors
        double *d_X_arr;
        const size_t d_X_arr_SIZE = X_ROWS * X_COLS * sizeof(*d_X_arr);
        double *d_X_n_arr;
        const size_t d_X_n_arr_SIZE = X_n_ROWS * X_n_COLS * sizeof(*d_X_n_arr);
        assert(d_X_arr_SIZE == d_X_n_arr_SIZE);
        double *d_gJ_arr;
        const size_t d_gJ_arr_SIZE = gJ_ROWS * gJ_COLS * sizeof(*d_gJ_arr);
        double *d_rand_g_arr;
        const size_t d_rand_g_arr_SIZE = rand_g_ROWS * rand_g_COLS * sizeof(*d_rand_g_arr);
        double *d_w1_arr;
        const size_t d_w1_arr_SIZE = nCells * nStates * sizeof(*d_w1_arr);
        int *d_V_idx_arr;
        const size_t d_V_idx_arr_SIZE = V_idx.size() * sizeof(*d_V_idx_arr);
        double *d_Cm;
        const size_t d_Cm_SIZE = Cm.size() * sizeof(*d_Cm);
        int *d_cell_type;
        const size_t d_cell_type_SIZE = cell_type.size() * sizeof(*d_cell_type);

        // Data allocation and copy
        cudaSafeCall( cudaMalloc( (void**)&d_X_arr, d_X_arr_SIZE ) );
        cudaSafeCall( cudaMemcpy(d_X_arr, X_arr, d_X_arr_SIZE, cudaMemcpyHostToDevice) );

        cudaSafeCall( cudaMalloc( (void**)&d_X_n_arr, d_X_n_arr_SIZE ) );
        // no copy needed

        cudaSafeCall( cudaMalloc( (void**)&d_gJ_arr, d_gJ_arr_SIZE ) );
        cudaSafeCall( cudaMemcpy(d_gJ_arr, gJ_arr, d_gJ_arr_SIZE, cudaMemcpyHostToDevice) );

        cudaSafeCall( cudaMalloc( (void**)&d_rand_g_arr, d_rand_g_arr_SIZE ) );
        cudaSafeCall( cudaMemcpy(d_rand_g_arr, rand_g_arr, d_rand_g_arr_SIZE, cudaMemcpyHostToDevice) );

        cudaSafeCall( cudaMalloc( (void**)&d_w1_arr, d_w1_arr_SIZE ) );
        // no copy needed

        cudaSafeCall( cudaMalloc( (void**)&d_V_idx_arr, d_V_idx_arr_SIZE ) );
        cudaSafeCall( cudaMemcpy(d_V_idx_arr, V_idx_arr, d_V_idx_arr_SIZE, cudaMemcpyHostToDevice) );

        cudaSafeCall( cudaMalloc( (void**)&d_Cm, d_Cm_SIZE ) );
        cudaSafeCall( cudaMemcpy( d_Cm, Cm.data(), d_Cm_SIZE, cudaMemcpyHostToDevice) );

        cudaSafeCall( cudaMalloc( (void**)&d_cell_type, d_cell_type_SIZE) );
        cudaSafeCall( cudaMemcpy( d_cell_type, cell_type.data(), d_cell_type_SIZE, cudaMemcpyHostToDevice) );

        const dim3 BLOCK(BLKDIM);
        const dim3 GRID((nCells + BLKDIM - 1)/BLKDIM);

	if (sim_steps > (int) 5/integration_step){
	//	cout << "\n\n" << "!!! OCCHIO !!!" << "\n\n";
		last_2s = (int) (sim_steps - (5/integration_step));
	} else { 
		last_2s = 0;
	}
	cudaCheckError();

        cout << "<< Sim start >>" << endl;

        // Cycle in time
        for(int steps = 0; steps < sim_steps; steps++) { //sim_steps

                f_sys<<< GRID, BLOCK >>>(nCells,
                                         d_X_arr,
                                         d_X_n_arr,
                                         X_COLS,
                                         d_w1_arr,
                                         nStates,
                                         d_rand_g_arr,
                                         rand_g_COLS,
                                         d_V_idx_arr,
                                         d_Cm,
                                         d_gJ_arr,
                                         gJ_COLS,
                                         WIDTH,
                                         LENGTH,
                                         t,
                                         integration_step,
                                         d_cell_type);
//		cout << "Size: " << d_w1_arr_SIZE << endl;
//		cout << "Cols: " << X_COLS << endl;
//		cout << "dCols: " << nStates << endl;
		cudaCheckError();


                //cudaSafeCall( cudaMemcpy(d_X_arr, d_X_n_arr, d_X_n_arr_SIZE, cudaMemcpyDeviceToDevice) );
                double *tmp = d_X_n_arr;
                d_X_n_arr = d_X_arr;
                d_X_arr = tmp;

                // save voltage vector only every "under_samp" steps
                if (steps > last_2s & steps % under_samp == 0) {
                        cudaSafeCall( cudaMemcpy(X_arr, d_X_arr, d_X_arr_SIZE, cudaMemcpyDeviceToHost) );
                        for (size_t j=0; j<X_ROWS; j++) {
                                fV << X_arr[j*X_COLS + V_idx_arr[j]] << endl;
                        }
			fT << t << endl; // save time vector

			
			// SAVE ALL STATE VARIABLES FOR CELL AT THE FRONTIERI
			for (size_t j = 0; j < front_cells.size(); j++) {
                         	for (size_t k = 0; k < X_COLS; k++) {
                                       	fStates << X_arr[j*X_COLS + k] << endl;
                                }
	                }
		
		}

		//cout << "step = " << steps << endl;
		//cout << setprecision(15) << X_arr[0*X_COLS+14] << "," << X_arr[19900*X_COLS+14] << "," << X_arr[23840*X_COLS+14] << "," << X_arr[22501+14] << endl;

                // save all states variables of all cells every N steps
                if (steps % 10000 == 0) {
                        cout << "Integrating ODEs... t = " << t << " s" << endl; // print progress
                        
                                /*
                        cudaSafeCall( cudaMemcpy(X_arr, d_X_arr, d_X_arr_SIZE, cudaMemcpyDeviceToHost) );
			fStates.close();
			fStates.open( output_states_string, ofstream::out | ofstream::trunc );			
			
			for (size_t j = 0; j < X_ROWS; j++) {
                                for (size_t k = 0; k < X_COLS; k++) {
                                        fStates << X_arr[j*X_COLS + k] << endl;
                                }
                        }
                */
                }
                t += integration_step; // advance in time
	}

        /* Copy data out */
        cudaSafeCall( cudaMemcpy(X_arr, d_X_arr, d_X_arr_SIZE, cudaMemcpyDeviceToHost) );
        cudaSafeCall( cudaMemcpy(X_n_arr, d_X_n_arr, d_X_n_arr_SIZE, cudaMemcpyDeviceToHost) );
        for (size_t i=0; i<X_ROWS; i++) {
                for (size_t j=0; j<X_COLS; j++) {
                        X[i][j] = X_arr[i*nStates + j];
                        X_n[i][j] = X_n_arr[i*nStates + j];
                }
        }

        delete X_arr;
        delete X_n_arr;
        delete gJ_arr;
        delete rand_g_arr;
        delete w1_arr;

        cudaFree( d_X_arr );
        cudaFree( d_X_n_arr );
        cudaFree( d_gJ_arr );
        cudaFree( d_rand_g_arr );
        cudaFree( d_w1_arr );
        cudaFree( d_V_idx_arr );
        cudaFree( d_Cm );
        cudaFree( d_cell_type );
}
