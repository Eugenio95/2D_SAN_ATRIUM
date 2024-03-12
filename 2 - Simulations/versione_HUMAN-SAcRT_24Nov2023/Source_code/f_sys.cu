// -*- mode: c++; c-basic-offset: 8; -*-

/*
    F_SYS kernel
    The kernel updates the state variables according to cell types and calculates the currents exchanged by neighbouring cells

*/

#include "koi_2011.cuh"
#include "fab_2017.cuh"
#include "mor_2016.cuh"
#include "update_Vgap.cuh"
#include "f_sys.hpp"
#include "const_def.hpp"

__global__
void f_sys(size_t nCells,       // number of cells
           const double *Ymat,  // nStates * nCells
           double *Y_n_mat,     // nStates * nCells
           size_t Ymat_COLS,    // n. of columns of Ymat
           double *dYmat,       // nStates * nCells
           size_t dYmat_COLS,   // n. of columns of dYmat
           const double *rand_g_mat, // nStates * num_g_rand
           size_t rand_g_COLS,  // n. of columns of rand_g_mat
           const int *V_idx,    // array containing index of membrane voltage state variable in X, different between cell types (size:nCells)
           const double *Cm,    // array containing membrane capacitances, different between cell types (size:nCells)
           const double *gJ,    // matrix containing gap junctional resistance values among cells (size: nCells x 4, every cells has Rgap with its 4 neighbours)
           size_t gJ_COLS,      // n. of columns of gJ
           const int WIDTH,     // width of the 2D tissue matrix (size: scalar)
           const int LENGTH,    // length of the 2D tissue matrix (size: scalar)
           const double time,   // time of simulation at current step
           const double h,      // integration step
           const int *cell_type) // array containing logicals expressing the type of the cell (size: nCells)
{
        const size_t j = threadIdx.x + blockIdx.x * blockDim.x;

        // Threads that are out of bound exit now
        if ( j >= nCells )
                return;

        const double *Y = Ymat + j*Ymat_COLS;
        double *dY = dYmat + j*dYmat_COLS;
        double *Y_n = Y_n_mat + j*Ymat_COLS;
        const double *rand_g = rand_g_mat + j*rand_g_COLS;
        
	for (size_t k = 0; k < Ymat_COLS; k++) {
                dY[k] = 0.0; // [MM] non sono certo che serva...
        }
        switch (cell_type[j]) {
        case idx_atr:
		koi_2011(Y, dY, time, rand_g, j, WIDTH, LENGTH);
                break;
        case idx_san:
		fab_2017(Y, dY, time, rand_g);
                break;
        case idx_fat:
                for (size_t s=0; s<nStates; s++) {
                        dY[s] = Y[s];
                }
                break;
	default: // case idx_fibro:
		mor_2016(Y, dY, time, rand_g);
                break;
                // default:
                // std::cerr << "Unknown cell type: " << cell_type[j] << std::endl;
                // std::abort();
	}
        
	/*if(j == 20000){
		printf("Vm = %f \n", Ymat[Ymat_COLS*j+V_idx[j]]);
	}*/
	
	
	dY[V_idx[j]] += update_Vgap(Ymat,
                                    Ymat_COLS,
                                    V_idx,
                                    Cm[j],
                                    gJ,
                                    gJ_COLS,
                                    j,
                                    WIDTH,
                                    LENGTH); // evaluate Igap

        /*if(j == 20000){
                printf("Vm_post = %f \n\n", Ymat[Ymat_COLS*j+V_idx[j]]);
        }*/


        for (size_t k = 0; k < Ymat_COLS; k++) {
                Y_n[k] = Y[k] + h * dY[k]; // update every temporary state
        }
}
