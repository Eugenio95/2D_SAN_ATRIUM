#ifndef F_SYS_H
#define F_SYS_H

__global__
void f_sys(size_t nCells,
           const double *Ymat,  // nStates * nCells
           double *Y_n_mat,     // nStates * nCells
           size_t Ymat_COLS,    // columns of Ymat
           double *dYmat,       // nSTates * nCells
           size_t dYmat_COLS,   // columns of dYmat
           const double *rand_g_mat, // nStates * num_g_rand
           size_t rand_g_COLS,  // columns of rand_g_mat
           const int *V_idx,    // array
           const double *Cm,    // array
           const double *gJ,    // matrix
           size_t gJ_COLS,
           const int WIDTH,
           const int LENGTH,
           const double time,
           const double h,
           const int *cell_type);

#endif
