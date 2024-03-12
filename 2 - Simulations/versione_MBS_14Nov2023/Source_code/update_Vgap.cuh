// -*- mode: c++; c-basic-offset: 8; -*-

/*
    update_Vgap function
    This function computes the gap junctional current exhancged among the current cell and its 4 neighbours

*/

#include <stdio.h>

__device__
double IDX(const double *data, size_t ncols, size_t i, size_t j)
{
        return data[i*ncols + j];
}

__device__
double update_Vgap(const double *V, 	 // membrane voltage matrix (size: nCells x nStates)
                   const size_t V_COLS,  // numero di colonne della matrice V[][]
                   const int *V_idx, 	 // array containing index of membrane voltage state variable in X, different between cell types (size:nCells)
                   const double Cm,	 // membrane capacitance (size: scalar)
                   const double *gJ,     // matrix containing gap junctional resistance values among cells (size: nCells x 4, every cells has Rgap with its 4 neighbours)
                   const size_t gJ_COLS, // numero di colonne della matrice gJ[][]
                   const size_t j, 	 // cell index
                   const int WIDTH,	 // width of the 2D tissue matrix (size: scalar)
                   const int LENGTH )	 // length of the 2D tissue matrix (size: scalar)
{
	//////////////////////////////////////////////////////////////

	/* 2D model function */

	//////////////////////////////////////////////////////////////

	// Cycle all the cells to compute Igap
	double Inet_U, Inet_D, Inet_L, Inet_R;

	if (j >= WIDTH) {
                Inet_L = (IDX(V, V_COLS, j - WIDTH, V_idx[j - WIDTH]) - IDX(V, V_COLS, j, V_idx[j])) * IDX(gJ, gJ_COLS, j, 0);
        } else {
                Inet_L = 0.0;
        }

	if (j < (LENGTH - 1) * WIDTH) {
                Inet_R = (IDX(V, V_COLS, j + WIDTH, V_idx[j + WIDTH]) - IDX(V, V_COLS, j, V_idx[j])) * IDX(gJ, gJ_COLS, j, 2);
        } else {
                Inet_R = 0.0;
        }

	if (j % WIDTH == 0) {
                Inet_U = 0.0;
        } else {
                Inet_U = (IDX(V, V_COLS, j - 1, V_idx[j - 1]) - IDX(V, V_COLS, j, V_idx[j])) * IDX(gJ, gJ_COLS, j, 1);
        }

	if (j % WIDTH < (WIDTH - 1)) {
                Inet_D = (IDX(V, V_COLS, j + 1, V_idx[j + 1]) - IDX(V, V_COLS, j, V_idx[j])) * IDX(gJ, gJ_COLS, j, 3);
        } else {
                Inet_D = 0.0;
        }

	const double Igap = (Inet_U + Inet_D + Inet_L + Inet_R) * (1 / Cm);

/*	if (j == 1900){
		printf("\n\n");
		printf("Vm = %f",  IDX(V, V_COLS, j, V_idx[j]));
		printf("||||| Igap = %.10e |||||", Igap );
		printf("gJ = %f", IDX(gJ, gJ_COLS, j, 3));
	//}
	
		printf("\n\n");
		printf("COLS = %lu, j = %lu", gJ_COLS, j);
		printf("V1 = %f - V2 = %f \n", IDX(V, V_COLS, j - WIDTH, V_idx[j - WIDTH]), IDX(V, V_COLS, j, V_idx[j]) );

		printf("I = %.10e \n\n",  (IDX(V, V_COLS, j - WIDTH, V_idx[j - WIDTH]) - IDX(V, V_COLS, j, V_idx[j])) *  IDX(gJ, gJ_COLS, j, 0) );
		printf("Inet_L = %.10e  \n\n ", Inet_L);

		printf("V1 = %f - V2 = %f \n", IDX(V, V_COLS, j + WIDTH, V_idx[j + WIDTH]), IDX(V, V_COLS, j, V_idx[j]) );
		printf("Inet_R = %.10e  \n\n ", Inet_U);

		printf("V1 = %f - V2 = %f \n", IDX(V, V_COLS, j - 1, V_idx[j - 1]),  IDX(V, V_COLS, j, V_idx[j]) );
		printf("Inet_U = %.10e  \n \n ", Inet_R);

		printf("V1 = %f - V2 = %f \n", IDX(V, V_COLS, j + 1, V_idx[j + 1]), IDX(V, V_COLS, j, V_idx[j]) );
		printf("Inet_D = %.10e \n \n", Inet_D);
	
		printf("||||| Igap = %.10e |||||", Igap );

	}
*/

	return Igap;
}
