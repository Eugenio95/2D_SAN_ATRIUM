
gJ_matrix_fib = gJ_matrix;
ind_fib = find(atrial_tissue == 3);

gJ_matrix_fib(ind_fib, :) = G_fib;
gJ_matrix_fib(ind_fib-200, 3) = G_fib;
gJ_matrix_fib(ind_fib-1, 4) = G_fib;
gJ_matrix_fib(ind_fib+200, 1) = G_fib;
gJ_matrix_fib(ind_fib+1, 2) = G_fib;

gJ_matrix_fib(gJ_matrix == 0) = 0;

gJ_matrix = gJ_matrix_fib;


