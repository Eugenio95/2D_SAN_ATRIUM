#ifndef CONST_DEF_HPP
#define CONST_DEF_HPP

// State variables
#define nStates 33

// Random conductances
#define num_g_rand 12 // Lind = 12, SEV = 9

// Membrane capacitances
#define Cm_SAN 57E-12
#define Cm_atrium 56E-12 // rescaled from 0.05E-12 to have current in [nA]
#define Cm_fibro 12.4E-12

// Cell types
#define idx_atr 0
#define idx_san 1
#define idx_fibro 3
#define idx_fat 9

// Keep only 1 timestep every 1000 -> to get sharper upstroke
#define integration_step 5E-6 
#define under_samp 40

// Length of simulation paramameters input file
#define sim_param_lgth 16

#endif
