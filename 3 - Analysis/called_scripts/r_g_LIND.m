function [g_rand] = r_g_LIND(sigma,num_cells, g_SEV)

% Valori nominali
P_Na = 0.0000014;   % nanolitre_per_second (in sodium_current)
g_Ca_L = 0.004;   % nanoS (in L_type_Ca_channel)
g_Ca_T = 0.006;   % nanoS (in T_type_Ca_channel)
g_to = 0.050002;   % nanoS (in Ca_independent_transient_outward_K_current)
g_Kr = 0.0035;   % nanoS (in delayed_rectifier_K_current)
g_Ks = 0.0025;   % nanoS (in delayed_rectifier_K_current)
g_K1 = 0.00508;   % nanoS (in inward_rectifier)
i_NaK_max = 0.06441;   % picoA (in sodium_potassium_pump)
i_CaP_max = 0.009509;   % picoA (in sarcolemmal_calcium_pump_current)
k_NaCa = 2.0e-5;   % picoA_per_millimolar_4 (in Na_Ca_ion_exchanger_current)
I_up_max = 2.8;   % picoA (in Ca_handling_by_the_SR)

%% NO P_up_Basal ??
nomCond = [P_Na, g_Ca_L, g_Ca_T, g_to, g_Kr, g_Ks, g_K1, i_NaK_max, i_CaP_max, k_NaCa, I_up_max]; %, I_up_max, 0]; 
nomCond = [nomCond, zeros(1, size(g_SEV, 2)-length(nomCond))];

%% Randomizzazione
for i = num_cells:-1:1
    g_rand(i,:) = nomCond .* exp( sigma * randn(1,length(nomCond)) );    
end

end

