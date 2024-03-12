function [g_rand] = r_g_MOR(sigma,num_cells)

   g_b_Na = 0.00607;   %% nanoS_per_picoF (in background_currents)
   g_K1 = 0.03;   %% nanoS_per_picoF (in inward_rectifier)
   g_ns = 0.018;   %% nanoS_per_picoF (in non_specific_current)
   i_NaK_max = 2.002;   %% picoA_per_picoF (in sodium_potassium_pump)
   g_to = 0.01652;   %% nanoS_per_picoF (in transient_outward_K_current)
   g_Kur = 0.6;   %% nanoS_per_picoF (in ultra_rapid_K_current)

%% All of them?
nomCond = [g_b_Na, g_K1, g_ns, i_NaK_max, g_to, g_Kur,  zeros(1, 12- 6)]; 

%% Randomizzazion
for i = num_cells:-1:1
    g_rand(i,:) = nomCond .* exp( sigma * randn(1,length(nomCond)) );    
end
   
end

