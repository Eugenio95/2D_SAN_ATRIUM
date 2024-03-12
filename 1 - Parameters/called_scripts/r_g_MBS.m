function [g_rand] = r_g_MBS(sigma,num_cells, g_SAN)

% Valori nominali
gNa     = 350.0;
gCaL    = 14.5;
gt      = 8.1750;
gKur    = 2.35;
gK1     = 3.44;
gKr     = 0.5;
gKs     = 1.0;
INaKmax = 70.8253;
ICaPmax = 2.0;
kNaCa   = 0.0083;
gf      = 1;
Pserca  = 7.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 70% gK1
% gK1 = gK1 * 0.75;
gNa = gNa * 2;
disp('!!! Double I_Na !!!')

nomCond = [gNa, gCaL, gt, gKur, gK1, gKr, gKs, INaKmax, ICaPmax, kNaCa, gf, Pserca]; %, I_up_max, 0]; 
nomCond = [nomCond, zeros(1, size(g_SAN, 2)-length(nomCond))];

%% Randomizzazione
for i = num_cells:-1:1
    g_rand(i,:) = nomCond .* exp( sigma * randn(1,length(nomCond)) );    
end

end
