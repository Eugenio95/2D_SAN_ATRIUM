function [g_rand] = r_g_KOIVU(sigma,num_cells, g_SAN)

% Valori nominali
PNa = 0.0018;   % m3_per_s_times_1e_minus_12 (in ina)
gCaL = 25.3125;   % nS (in ical)
gIf = 1.0;   % nS (in if)
gK1 = 3.825*0.9; % nS (in ik1)
gKr = 0.5;   % nS (in ikr)
gKs = 1.0;   % nS (in iks)
INaKmax = 70.8253;   % pA (in inak)
ICaPmax = 2.0;   % pA (in icap)
kNaCa   = 0.0084;   % m12_A_per_mol4_times_1e_minus_12 (in inaca)
cpumps  = 0.04;   % mM (in serca)
gto     = 1.09*7.5;
gKur    = 0.89*2.75;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 70% gK1
% gK1 = gK1 * 0.75;
PNa = PNa * 2;
disp('!!! Double I_Na !!!')

nomCond = [PNa, gCaL, gIf, gK1, gKr, gKs, INaKmax, ICaPmax, kNaCa, cpumps, gto, gKur]; %, I_up_max, 0]; 
nomCond = [nomCond, zeros(1, size(g_SAN, 2)-length(nomCond))];

%% Randomizzazione
for i = num_cells:-1:1
    g_rand(i,:) = nomCond .* exp( sigma * randn(1,length(nomCond)) );    
end

end
