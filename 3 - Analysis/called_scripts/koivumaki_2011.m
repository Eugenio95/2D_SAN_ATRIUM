%===============================================================================
% CellML file:   C:\mia_roba\Dottorato\Archivio Modelli\Koivumaki\human_atrial\koivumaki-2011-pmr.cellml
% CellML model:  koivumaki_2011
% Date and time: 6.10.2021 at 16.16.31
%-------------------------------------------------------------------------------
% Conversion from CellML 1.0 to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2021 Dr Alan Garny
%    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%-------------------------------------------------------------------------------
% http://www.cellml.org/
%===============================================================================

function dY = koivumaki_2011(time, Y, durStim, iStim)

global numDiffEqEval
numDiffEqEval = numDiffEqEval+1;
%-------------------------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------------------------

% Y = [0.6189225, 0.6076289, 0.5905266, 0.5738108, 1.35496500000000013e-4, 1.38142100000000014e-4, 1.44208699999999994e-4, 1.56184399999999995e-4, 1.61937700000000013e-4, 1.06091699999999996e-5, 0.9988566, 0.9988624, 0.9744374, 5.62066499999999969e-2, 4.18941700000000008e-5, 4.10975100000000003e-3, 3.11170299999999984e-4, 0.9751094, 0.90391, 0.9039673, 2.7758119999999999e-3, 9.59425800000000026e-4, 0.954338, -75.42786, 134.6313, 0.1925362, 0.2010345, 0.2163122, 0.2455297, 0.9993722, 0.9995086, 0.9995604, 0.9999717, 9.47851400000000044e-5, 7.76550300000000031e-5, 5.67494700000000006e-5, 3.97509699999999973e-5, 4.63856499999999988e-3, 4.5120780000000001e-3, 4.32640899999999981e-3, 4.25044500000000026e-3, 9.28686, 8.691504];

% YNames = {'CaSR1', 'CaSR2', 'CaSR3', 'CaSR4', 'Cai1', 'Cai2', 'Cai3', 'Cai4', 'Cass', 'd', 'f1', 'f2', 'fca', 'y', 'pa', 'n', 'r1', 's1', 'h1', 'h2', 'm', 'r2', 's2', 'V', 'Ki', 'a11', 'a21', 'a31', 'ass1', 'c1', 'c2', 'c3', 'css', 'o1', 'o2', 'o3', 'oss', 'a12', 'a22', 'a32', 'ass2', 'Nai', 'Nass'};
% YUnits = {'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'mV', 'mM', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'hertz', 'mM', 'mM', 'mM', 'mM', 'mM', 'mM'};
% YComponents = {'calcium', 'calcium', 'calcium', 'calcium', 'calcium', 'calcium', 'calcium', 'calcium', 'calcium', 'ical', 'ical', 'ical', 'ical', 'if', 'ikr', 'iks', 'ikur', 'ikur', 'ina', 'ina', 'ina', 'it', 'it', 'membrane', 'potassium', 'ryr', 'ryr', 'ryr', 'ryr', 'ryr', 'ryr', 'ryr', 'ryr', 'ryr', 'ryr', 'ryr', 'ryr', 'serca', 'serca', 'serca', 'serca', 'sodium', 'sodium'};

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: CaSR1 (mM) (in calcium)
% 2: CaSR2 (mM) (in calcium)
% 3: CaSR3 (mM) (in calcium)
% 4: CaSR4 (mM) (in calcium)
% 5: Cai1 (mM) (in calcium)
% 6: Cai2 (mM) (in calcium)
% 7: Cai3 (mM) (in calcium)
% 8: Cai4 (mM) (in calcium)
% 9: Cass (mM) (in calcium)
% 10: d (hertz) (in ical)
% 11: f1 (hertz) (in ical)
% 12: f2 (hertz) (in ical)
% 13: fca (hertz) (in ical)
% 14: y (hertz) (in if)
% 15: pa (hertz) (in ikr)
% 16: n (hertz) (in iks)
% 17: r1 (hertz) (r in ikur)
% 18: s1 (hertz) (s in ikur)
% 19: h1 (hertz) (in ina)
% 20: h2 (hertz) (in ina)
% 21: m (hertz) (in ina)
% 22: r2 (hertz) (r in it)
% 23: s2 (hertz) (s in it)
% 24: V (mV) (in membrane)
% 25: Ki (mM) (in potassium)
% 26: a11 (hertz) (a1 in ryr)
% 27: a21 (hertz) (a2 in ryr)
% 28: a31 (hertz) (a3 in ryr)
% 29: ass1 (hertz) (ass in ryr)
% 30: c1 (hertz) (in ryr)
% 31: c2 (hertz) (in ryr)
% 32: c3 (hertz) (in ryr)
% 33: css (hertz) (in ryr)
% 34: o1 (hertz) (in ryr)
% 35: o2 (hertz) (in ryr)
% 36: o3 (hertz) (in ryr)
% 37: oss (hertz) (in ryr)
% 38: a12 (mM) (a1 in serca)
% 39: a22 (mM) (a2 in serca)
% 40: a32 (mM) (a3 in serca)
% 41: ass2 (mM) (ass in serca)
% 42: Nai (mM) (in sodium)
% 43: Nass (mM) (in sodium)

%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------

BCa = 0.024;   % mM (in calcium)
CSQN = 6.7;   % mM (in calcium)
DCa = 780.0;   % m2_per_s_times_1e_minus_12 (in calcium)
DCaBm = 25.0;   % m2_per_s_times_1e_minus_12 (in calcium)
DCaSR = 44.0;   % m2_per_s_times_1e_minus_12 (in calcium)
KdBCa = 0.00238;   % mM (in calcium)
KdCSQN = 0.8;   % mM (in calcium)
KdSLhigh = 0.013;   % mM (in calcium)
KdSLlow = 1.1;   % mM (in calcium)
SLhigh = 13.0;   % mM (in calcium)
SLlow = 165.0;   % mM (in calcium)
kSRleak = 0.006;   % hertz (in calcium)
Cm = 0.05;   % nF (in cell)
Vss = 4.99231999999999966e-5;   % nL (in cell)
dx = 1.625;   % um (in cell)
lcell = 122.051;   % um (in cell)
piGreco = 3.14159265358979312e0;   % dimensionless (in cell)
rjunct = 6.5;   % um (in cell)
Cao = 1.8;   % mM (in extra)
Ko = 5.4;   % mM (in extra)
Nao = 130.0;   % mM (in extra)
gCab = 0.0952;   % nS (in icab)
ECa_app = 60.0;   % mV (in ical)
gCaL = 25.3125;   % nS (in ical)
ical_fca_tau = 0.002;   % second (in ical)
kCa = 0.001;   % mM (in ical)
kCan = 2.0;   % dimensionless (in ical)
ICaPmax = 2.0;   % pA (in icap)
kCaP = 0.0005;   % mM (in icap)
gIf = 1.0;   % nS (in if)
gKr = 0.5;   % nS (in ikr)
gKs = 1.0;   % nS (in iks)
gNab = 0.060599;   % nS (in inab)
dNaCa = 0.0003;   % m12_per_mol4 (in inaca)
fCaNCX = 1.0;   % dimensionless (in inaca)
gam = 0.45;   % dimensionless (in inaca)
kNaCa = 0.0084;   % m12_A_per_mol4_times_1e_minus_12 (in inaca)
INaKmax = 70.8253;   % pA (in inak)
kNaKK = 1.0;   % mM (in inak)
kNaKNa = 11.0;   % dimensionless (in inak)
PNa = 0.0018;   % m3_per_s_times_1e_minus_12 (in ina)
F = 96487.0;   % C_per_mol (in phys)
R = 8314.0;   % mJ_per_mol_per_K (in phys)
T = 306.15;   % kelvin (in phys)
tau_act = 0.01875;   % second (in ryr)
tau_actss = 0.005;   % second (in ryr)
tau_adapt = 1.0;   % second (in ryr)
tau_inact = 0.0875;   % second (in ryr)
tau_inactss = 0.015;   % second (in ryr)
SERCAKmf = 0.00025;   % mM (in serca)
SERCAKmr = 1.8;   % mM (in serca)
cpumps = 0.04;   % mM (in serca)
k4 = 7.5;   % hertz (in serca)
DNa = 0.12;   % m2_per_s_times_1e_minus_12 (in sodium)
KdBNa = 10.0;   % mM (in sodium)

duration = durStim;   % second (in stimulus)
offset = 0.01;   % second (in stimulus)
period = 1;   % second (in stimulus)

% PNa = PNa * 2;
% gKs = gKs * 0.5;
% gKr = gKr * 0.5;

%-------------------------------------------------------------------------------
% Computed variables
%-------------------------------------------------------------------------------

% JCa1 (kat_times_1e_minus_12) (in calcium)
% JCa2 (kat_times_1e_minus_12) (in calcium)
% JCa3 (kat_times_1e_minus_12) (in calcium)
% JCa4 (kat_times_1e_minus_12) (in calcium)
% JCass (kat_times_1e_minus_12) (in calcium)
% JSRCa1 (kat_times_1e_minus_12) (in calcium)
% JSRCa2 (kat_times_1e_minus_12) (in calcium)
% JSRCa3 (kat_times_1e_minus_12) (in calcium)
% JSRCa4 (kat_times_1e_minus_12) (in calcium)
% JSRCaleak1 (kat_times_1e_minus_12) (in calcium)
% JSRCaleak2 (kat_times_1e_minus_12) (in calcium)
% JSRCaleak3 (kat_times_1e_minus_12) (in calcium)
% JSRCaleakss (kat_times_1e_minus_12) (in calcium)
% Jj_nj (kat_times_1e_minus_12) (in calcium)
% calcium_CaSR1_beta (dimensionless) (in calcium)
% calcium_CaSR2_beta (dimensionless) (in calcium)
% calcium_CaSR3_beta (dimensionless) (in calcium)
% calcium_CaSR4_beta (dimensionless) (in calcium)
% calcium_Cai1_beta (dimensionless) (in calcium)
% calcium_Cai1_gamma (dimensionless) (in calcium)
% calcium_Cai2_beta (dimensionless) (in calcium)
% calcium_Cai2_gamma (dimensionless) (in calcium)
% calcium_Cai3_beta (dimensionless) (in calcium)
% calcium_Cai3_gamma (dimensionless) (in calcium)
% calcium_Cai4_beta (dimensionless) (in calcium)
% calcium_Cai4_gamma (dimensionless) (in calcium)
% calcium_Cass_beta (dimensionless) (in calcium)
% calcium_Cass_i_tot (pA) (in calcium)
% Aj_nj (um2) (in cell)
% VSR1 (nL) (in cell)
% VSR2 (nL) (in cell)
% VSR3 (nL) (in cell)
% VSR4 (nL) (in cell)
% Vcytosol (nL) (in cell)
% Vnonjunct1 (nL) (in cell)
% Vnonjunct2 (nL) (in cell)
% Vnonjunct3 (nL) (in cell)
% Vnonjunct4 (nL) (in cell)
% Vnonjunct_Nai (nL) (in cell)
% xj_nj (um) (in cell)
% xj_nj_Nai (um) (in cell)
% ICab (pA) (in icab)
% ICaL (pA) (in ical)
% f_inf (dimensionless) (in ical)
% ical_d_inf (dimensionless) (in ical)
% ical_d_tau (second) (in ical)
% ical_f1_tau (second) (in ical)
% ical_f2_tau (second) (in ical)
% ical_fca_inf (dimensionless) (in ical)
% ICaP (pA) (in icap)
% If (pA) (in if)
% IfK (pA) (in if)
% IfNa (pA) (in if)
% if_y_inf (dimensionless) (in if)
% if_y_tau (second) (in if)
% IK1 (pA) (in ik1)
% gK1 (nS) (in ik1)
% IKr (pA) (in ikr)
% ikr_pa_inf (dimensionless) (in ikr)
% ikr_pa_tau (second) (in ikr)
% piIKr (dimensionless) (in ikr)
% IKs (pA) (in iks)
% iks_n_inf (dimensionless) (in iks)
% iks_n_tau (second) (in iks)
% IKur (pA) (in ikur)
% gKur (nS) (in ikur)
% ikur_r_inf (dimensionless) (in ikur)
% ikur_r_tau (second) (in ikur)
% ikur_s_inf (dimensionless) (in ikur)
% ikur_s_tau (second) (in ikur)
% INab (pA) (in inab)
% INaCa (pA) (in inaca)
% INaK (pA) (in inak)
% Nass15 (dimensionless) (in inak)
% INa (pA) (in ina)
% h_inf (dimensionless) (in ina)
% ina_h1_tau (second) (in ina)
% ina_h2_tau (second) (in ina)
% ina_m_inf (dimensionless) (in ina)
% ina_m_tau (second) (in ina)
% It (pA) (in it)
% gt (nS) (in it)
% it_r_inf (dimensionless) (in it)
% it_r_tau (second) (in it)
% it_s_inf (dimensionless) (in it)
% it_s_tau (second) (in it)
% i_ion (pA) (in membrane)
% ECa (mV) (in nernst)
% EK (mV) (in nernst)
% ENa (mV) (in nernst)
% FRT (per_mV) (in phys)
% RTF (mV) (in phys)
% i_tot (pA) (in potassium)
% Jrel1 (kat_times_1e_minus_12) (in ryr)
% Jrel2 (kat_times_1e_minus_12) (in ryr)
% Jrel3 (kat_times_1e_minus_12) (in ryr)
% Jrelss (kat_times_1e_minus_12) (in ryr)
% SRCa1 (dimensionless) (in ryr)
% SRCa2 (dimensionless) (in ryr)
% SRCa3 (dimensionless) (in ryr)
% SRCass (dimensionless) (in ryr)
% ainf1 (dimensionless) (in ryr)
% ainf2 (dimensionless) (in ryr)
% ainf3 (dimensionless) (in ryr)
% ainfss (dimensionless) (in ryr)
% cinf1 (dimensionless) (in ryr)
% cinf2 (dimensionless) (in ryr)
% cinf3 (dimensionless) (in ryr)
% cinfss (dimensionless) (in ryr)
% nu1 (m3_per_s_times_1e_minus_12) (in ryr)
% nu2 (m3_per_s_times_1e_minus_12) (in ryr)
% nu3 (m3_per_s_times_1e_minus_12) (in ryr)
% nuss (m3_per_s_times_1e_minus_12) (in ryr)
% oinf1 (dimensionless) (in ryr)
% oinf2 (dimensionless) (in ryr)
% oinf3 (dimensionless) (in ryr)
% oinfss (dimensionless) (in ryr)
% J_SERCASR1 (kat_times_1e_minus_12) (in serca)
% J_SERCASR2 (kat_times_1e_minus_12) (in serca)
% J_SERCASR3 (kat_times_1e_minus_12) (in serca)
% J_SERCASRss (kat_times_1e_minus_12) (in serca)
% J_bulkSERCA1 (kat_times_1e_minus_12) (in serca)
% J_bulkSERCA2 (kat_times_1e_minus_12) (in serca)
% J_bulkSERCA3 (kat_times_1e_minus_12) (in serca)
% J_bulkSERCAss (kat_times_1e_minus_12) (in serca)
% k1 (m6_per_s_per_mol2) (in serca)
% k2 (hertz) (in serca)
% k3 (m6_per_s_per_mol2) (in serca)
% BNa (mM) (in sodium)
% JNa (kat_times_1e_minus_12) (in sodium)
% betaNass (dimensionless) (in sodium)
% i_ss (pA) (in sodium)
% amplitude (pA) (in stimulus)
% i_stim (pA) (in stimulus)

%-------------------------------------------------------------------------------
% Computation
%-------------------------------------------------------------------------------

% time (second)

calcium_CaSR1_beta = 1.0/(1.0+CSQN*KdCSQN/(Y(1)+KdCSQN)^2.0);
k3 = k4/SERCAKmr^2.0;
Vnonjunct1 = ((1.0*dx)^2.0-(0.0*dx)^2.0)*piGreco*lcell*0.5*1.0e-6;
J_SERCASR1 = (-k3*Y(1)^2.0*(cpumps-Y(38))+k4*Y(38))*Vnonjunct1*2.0;
JSRCaleak1 = kSRleak*(Y(1)-Y(5))*Vnonjunct1;
nu1 = 1.0*Vnonjunct1;
SRCa1 = 1.0-1.0/(1.0+exp((Y(1)-0.3)/0.1));
Jrel1 = nu1*Y(34)*Y(30)*SRCa1*(Y(1)-Y(5));
JSRCa1 = J_SERCASR1-JSRCaleak1-Jrel1;
VSR1 = 0.05*Vnonjunct1/2.0*0.9;
dY(1, 1) = calcium_CaSR1_beta*DCaSR*((Y(2)-2.0*Y(1)+Y(1))/dx^2.0+(Y(2)-Y(1))/(2.0*1.0*dx^2.0))+JSRCa1/VSR1*calcium_CaSR1_beta;
calcium_CaSR2_beta = 1.0/(1.0+CSQN*KdCSQN/(Y(2)+KdCSQN)^2.0);
Vnonjunct2 = ((2.0*dx)^2.0-(1.0*dx)^2.0)*piGreco*lcell*0.5*1.0e-6;
J_SERCASR2 = (-k3*Y(2)^2.0*(cpumps-Y(39))+k4*Y(39))*Vnonjunct2*2.0;
JSRCaleak2 = kSRleak*(Y(2)-Y(6))*Vnonjunct2;
nu2 = 1.0*Vnonjunct2;
SRCa2 = 1.0-1.0/(1.0+exp((Y(2)-0.3)/0.1));
Jrel2 = nu2*Y(35)*Y(31)*SRCa2*(Y(2)-Y(6));
JSRCa2 = J_SERCASR2-JSRCaleak2-Jrel2;
VSR2 = 0.05*Vnonjunct2/2.0*0.9;
dY(2, 1) = calcium_CaSR2_beta*DCaSR*((Y(3)-2.0*Y(2)+Y(1))/dx^2.0+(Y(3)-Y(1))/(2.0*2.0*dx^2.0))+JSRCa2/VSR2*calcium_CaSR2_beta;
calcium_CaSR3_beta = 1.0/(1.0+CSQN*KdCSQN/(Y(3)+KdCSQN)^2.0);
Vnonjunct3 = ((3.0*dx)^2.0-(2.0*dx)^2.0)*piGreco*lcell*0.5*1.0e-6;
J_SERCASR3 = (-k3*Y(3)^2.0*(cpumps-Y(40))+k4*Y(40))*Vnonjunct3*2.0;
JSRCaleak3 = kSRleak*(Y(3)-Y(7))*Vnonjunct3;
nu3 = 1.0*Vnonjunct3;
SRCa3 = 1.0-1.0/(1.0+exp((Y(3)-0.3)/0.1));
Jrel3 = nu3*Y(36)*Y(32)*SRCa3*(Y(3)-Y(7));
JSRCa3 = J_SERCASR3-JSRCaleak3-Jrel3;
VSR3 = 0.05*Vnonjunct3/2.0*0.9;
dY(3, 1) = calcium_CaSR3_beta*DCaSR*((Y(4)-2.0*Y(3)+Y(2))/dx^2.0+(Y(4)-Y(2))/(2.0*3.0*dx^2.0))+JSRCa3/VSR3*calcium_CaSR3_beta;
calcium_CaSR4_beta = 1.0/(1.0+CSQN*KdCSQN/(Y(4)+KdCSQN)^2.0);
J_SERCASRss = (-k3*Y(4)^2.0*(cpumps-Y(41))+k4*Y(41))*Vss*2.0;
JSRCaleakss = kSRleak*(Y(4)-Y(9))*Vss;
nuss = 625.0*Vss;
SRCass = 1.0-1.0/(1.0+exp((Y(4)-0.3)/0.1));
Jrelss = nuss*Y(37)*Y(33)*SRCass*(Y(4)-Y(9));
JSRCa4 = J_SERCASRss-JSRCaleakss-Jrelss;
Vnonjunct4 = ((4.0*dx)^2.0-(3.0*dx)^2.0)*piGreco*lcell*0.5*1.0e-6;
VSR4 = 0.05*Vnonjunct4/2.0*0.9;
dY(4, 1) = calcium_CaSR4_beta*DCaSR*((Y(4)-2.0*Y(4)+Y(3))/dx^2.0+(Y(4)-Y(3))/(2.0*4.0*dx^2.0))+JSRCa4/VSR4*calcium_CaSR4_beta;
calcium_Cai1_gamma = BCa*KdBCa/(Y(5)+KdBCa)^2.0;
calcium_Cai1_beta = 1.0/(1.0+calcium_Cai1_gamma);
k1 = 1000.0^2.0*k4;
k2 = k1*SERCAKmf^2.0;
J_bulkSERCA1 = (k1*Y(5)^2.0*(cpumps-Y(38))-k2*Y(38))*Vnonjunct1*2.0;
JCa1 = -J_bulkSERCA1+JSRCaleak1+Jrel1;
dY(5, 1) = calcium_Cai1_beta*(DCa+calcium_Cai1_gamma*DCaBm)*((Y(6)-2.0*Y(5)+Y(5))/dx^2.0+(Y(6)-Y(5))/(2.0*1.0*dx^2.0))-2.0*calcium_Cai1_beta*calcium_Cai1_gamma*DCaBm/(KdBCa+Y(5))*((Y(6)-Y(5))/(2.0*dx))^2.0+JCa1/Vnonjunct1*calcium_Cai1_beta;
calcium_Cai2_gamma = BCa*KdBCa/(Y(6)+KdBCa)^2.0;
calcium_Cai2_beta = 1.0/(1.0+calcium_Cai2_gamma);
J_bulkSERCA2 = (k1*Y(6)^2.0*(cpumps-Y(39))-k2*Y(39))*Vnonjunct2*2.0;
JCa2 = -J_bulkSERCA2+JSRCaleak2+Jrel2;
dY(6, 1) = calcium_Cai2_beta*(DCa+calcium_Cai2_gamma*DCaBm)*((Y(7)-2.0*Y(6)+Y(5))/dx^2.0+(Y(7)-Y(5))/(2.0*2.0*dx^2.0))-2.0*calcium_Cai2_beta*calcium_Cai2_gamma*DCaBm/(KdBCa+Y(6))*((Y(7)-Y(5))/(2.0*dx))^2.0+JCa2/Vnonjunct2*calcium_Cai2_beta;
calcium_Cai3_gamma = BCa*KdBCa/(Y(7)+KdBCa)^2.0;
calcium_Cai3_beta = 1.0/(1.0+calcium_Cai3_gamma);
J_bulkSERCA3 = (k1*Y(7)^2.0*(cpumps-Y(40))-k2*Y(40))*Vnonjunct3*2.0;
JCa3 = -J_bulkSERCA3+JSRCaleak3+Jrel3;
dY(7, 1) = calcium_Cai3_beta*(DCa+calcium_Cai3_gamma*DCaBm)*((Y(8)-2.0*Y(7)+Y(6))/dx^2.0+(Y(8)-Y(6))/(2.0*3.0*dx^2.0))-2.0*calcium_Cai3_beta*calcium_Cai3_gamma*DCaBm/(KdBCa+Y(7))*((Y(8)-Y(6))/(2.0*dx))^2.0+JCa3/Vnonjunct3*calcium_Cai3_beta;
calcium_Cai4_gamma = BCa*KdBCa/(Y(8)+KdBCa)^2.0;
calcium_Cai4_beta = 1.0/(1.0+calcium_Cai4_gamma);
Aj_nj = piGreco*rjunct*2.0*lcell*0.5;
xj_nj = 0.02/2.0+dx/2.0;
Jj_nj = DCa*Aj_nj/xj_nj*(Y(9)-Y(8))*1.0e-6;
JCa4 = Jj_nj;
dY(8, 1) = calcium_Cai4_beta*(DCa+calcium_Cai4_gamma*DCaBm)*((Y(8)-2.0*Y(8)+Y(7))/dx^2.0+(Y(8)-Y(7))/(2.0*4.0*dx^2.0))-2.0*calcium_Cai4_beta*calcium_Cai4_gamma*DCaBm/(KdBCa+Y(8))*((Y(8)-Y(7))/(2.0*dx))^2.0+JCa4/Vnonjunct4*calcium_Cai4_beta;
calcium_Cass_beta = 1.0/(1.0+SLlow*KdSLlow/(Y(9)+KdSLlow)^2.0+SLhigh*KdSLhigh/(Y(9)+KdSLhigh)^2.0+BCa*KdBCa/(Y(9)+KdBCa)^2.0);
J_bulkSERCAss = (k1*Y(9)^2.0*(cpumps-Y(41))-k2*Y(41))*Vss*2.0;
JCass = -Jj_nj+JSRCaleakss-J_bulkSERCAss+Jrelss;
ICaL = gCaL*Y(10)*Y(13)*Y(11)*Y(12)*(Y(24)-ECa_app);
RTF = R*T/F;
ECa = RTF*log(Cao/Y(9))/2.0;
ICab = gCab*(Y(24)-ECa);
ICaP = ICaPmax*Y(9)/(kCaP+Y(9));
FRT = 1.0/RTF;
INaCa = kNaCa*(exp(gam*Y(24)*FRT)*Y(43)^3.0*Cao-exp((gam-1.0)*Y(24)*FRT)*Nao^3.0*Y(9)*fCaNCX)/(1.0+dNaCa*(Nao^3.0*Y(9)*fCaNCX+Y(43)^3.0*Cao));
calcium_Cass_i_tot = -ICaL-ICab-ICaP+2.0*INaCa;
dY(9, 1) = calcium_Cass_beta*(JCass/Vss+calcium_Cass_i_tot/(2.0*Vss*F));
Vnonjunct_Nai = Vnonjunct1+Vnonjunct2+Vnonjunct3+Vnonjunct4;
Vcytosol = Vnonjunct_Nai+Vss;
xj_nj_Nai = 0.02/2.0+2.0*dx;
ical_d_inf = 1.0/(1.0+exp((Y(24)+9.0)/-5.8));
ical_d_tau = 0.0027*exp(-((Y(24)+35.0)/30.0)^2.0)+0.002;
dY(10, 1) = (ical_d_inf-Y(10))/ical_d_tau;
f_inf = 1.0/(1.0+exp((Y(24)+27.4)/7.1));
ical_f1_tau = 0.98698*exp(-((Y(24)+30.16047)/7.09396)^2.0)+0.04275/(1.0+exp((Y(24)-51.61555)/-80.61331))+0.03576/(1.0+exp((Y(24)+29.57272)/13.21758))-0.00821;
dY(11, 1) = (f_inf-Y(11))/ical_f1_tau;
ical_f2_tau = 1.3323*exp(-((Y(24)+40.0)/14.2)^2.0)+0.0626;
dY(12, 1) = (f_inf-Y(12))/ical_f2_tau;
ical_fca_inf = 1.0-1.0/(1.0+(kCa/Y(9))^kCan);
dY(13, 1) = (ical_fca_inf-Y(13))/ical_fca_tau;
EK = RTF*log(Ko/Y(25));
IfK = gIf*Y(14)*(1.0-0.2677)*(Y(24)-EK);
ENa = RTF*log(Nao/Y(43));
IfNa = gIf*Y(14)*0.2677*(Y(24)-ENa);
If = IfK+IfNa;
if_y_inf = 1.0/(1.0+exp((Y(24)+97.82874)/12.48025));
if_y_tau = 1.0/(0.00332*exp(-Y(24)/16.54103)+23.71839*exp(Y(24)/16.54103));
dY(14, 1) = (if_y_inf-Y(14))/if_y_tau;
gK1 = 3.825*0.9;
IK1 = gK1*(Ko*1.0)^0.4457*(Y(24)-EK)/(1.0+exp(1.5*(Y(24)-EK+3.6)*FRT));
piIKr = 1.0/(1.0+exp((Y(24)+55.0)/24.0));
IKr = gKr*Y(15)*piIKr*(Y(24)-EK);
ikr_pa_inf = 1.0/(1.0+exp((Y(24)+15.0)/-6.0));
ikr_pa_tau = 0.21718*exp(-((Y(24)+20.1376)/22.1996)^2.0)+0.03118;
dY(15, 1) = (ikr_pa_inf-Y(15))/ikr_pa_tau;
IKs = gKs*Y(16)*(Y(24)-EK);
iks_n_inf = 1.0/(1.0+exp((Y(24)-19.9)/-12.7));
iks_n_tau = 0.4*exp(-((Y(24)-20.0)/20.0)^2.0)+0.7;
dY(16, 1) = (iks_n_inf-Y(16))/iks_n_tau;
gKur = 0.89*2.75;
IKur = gKur*Y(17)*Y(18)*(Y(24)-EK);
ikur_r_inf = 1.0/(1.0+exp((Y(24)+6.0)/-8.6));
ikur_r_tau = 0.009/(1.0+exp((Y(24)+5.0)/12.0))+0.0005;
ikur_s_inf = 1.0/(1.0+exp((Y(24)+7.5)/10.0));
ikur_s_tau = 0.59/(1.0+exp((Y(24)+60.0)/10.0))+3.05;
dY(17, 1) = (ikur_r_inf-Y(17))/ikur_r_tau;
dY(18, 1) = (ikur_s_inf-Y(18))/ikur_s_tau;
INa = PNa*Y(21)^3.0*(0.9*Y(19)+0.1*Y(20))*Nao*Y(24)*F*FRT*(exp((Y(24)-ENa)*FRT)-1.0)/(exp(Y(24)*FRT)-1.0);
h_inf = 1.0/(1.0+exp((Y(24)+63.6)/5.3));
ina_h1_tau = 0.03/(1.0+exp((Y(24)+35.1)/3.2))+0.0003;
dY(19, 1) = (h_inf-Y(19))/ina_h1_tau;
ina_h2_tau = 0.12/(1.0+exp((Y(24)+35.1)/3.2))+0.003;
dY(20, 1) = (h_inf-Y(20))/ina_h2_tau;
ina_m_inf = 1.0/(1.0+exp((Y(24)+27.12)/-8.21));
ina_m_tau = 4.2e-5*exp(-((Y(24)+25.57)/28.8)^2.0)+2.4e-5;
dY(21, 1) = (ina_m_inf-Y(21))/ina_m_tau;
INab = gNab*(Y(24)-ENa);
Nass15 = (Y(43)*1.0)^1.5;
INaK = INaKmax*Ko/(Ko+kNaKK)*Nass15/(Nass15+kNaKNa^1.5)*(Y(24)+150.0)/(Y(24)+200.0);
gt = 1.09*7.5;
It = gt*Y(22)*Y(23)*(Y(24)-EK);
it_r_inf = 1.0/(1.0+exp((Y(24)-1.0)/-11.0));
it_r_tau = 0.0035*exp(-((Y(24)+0.0)/30.0)^2.0)+0.0015;
it_s_inf = 1.0/(1.0+exp((Y(24)+40.5)/11.5));
it_s_tau = 0.025635*exp(-((Y(24)+52.45)/15.8827)^2.0)+0.01414;
dY(22, 1) = (it_r_inf-Y(22))/it_r_tau;
dY(23, 1) = (it_s_inf-Y(23))/it_s_tau;
i_ion = INa+ICaL+It+IKur+IK1+IKr+IKs+INab+ICab+INaK+ICaP+INaCa+If;
amplitude = iStim;

if (time-offset-period*floor((time-offset)/period) < duration)
   i_stim = amplitude;
else
   i_stim = 0.0;
end

dY(24, 1) = -(i_ion+i_stim)/Cm;
i_tot = It+IKur+IK1+IKr+IKs-2.0*INaK+IfK+i_stim;
dY(25, 1) = -i_tot/(Vcytosol*F);
ainf1 = 0.505-0.427/(1.0+exp((Y(5)*1000.0-0.29)/0.082));
dY(26, 1) = (ainf1-Y(26))/tau_adapt;
ainf2 = 0.505-0.427/(1.0+exp((Y(6)*1000.0-0.29)/0.082));
dY(27, 1) = (ainf2-Y(27))/tau_adapt;
ainf3 = 0.505-0.427/(1.0+exp((Y(7)*1000.0-0.29)/0.082));
dY(28, 1) = (ainf3-Y(28))/tau_adapt;
ainfss = 0.505-0.427/(1.0+exp((Y(9)*1000.0-0.29)/0.082));
dY(29, 1) = (ainfss-Y(29))/tau_adapt;
cinf1 = 1.0/(1.0+exp((Y(5)*1000.0-(Y(26)+0.02))/0.01));
dY(30, 1) = (cinf1-Y(30))/tau_inact;
cinf2 = 1.0/(1.0+exp((Y(6)*1000.0-(Y(27)+0.02))/0.01));
dY(31, 1) = (cinf2-Y(31))/tau_inact;
cinf3 = 1.0/(1.0+exp((Y(7)*1000.0-(Y(28)+0.02))/0.01));
dY(32, 1) = (cinf3-Y(32))/tau_inact;
cinfss = 1.0/(1.0+exp((Y(9)*1000.0-(Y(29)+0.02))/0.01));
dY(33, 1) = (cinfss-Y(33))/tau_inactss;
oinf1 = 1.0-1.0/(1.0+exp((Y(5)*1000.0-(Y(26)+0.22))/0.03));
dY(34, 1) = (oinf1-Y(34))/tau_act;
oinf2 = 1.0-1.0/(1.0+exp((Y(6)*1000.0-(Y(27)+0.22))/0.03));
dY(35, 1) = (oinf2-Y(35))/tau_act;
oinf3 = 1.0-1.0/(1.0+exp((Y(7)*1000.0-(Y(28)+0.22))/0.03));
dY(36, 1) = (oinf3-Y(36))/tau_act;
oinfss = 1.0-1.0/(1.0+exp((Y(9)*1000.0-(Y(29)+0.22))/0.03));
dY(37, 1) = (oinfss-Y(37))/tau_actss;
dY(38, 1) = 0.5*(-J_SERCASR1+J_bulkSERCA1)/Vnonjunct1;
dY(39, 1) = 0.5*(-J_SERCASR2+J_bulkSERCA2)/Vnonjunct2;
dY(40, 1) = 0.5*(-J_SERCASR3+J_bulkSERCA3)/Vnonjunct3;
dY(41, 1) = 0.5*(-J_SERCASRss+J_bulkSERCAss)/Vss;
BNa = 0.49*2.31;
JNa = DNa*Aj_nj/xj_nj_Nai*(Y(43)-Y(42))*1.0e-6;
dY(42, 1) = JNa/Vnonjunct_Nai;
betaNass = 1.0/(1.0+BNa*KdBNa/(Y(43)+KdBNa)^2.0);
i_ss = INa+INab+3.0*INaK+3.0*INaCa+IfNa;
dY(43, 1) = betaNass*(-JNa/Vss-i_ss/(Vss*F));

%===============================================================================
% End of file
%===============================================================================
