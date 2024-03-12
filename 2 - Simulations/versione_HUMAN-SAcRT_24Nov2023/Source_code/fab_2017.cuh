//------------------------------------------------------------------------------
// State variables
//------------------------------------------------------------------------------

// 0: I (dimensionless) (in Ca_SR_release)
// 1: O (dimensionless) (in Ca_SR_release)
// 2: R1 (dimensionless) (R in Ca_SR_release)
// 3: RI (dimensionless) (in Ca_SR_release)
// 4: fCMi (dimensionless) (in Ca_buffering)
// 5: fCMs (dimensionless) (in Ca_buffering)
// 6: fCQ (dimensionless) (in Ca_buffering)
// 7: fTC (dimensionless) (in Ca_buffering)
// 8: fTMC (dimensionless) (in Ca_buffering)
// 9: fTMM (dimensionless) (in Ca_buffering)
// 10: Ca_jsr (millimolar) (in Ca_dynamics)
// 11: Ca_nsr (millimolar) (in Ca_dynamics)
// 12: Ca_sub (millimolar) (in Ca_dynamics)
// 13: Cai (millimolar) (in Ca_dynamics)
// 14: V_ode (millivolt) (in Membrane)
// 15: Nai_ (millimolar) (in Nai_concentration)
// 16: dL (dimensionless) (in i_CaL_dL_gate)
// 17: fCa (dimensionless) (in i_CaL_fCa_gate)
// 18: fL (dimensionless) (in i_CaL_fL_gate)
// 19: dT (dimensionless) (in i_CaT_dT_gate)
// 20: fT (dimensionless) (in i_CaT_fT_gate)
// 21: a (dimensionless) (in i_KACh_a_gate)
// 22: paF (dimensionless) (in i_Kr_pa_gate)
// 23: paS (dimensionless) (in i_Kr_pa_gate)
// 24: piy (dimensionless) (in i_Kr_pi_gate)
// 25: n (dimensionless) (in i_Ks_n_gate)
// 26: r_Kur (dimensionless) (in i_Kur_rKur_gate)
// 27: s_Kur (dimensionless) (in i_Kur_sKur_gate)
// 28: h (dimensionless) (in i_Na_h_gate)
// 29: m (dimensionless) (in i_Na_m_gate)
// 30: y (dimensionless) (in i_f_y_gate)
// 31: q (dimensionless) (in i_to_q_gate)
// 32: r (dimensionless) (in i_to_r_gate)

#include <math.h>

__device__
void fab_2017(const double *Y,
              double *dY,
              const double time,
              const double *rand_g)
{

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

double EC50_SR;   // millimolar (in Ca_SR_release)
double HSR;   // dimensionless (in Ca_SR_release)
double MaxSR;   // dimensionless (in Ca_SR_release)
double MinSR;   // dimensionless (in Ca_SR_release)
double kiCa;   // per_millimolar_second (in Ca_SR_release)
double kim;   // per_second (in Ca_SR_release)
double koCa;   // per_millimolar2_second (in Ca_SR_release)
double kom;   // per_second (in Ca_SR_release)
double ks;   // per_second (in Ca_SR_release)
double CM_tot;   // millimolar (in Ca_buffering)
double CQ_tot;   // millimolar (in Ca_buffering)
double Mgi;   // millimolar (in Ca_buffering)
double TC_tot;   // millimolar (in Ca_buffering)
double TMC_tot;   // millimolar (in Ca_buffering)
double kb_CM;   // per_second (in Ca_buffering)
double kb_CQ;   // per_second (in Ca_buffering)
double kb_TC;   // per_second (in Ca_buffering)
double kb_TMC;   // per_second (in Ca_buffering)
double kb_TMM;   // per_second (in Ca_buffering)
double kf_CM;   // per_millimolar_second (in Ca_buffering)
double kf_CQ;   // per_millimolar_second (in Ca_buffering)
double kf_TC;   // per_millimolar_second (in Ca_buffering)
double kf_TMC;   // per_millimolar_second (in Ca_buffering)
double kf_TMM;   // per_millimolar_second (in Ca_buffering)
double K_up;   // millimolar (in Ca_intracellular_fluxes)
double P_up_basal;   // millimolar_per_second (in Ca_intracellular_fluxes)
double slope_up;   // millimolar (in Ca_intracellular_fluxes)
double tau_dif_Ca;   // second (in Ca_intracellular_fluxes)
double tau_tr;   // second (in Ca_intracellular_fluxes)
double L_cell;   // micrometre (in Cell_parameters)
double L_sub;   // micrometre (in Cell_parameters)
double R_cell;   // micrometre (in Cell_parameters)
double V_i_part;   // dimensionless (in Cell_parameters)
double V_jsr_part;   // dimensionless (in Cell_parameters)
double V_nsr_part;   // dimensionless (in Cell_parameters)
double Cao;   // millimolar (in Ionic_values)
double Ki;   // millimolar (in Ionic_values)
double Ko;   // millimolar (in Ionic_values)
double Nao;   // millimolar (in Ionic_values)
double C;   // microF (in Membrane)
double F;   // coulomb_per_mole (in Membrane)
double R2;   // joule_per_kilomole_kelvin (R in Membrane)
double T;   // kelvin (in Membrane)
double clamp_mode;   // dimensionless (in Membrane)
double Nai_clamp;   // dimensionless (in Nai_concentration)
double ACh;   // millimolar (in Rate_modulation_experiments)
double Iso_1_uM;   // dimensionless (in Rate_modulation_experiments)
double V_holding;   // millivolt (in Voltage_clamp)
double V_test;   // millivolt (in Voltage_clamp)
double t_holding;   // second (in Voltage_clamp)
double t_test;   // second (in Voltage_clamp)
double V_dL;   // millivolt (in i_CaL_dL_gate)
double k_dL;   // millivolt (in i_CaL_dL_gate)
double Km_fCa;   // millimolar (in i_CaL_fCa_gate)
double alpha_fCa;   // per_second (in i_CaL_fCa_gate)
double k_fL;   // millivolt (in i_CaL_fL_gate)
double shift_fL;   // millivolt (in i_CaL_fL_gate)
double P_CaL;   // nanoA_per_millimolar (in i_CaL)
double offset_fT;   // second (in i_CaT_fT_gate)
double P_CaT;   // nanoA_per_millimolar (in i_CaT)
double ACh_on;   // dimensionless (in i_KACh)
double g_KACh;   // microS (in i_KACh)
double g_Kr;   // microS (in i_Kr)
double g_Ks;   // microS (in i_Ks)
double g_Kur;   // microS (in i_Kur)
double K1ni;   // millimolar (in i_NaCa)
double K1no;   // millimolar (in i_NaCa)
double K2ni;   // millimolar (in i_NaCa)
double K2no;   // millimolar (in i_NaCa)
double K3ni;   // millimolar (in i_NaCa)
double K3no;   // millimolar (in i_NaCa)
double K_NaCa;   // nanoA (in i_NaCa)
double Kci;   // millimolar (in i_NaCa)
double Kcni;   // millimolar (in i_NaCa)
double Kco;   // millimolar (in i_NaCa)
double Qci;   // dimensionless (in i_NaCa)
double Qco;   // dimensionless (in i_NaCa)
double Qn;   // dimensionless (in i_NaCa)
double blockade_NaCa;   // dimensionless (in i_NaCa)
double Km_Kp;   // millimolar (in i_NaK)
double Km_Nap;   // millimolar (in i_NaK)
double i_NaK_max;   // nanoA (in i_NaK)
double delta_m;   // millivolt (in i_Na_m_gate)
double g_Na;   // microS (in i_Na)
double g_Na_L;   // microS (in i_Na)
double y_shift;   // millivolt (in i_f_y_gate)
double Km_f;   // millimolar (in i_f)
double alpha;   // dimensionless (in i_f)
double blockade;   // dimensionless (in i_f)
double g_f;   // microS (in i_f)
double g_to;   // microS (in i_to)

//------------------------------------------------------------------------------
// Computed variables
//------------------------------------------------------------------------------

//double P_tot;   // dimensionless (in Ca_SR_release)
//double diff;   // millimolar (in Ca_SR_release)
double j_SRCarel;   // millimolar_per_second (in Ca_SR_release)
double kCaSR;   // dimensionless (in Ca_SR_release)
double kiSRCa;   // per_millimolar_second (in Ca_SR_release)
double koSRCa;   // per_millimolar2_second (in Ca_SR_release)
double delta_fCMi;   // per_second (in Ca_buffering)
double delta_fCMs;   // per_second (in Ca_buffering)
double delta_fCQ;   // per_second (in Ca_buffering)
double delta_fTC;   // per_second (in Ca_buffering)
double delta_fTMC;   // per_second (in Ca_buffering)
double delta_fTMM;   // per_second (in Ca_buffering)
double P_up;   // millimolar_per_second (in Ca_intracellular_fluxes)
double b_up;   // dimensionless (in Ca_intracellular_fluxes)
double j_Ca_dif;   // millimolar_per_second (in Ca_intracellular_fluxes)
double j_tr;   // millimolar_per_second (in Ca_intracellular_fluxes)
double j_up;   // millimolar_per_second (in Ca_intracellular_fluxes)
double V_cell;   // millimetre3 (in Cell_parameters)
double V_i;   // millimetre3 (in Cell_parameters)
double V_jsr;   // millimetre3 (in Cell_parameters)
double V_nsr;   // millimetre3 (in Cell_parameters)
double V_sub;   // millimetre3 (in Cell_parameters)
//double E_Ca;   // millivolt (in Ionic_values)
double E_K;   // millivolt (in Ionic_values)
double E_Na;   // millivolt (in Ionic_values)
double RTONF;   // millivolt (in Membrane)
double V;   // millivolt (in Membrane)
double i_tot;   // nanoA (in Membrane)
double Nai;   // millimolar (in Nai_concentration)
double V_clamp;   // millivolt (in Voltage_clamp)
double Iso_shift_dL;   // millivolt (in i_CaL_dL_gate)
double Iso_slope_dL;   // dimensionless (in i_CaL_dL_gate)
double adVm;   // millivolt (in i_CaL_dL_gate)
double alpha_dL;   // per_second (in i_CaL_dL_gate)
double bdVm;   // millivolt (in i_CaL_dL_gate)
double beta_dL;   // per_second (in i_CaL_dL_gate)
double dL_infinity;   // dimensionless (in i_CaL_dL_gate)
double tau_dL;   // second (in i_CaL_dL_gate)
double fCa_infinity;   // dimensionless (in i_CaL_fCa_gate)
double tau_fCa;   // second (in i_CaL_fCa_gate)
double fL_infinity;   // dimensionless (in i_CaL_fL_gate)
double tau_fL;   // second (in i_CaL_fL_gate)
double ACh_block;   // dimensionless (in i_CaL)
double Iso_increase_1;   // dimensionless (Iso_increase in i_CaL)
double i_CaL;   // nanoA (in i_CaL)
double i_siCa;   // nanoA (in i_CaL)
double i_siK;   // nanoA (in i_CaL)
double i_siNa;   // nanoA (in i_CaL)
double dT_infinity;   // dimensionless (in i_CaT_dT_gate)
double tau_dT;   // second (in i_CaT_dT_gate)
double fT_infinity;   // dimensionless (in i_CaT_fT_gate)
double tau_fT;   // second (in i_CaT_fT_gate)
double i_CaT;   // nanoA (in i_CaT)
double a_infinity;   // dimensionless (in i_KACh_a_gate)
double alpha_a;   // per_second (in i_KACh_a_gate)
double beta_a;   // per_second (in i_KACh_a_gate)
double tau_a;   // second (in i_KACh_a_gate)
double i_KACh;   // nanoA (in i_KACh)
//double alfapaF;   // per_second (in i_Kr_pa_gate)
//double betapaF;   // per_second (in i_Kr_pa_gate)
double pa_infinity;   // dimensionless (in i_Kr_pa_gate)
double tau_paF;   // second (in i_Kr_pa_gate)
double tau_paS;   // second (in i_Kr_pa_gate)
double pi_infinity;   // dimensionless (in i_Kr_pi_gate)
double tau_pi;   // second (in i_Kr_pi_gate)
double i_Kr;   // nanoA (in i_Kr)
double Iso_shift_1;   // millivolt (Iso_shift in i_Ks_n_gate)
double alpha_n;   // per_second (in i_Ks_n_gate)
double beta_n;   // per_second (in i_Ks_n_gate)
double n_infinity;   // dimensionless (in i_Ks_n_gate)
double tau_n;   // second (in i_Ks_n_gate)
double E_Ks;   // millivolt (in i_Ks)
double i_Ks;   // nanoA (in i_Ks)
double r_Kur_infinity;   // dimensionless (in i_Kur_rKur_gate)
double tau_r_Kur;   // second (in i_Kur_rKur_gate)
double s_Kur_infinity;   // dimensionless (in i_Kur_sKur_gate)
double tau_s_Kur;   // second (in i_Kur_sKur_gate)
double i_Kur;   // nanoA (in i_Kur)
double d_i;   // dimensionless (in i_NaCa)
double d_o;   // dimensionless (in i_NaCa)
double i_NaCa;   // nanoA (in i_NaCa)
double k12;   // dimensionless (in i_NaCa)
double k14;   // dimensionless (in i_NaCa)
double k21;   // dimensionless (in i_NaCa)
double k23;   // dimensionless (in i_NaCa)
double k32;   // dimensionless (in i_NaCa)
double k34;   // dimensionless (in i_NaCa)
double k41;   // dimensionless (in i_NaCa)
double k43;   // dimensionless (in i_NaCa)
double x1;   // dimensionless (in i_NaCa)
double x2;   // dimensionless (in i_NaCa)
double x3;   // dimensionless (in i_NaCa)
double x4;   // dimensionless (in i_NaCa)
double Iso_increase_2;   // dimensionless (Iso_increase in i_NaK)
double i_NaK;   // nanoA (in i_NaK)
double alpha_h;   // per_second (in i_Na_h_gate)
double beta_h;   // per_second (in i_Na_h_gate)
double h_infinity;   // dimensionless (in i_Na_h_gate)
double tau_h;   // second (in i_Na_h_gate)
double E0_m;   // millivolt (in i_Na_m_gate)
double alpha_m;   // per_second (in i_Na_m_gate)
double beta_m;   // per_second (in i_Na_m_gate)
double m_infinity;   // dimensionless (in i_Na_m_gate)
double tau_m;   // second (in i_Na_m_gate)
double E_mh;   // millivolt (in i_Na)
double i_Na;   // nanoA (in i_Na)
double i_Na_;   // nanoA (in i_Na)
double i_Na_L;   // nanoA (in i_Na)
double ACh_shift;   // millivolt (in i_f_y_gate)
double Iso_shift_2;   // millivolt (Iso_shift in i_f_y_gate)
double tau_y;   // second (in i_f_y_gate)
double y_infinity;   // dimensionless (in i_f_y_gate)
double G_f;   // microS (in i_f)
double G_f_K;   // microS (in i_f)
double G_f_Na;   // microS (in i_f)
double g_f_K;   // microS (in i_f)
double g_f_Na;   // microS (in i_f)
double i_f;   // nanoA (in i_f)
double i_fK;   // nanoA (in i_f)
double i_fNa;   // nanoA (in i_f)
double q_infinity;   // dimensionless (in i_to_q_gate)
double tau_q;   // second (in i_to_q_gate)
double r_infinity;   // dimensionless (in i_to_r_gate)
double tau_r;   // second (in i_to_r_gate)
double i_to;   // nanoA (in i_to)
	

   //---------------------------------------------------------------------------
   // Constants
   //---------------------------------------------------------------------------

   EC50_SR = 0.45;   // millimolar (in Ca_SR_release)
   HSR = 2.5;   // dimensionless (in Ca_SR_release)
   MaxSR = 15.0;   // dimensionless (in Ca_SR_release)
   MinSR = 1.0;   // dimensionless (in Ca_SR_release)
   kiCa = 500.0;   // per_millimolar_second (in Ca_SR_release)
   kim = 5.0;   // per_second (in Ca_SR_release)
   koCa = 10000.0;   // per_millimolar2_second (in Ca_SR_release)
   kom = 660.0;   // per_second (in Ca_SR_release)
   ks = 148041085.1;   // per_second (in Ca_SR_release)
   CM_tot = 0.045;   // millimolar (in Ca_buffering)
   CQ_tot = 10.0;   // millimolar (in Ca_buffering)
   Mgi = 2.5;   // millimolar (in Ca_buffering)
   TC_tot = 0.031;   // millimolar (in Ca_buffering)
   TMC_tot = 0.062;   // millimolar (in Ca_buffering)
   kb_CM = 542.0;   // per_second (in Ca_buffering)
   kb_CQ = 445.0;   // per_second (in Ca_buffering)
   kb_TC = 446.0;   // per_second (in Ca_buffering)
   kb_TMC = 7.51;   // per_second (in Ca_buffering)
   kb_TMM = 751.0;   // per_second (in Ca_buffering)
   kf_CM = 1.642e6;   // per_millimolar_second (in Ca_buffering)
   kf_CQ = 175.4;   // per_millimolar_second (in Ca_buffering)
   kf_TC = 88800.0;   // per_millimolar_second (in Ca_buffering)
   kf_TMC = 227700.0;   // per_millimolar_second (in Ca_buffering)
   kf_TMM = 2277.0;   // per_millimolar_second (in Ca_buffering)
   K_up = 0.000286113;   // millimolar (in Ca_intracellular_fluxes)
   P_up_basal = 5.0;   // millimolar_per_second (in Ca_intracellular_fluxes)
   slope_up = 5.0e-5;   // millimolar (in Ca_intracellular_fluxes)
   tau_dif_Ca = 5.469e-5;   // second (in Ca_intracellular_fluxes)
   tau_tr = 0.04;   // second (in Ca_intracellular_fluxes)
   L_cell = 67.0;   // micrometre (in Cell_parameters)
   L_sub = 0.02;   // micrometre (in Cell_parameters)
   R_cell = 3.9;   // micrometre (in Cell_parameters)
   V_i_part = 0.46;   // dimensionless (in Cell_parameters)
   V_jsr_part = 0.0012;   // dimensionless (in Cell_parameters)
   V_nsr_part = 0.0116;   // dimensionless (in Cell_parameters)
   Cao = 1.8;   // millimolar (in Ionic_values)
   Ki = 140.0;   // millimolar (in Ionic_values)
   Ko = 5.4;   // millimolar (in Ionic_values)
   Nao = 140.0;   // millimolar (in Ionic_values)
   C = 5.7e-5;   // microF (in Membrane)
   F = 96485.3415;   // coulomb_per_mole (in Membrane)
   R2 = 8314.472;   // joule_per_kilomole_kelvin (R in Membrane)
   T = 310.0;   // kelvin (in Membrane)
   clamp_mode = 0.0;   // dimensionless (in Membrane)
   Nai_clamp = 1.0;   // dimensionless (in Nai_concentration)
   ACh = 0.0;   // millimolar (in Rate_modulation_experiments)
   Iso_1_uM = 0.0;   // dimensionless (in Rate_modulation_experiments)
   V_holding = -45.0;   // millivolt (in Voltage_clamp)
   V_test = -35.0;   // millivolt (in Voltage_clamp)
   t_holding = 0.5;   // second (in Voltage_clamp)
   t_test = 0.5;   // second (in Voltage_clamp)
   V_dL = -16.4508;   // millivolt (in i_CaL_dL_gate)
   k_dL = 4.3371;   // millivolt (in i_CaL_dL_gate)
   Km_fCa = 0.000338;   // millimolar (in i_CaL_fCa_gate)
   alpha_fCa = 0.0075;   // per_second (in i_CaL_fCa_gate)
   k_fL = 0.0;   // millivolt (in i_CaL_fL_gate)
   shift_fL = 0.0;   // millivolt (in i_CaL_fL_gate)
   offset_fT = 0.0;   // second (in i_CaT_fT_gate)
   ACh_on = 1.0;   // dimensionless (in i_KACh)
   g_KACh = 0.00345;   // microS (in i_KACh)
   //g_Kur = 0.1539e-3;   // microS (in i_Kur)
   K1ni = 395.3;   // millimolar (in i_NaCa)
   K1no = 1628.0;   // millimolar (in i_NaCa)
   K2ni = 2.289;   // millimolar (in i_NaCa)
   K2no = 561.4;   // millimolar (in i_NaCa)
   K3ni = 26.44;   // millimolar (in i_NaCa)
   K3no = 4.663;   // millimolar (in i_NaCa)
   Kci = 0.0207;   // millimolar (in i_NaCa)
   Kcni = 26.44;   // millimolar (in i_NaCa)
   Kco = 3.663;   // millimolar (in i_NaCa)
   Qci = 0.1369;   // dimensionless (in i_NaCa)
   Qco = 0.0;   // dimensionless (in i_NaCa)
   Qn = 0.4315;   // dimensionless (in i_NaCa)
   blockade_NaCa = 0.0;   // dimensionless (in i_NaCa)
   Km_Kp = 1.4;   // millimolar (in i_NaK)
   Km_Nap = 14.0;   // millimolar (in i_NaK)
   delta_m = 1.0e-5;   // millivolt (in i_Na_m_gate)
   g_Na_L = 0.0;   // microS (in i_Na)
   y_shift = 0.0;   // millivolt (in i_f_y_gate)
   Km_f = 45.0;   // millimolar (in i_f)
   alpha = 0.5927;   // dimensionless (in i_f)
   blockade = 0.0;   // dimensionless (in i_f)
   //g_to = 3.5e-3;   // microS (in i_to)

   P_CaL      = rand_g[0]; // 0.4578;   // nanoA_per_millimolar (in i_CaL)
   P_CaT      = rand_g[1]; // 0.04132;   // nanoA_per_millimolar (in i_CaT)
   g_Kr       = rand_g[2]; // 0.00424;   // microS (in i_Kr)
   K_NaCa     = rand_g[3]; // 3.343;   // nanoA (in i_NaCa)
   i_NaK_max  = rand_g[4]; // 0.08105;   // nanoA (in i_NaK)
   g_Na       = rand_g[5]; // 0.0223;   // microS (in i_Na)
   g_Ks       = rand_g[6]; // 0.00065;   // microS (in i_Ks)
   g_f        = rand_g[7]; // 0.00427;   // microS (in i_f)
   g_to       = rand_g[8];
   g_Kur      = rand_g[9];
   P_up_basal = rand_g[10];
   
   //---------------------------------------------------------------------------
   // Computed variables
   //---------------------------------------------------------------------------

   V_sub = 0.000000001*2.0*3.14159265358979*L_sub*(R_cell-L_sub/2.0)*L_cell;

   if (Iso_1_uM > 0.0)
      b_up = -0.25;
   else if (ACh > 0.0)
      b_up = 0.7*ACh/(0.00009+ACh);
   else
      b_up = 0.0;

   P_up = P_up_basal*(1.0-b_up);
   V_cell = 0.000000001*3.14159265358979*pow(R_cell, 2.0)*L_cell;
   V_nsr = V_nsr_part*V_cell;
   V_i = V_i_part*V_cell-V_sub;
   V_jsr = V_jsr_part*V_cell;
   RTONF = R2*T/F;
   k34 = Nao/(K3no+Nao);
   E_K = RTONF*log(Ko/Ki);
   G_f = g_f/(Ko/(Ko+Km_f));
   G_f_K = G_f/(alpha+1.0);
   G_f_Na = alpha*G_f_K;
   g_f_Na = G_f_Na*Ko/(Ko+Km_f);
   g_f_K = G_f_K*Ko/(Ko+Km_f);

   if (Iso_1_uM > 0.0)
      g_Ks = 1.2*g_Ks;
   else
      g_Ks = g_Ks;

   if (Iso_1_uM > 0.0)
      Iso_increase_2 = 1.2;
   else
      Iso_increase_2 = 1.0;

   ACh_block = 0.31*ACh/(ACh+0.00009);

   if (Iso_1_uM > 0.0)
      Iso_increase_1 = 1.23;
   else
      Iso_increase_1 = 1.0;

   if (Iso_1_uM > 0.0)
      Iso_shift_dL = -8.0;
   else
      Iso_shift_dL = 0.0;

   if (Iso_1_uM > 0.0)
      Iso_slope_dL = -27.0;
   else
      Iso_slope_dL = 0.0;

   alpha_a = (3.5988-0.025641)/(1.0+0.0000012155/pow(1.0*ACh, 1.6951))+0.025641;

   if (Iso_1_uM > 0.0)
      Iso_shift_1 = -14.0;
   else
      Iso_shift_1 = 0.0;

   if (ACh > 0.0)
      ACh_shift = -1.0-9.898*pow(1.0*ACh, 0.618)/(pow(1.0*ACh, 0.618)+0.00122423);
   else
      ACh_shift = 0.0;

   if (Iso_1_uM > 0.0)
      Iso_shift_2 = 7.5;
   else
      Iso_shift_2 = 0.0;

//------------------------------------------------------------------------------
// Computation
//------------------------------------------------------------------------------

   // time: time (second)

   j_SRCarel = ks*Y[1]*(Y[10]-Y[12]);
   //diff = Y[10]-Y[12];
   kCaSR = MaxSR-(MaxSR-MinSR)/(1.0+pow(EC50_SR/Y[10], HSR));
   koSRCa = koCa/kCaSR;
   kiSRCa = kiCa*kCaSR;
   dY[2] = kim*Y[3]-kiSRCa*Y[12]*Y[2]-(koSRCa*pow(Y[12], 2.0)*Y[2]-kom*Y[1]);
   dY[1] = koSRCa*pow(Y[12], 2.0)*Y[2]-kom*Y[1]-(kiSRCa*Y[12]*Y[1]-kim*Y[0]);
   dY[0] = kiSRCa*Y[12]*Y[1]-kim*Y[0]-(kom*Y[0]-koSRCa*pow(Y[12], 2.0)*Y[3]);
   dY[3] = kom*Y[0]-koSRCa*pow(Y[12], 2.0)*Y[3]-(kim*Y[3]-kiSRCa*Y[12]*Y[2]);
   //P_tot = Y[2]+Y[1]+Y[0]+Y[3];
   delta_fTC = kf_TC*Y[13]*(1.0-Y[7])-kb_TC*Y[7];
   dY[7] = delta_fTC;
   delta_fTMC = kf_TMC*Y[13]*(1.0-(Y[8]+Y[9]))-kb_TMC*Y[8];
   dY[8] = delta_fTMC;
   delta_fTMM = kf_TMM*Mgi*(1.0-(Y[8]+Y[9]))-kb_TMM*Y[9];
   dY[9] = delta_fTMM;
   delta_fCMi = kf_CM*Y[13]*(1.0-Y[4])-kb_CM*Y[4];
   dY[4] = delta_fCMi;
   delta_fCMs = kf_CM*Y[12]*(1.0-Y[5])-kb_CM*Y[5];
   dY[5] = delta_fCMs;
   delta_fCQ = kf_CQ*Y[10]*(1.0-Y[6])-kb_CQ*Y[6];
   dY[6] = delta_fCQ;
   j_Ca_dif = (Y[12]-Y[13])/tau_dif_Ca;
   j_up = P_up/(1.0+exp((-Y[13]+K_up)/slope_up));
   dY[13] = 1.0*(j_Ca_dif*V_sub-j_up*V_nsr)/V_i-(CM_tot*delta_fCMi+TC_tot*delta_fTC+TMC_tot*delta_fTMC);

   if ((time > t_holding) && (time < t_holding+t_test))
      V_clamp = V_test;
   else
      V_clamp = V_holding;

   if (clamp_mode >= 1.0)
      V = V_clamp;
   else
      V = Y[14];

   i_siCa = 2.0*P_CaL*(V-0.0)/(RTONF*(1.0-exp(-1.0*(V-0.0)*2.0/RTONF)))*(Y[12]-Cao*exp(-2.0*(V-0.0)/RTONF))*Y[16]*Y[18]*Y[17];
   i_CaT = 2.0*P_CaT*V/(RTONF*(1.0-exp(-1.0*V*2.0/RTONF)))*(Y[12]-Cao*exp(-2.0*V/RTONF))*Y[19]*Y[20];
   k32 = exp(Qn*V/(2.0*RTONF));
   Nai = Y[15];
   k43 = Nai/(K3ni+Nai);
   d_i = 1.0+Y[12]/Kci*(1.0+exp(-Qci*V/RTONF)+Nai/Kcni)+Nai/K1ni*(1.0+Nai/K2ni*(1.0+Nai/K3ni));
   k14 = Nai/K1ni*Nai/K2ni*(1.0+Nai/K3ni)*exp(Qn*V/(2.0*RTONF))/d_i;
   k12 = Y[12]/Kci*exp(-Qci*V/RTONF)/d_i;
   k41 = exp(-Qn*V/(2.0*RTONF));
   x2 = k32*k43*(k14+k12)+k41*k12*(k34+k32);
   d_o = 1.0+Cao/Kco*(1.0+exp(Qco*V/RTONF))+Nao/K1no*(1.0+Nao/K2no*(1.0+Nao/K3no));
   k21 = Cao/Kco*exp(Qco*V/RTONF)/d_o;
   k23 = Nao/K1no*Nao/K2no*(1.0+Nao/K3no)*exp(-Qn*V/(2.0*RTONF))/d_o;
   x1 = k41*k34*(k23+k21)+k21*k32*(k43+k41);
   x3 = k14*k43*(k23+k21)+k12*k23*(k43+k41);
   x4 = k23*k34*(k14+k12)+k14*k21*(k34+k32);
   i_NaCa = (1.0-blockade_NaCa)*K_NaCa*(x2*k21-x1*k12)/(x1+x2+x3+x4);
   dY[12] = j_SRCarel*V_jsr/V_sub-((i_siCa+i_CaT-2.0*i_NaCa)/(2.0*F*V_sub)+j_Ca_dif+CM_tot*delta_fCMs);
   j_tr = (Y[11]-Y[10])/tau_tr;
   dY[11] = j_up-j_tr*V_jsr/V_nsr;
   dY[10] = j_tr-(j_SRCarel+CQ_tot*delta_fCQ);
   E_Na = RTONF*log(Nao/Nai);
   //E_Ca = 0.5*RTONF*log(Cao/Y[12]);
   i_fNa = Y[30]*g_f_Na*(V-E_Na)*(1.0-blockade);
   i_fK = Y[30]*g_f_K*(V-E_K)*(1.0-blockade);
   i_f = i_fNa+i_fK;
   i_Kr = g_Kr*(V-E_K)*(0.9*Y[22]+0.1*Y[23])*Y[24];
   E_Ks = RTONF*log((Ko+0.12*Nao)/(Ki+0.12*Nai));
   i_Ks = g_Ks*(V-E_Ks)*pow(Y[25], 2.0);
   i_to = g_to*(V-E_K)*Y[31]*Y[32];
   i_NaK = Iso_increase_2*i_NaK_max*pow(1.0+pow(Km_Kp/Ko, 1.2), -1.0)*pow(1.0+pow(Km_Nap/Nai, 1.3), -1.0)*pow(1.0+exp(-(V-E_Na+110.0)/20.0), -1.0);
   E_mh = RTONF*log((Nao+0.12*Ko)/(Nai+0.12*Ki));
   i_Na_ = g_Na*pow(Y[29], 3.0)*Y[28]*(V-E_mh);
   i_Na_L = g_Na_L*pow(Y[29], 3.0)*(V-E_mh);
   i_Na = i_Na_+i_Na_L;
   i_siK = 0.000365*P_CaL*(V-0.0)/(RTONF*(1.0-exp(-1.0*(V-0.0)/RTONF)))*(Ki-Ko*exp(-1.0*(V-0.0)/RTONF))*Y[16]*Y[18]*Y[17];
   i_siNa = 0.0000185*P_CaL*(V-0.0)/(RTONF*(1.0-exp(-1.0*(V-0.0)/RTONF)))*(Nai-Nao*exp(-1.0*(V-0.0)/RTONF))*Y[16]*Y[18]*Y[17];
   i_CaL = (i_siCa+i_siK+i_siNa)*(1.0-ACh_block)*1.0*Iso_increase_1;

   if (ACh > 0.0)
      i_KACh = ACh_on*g_KACh*(V-E_K)*(1.0+exp((V+20.0)/20.0))*Y[21];
   else
      i_KACh = 0.0;

   i_Kur = g_Kur*Y[26]*Y[27]*(V-E_K);
   i_tot = i_f+i_Kr+i_Ks+i_to+i_NaK+i_NaCa+i_Na+i_CaL+i_CaT+i_KACh+i_Kur;
   dY[14] = -i_tot/C;
   dY[15] = (1.0-Nai_clamp)*-1.0*(i_Na+i_fNa+i_siNa+3.0*i_NaK+3.0*i_NaCa)/(1.0*(V_i+V_sub)*F);
   dL_infinity = 1.0/(1.0+exp(-(V-V_dL-Iso_shift_dL)/(k_dL*(1.0+Iso_slope_dL/100.0))));

   if (V == -41.8)
      adVm = -41.80001;
   else if (V == 0.0)
      adVm = 0.0;
   else if (V == -6.8)
      adVm = -6.80001;
   else
      adVm = V;

   alpha_dL = -0.02839*(adVm+41.8)/(exp(-(adVm+41.8)/2.5)-1.0)-0.0849*(adVm+6.8)/(exp(-(adVm+6.8)/4.8)-1.0);

   if (V == -1.8)
      bdVm = -1.80001;
   else
      bdVm = V;

   beta_dL = 0.01143*(bdVm+1.8)/(exp((bdVm+1.8)/2.5)-1.0);
   tau_dL = 0.001/(alpha_dL+beta_dL);
   dY[16] = (dL_infinity-Y[16])/tau_dL;
   fCa_infinity = Km_fCa/(Km_fCa+Y[12]);
   tau_fCa = 0.001*fCa_infinity/alpha_fCa;
   dY[17] = (fCa_infinity-Y[17])/tau_fCa;
   fL_infinity = 1.0/(1.0+exp((V+37.4+shift_fL)/(5.3+k_fL)));
   tau_fL = 0.001*(44.3+230.0*exp(-pow((V+36.0)/10.0, 2.0)));
   dY[18] = (fL_infinity-Y[18])/tau_fL;
   dT_infinity = 1.0/(1.0+exp(-(V+38.3)/5.5));
   tau_dT = 0.001/(1.068*exp((V+38.3)/30.0)+1.068*exp(-(V+38.3)/30.0));
   dY[19] = (dT_infinity-Y[19])/tau_dT;
   fT_infinity = 1.0/(1.0+exp((V+58.7)/3.8));
   tau_fT = 1.0/(16.67*exp(-(V+75.0)/83.3)+16.67*exp((V+75.0)/15.38))+offset_fT;
   dY[20] = (fT_infinity-Y[20])/tau_fT;
   beta_a = 10.0*exp(0.0133*(V+40.0));
   a_infinity = alpha_a/(alpha_a+beta_a);
   tau_a = 1.0/(alpha_a+beta_a);
   dY[21] = (a_infinity-Y[21])/tau_a;
   //alfapaF = 1.0/(1.0+exp(-(V+23.2)/6.6))/(0.84655354/(37.2*exp(V/11.9)+0.96*exp(-V/18.5)));
   //betapaF = 4.0*((37.2*exp(V/15.9)+0.96*exp(-V/22.5))/0.84655354-1.0/(1.0+exp(-(V+23.2)/10.6))/(0.84655354/(37.2*exp(V/15.9)+0.96*exp(-V/22.5))));
   pa_infinity = 1.0/(1.0+exp(-(V+10.0144)/7.6607));
   tau_paS = 0.84655354/(4.2*exp(V/17.0)+0.15*exp(-V/21.6));
   tau_paF = 1.0/(30.0*exp(V/10.0)+exp(-V/12.0));
   dY[23] = (pa_infinity-Y[23])/tau_paS;
   dY[22] = (pa_infinity-Y[22])/tau_paF;
   tau_pi = 1.0/(100.0*exp(-V/54.645)+656.0*exp(V/106.157));
   pi_infinity = 1.0/(1.0+exp((V+28.6)/17.1));
   dY[24] = (pi_infinity-Y[24])/tau_pi;
   n_infinity = sqrt(1.0/(1.0+exp(-(V+0.6383-Iso_shift_1)/10.7071)));
   alpha_n = 28.0/(1.0+exp(-(V-40.0-Iso_shift_1)/3.0));
   beta_n = 1.0*exp(-(V-Iso_shift_1-5.0)/25.0);
   tau_n = 1.0/(alpha_n+beta_n);
   dY[25] = (n_infinity-Y[25])/tau_n;
   r_Kur_infinity = 1.0/(1.0+exp((V+6.0)/-8.6));
   tau_r_Kur = 0.009/(1.0+exp((V+5.0)/12.0))+0.0005;
   dY[26] = (r_Kur_infinity-Y[26])/tau_r_Kur;
   s_Kur_infinity = 1.0/(1.0+exp((V+7.5)/10.0));
   tau_s_Kur = 0.59/(1.0+exp((V+60.0)/10.0))+3.05;
   dY[27] = (s_Kur_infinity-Y[27])/tau_s_Kur;
   h_infinity = 1.0/(1.0+exp((V+69.804)/4.4565));
   alpha_h = 20.0*exp(-0.125*(V+75.0));
   beta_h = 2000.0/(320.0*exp(-0.1*(V+75.0))+1.0);
   tau_h = 1.0/(alpha_h+beta_h);
   dY[28] = (h_infinity-Y[28])/tau_h;
   m_infinity = 1.0/(1.0+exp(-(V+42.0504)/8.3106));
   E0_m = V+41.0;

   if (fabs(E0_m) < delta_m)
      alpha_m = 2000.0;
   else
      alpha_m = 200.0*E0_m/(1.0-exp(-0.1*E0_m));

   beta_m = 8000.0*exp(-0.056*(V+66.0));
   tau_m = 1.0/(alpha_m+beta_m);
   dY[29] = (m_infinity-Y[29])/tau_m;
   tau_y = 1.0/(0.36*(V+148.8-ACh_shift-Iso_shift_2)/(exp(0.066*(V+148.8-ACh_shift-Iso_shift_2))-1.0)+0.1*(V+87.3-ACh_shift-Iso_shift_2)/(1.0-exp(-0.2*(V+87.3-ACh_shift-Iso_shift_2))))-0.054;

   if (V < -(80.0-ACh_shift-Iso_shift_2-y_shift))
      y_infinity = 0.01329+0.99921/(1.0+exp((V+97.134-ACh_shift-Iso_shift_2-y_shift)/8.1752));
   else
      y_infinity = 0.0002501*exp(-(V-ACh_shift-Iso_shift_2-y_shift)/12.861);

   dY[30] = (y_infinity-Y[30])/tau_y;
   q_infinity = 1.0/(1.0+exp((V+49.0)/13.0));
   tau_q = 0.001*0.6*(65.17/(0.57*exp(-0.08*(V+44.0))+0.065*exp(0.1*(V+45.93)))+10.1);
   dY[31] = (q_infinity-Y[31])/tau_q;
   r_infinity = 1.0/(1.0+exp(-(V-19.3)/15.0));
   tau_r = 0.001*0.66*1.4*(15.59/(1.037*exp(0.09*(V+30.61))+0.369*exp(-0.12*(V+23.84)))+2.98);
   dY[32] = (r_infinity-Y[32])/tau_r;
}

//==============================================================================
// End of file
//==============================================================================
