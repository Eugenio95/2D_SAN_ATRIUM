
#include <math.h>
       
__device__
void koi_2011(const double *Y,
              double *dY,
              const double time,
              const double *rand_g,
	          const int j, 
              const int WIDTH,
              const int LENGTH)
{

//------------------------------------------------------------------------------
// State variables
//------------------------------------------------------------------------------

// 0: CaSR1 (mM) (in calcium)
// 1: CaSR2 (mM) (in calcium)
// 2: CaSR3 (mM) (in calcium)
// 3: CaSR4 (mM) (in calcium)
// 4: Cai1 (mM) (in calcium)
// 5: Cai2 (mM) (in calcium)
// 6: Cai3 (mM) (in calcium)
// 7: Cai4 (mM) (in calcium)
// 8: Cass (mM) (in calcium)
// 9: d (hertz) (in ical)
// 10: f1 (hertz) (in ical)
// 11: f2 (hertz) (in ical)
// 12: fca (hertz) (in ical)
// 13: y (hertz) (in if)
// 14: pa (hertz) (in ikr)
// 15: n (hertz) (in iks)
// 16: r1 (hertz) (r in ikur)
// 17: s1 (hertz) (s in ikur)
// 18: h1 (hertz) (in ina)
// 19: h2 (hertz) (in ina)
// 20: m (hertz) (in ina)
// 21: r2 (hertz) (r in it)
// 22: s2 (hertz) (s in it)
// 23: V (mV) (in membrane)
// 24: Ki (mM) (in potassium)
// 25: a11 (hertz) (a1 in ryr)
// 26: a21 (hertz) (a2 in ryr)
// 27: a31 (hertz) (a3 in ryr)
// 28: ass1 (hertz) (ass in ryr)
// 29: c1 (hertz) (in ryr)
// 30: c2 (hertz) (in ryr)
// 31: c3 (hertz) (in ryr)
// 32: css (hertz) (in ryr)
// 33: o1 (hertz) (in ryr)
// 34: o2 (hertz) (in ryr)
// 35: o3 (hertz) (in ryr)
// 36: oss (hertz) (in ryr)
// 37: a12 (mM) (a1 in serca)
// 38: a22 (mM) (a2 in serca)
// 39: a32 (mM) (a3 in serca)
// 40: ass2 (mM) (ass in serca)
// 41: Nai (mM) (in sodium)
// 42: Nass (mM) (in sodium)

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

double BCa;   // mM (in calcium)
double CSQN;   // mM (in calcium)
double DCa;   // m2_per_s_times_1e_minus_12 (in calcium)
double DCaBm;   // m2_per_s_times_1e_minus_12 (in calcium)
double DCaSR;   // m2_per_s_times_1e_minus_12 (in calcium)
double KdBCa;   // mM (in calcium)
double KdCSQN;   // mM (in calcium)
double KdSLhigh;   // mM (in calcium)
double KdSLlow;   // mM (in calcium)
double SLhigh;   // mM (in calcium)
double SLlow;   // mM (in calcium)
double kSRleak;   // hertz (in calcium)
double Cm;   // nF (in cell)
double Vss;   // nL (in cell)
double dx;   // um (in cell)
double lcell;   // um (in cell)
double piGreco;   // dimensionless (in cell)
double rjunct;   // um (in cell)
double Cao;   // mM (in extra)
double Ko;   // mM (in extra)
double Nao;   // mM (in extra)
double gCab;   // nS (in icab)
double ECa_app;   // mV (in ical)
double gCaL;   // nS (in ical)
double ical_fca_tau;   // second (in ical)
double kCa;   // mM (in ical)
double kCan;   // dimensionless (in ical)
double ICaPmax;   // pA (in icap)
double kCaP;   // mM (in icap)
double gIf;   // nS (in if)
double gKr;   // nS (in ikr)
double gKs;   // nS (in iks)
double gNab;   // nS (in inab)
double dNaCa;   // m12_per_mol4 (in inaca)
double fCaNCX;   // dimensionless (in inaca)
double gam;   // dimensionless (in inaca)
double kNaCa;   // m12_A_per_mol4_times_1e_minus_12 (in inaca)
double INaKmax;   // pA (in inak)
double kNaKK;   // mM (in inak)
double kNaKNa;   // dimensionless (in inak)
double PNa;   // m3_per_s_times_1e_minus_12 (in ina)
double F;   // C_per_mol (in phys)
double R;   // mJ_per_mol_per_K (in phys)
double T;   // kelvin (in phys)
double tau_act;   // second (in ryr)
double tau_actss;   // second (in ryr)
double tau_adapt;   // second (in ryr)
double tau_inact;   // second (in ryr)
double tau_inactss;   // second (in ryr)
double SERCAKmf;   // mM (in serca)
double SERCAKmr;   // mM (in serca)
double cpumps;   // mM (in serca)
double k4;   // hertz (in serca)
double DNa;   // m2_per_s_times_1e_minus_12 (in sodium)
double KdBNa;   // mM (in sodium)
double duration;   // second (in stimulus)
double stim_start;   // second (in stimulus)
double stim_end;   // second (in stimulus)
double period;   // second (in stimulus)

//------------------------------------------------------------------------------
// Computed variables
//------------------------------------------------------------------------------

double JCa1;   // kat_times_1e_minus_12 (in calcium)
double JCa2;   // kat_times_1e_minus_12 (in calcium)
double JCa3;   // kat_times_1e_minus_12 (in calcium)
double JCa4;   // kat_times_1e_minus_12 (in calcium)
double JCass;   // kat_times_1e_minus_12 (in calcium)
double JSRCa1;   // kat_times_1e_minus_12 (in calcium)
double JSRCa2;   // kat_times_1e_minus_12 (in calcium)
double JSRCa3;   // kat_times_1e_minus_12 (in calcium)
double JSRCa4;   // kat_times_1e_minus_12 (in calcium)
double JSRCaleak1;   // kat_times_1e_minus_12 (in calcium)
double JSRCaleak2;   // kat_times_1e_minus_12 (in calcium)
double JSRCaleak3;   // kat_times_1e_minus_12 (in calcium)
double JSRCaleakss;   // kat_times_1e_minus_12 (in calcium)
double Jj_nj;   // kat_times_1e_minus_12 (in calcium)
double calcium_CaSR1_beta;   // dimensionless (in calcium)
double calcium_CaSR2_beta;   // dimensionless (in calcium)
double calcium_CaSR3_beta;   // dimensionless (in calcium)
double calcium_CaSR4_beta;   // dimensionless (in calcium)
double calcium_Cai1_beta;   // dimensionless (in calcium)
double calcium_Cai1_gamma;   // dimensionless (in calcium)
double calcium_Cai2_beta;   // dimensionless (in calcium)
double calcium_Cai2_gamma;   // dimensionless (in calcium)
double calcium_Cai3_beta;   // dimensionless (in calcium)
double calcium_Cai3_gamma;   // dimensionless (in calcium)
double calcium_Cai4_beta;   // dimensionless (in calcium)
double calcium_Cai4_gamma;   // dimensionless (in calcium)
double calcium_Cass_beta;   // dimensionless (in calcium)
double calcium_Cass_i_tot;   // pA (in calcium)
double Aj_nj;   // um2 (in cell)
double VSR1;   // nL (in cell)
double VSR2;   // nL (in cell)
double VSR3;   // nL (in cell)
double VSR4;   // nL (in cell)
double Vcytosol;   // nL (in cell)
double Vnonjunct1;   // nL (in cell)
double Vnonjunct2;   // nL (in cell)
double Vnonjunct3;   // nL (in cell)
double Vnonjunct4;   // nL (in cell)
double Vnonjunct_Nai;   // nL (in cell)
double xj_nj;   // um (in cell)
double xj_nj_Nai;   // um (in cell)
double ICab;   // pA (in icab)
double ICaL;   // pA (in ical)
double f_inf;   // dimensionless (in ical)
double ical_d_inf;   // dimensionless (in ical)
double ical_d_tau;   // second (in ical)
double ical_f1_tau;   // second (in ical)
double ical_f2_tau;   // second (in ical)
double ical_fca_inf;   // dimensionless (in ical)
double ICaP;   // pA (in icap)
double If;   // pA (in if)
double IfK;   // pA (in if)
double IfNa;   // pA (in if)
double if_y_inf;   // dimensionless (in if)
double if_y_tau;   // second (in if)
double IK1;   // pA (in ik1)
double gK1;   // nS (in ik1)
double IKr;   // pA (in ikr)
double ikr_pa_inf;   // dimensionless (in ikr)
double ikr_pa_tau;   // second (in ikr)
double piIKr;   // dimensionless (in ikr)
double IKs;   // pA (in iks)
double iks_n_inf;   // dimensionless (in iks)
double iks_n_tau;   // second (in iks)
double IKur;   // pA (in ikur)
double gKur;   // nS (in ikur)
double ikur_r_inf;   // dimensionless (in ikur)
double ikur_r_tau;   // second (in ikur)
double ikur_s_inf;   // dimensionless (in ikur)
double ikur_s_tau;   // second (in ikur)
double INab;   // pA (in inab)
double INaCa;   // pA (in inaca)
double INaK;   // pA (in inak)
double Nass15;   // dimensionless (in inak)
double INa;   // pA (in ina)
double h_inf;   // dimensionless (in ina)
double ina_h1_tau;   // second (in ina)
double ina_h2_tau;   // second (in ina)
double ina_m_inf;   // dimensionless (in ina)
double ina_m_tau;   // second (in ina)
double It;   // pA (in it)
double gt;   // nS (in it)
double it_r_inf;   // dimensionless (in it)
double it_r_tau;   // second (in it)
double it_s_inf;   // dimensionless (in it)
double it_s_tau;   // second (in it)
double i_ion;   // pA (in membrane)
double ECa;   // mV (in nernst)
double EK;   // mV (in nernst)
double ENa;   // mV (in nernst)
double FRT;   // per_mV (in phys)
double RTF;   // mV (in phys)
double i_tot;   // pA (in potassium)
double Jrel1;   // kat_times_1e_minus_12 (in ryr)
double Jrel2;   // kat_times_1e_minus_12 (in ryr)
double Jrel3;   // kat_times_1e_minus_12 (in ryr)
double Jrelss;   // kat_times_1e_minus_12 (in ryr)
double SRCa1;   // dimensionless (in ryr)
double SRCa2;   // dimensionless (in ryr)
double SRCa3;   // dimensionless (in ryr)
double SRCass;   // dimensionless (in ryr)
double ainf1;   // dimensionless (in ryr)
double ainf2;   // dimensionless (in ryr)
double ainf3;   // dimensionless (in ryr)
double ainfss;   // dimensionless (in ryr)
double cinf1;   // dimensionless (in ryr)
double cinf2;   // dimensionless (in ryr)
double cinf3;   // dimensionless (in ryr)
double cinfss;   // dimensionless (in ryr)
double nu1;   // m3_per_s_times_1e_minus_12 (in ryr)
double nu2;   // m3_per_s_times_1e_minus_12 (in ryr)
double nu3;   // m3_per_s_times_1e_minus_12 (in ryr)
double nuss;   // m3_per_s_times_1e_minus_12 (in ryr)
double oinf1;   // dimensionless (in ryr)
double oinf2;   // dimensionless (in ryr)
double oinf3;   // dimensionless (in ryr)
double oinfss;   // dimensionless (in ryr)
double J_SERCASR1;   // kat_times_1e_minus_12 (in serca)
double J_SERCASR2;   // kat_times_1e_minus_12 (in serca)
double J_SERCASR3;   // kat_times_1e_minus_12 (in serca)
double J_SERCASRss;   // kat_times_1e_minus_12 (in serca)
double J_bulkSERCA1;   // kat_times_1e_minus_12 (in serca)
double J_bulkSERCA2;   // kat_times_1e_minus_12 (in serca)
double J_bulkSERCA3;   // kat_times_1e_minus_12 (in serca)
double J_bulkSERCAss;   // kat_times_1e_minus_12 (in serca)
double k1;   // m6_per_s_per_mol2 (in serca)
double k2;   // hertz (in serca)
double k3;   // m6_per_s_per_mol2 (in serca)
double BNa;   // mM (in sodium)
double JNa;   // kat_times_1e_minus_12 (in sodium)
double betaNass;   // dimensionless (in sodium)
double i_ss;   // pA (in sodium)
double amplitude;   // pA (in stimulus)
double i_stim;   // pA (in stimulus)

//------------------------------------------------------------------------------
// Initialisation
//------------------------------------------------------------------------------

   //---------------------------------------------------------------------------
   // Constants
   //---------------------------------------------------------------------------

   BCa = 0.024;   // mM (in calcium)
   CSQN = 6.7;   // mM (in calcium)
   DCa = 780.0;   // m2_per_s_times_1e_minus_12 (in calcium)
   DCaBm = 25.0;   // m2_per_s_times_1e_minus_12 (in calcium)
   DCaSR = 44.0;   // m2_per_s_times_1e_minus_12 (in calcium)
   KdBCa = 0.00238;   // mM (in calcium)
   KdCSQN = 0.8;   // mM (in calcium)
   KdSLhigh = 0.013;   // mM (in calcium)
   KdSLlow = 1.1;   // mM (in calcium)
   SLhigh = 13.0;   // mM (in calcium)
   SLlow = 165.0;   // mM (in calcium)
   kSRleak = 0.006;   // hertz (in calcium)
   Cm = 0.05;   // nF (in cell)
   Vss = 4.99231999999999966e-5;   // nL (in cell)
   dx = 1.625;   // um (in cell)
   lcell = 122.051;   // um (in cell)
   piGreco = 3.14159265358979312e0;   // dimensionless (in cell)
   rjunct = 6.5;   // um (in cell)
   Cao = 1.8;   // mM (in extra)
   Ko = 5.4;   // mM (in extra)
   Nao = 130.0;   // mM (in extra)
   gCab = 0.0952;   // nS (in icab)
   ECa_app = 60.0;   // mV (in ical)
   ical_fca_tau = 0.002;   // second (in ical)
   kCa = 0.001;   // mM (in ical)
   kCan = 2.0;   // dimensionless (in ical)
   kCaP = 0.0005;   // mM (in icap)
   gNab = 0.060599;   // nS (in inab)
   dNaCa = 0.0003;   // m12_per_mol4 (in inaca)
   fCaNCX = 1.0;   // dimensionless (in inaca)
   gam = 0.45;   // dimensionless (in inaca)
   kNaKK = 1.0;   // mM (in inak)
   kNaKNa = 11.0;   // dimensionless (in inak)
   F = 96487.0;   // C_per_mol (in phys)
   R = 8314.0;   // mJ_per_mol_per_K (in phys)
   T = 306.15;   // kelvin (in phys)
   tau_act = 0.01875;   // second (in ryr)
   tau_actss = 0.005;   // second (in ryr)
   tau_adapt = 1.0;   // second (in ryr)
   tau_inact = 0.0875;   // second (in ryr)
   tau_inactss = 0.015;   // second (in ryr)
   SERCAKmf = 0.00025;   // mM (in serca)
   SERCAKmr = 1.8;   // mM (in serca)
   k4 = 7.5;   // hertz (in serca)
   DNa = 0.12;   // m2_per_s_times_1e_minus_12 (in sodium)
   KdBNa = 10.0;   // mM (in sodium)

   duration   = 0.001;   // second (in stimulus)
   stim_start = 0.01;   // second (in stimulus)
   stim_end   = 2.0;
   period     = 0.5;   // second (in stimulus)

   PNa     = rand_g[0]; // *  2; // 0.0018;   // m3_per_s_times_1e_minus_12 (in ina)
   gCaL    = rand_g[1]; //  25.3125;   // nS (in ical)
   gIf     = rand_g[2]; //  1.0;   // nS (in if)
   gK1     = rand_g[3]; // * 0.75; 
   gKr     = rand_g[4]; // * 0.5; //  0.5;   // nS (in ikr)
   gKs     = rand_g[5]; //  1.0;   // nS (in iks)
   INaKmax = rand_g[6]; //  70.8253;   // pA (in inak)
   ICaPmax = rand_g[7]; //  2.0;   // pA (in icap)
   kNaCa   = rand_g[8]; //  0.0084;   // m12_A_per_mol4_times_1e_minus_12 (in inaca)
   cpumps  = rand_g[9]; //  0.04;   // mM (in serca)
   gt      = rand_g[10]; // * 0.5;
   gKur    = rand_g[11];

   const int idx_x = (int) j/WIDTH + 1;
   const int idx_y = j%LENGTH;

   //---------------------------------------------------------------------------
   // Computed variables
   //---------------------------------------------------------------------------

   k3 = k4/pow(SERCAKmr, 2.0);
   Vnonjunct1 = (pow(1.0*dx, 2.0)-pow(0.0*dx, 2.0))*piGreco*lcell*0.5*1.0e-6;
   nu1 = 1.0*Vnonjunct1;
   VSR1 = 0.05*Vnonjunct1/2.0*0.9;
   Vnonjunct2 = (pow(2.0*dx, 2.0)-pow(1.0*dx, 2.0))*piGreco*lcell*0.5*1.0e-6;
   nu2 = 1.0*Vnonjunct2;
   VSR2 = 0.05*Vnonjunct2/2.0*0.9;
   Vnonjunct3 = (pow(3.0*dx, 2.0)-pow(2.0*dx, 2.0))*piGreco*lcell*0.5*1.0e-6;
   nu3 = 1.0*Vnonjunct3;
   VSR3 = 0.05*Vnonjunct3/2.0*0.9;
   nuss = 625.0*Vss;
   Vnonjunct4 = (pow(4.0*dx, 2.0)-pow(3.0*dx, 2.0))*piGreco*lcell*0.5*1.0e-6;
   VSR4 = 0.05*Vnonjunct4/2.0*0.9;
   k1 = pow(1000.0, 2.0)*k4;
   k2 = k1*pow(SERCAKmf, 2.0);
   Aj_nj = piGreco*rjunct*2.0*lcell*0.5;
   xj_nj = 0.02/2.0+dx/2.0;
   RTF = R*T/F;
   FRT = 1.0/RTF;
   Vnonjunct_Nai = Vnonjunct1+Vnonjunct2+Vnonjunct3+Vnonjunct4;
   Vcytosol = Vnonjunct_Nai+Vss;
   xj_nj_Nai = 0.02/2.0+2.0*dx;
   //gK1 = 3.825*0.9;
   //gKur = 0.89*2.75;
   //gt = 1.09*7.5;
   amplitude = -2500.0*2;
   BNa = 0.49*2.31;


//------------------------------------------------------------------------------
// Computation
//------------------------------------------------------------------------------

   // time: time (second)

   calcium_CaSR1_beta = 1.0/(1.0+CSQN*KdCSQN/pow(Y[0]+KdCSQN, 2.0));
   J_SERCASR1 = (-k3*pow(Y[0], 2.0)*(cpumps-Y[37])+k4*Y[37])*Vnonjunct1*2.0;
   JSRCaleak1 = kSRleak*(Y[0]-Y[4])*Vnonjunct1;
   SRCa1 = 1.0-1.0/(1.0+exp((Y[0]-0.3)/0.1));
   Jrel1 = nu1*Y[33]*Y[29]*SRCa1*(Y[0]-Y[4]);
   JSRCa1 = J_SERCASR1-JSRCaleak1-Jrel1;
   dY[0] = calcium_CaSR1_beta*DCaSR*((Y[1]-2.0*Y[0]+Y[0])/pow(dx, 2.0)+(Y[1]-Y[0])/(2.0*1.0*pow(dx, 2.0)))+JSRCa1/VSR1*calcium_CaSR1_beta;
   calcium_CaSR2_beta = 1.0/(1.0+CSQN*KdCSQN/pow(Y[1]+KdCSQN, 2.0));
   J_SERCASR2 = (-k3*pow(Y[1], 2.0)*(cpumps-Y[38])+k4*Y[38])*Vnonjunct2*2.0;
   JSRCaleak2 = kSRleak*(Y[1]-Y[5])*Vnonjunct2;
   SRCa2 = 1.0-1.0/(1.0+exp((Y[1]-0.3)/0.1));
   Jrel2 = nu2*Y[34]*Y[30]*SRCa2*(Y[1]-Y[5]);
   JSRCa2 = J_SERCASR2-JSRCaleak2-Jrel2;
   dY[1] = calcium_CaSR2_beta*DCaSR*((Y[2]-2.0*Y[1]+Y[0])/pow(dx, 2.0)+(Y[2]-Y[0])/(2.0*2.0*pow(dx, 2.0)))+JSRCa2/VSR2*calcium_CaSR2_beta;
   calcium_CaSR3_beta = 1.0/(1.0+CSQN*KdCSQN/pow(Y[2]+KdCSQN, 2.0));
   J_SERCASR3 = (-k3*pow(Y[2], 2.0)*(cpumps-Y[39])+k4*Y[39])*Vnonjunct3*2.0;
   JSRCaleak3 = kSRleak*(Y[2]-Y[6])*Vnonjunct3;
   SRCa3 = 1.0-1.0/(1.0+exp((Y[2]-0.3)/0.1));
   Jrel3 = nu3*Y[35]*Y[31]*SRCa3*(Y[2]-Y[6]);
   JSRCa3 = J_SERCASR3-JSRCaleak3-Jrel3;
   dY[2] = calcium_CaSR3_beta*DCaSR*((Y[3]-2.0*Y[2]+Y[1])/pow(dx, 2.0)+(Y[3]-Y[1])/(2.0*3.0*pow(dx, 2.0)))+JSRCa3/VSR3*calcium_CaSR3_beta;
   calcium_CaSR4_beta = 1.0/(1.0+CSQN*KdCSQN/pow(Y[3]+KdCSQN, 2.0));
   J_SERCASRss = (-k3*pow(Y[3], 2.0)*(cpumps-Y[40])+k4*Y[40])*Vss*2.0;
   JSRCaleakss = kSRleak*(Y[3]-Y[8])*Vss;
   SRCass = 1.0-1.0/(1.0+exp((Y[3]-0.3)/0.1));
   Jrelss = nuss*Y[36]*Y[32]*SRCass*(Y[3]-Y[8]);
   JSRCa4 = J_SERCASRss-JSRCaleakss-Jrelss;
   dY[3] = calcium_CaSR4_beta*DCaSR*((Y[3]-2.0*Y[3]+Y[2])/pow(dx, 2.0)+(Y[3]-Y[2])/(2.0*4.0*pow(dx, 2.0)))+JSRCa4/VSR4*calcium_CaSR4_beta;
   calcium_Cai1_gamma = BCa*KdBCa/pow(Y[4]+KdBCa, 2.0);
   calcium_Cai1_beta = 1.0/(1.0+calcium_Cai1_gamma);
   J_bulkSERCA1 = (k1*pow(Y[4], 2.0)*(cpumps-Y[37])-k2*Y[37])*Vnonjunct1*2.0;
   JCa1 = -J_bulkSERCA1+JSRCaleak1+Jrel1;
   dY[4] = calcium_Cai1_beta*(DCa+calcium_Cai1_gamma*DCaBm)*((Y[5]-2.0*Y[4]+Y[4])/pow(dx, 2.0)+(Y[5]-Y[4])/(2.0*1.0*pow(dx, 2.0)))-2.0*calcium_Cai1_beta*calcium_Cai1_gamma*DCaBm/(KdBCa+Y[4])*pow((Y[5]-Y[4])/(2.0*dx), 2.0)+JCa1/Vnonjunct1*calcium_Cai1_beta;
   calcium_Cai2_gamma = BCa*KdBCa/pow(Y[5]+KdBCa, 2.0);
   calcium_Cai2_beta = 1.0/(1.0+calcium_Cai2_gamma);
   J_bulkSERCA2 = (k1*pow(Y[5], 2.0)*(cpumps-Y[38])-k2*Y[38])*Vnonjunct2*2.0;
   JCa2 = -J_bulkSERCA2+JSRCaleak2+Jrel2;
   dY[5] = calcium_Cai2_beta*(DCa+calcium_Cai2_gamma*DCaBm)*((Y[6]-2.0*Y[5]+Y[4])/pow(dx, 2.0)+(Y[6]-Y[4])/(2.0*2.0*pow(dx, 2.0)))-2.0*calcium_Cai2_beta*calcium_Cai2_gamma*DCaBm/(KdBCa+Y[5])*pow((Y[6]-Y[4])/(2.0*dx), 2.0)+JCa2/Vnonjunct2*calcium_Cai2_beta;
   calcium_Cai3_gamma = BCa*KdBCa/pow(Y[6]+KdBCa, 2.0);
   calcium_Cai3_beta = 1.0/(1.0+calcium_Cai3_gamma);
   J_bulkSERCA3 = (k1*pow(Y[6], 2.0)*(cpumps-Y[39])-k2*Y[39])*Vnonjunct3*2.0;
   JCa3 = -J_bulkSERCA3+JSRCaleak3+Jrel3;
   dY[6] = calcium_Cai3_beta*(DCa+calcium_Cai3_gamma*DCaBm)*((Y[7]-2.0*Y[6]+Y[5])/pow(dx, 2.0)+(Y[7]-Y[5])/(2.0*3.0*pow(dx, 2.0)))-2.0*calcium_Cai3_beta*calcium_Cai3_gamma*DCaBm/(KdBCa+Y[6])*pow((Y[7]-Y[5])/(2.0*dx), 2.0)+JCa3/Vnonjunct3*calcium_Cai3_beta;
   calcium_Cai4_gamma = BCa*KdBCa/pow(Y[7]+KdBCa, 2.0);
   calcium_Cai4_beta = 1.0/(1.0+calcium_Cai4_gamma);
   Jj_nj = DCa*Aj_nj/xj_nj*(Y[8]-Y[7])*1.0e-6;
   JCa4 = Jj_nj;
   dY[7] = calcium_Cai4_beta*(DCa+calcium_Cai4_gamma*DCaBm)*((Y[7]-2.0*Y[7]+Y[6])/pow(dx, 2.0)+(Y[7]-Y[6])/(2.0*4.0*pow(dx, 2.0)))-2.0*calcium_Cai4_beta*calcium_Cai4_gamma*DCaBm/(KdBCa+Y[7])*pow((Y[7]-Y[6])/(2.0*dx), 2.0)+JCa4/Vnonjunct4*calcium_Cai4_beta;
   calcium_Cass_beta = 1.0/(1.0+SLlow*KdSLlow/pow(Y[8]+KdSLlow, 2.0)+SLhigh*KdSLhigh/pow(Y[8]+KdSLhigh, 2.0)+BCa*KdBCa/pow(Y[8]+KdBCa, 2.0));
   J_bulkSERCAss = (k1*pow(Y[8], 2.0)*(cpumps-Y[40])-k2*Y[40])*Vss*2.0;
   JCass = -Jj_nj+JSRCaleakss-J_bulkSERCAss+Jrelss;
   ICaL = gCaL*Y[9]*Y[12]*Y[10]*Y[11]*(Y[23]-ECa_app);
   ECa = RTF*log(Cao/Y[8])/2.0;
   ICab = gCab*(Y[23]-ECa);
   ICaP = ICaPmax*Y[8]/(kCaP+Y[8]);
   INaCa = kNaCa*(exp(gam*Y[23]*FRT)*pow(Y[42], 3.0)*Cao-exp((gam-1.0)*Y[23]*FRT)*pow(Nao, 3.0)*Y[8]*fCaNCX)/(1.0+dNaCa*(pow(Nao, 3.0)*Y[8]*fCaNCX+pow(Y[42], 3.0)*Cao));
   calcium_Cass_i_tot = -ICaL-ICab-ICaP+2.0*INaCa;
   dY[8] = calcium_Cass_beta*(JCass/Vss+calcium_Cass_i_tot/(2.0*Vss*F));
   ical_d_inf = 1.0/(1.0+exp((Y[23]+9.0)/-5.8));
   ical_d_tau = 0.0027*exp(-pow((Y[23]+35.0)/30.0, 2.0))+0.002;
   dY[9] = (ical_d_inf-Y[9])/ical_d_tau;
   f_inf = 1.0/(1.0+exp((Y[23]+27.4)/7.1));
   ical_f1_tau = 0.98698*exp(-pow((Y[23]+30.16047)/7.09396, 2.0))+0.04275/(1.0+exp((Y[23]-51.61555)/-80.61331))+0.03576/(1.0+exp((Y[23]+29.57272)/13.21758))-0.00821;
   dY[10] = (f_inf-Y[10])/ical_f1_tau;
   ical_f2_tau = 1.3323*exp(-pow((Y[23]+40.0)/14.2, 2.0))+0.0626;
   dY[11] = (f_inf-Y[11])/ical_f2_tau;
   ical_fca_inf = 1.0-1.0/(1.0+pow(kCa/Y[8], kCan));
   dY[12] = (ical_fca_inf-Y[12])/ical_fca_tau;
   EK = RTF*log(Ko/Y[24]);
   IfK = gIf*Y[13]*(1.0-0.2677)*(Y[23]-EK);
   ENa = RTF*log(Nao/Y[42]);
   IfNa = gIf*Y[13]*0.2677*(Y[23]-ENa);
   If = IfK+IfNa;
   if_y_inf = 1.0/(1.0+exp((Y[23]+97.82874)/12.48025));
   if_y_tau = 1.0/(0.00332*exp(-Y[23]/16.54103)+23.71839*exp(Y[23]/16.54103));
   dY[13] = (if_y_inf-Y[13])/if_y_tau;
   IK1 = gK1*pow(Ko*1.0, 0.4457)*(Y[23]-EK)/(1.0+exp(1.5*(Y[23]-EK+3.6)*FRT));
   piIKr = 1.0/(1.0+exp((Y[23]+55.0)/24.0));
   IKr = gKr*Y[14]*piIKr*(Y[23]-EK);
   ikr_pa_inf = 1.0/(1.0+exp((Y[23]+15.0)/-6.0));
   ikr_pa_tau = 0.21718*exp(-pow((Y[23]+20.1376)/22.1996, 2.0))+0.03118;
   dY[14] = (ikr_pa_inf-Y[14])/ikr_pa_tau;
   IKs = gKs*Y[15]*(Y[23]-EK);
   iks_n_inf = 1.0/(1.0+exp((Y[23]-19.9)/-12.7));
   iks_n_tau = 0.4*exp(-pow((Y[23]-20.0)/20.0, 2.0))+0.7;
   dY[15] = (iks_n_inf-Y[15])/iks_n_tau;
   IKur = gKur*Y[16]*Y[17]*(Y[23]-EK);
   ikur_r_inf = 1.0/(1.0+exp((Y[23]+6.0)/-8.6));
   ikur_r_tau = 0.009/(1.0+exp((Y[23]+5.0)/12.0))+0.0005;
   ikur_s_inf = 1.0/(1.0+exp((Y[23]+7.5)/10.0));
   ikur_s_tau = 0.59/(1.0+exp((Y[23]+60.0)/10.0))+3.05;
   dY[16] = (ikur_r_inf-Y[16])/ikur_r_tau;
   dY[17] = (ikur_s_inf-Y[17])/ikur_s_tau;
   INa = PNa*pow(Y[20], 3.0)*(0.9*Y[18]+0.1*Y[19])*Nao*Y[23]*F*FRT*(exp((Y[23]-ENa)*FRT)-1.0)/(exp(Y[23]*FRT)-1.0);
   h_inf = 1.0/(1.0+exp((Y[23]+63.6)/5.3));
   ina_h1_tau = 0.03/(1.0+exp((Y[23]+35.1)/3.2))+0.0003;
   dY[18] = (h_inf-Y[18])/ina_h1_tau;
   ina_h2_tau = 0.12/(1.0+exp((Y[23]+35.1)/3.2))+0.003;
   dY[19] = (h_inf-Y[19])/ina_h2_tau;
   ina_m_inf = 1.0/(1.0+exp((Y[23]+27.12)/-8.21));
   ina_m_tau = 4.2e-5*exp(-pow((Y[23]+25.57)/28.8, 2.0))+2.4e-5;
   dY[20] = (ina_m_inf-Y[20])/ina_m_tau;
   INab = gNab*(Y[23]-ENa);
   Nass15 = pow(Y[42]*1.0, 1.5);
   INaK = INaKmax*Ko/(Ko+kNaKK)*Nass15/(Nass15+pow(kNaKNa, 1.5))*(Y[23]+150.0)/(Y[23]+200.0);
   It = gt*Y[21]*Y[22]*(Y[23]-EK);
   it_r_inf = 1.0/(1.0+exp((Y[23]-1.0)/-11.0));
   it_r_tau = 0.0035*exp(-pow((Y[23]+0.0)/30.0, 2.0))+0.0015;
   it_s_inf = 1.0/(1.0+exp((Y[23]+40.5)/11.5));
   it_s_tau = 0.025635*exp(-pow((Y[23]+52.45)/15.8827, 2.0))+0.01414;
   dY[21] = (it_r_inf-Y[21])/it_r_tau;
   dY[22] = (it_s_inf-Y[22])/it_s_tau;
   i_ion = INa+ICaL+It+IKur+IK1+IKr+IKs+INab+ICab+INaK+ICaP+INaCa+If;

   if ((idx_x > 157) && (idx_x < 173) && (idx_y > 37) && (idx_y < 43) && (time >= stim_start) && (time <= stim_end)  && time-stim_start-period*floor((time-stim_start)/period) < duration)
      i_stim = amplitude;
   else
      i_stim = 0.0;

   dY[23] = -(i_ion+i_stim)/Cm;
   i_tot = It+IKur+IK1+IKr+IKs-2.0*INaK+IfK+i_stim;
   dY[24] = -i_tot/(Vcytosol*F);
   ainf1 = 0.505-0.427/(1.0+exp((Y[4]*1000.0-0.29)/0.082));
   dY[25] = (ainf1-Y[25])/tau_adapt;
   ainf2 = 0.505-0.427/(1.0+exp((Y[5]*1000.0-0.29)/0.082));
   dY[26] = (ainf2-Y[26])/tau_adapt;
   ainf3 = 0.505-0.427/(1.0+exp((Y[6]*1000.0-0.29)/0.082));
   dY[27] = (ainf3-Y[27])/tau_adapt;
   ainfss = 0.505-0.427/(1.0+exp((Y[8]*1000.0-0.29)/0.082));
   dY[28] = (ainfss-Y[28])/tau_adapt;
   cinf1 = 1.0/(1.0+exp((Y[4]*1000.0-(Y[25]+0.02))/0.01));
   dY[29] = (cinf1-Y[29])/tau_inact;
   cinf2 = 1.0/(1.0+exp((Y[5]*1000.0-(Y[26]+0.02))/0.01));
   dY[30] = (cinf2-Y[30])/tau_inact;
   cinf3 = 1.0/(1.0+exp((Y[6]*1000.0-(Y[27]+0.02))/0.01));
   dY[31] = (cinf3-Y[31])/tau_inact;
   cinfss = 1.0/(1.0+exp((Y[8]*1000.0-(Y[28]+0.02))/0.01));
   dY[32] = (cinfss-Y[32])/tau_inactss;
   oinf1 = 1.0-1.0/(1.0+exp((Y[4]*1000.0-(Y[25]+0.22))/0.03));
   dY[33] = (oinf1-Y[33])/tau_act;
   oinf2 = 1.0-1.0/(1.0+exp((Y[5]*1000.0-(Y[26]+0.22))/0.03));
   dY[34] = (oinf2-Y[34])/tau_act;
   oinf3 = 1.0-1.0/(1.0+exp((Y[6]*1000.0-(Y[27]+0.22))/0.03));
   dY[35] = (oinf3-Y[35])/tau_act;
   oinfss = 1.0-1.0/(1.0+exp((Y[8]*1000.0-(Y[28]+0.22))/0.03));
   dY[36] = (oinfss-Y[36])/tau_actss;
   dY[37] = 0.5*(-J_SERCASR1+J_bulkSERCA1)/Vnonjunct1;
   dY[38] = 0.5*(-J_SERCASR2+J_bulkSERCA2)/Vnonjunct2;
   dY[39] = 0.5*(-J_SERCASR3+J_bulkSERCA3)/Vnonjunct3;
   dY[40] = 0.5*(-J_SERCASRss+J_bulkSERCAss)/Vss;
   JNa = DNa*Aj_nj/xj_nj_Nai*(Y[42]-Y[41])*1.0e-6;
   dY[41] = JNa/Vnonjunct_Nai;
   betaNass = 1.0/(1.0+BNa*KdBNa/pow(Y[42]+KdBNa, 2.0));
   i_ss = INa+INab+3.0*INaK+3.0*INaCa+IfNa;
   dY[42] = betaNass*(-JNa/Vss-i_ss/(Vss*F));
}

//==============================================================================
// End of file
//==============================================================================
