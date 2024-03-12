
#include <math.h>

__device__
        void mbs_2023(const double *Y,
                      double *dY,
                      const double time,
                      const double *rand_g,
	                  const int j, 
                      const int WIDTH,
                      const int LENGTH)
{

  double ICaL;
  double ICaP;
  double ICab;
  double INa;
  double INaCa;
  double INaK;
  double INaK_tmp;
  double INa_tmp;
  double INab;
  double INahinf;
  double IfNa;
  double JNa;
  double JSRCaleak1bc;
  double JSRCaleakss;
  double J_SERCASRbc;
  double J_SERCASRss;
  double J_bulkSERCAbc;
  double J_bulkSERCAss;
  double Jrelbc;
  double Jrelss;
  double Jss_bc;
  double RyRcinfbc;
  double RyRcinfss;
  double RyRtauinactss;
  double a;
  double b_a;
  double betaSRbc;
  double betaSRss;
  double c_a;
  double d_a;
  double e_a;
  double f_a;
  double g_a;
  double h_a;
  double Istim;

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*                state variables                         % */
  /* EC index                                             EP only index     */
  /* 1-4:  y_XB              */
  /* 5-20: y_RU             T_{i-1},T_i,T_{i+1},C_i */
  /* 21:   y_B_nso */
  /* 22: V                                                  1 */
  /* 23: m     INa                                          2    */
  /* 24: h1    INa                                          3 */
  /* 25: h2    INa                                          4 */
  /* 26: d     ICaL                                         5    */
  /* 27: f1    ICaL not in use */
  /* 28: f2    ICaL                                         6 */
  /* 29: fca   ICaL                                         7    */
  /* 30: r     Ito                                          8        */
  /* 31: s     Ito                                          9 */
  /* 32: susr  IKur                                         10 */
  /* 33: suss  IKur                                         11 */
  /* 34: n     IKs                                          12   */
  /* 35: pa    IKr                                          13 */
  /* 36: y     If                                           14   */
  /* 37: o     RyRss                                        15 */
  /* 38: c     RyRss                                        16         */
  /* 39: a     RyRss                                        17             */
  /* 40: o     RyRbc                                        18   */
  /* 41: c     RyRbc                                        19 */
  /* 42: a     RyRbc                                        20 */
  /* 43: SERCAbc                                            21 */
  /* 44: SERCAss                                            22 */
  /* 45: Nass                                               23 */
  /* 46: Nai                                                24     */
  /* 47: Ki   not in use                                     */
  /* 48: Cass                                               25 */
  /* 49: Cabc                                               26     */
  /* 50: CaSRbc                                             27     */
  /* 51: CaSRss                                             28     */
  /*  input */
          
  /* ********************************************** */
  /*  Stimulus */
          /* ********************************************** */
  const int idx_x = (int) j/WIDTH + 1;
  const int idx_y = j%LENGTH;
  double duration  = 0.005;
  double stim_start = 0.5;
  double stim_end  = 1.5;
  double period    = 0.5;
  double amplitude = 920.0*2.0;
  if ((idx_x > 157) && (idx_x < 173) && (idx_y > 37) && (idx_y < 53) && (time >= stim_start) && (time <= stim_end) && (time-stim_start-period*floor((time-stim_start)/period) < duration) ) {
    Istim = -amplitude;
  } else {
    Istim = 0.0;
  }

  /* ********************************************** */
  /*  Maximal conductances */
  /* ********************************************** */
 double gNa     = rand_g[0];
 double gCaL    = rand_g[1];
 double gt      = rand_g[2];
 double gKur    = rand_g[3];
 double gK1     = rand_g[4];
 double gKr     = rand_g[5];
 double gKs     = rand_g[6];
 double INaKmax = rand_g[7];
 double gCap    = rand_g[8];
 double kNaCa   = rand_g[9];
 double gf      = rand_g[10];
 double Pserca  = rand_g[11];

 double gNab = 0.0606;
 double gCab = 0.0850;

  /* ********************************************** */
  /*  Analytical equations */
  /* ********************************************** */

  /*  INa Skibsbye2016 */
  INa_tmp = Y[0] - 26.724710064568281 * log(130.0 / Y[22]);
  INa = gNa * pow(Y[1], 3.0) * (Y[2] * Y[3]) * INa_tmp;
  INahinf = 1.0 / (exp((Y[0] + 67.0) / 5.6) + 1.0);
  a = (Y[0] + 48.0) / 15.0;
  /*  ICaL  */
  b_a = Y[24] / 0.0006;
  ICaL = gCaL* Y[4] * Y[6] * Y[5] * (Y[0] - 60.0);
  c_a = (Y[0] + 35.0) / 30.0;
  d_a = (Y[0] + 40.0) / 14.2;
  /*  It */
  e_a = Y[0] / 30.0;
  f_a = (Y[0] + 52.45) / 15.8827;
  /*  Maleckar et al. */
  /*  IKur */
  /*  Maleckar et al. */
  /*  Maleckar et al. */
  /*  Maleckar et al. */
  /*  Maleckar et al. */
  /*  IKs */
  g_a = (Y[0] - 20.0) / 20.0;
  /*  IKr */
  h_a = (Y[0] + 20.1376) / 22.1996;
  /*  IK1 */
  /*  Background leaks */
  INab = gNab * INa_tmp;
  ICab = gCab * (Y[0] - 13.362355032284141 * log(1.8 / Y[24]));
  /*  INaK */
  INaK_tmp = pow(Y[22], 1.5);
  INaK = INaKmax * 0.8438 * INaK_tmp / (INaK_tmp + 36.4828726939094) *
         (Y[0] + 150.0) / (Y[0] + 200.0);
  /*  INaCa\ */
  INaK_tmp = pow(Y[22], 3.0);
  INaCa = kNaCa *
          ((exp(0.45 * Y[0] * 96487.0 / 8314.0 / 310.15) * INaK_tmp * 1.8 -
            exp(-0.55 * Y[0] * 96487.0 / 8314.0 / 310.15) * 2.197E+6 * Y[24]) /
           (0.0003 * (2.197E+6 * Y[24] + INaK_tmp * 1.8) + 1.0));
  /*  ICaP */
  ICaP =gCap * Y[24] / (Y[24] + 0.0005);
  /*  If, Zorn-Pauly LAW fit */
  IfNa = gf * Y[13] * (0.2677 * INa_tmp);
  /*  Ca buffers */
  betaSRbc = 1.0 / (5.36 / ((Y[26] + 0.8) * (Y[26] + 0.8)) + 1.0);
  betaSRss = 1.0 / (5.36 / ((Y[27] + 0.8) * (Y[27] + 0.8)) + 1.0);
  /*   Diffusion from junct to non-junct */
  Jss_bc = 2.1748606163574806E+6 * (Y[24] - Y[25]) * 1.0E-6;
  /*  Naflux in 1 nl volume */
  JNa = 334.593940978074 * (Y[22] - Y[23]) * 1.0E-6;
  /*  SERCA fluxes */
  J_SERCASRbc =
      (-4.7391448307935731 * (Y[26] * Y[26]) * (0.04 - Y[20]) + Pserca * Y[20]) *
      0.0066366144807084373 * 2.0;
  /*  in 1 nl volume */
  J_bulkSERCAbc =
      (7.5E+6 * (Y[25] * Y[25]) * (0.04 - Y[20]) - 0.7728075 * Y[20]) *
      0.0066366144807084373 * 2.0;
  /*  in 1 nl volume */
  J_SERCASRss =
      (-4.7391448307935731 * (Y[27] * Y[27]) * (0.04 - Y[21]) + Pserca * Y[21]) *
      4.99232E-5 * 2.0;
  /*  in 1 nl volume */
  J_bulkSERCAss =
      (7.5E+6 * (Y[24] * Y[24]) * (0.04 - Y[21]) - 0.7728075 * Y[21]) *
      4.99232E-5 * 2.0;
  /*  in 1 nl volume */
  /*  RyR */
  RyRcinfss = 1.0 / (exp((1000.0 * Y[24] - (Y[16] + 0.02)) / 0.01) + 1.0);
  INa_tmp = Y[27] - Y[24];
  Jrelss = 1.22 * (0.031201999999999997 * Y[14] * Y[15] *
                   (1.0 - 1.0 / (exp((Y[27] - 0.8) / 0.1) + 1.0)) * INa_tmp);
  /* 22% increase by Skibsbye 2016  */
  /* ------------------------------------------------------------------- */
  /*  RyR gates time constant */
  if (RyRcinfss >= Y[15]) {
    RyRtauinactss = 0.48;
  } else {
    RyRtauinactss = 0.06;
    /* s */
  }
  /* ------------------------------------------------------------- */
  RyRcinfbc = 1.0 / (exp((1000.0 * Y[25] - (Y[19] + 0.02)) / 0.01) + 1.0);
  INaK_tmp = Y[26] - Y[25];
  Jrelbc = 1.36 * (0.0066366144807084373 * Y[17] * Y[18] *
                   (1.0 - 1.0 / (exp((Y[26] - 0.8) / 0.1) + 1.0)) * INaK_tmp);
  /*  SR leak fluxes */
  JSRCaleak1bc = 0.006 * INaK_tmp * 0.0066366144807084373;
  JSRCaleakss = 0.006 * INa_tmp * 4.99232E-5;
  /*  Cafluxes in 1 nl volume */
  /* ********************************************** */
  /*  Differential equations */
  /* ********************************************** */
  /* 1: INa 4: It 5: IK1 6: IKr 7:IKs 8: INaK 9: INaCa 10: ICaL 11: IKur */
  /*  1 V */
  dY[0] = (((((((((((((INa + ICaL) +
                      gt * Y[7] * Y[8] * (Y[0] - -85.8248255090355)) +
                     gKur * Y[9] * Y[10] * (Y[0] - -85.8248255090355)) +
                    gK1 * pow(5.4, 0.4457) * (Y[0] - -85.8248255090355) /
                        (exp(1.5 * ((Y[0] - -85.8248255090355) + 3.6) *
                             96487.0 / 8314.0 / 310.15) +
                         1.0)) +
                   gKr * Y[12] * (1.0 / (exp((Y[0] + 55.0) / 24.0) + 1.0)) *
                       (Y[0] - -85.8248255090355)) +
                  Y[11] * (Y[0] - -85.8248255090355)) +
                 INab) +
                ICab) +
               INaK) +
              ICaP) +
             INaCa) +
            (gf * Y[13] * (0.7323 * (Y[0] - -85.8248255090355)) + IfNa)) +
           (double)Istim) /
          -0.056;
  /*  currents are in (pA) */
  /*  2 3 4 INa */
  dY[1] = (1.0 / (exp((Y[0] + 36.3) / -7.8) + 1.0) - Y[1]) /
          ((0.00013 * exp(-(a * a)) + 1.0E-5) +
           4.5E-5 / (exp((Y[0] + 42.0) / -5.0) + 1.0));
  INaK_tmp = (exp((Y[0] + 41.0) / 5.5) + 1.0) + exp(-(Y[0] + 41.0) / 14.0);
  INa_tmp = exp(-(Y[0] + 79.0) / 14.0) + 1.0;
  dY[2] = (INahinf - Y[2]) / ((0.034 / INaK_tmp + 7.0E-5) + 0.0002 / INa_tmp);
  dY[3] = (INahinf - Y[3]) / ((0.15 / INaK_tmp + 0.0007) + 0.002 / INa_tmp);
  /*  5 6 7 8 ICaL */
  dY[4] = (1.0 / (exp((Y[0] + 9.0) / -5.8) + 1.0) - Y[4]) /
          (0.0018000000000000002 * exp(-(c_a * c_a)) + 0.0005);
  /* dY(i_ICaLf1) = (ICaLfinf - Y(i_ICaLf1))/ICaLf1tau; */
  dY[5] = (1.0 / (exp((Y[0] + 27.4) / 7.1) + 1.0) - Y[5]) /
          (1.34 * exp(-(d_a * d_a)) + 0.04471428571428572);
  dY[6] = (1.0 / (b_a * b_a + 1.0) - Y[6]) / 0.002;
  /*  9 10 It */
  dY[7] = (1.0 / (exp((Y[0] - 1.0) / -11.0) + 1.0) - Y[7]) /
          (0.0035 * exp(-(e_a * e_a)) + 0.0015);
  dY[8] = (1.0 / (exp((Y[0] + 40.5) / 11.5) + 1.0) - Y[8]) /
          (0.025635 * exp(-(f_a * f_a)) + 0.01414);
  /*  11 12 IKur */
  dY[9] = (1.0 / (exp((Y[0] + 6.0) / -8.6) + 1.0) - Y[9]) /
          (0.009 / (exp((Y[0] + 5.0) / 12.0) + 1.0) + 0.0005);
  dY[10] = (1.0 / (exp((Y[0] + 7.5) / 10.0) + 1.0) - Y[10]) /
           (0.59 / (exp((Y[0] + 60.0) / 10.0) + 1.0) + 3.05);
  /*  13 IKs */
  dY[11] = (gKs / (exp((Y[0] - 19.9) / -12.7) + 1.0) - Y[11]) /
           (0.4 * exp(-(g_a * g_a)) + 0.7);
  /*  14 IKr */
  dY[12] = (1.0 / (exp((Y[0] + 15.0) / -6.0) + 1.0) - Y[12]) /
           (0.21718 * exp(-(h_a * h_a)) + 0.03118);
  /*  15 If */
  dY[13] = (1.0 / (exp((Y[0] + 97.82874) / 12.48025) + 1.0) - Y[13]) /
           (1.0 / (0.00332 * exp(-Y[0] / 16.54103) +
                   23.71839 * exp(Y[0] / 16.54103)));
  /*  16 17 18 19 20 21 RyR */
  dY[14] =
      ((1.0 - 1.0 / (exp((1000.0 * Y[24] - (Y[16] + 0.22)) / 0.03) + 1.0)) -
       Y[14]) /
      0.005;
  dY[15] = (RyRcinfss - Y[15]) / RyRtauinactss;
  dY[16] = (0.4 - 0.318 / (exp((1000.0 * Y[24] - 0.29) / 0.082) + 1.0)) - Y[16];
  dY[17] =
      ((1.0 - 1.0 / (exp((1000.0 * Y[25] - (Y[19] + 0.22)) / 0.03) + 1.0)) -
       Y[17]) /
      0.01875;
  if (RyRcinfbc >= Y[18]) {
    INa_tmp = 0.7;
  } else {
    INa_tmp = 0.0875;
  }
  dY[18] = (RyRcinfbc - Y[18]) / INa_tmp;
  dY[19] =
      (0.33 - 0.236 / (exp((1000.0 * Y[25] - 0.34) / 0.082) + 1.0)) - Y[19];
  /*  22 23 SERCACa */
  dY[20] = 0.5 * (-J_SERCASRbc + J_bulkSERCAbc) / 0.0066366144807084373;
  dY[21] = 0.5 * (-J_SERCASRss + J_bulkSERCAss) / 4.99232E-5;
  /*  24 25 26  Nai & Ki */
  dY[22] =
      1.0 / (11.318999999999999 / ((Y[22] + 10.0) * (Y[22] + 10.0)) + 1.0) *
      (-JNa / 4.99232E-5 -
       ((((INa + INab) + 3.0 * INaK) + 3.0 * INaCa) + IfNa) / 4.8169397984);
  dY[23] = JNa / 0.0066366144807084373;
  /* dY(i_Ki) =0; */
  /* -(It + IKur + IK1 + IKr + IKs - 2*INaK + IfK + Istim) / (Vcytosol*FMStr.F);
   */
  /*  27  28 Ca */
  dY[24] = 1.0 /
           (((181.50000000000003 / ((Y[24] + 1.1) * (Y[24] + 1.1)) + 1.0) +
             0.16899999999999998 / ((Y[24] + 0.013) * (Y[24] + 0.013))) +
            5.712E-5 / ((Y[24] + 0.00238) * (Y[24] + 0.00238))) *
           ((((-Jss_bc + JSRCaleakss) - J_bulkSERCAss) + Jrelss) / 4.99232E-5 +
            (((-ICaL - ICab) - ICaP) + 2.0 * INaCa) / 9.6338795968);
  dY[25] = (((Jss_bc - J_bulkSERCAbc) + JSRCaleak1bc) + Jrelbc) /
           0.0066366144807084373 *
           (1.0 / (5.712E-5 / ((Y[25] + 0.00238) * (Y[25] + 0.00238)) + 1.0));
  INaK_tmp = Y[27] - Y[26];
  dY[26] = betaSRbc * 44.0 * (INaK_tmp / 2.640625 + INaK_tmp / 15.84375) +
           ((J_SERCASRbc - JSRCaleak1bc) - Jrelbc) / 8.3994652021466166E-5 *
               betaSRbc;
  dY[27] =
      betaSRss * 44.0 * ((-Y[27] + Y[26]) / 2.640625 + INaK_tmp / 21.125) +
      ((J_SERCASRss - JSRCaleakss) - Jrelss) / 6.6452445794473717E-5 * betaSRss;
}

