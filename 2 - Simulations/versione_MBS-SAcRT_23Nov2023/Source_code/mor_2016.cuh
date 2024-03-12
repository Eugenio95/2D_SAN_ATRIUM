// -*- mode: c++; c-basic-offset: 8; -*-

/*
   Atrial fibroblast function (Morgan et al. 2016)

   Y[0] -> V (millivolt) (in membrane)
   Y[1] -> oa (dimensionless) (in transient_outward_K_current_oa_gate)
   Y[2] -> oi (dimensionless) (in transient_outward_K_current_oi_gate)
   Y[3] -> ua (dimensionless) (in ultra_rapid_K_current_ua_gate)
   Y[4] -> ui (dimensionless) (in ultra_rapid_K_current_ui_gate)

*/

#include <math.h>

__device__
void mor_2016(const double *Y,
              double *dY,
              const double time,
              const double *rand_g)
{
    const double K_c = 5.4;   // millimolar (in membrane)
    const double K_i = 139.0;   // millimolar (in membrane)
    const double Na_c = 140.0;   // millimolar (in membrane)
    const double Na_i = 11.2;   // millimolar (in membrane)
    const double RToF = 26.54;   // millivolt (in membrane)
    const double stim_amplitude = -100.0;   // picoA_per_picoF (in membrane)
    const double stim_duration = 2.0;   // millisecond (in membrane)
    const double stim_end = 10000.0;   // millisecond (in membrane)
    const double stim_period = 5000.0;   // millisecond (in membrane)
    const double stim_start = 100.0;   // millisecond (in membrane)
    const double B = -200.0;   // millivolt (in sodium_potassium_pump)
    const double V_rev = -150.0;   // millivolt (in sodium_potassium_pump)
    const double k_NaK_K = 1.0;   // millimolar (in sodium_potassium_pump)
    const double k_NaK_Na = 11.0;   // millimolar (in sodium_potassium_pump)

    const double g_b_Na = 0.00607;   // nanoS_per_picoF (in background_currents)
    const double g_K1 = 0.03;   // nanoS_per_picoF (in inward_rectifier)
    const double g_ns = 0.018;   // nanoS_per_picoF (in non_specific_current)
    const double i_NaK_max = 2.002;   // picoA_per_picoF (in sodium_potassium_pump)
    const double g_to = 0.01652;   // nanoS_per_picoF (in transient_outward_K_current)
    const double g_Kur = 0.6;   // nanoS_per_picoF (in ultra_rapid_K_current)
    const double E_Na = RToF * log(Na_c / Na_i);
    const double E_K = RToF * log(K_c / K_i);

    const double i_b_Na = g_b_Na*(Y[0]-E_Na);
    const double i_K1 = g_K1*(Y[0]+86.75)/(1.0+exp((Y[0]+20.0)/20.0));

    double i_Stim;

    if ((time >= stim_start) && (time <= stim_end) && (time-stim_start-floor((time-stim_start)/stim_period)*stim_period <= stim_duration))
        i_Stim = stim_amplitude;
    else
        i_Stim = 0.0;

    const double i_ns = g_ns*Y[0];
    const double i_Kur = (g_Kur*0.005+0.05/(1.0+exp(-(Y[0]-15.0)/13.0)))*pow(Y[3], 3.0)*Y[4]*(Y[0]-E_K);
    const double i_to = g_to*pow(Y[1], 3.0)*Y[2]*(Y[0]-E_K);
    const double i_NaK = i_NaK_max*K_c/(K_c+k_NaK_K)*pow(Na_i, 1.5)/(pow(Na_i, 1.5)+pow(k_NaK_Na, 1.5))*(Y[0]-V_rev)/(Y[0]-B);
    dY[0] = -(i_ns+i_Kur+i_to+i_K1+i_b_Na+i_NaK+i_Stim);
    const double alpha_oa = 0.65/(exp((Y[0]+10.0)/-8.5)+exp((Y[0]-30.0)/-59.0));
    const double beta_oa = 0.65/(2.5+exp((Y[0]+82.0)/17.0));
    const double oa_inf = 1.0/(1.0+exp((Y[0]+20.47)/-17.54));
    const double tau_oa = 15.0/(alpha_oa+beta_oa);
    dY[1] = (oa_inf-Y[1])/tau_oa;
    const double alpha_oi = 1.0/(18.53+exp((Y[0]+113.7)/10.95));
    const double beta_oi = 1.0/(35.56+exp((Y[0]+1.26)/-7.44));
    const double oi_inf = 1.0/(1.0+exp((Y[0]+43.1)/5.3));
    const double tau_oi = 15.0/(alpha_oi+beta_oi);
    dY[2] = (oi_inf-Y[2])/tau_oi;
    const double alpha_ua = 0.65/(exp((Y[0]+10.0)/-8.5)+exp((Y[0]-30.0)/-59.0));
    const double beta_ua = 0.65/(2.5+exp((Y[0]+82.0)/17.0));
    const double ua_inf = 1.0/(1.0+exp((Y[0]+33.3)/-9.6));
    const double tau_ua = 1.0/(alpha_ua+beta_ua);
    dY[3] = (ua_inf-Y[3])/tau_ua;
    const double alpha_ui = 1.0/(21.0+exp((Y[0]-185.0)/28.0));
    const double beta_ui = 1.0/exp((Y[0]-158.0)/-16.0);
    const double ui_inf = 1.0/(1.0+exp((Y[0]-99.45)/27.48));
    const double tau_ui = 5.0/(alpha_ui+beta_ui);
    dY[4] = (ui_inf-Y[4])/tau_ui;
}
