clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% case 1:
Kpv_VSI		= 0.0260;       % VSI voltage loop p controller
Kiv_VSI		= 41.89;        % VSI voltage loop i controller
Kpi_VSI		= 0.0102*270;   % VSI current loop p controller
Kii_VSI		= 5.585*270;    % VSI current loop i controller
Kpv_AFE		= 0.011;        % AFE voltage loop p controller
Kiv_AFE		= 3.49;         % AFE voltage loop i controller
Kpi_AFE		= 0.0119*270;   % AFE current loop p contrller
Kii_AFE		= 2.5598*270;   % AFE current loop i controller
FWPI_a1     = 26.7622;      % AFE Phase-Locked Loop pi controller 
FWPI_a2     = -24.9784;     % AFE Phase-Locked Loop pi controller
%% load resistor
R=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AFE parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vse=57.5;
Vsm=Vse*sqrt(2);
%% boost inductor parameters
Lboost_a    = 474e-6;           % inductor No. 2
RLboost_a   = 0.080;
Lboost_b    = 471e-6;           % inductor No. 1
RLboost_b   = 0.084;
Lboost_c    = 463e-6;           % inductor No. 3
RLboost_c   = 0.110;
%% dc link cap parameters
Cdc_afe     = 35e-6;
RCdc_afe    = 0.049;
LCdc_afe    = 0.65e-6;
%% load resistor
Rdc         = 90;
%% dead time and PWM parameters
DeadTime    = 0.5e-6;
PWM_Cycle   = 4000;
%% system parameters
f=400; w=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
m           =0.4;
%parameters for IGBT 6MBP30RH060-50
Vfs         =1.0;
Vfd         =0.75;
Rsd         =0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VSI parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% output filter inductor parameters
La          = 962e-6;       % inductor No. 1
RLa         = 0.123;        
Lb          = 985e-6;       % inductor No. 2
RLb         = 0.132;
Lc          = 972e-6;       % inductor No. 3
RLc         = 0.089;
%% output filter cap parameters
Ca          = 41.308e-6;    % cap No. 1
LCa         = 0.61e-6;
RCa         = 0.051;
Cb          = 41.308e-6;    % cap No. 2
LCb         = 0.61e-6;
RCb         = 0.054;
Cc          = 41.308e-6;    % cap No. 3
LCc         = 0.63e-6;
RCc         = 0.052;
%% dc link cap parameters
Cdc_vsi     = 10e-6;
RCdc_vsi    = 0.046;

%% input dc voltage
Vdc         = 270;
%% signal conditioning filter
num_i       = [0.214088836997569 -4.785718886316125e+04 6.011682697325629e+09 ...
    3.406392062546713e+13 8.240536716699668e+14];
den_i       = [1 1.537980654335390e+05 7.240463669538736e+09 ...
    3.437748696152356e+13 7.864394660289302e+14];
tf_i_filter = tf(num_i,den_i);
num_v       = [0.027689282410801 -6.160200000469789e+03 5.528082636005698e+08];
den_v       = [1 3.703396738686508e+04 5.528375486594036e+08];
tf_v_filter = tf(num_v,den_v);

fprintf('\ndone\n')