clc
clear all;
close all

Bode_O=bodeoptions;
Bode_O.XLabel.FontSize=14;
Bode_O.YLabel.FontSize=14;
Bode_O.TickLabel.FontSize=14;
Bode_O.Title.FontSize=14;
Bode_O.title.String=' ';
Bode_O.Grid='on';
Bode_O.XLim={[1e-4 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='off';

% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
Tstart=0.015;
Rs=1e-6;
Vdcref=270;
alpha=1e10;         %% 0.5
Cdc=150e-6; 
R=90;
RCdc=0.070;
P=Vdcref^2/R;
Vse=57.5;
Vse=Vse+P/3/Vse*Rs*1;
Vsm=Vse*sqrt(2);
Vsm_o=Vse*sqrt(2);
L=500e-6; 
RL=100e-3;
f=60;
omega=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
s=tf([1 0],[0 1]);
%% on board signal condition filter
%% signal conditioning filter
C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));


Tdelay = 1.5/fsw;
Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);%exp(-0.5/fsw*s);
Gsvmdq = -(1-.5*Tdelay*s)*0/(1+.5*Tdelay*s);
% Gsvm = [Gsvmdd*cos(omega*Tdelay) 0*Gsvmdd*sin(omega*Tdelay);-0*Gsvmdd*sin(omega*Tdelay) Gsvmdd*cos(omega*Tdelay)];
Gsvm = [Gsvmdd 1*Gsvmdq;-1*Gsvmdq Gsvmdd];

I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
    Gdel = I;
    Ki = tf_filter_dq;
    Kv = tf_filter_dq;%Hv;

% %parameters for IGBT
% Vfs= 0.5;
% Vfd=0.7;
% Rsd=0.05;
% DT=1e-6;
% DeadTime=1e-6;
% PWM_Cycle=4000;
%% PLL
DEF_pll=10;							%control loop for PLL (Hz)
DEF_pll_damp=0.707;					%damping factor for the PLL controller
DEF_Vin=57.5;						%input phase to neutral rms voltage
DEF_Tsw=0.00005;						%switching period
FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
tf_pll=d2c(tf_pll_z);

% fv=100; fi=1000;
% Lboost_con=500e-6;
% RLboost_con=100e-3;
% Cdc_con=150e-6;
% Rdc_con=90;
% kpv=2*pi*fv*Cdc_con; kiv=2*pi*fv/Rdc_con;
% kpi=2*pi*fi*Lboost_con/Vdcref; kii=2*pi*fi*RLboost_con/Vdcref;

Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd=Vsdq(1);
Vsq=Vsdq(2);
Dd=0.3777;
Dq=0.0358;
Id=-7.6907;
Iq=0;
fprintf('initialization is done!\n')


% %% Gvd model linearization:
% model = 'AFE_avg_Gvd';
% 
% %% Create the linearization I/O as specified in AFE_avg_Gvd
% ios(3) = linio('AFE_avg_Gvd/vdc',1,'out');
% ios(2) = linio('AFE_avg_Gvd/Dq',1,'in');
% ios(1) = linio('AFE_avg_Gvd/Dd',1,'in');
% 
% %% Linearize the model
% Gvd_avg_sim = linearize(model,0.5,ios);
% 
% %% Gve model linearization:
% model = 'AFE_avg_Gve';
% 
% %% Create the linearization I/O as specified in AFE_avg_Gve
% ios(3) = linio('AFE_avg_Gve/vdc',1,'out');
% ios(2) = linio('AFE_avg_Gve/vpq',1,'in');
% ios(1) = linio('AFE_avg_Gve/vpd',1,'in');
% 
% %% Linearize the model
% Gve_avg_sim = linearize(model,0.5,ios);

%% Gid model linearization:
model = 'VSI_GT_avg_Gid';

%% Create the linearization I/O as specified in AFE_avg_Gid
ios(4) = linio('VSI_GT_avg_Gid/ilq',1,'out');
ios(3) = linio('VSI_GT_avg_Gid/ild',1,'out');
ios(2) = linio('VSI_GT_avg_Gid/Dq',1,'in');
ios(1) = linio('VSI_GT_avg_Gid/Dd',1,'in');

%% Linearize the model
Gid_avg_sim = linearize(model,0.5,ios);

%% Zin model linearization:
model = 'VSI_GT_avg_Zol';

sys = linearize('VSI_GT_avg_Zol',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];
Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];
Zin_avg_sim=Vavg1/Ilavg1;

%% Zin_pll model linearization:
model = 'VSI_GT_avg_Zol_pll';

sys = linearize('VSI_GT_avg_Zol_pll',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];
Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];
Zin_pll_avg_sim=Vavg1/Ilavg1;

%% Calculation based on model
% s=tf([1 0],[0 1]);
% Zdc=R/(R*Cdc*s+1);
GL=1*L*s+RL;
omega=omega;
Vdc=270;
Zo_cal=[GL -1*omega*L; 1*omega*L GL];
Yin_cal=1/Zo_cal;
Dt=[Dd Dq;0 0];
It=[Id Iq;0 0];

den1=(GL)*(GL)-(-1*omega*L)*(1*omega*L);
num11=(Vdc)*(GL);
num12=(Vdc)*(1*omega*L);
num22=num11;
num21=-num12;
Gid_cal=[-num11/den1 -num12/den1;-num21/den1 -num22/den1];




% Vse=57.5;
% Vsm=Vse*sqrt(2);
% Vsdq=[sqrt(3/2)*Vsm; 0];
% Vsd=Vsdq(1)
% Vsq=Vsdq(2)
E0=Vsd;
Gpll=tf_pll/(s+E0*tf_pll);
Gipll=[0 Iq*Gpll;0 -Id*Gpll];
Gdpll=[0 -Dq*Gpll;0 +Dd*Gpll];

% Zin_ol_pll=1/(Gppll+Yin_ol*Gspll);
Zin_pll_cal=1/(1/Zo_cal+Gid_cal*Gdpll*Kv);
Yin_cal=1/Zo_cal;
Ypll=Gid_cal*Gdpll;
Yin_pll_cal=Yin_cal+Ypll;
%% calculation and simulation comparison
figure(1)
bode(Zin_avg_sim,Zo_cal,Bode_O)
legend('Zin\_avg\_sim','Zin\_cal')
Bode_Darklines(3)

% figure(2)
% bode(Gve_avg_sim,Gve_cal(1,:),Bode_O)
% legend('Gve\_avg\_sim','Gve\_cal')
% Bode_Darklines(3)

figure(3)
bode(Gid_avg_sim,Gid_cal,Bode_O)
legend('Gid\_avg\_sim','Gid\_cal')
Bode_Darklines(3)

% figure(4)
% bode(Gvd_avg_sim,Gvd_cal(1,:),Bode_O)
% legend('Gvd\_avg\_sim','Gvd\_cal')
% Bode_Darklines(3)

figure(5)
bode(Zin_pll_avg_sim,Zin_pll_cal,Bode_O)
legend('Zin\_pll\_avg\_sim','Zin\_pll\_cal')
Bode_Darklines(3)
Zin_pll_1_avg_sim=Zin_pll_avg_sim;
Zin_pll_1_cal=Zin_pll_cal;
% save('ZAFE_in_pll_800.mat','Zin_pll_1_avg_sim','Zin_pll_1_cal');

figure(6)
bode(Yin_cal,Ypll,Yin_pll_cal,Bode_O)
legend('Yin','Y\_PLL\_o','Yin\_ol\_PLL\_cal')
Bode_Darklines(3)

figure(7)
bode(Gid_cal,Ypll,Yin_cal,Bode_O)
Bode_Darklines(3)
