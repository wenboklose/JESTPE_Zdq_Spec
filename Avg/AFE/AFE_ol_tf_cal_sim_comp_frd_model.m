clc
clear all;
close all
%% voltage and current control loop as well as PLL parameters for AFE which is fixed for all the cases
Kpv_AFE      = 0.0462;          % AFE voltage loop p controller
Kiv_AFE      = 4.5815;          % AFE voltage loop i controller
Kpi_AFE      = 0.0116*270;      % AFE current loop p controller
Kii_AFE      = 46.5421*270;     % AFE current loop i controller
FWPI_a1      = 8.9207;          % AFE Phase-Locked Loop pi controller
FWPI_a2      = -8.7225;         % AFE Phase-Locked Loop pi controller
%% opertions for Bode plot display
Bode_O=bodeoptions;
Bode_O.XLabel.FontSize=14;
Bode_O.YLabel.FontSize=14;
Bode_O.TickLabel.FontSize=14;
Bode_O.Title.FontSize=14;
Bode_O.title.String=' ';
Bode_O.Grid='on';
Bode_O.XLim={[1 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='off';

% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
Tstart=0.015;
Rs=1e-6;
Vdcref=270;
alpha=1e10;         %% 0.5
Cdc=100e-6; 
R=90;
RCdc=0.070;
P=Vdcref^2/R;
Vse=57.5;
Vse=Vse+P/3/Vse*Rs*1;
Vsm=Vse*sqrt(2);
Vsm_o=Vse*sqrt(2);
L=500e-6; 
RL=100e-3;
f=400; w=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
s=tf([1 0],[0 1]);
Tdelay = 1.5/fsw;
Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);%exp(-0.5/fsw*s);
Gsvmdq = -(1-.5*Tdelay*s)*1/(1+.5*Tdelay*s);
Gsvm = [Gsvmdd 0*Gsvmdq;-0*Gsvmdq Gsvmdd];
Gsvm = tf(JF_DQFromABC(Gsvmdd,w));
Gdel = Gsvm;
%% assign control and PLL parameters for models
kpi=Kpi_AFE/270;
kii=Kii_AFE/270;
kpv=Kpv_AFE;
kiv=Kiv_AFE;
tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],Tsw);
tf_pll=d2c(tf_pll_z);

Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd=Vsdq(1);
Vsq=Vsdq(2);
Dd=0.3657;
Dq=-0.0382;
Id=8.2109;
Iq=0.0227;
fprintf('initialization is done!\n')

%% Gvd model linearization:
model = 'AFE_avg_Gvd';

%% Create the linearization I/O as specified in AFE_avg_Gvd
ios(3) = linio('AFE_avg_Gvd/vdc',1,'out');
ios(2) = linio('AFE_avg_Gvd/Dq',1,'in');
ios(1) = linio('AFE_avg_Gvd/Dd',1,'in');

%% Linearize the model
Gvd_avg_sim = linearize(model,0.5,ios);

%% Gve model linearization:
model = 'AFE_avg_Gve';

%% Create the linearization I/O as specified in AFE_avg_Gve
ios(3) = linio('AFE_avg_Gve/vdc',1,'out');
ios(2) = linio('AFE_avg_Gve/vpq',1,'in');
ios(1) = linio('AFE_avg_Gve/vpd',1,'in');

%% Linearize the model
Gve_avg_sim = linearize(model,0.5,ios);

%% Gid model linearization:
model = 'AFE_avg_Gid';

%% Create the linearization I/O as specified in AFE_avg_Gid
ios(4) = linio('AFE_avg_Gid/ilq',1,'out');
ios(3) = linio('AFE_avg_Gid/ild',1,'out');
ios(2) = linio('AFE_avg_Gid/Dq',1,'in');
ios(1) = linio('AFE_avg_Gid/Dd',1,'in');

%% Linearize the model
Gid_avg_sim = linearize(model,0.5,ios);

%% Zin model linearization:
model = 'AFE_avg_Zol';

sys = linearize('AFE_avg_Zol',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];
Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];
Zin_avg_sim=Vavg1/Ilavg1;

%% calculation and simulation comparison
figure(1)
bode(Zin_avg_sim,Bode_O)
legend('Zin\_avg\_sim')
Bode_Darklines(3)

figure(2)
bode(Gve_avg_sim,Bode_O)
legend('Gve\_avg\_sim')
Bode_Darklines(3)

figure(3)
bode(Gid_avg_sim,Bode_O)
legend('Gid\_avg\_sim')
Bode_Darklines(3)

figure(4)
bode(Gvd_avg_sim,Bode_O)
legend('Gvd\_avg\_sim')
Bode_Darklines(3)