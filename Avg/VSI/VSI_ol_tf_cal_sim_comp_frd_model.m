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
Bode_O.XLim={[1 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='off';

% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
%% for phase to netural derivation
s=tf([1,0],[0,1]);
I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
L=1.0e-3;
RL=0.15;
Vdc=270;
Vod=99.6;
Voq=0;
R=12000;
C=31.3e-6;
Cdc=150e-6;
f=60;%400;
omega=2*pi*f;
w=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;

%% on board signal condition filter
%% signal conditioning filter
C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));

%% PWM delay
Gdeldd = exp((-1/fsw)*s)*cos(omega*(1/fsw));
Gdeldq = exp((-1/fsw)*s)*sin(omega*(1/fsw));
Tdelay = 1.5/fsw;
Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);
Gsvmdq = -(1-.5*Tdelay*s)*0/(1+.5*Tdelay*s);
Gsvm = [Gsvmdd 0*Gsvmdq;-0*Gsvmdq Gsvmdd];
Gsvm = tf(JF_DQFromABC(Gsvmdd,omega));
Gdel = Gsvm;


Vdc=270; Vdcref=270;
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsdq=[sqrt(3/2)*Vsm; 0];
Vdref=Vsdq(1); Vqref=Vsdq(2);
omega = f*2*pi;
Rac =  R;
Vod = Vdref;
Voq = Vqref;
Rl=RL;
Dd = -(omega*L*Voq-Rac*Vod-Rl*Vod+Rl*omega*C*Voq*Rac+omega^2*L*C*Vod*Rac)/(Rac*Vdc);
Dq = -(-Rl*Voq-Rac*Voq-omega*L*Vod+omega^2*L*C*Voq*Rac-Rl*omega*C*Vod*Rac)/(Rac*Vdc);
Idref = -(-Vod+omega*C*Voq*Rac)/Rac;
Iqref = (Voq+omega*C*Vod*Rac)/Rac;

%% Zo_VSI_sim model linearization:
model = 'VSI_o_loop_Z_ac';

sys = linearize('VSI_o_loop_Z_ac',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];

Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];

Isavg1=[H(5,1) H(5,2);
    H(6,1) H(6,2);];
Zin_avg_sim=Vavg1/Isavg1;
% w=2*pi*f;
%% Gii_VSI_sim model linearization:
model = 'VSI_o_loop_Gii';

sys = linearize('VSI_o_loop_Gii',0.5);
%% Gii_VSI model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];

Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];

Isavg1=[H(5,1) H(5,2);
    H(6,1) H(6,2);];
Gii_avg_sim=Vavg1/Isavg1;
%% Gid model linearization:
model = 'VSI_o_loop_Gid';

%% Create the linearization I/O as specified
ios(4) = linio('VSI_o_loop_Gid/ilq',1,'out');
ios(3) = linio('VSI_o_loop_Gid/ild',1,'out');
ios(2) = linio('VSI_o_loop_Gid/Dq',1,'in');
ios(1) = linio('VSI_o_loop_Gid/Dd',1,'in');

%% Linearize the model
Gid_avg_sim = linearize(model,0.5,ios);

%% Gvd model linearization:
model = 'VSI_o_loop_Gvd';

%% Create the linearization I/O as specified in
ios(4) = linio('VSI_o_loop_Gvd/vq',1,'out');
ios(3) = linio('VSI_o_loop_Gvd/vd',1,'out');
ios(2) = linio('VSI_o_loop_Gvd/Dq',1,'in');
ios(1) = linio('VSI_o_loop_Gvd/Dd',1,'in');

%% Linearize the model
Gvd_avg_sim = linearize(model,0.5,ios);

%% Gig_VSI_sim model linearization:
sys = linearize('VSI_o_loop_Gig',0.5);
%% Gig_VSI model linearization
H=ss(sys);

Vavg=H(1,1);

Ilavg=H(2,1);

Isavg=H(3,1);
Zsavg=Isavg/Vavg; Zlavg=Ilavg/Vavg;

Gig_avg_sim = [Zlavg 0;Zsavg 0];

%% Gig_VSI_sim model linearization:
sys = linearize('VSI_o_loop_Gvg',0.5);
%% Gii_VSI model linearization
H=ss(sys);

Vavg=H(1,1);

Ilavg=H(2,1);

Isavg=H(3,1);
Zsavg=Isavg/Vavg; Zlavg=Ilavg/Vavg;

Gvg_avg_sim = [Zlavg 0;Zsavg 0];

%% calculation and simulation comparison
figure(1)
bode(Zin_avg_sim,Bode_O)
legend('Zo\_avg\_sim')
Bode_Darklines(3)

figure(2)
bode(Gii_avg_sim,Bode_O)
legend('Gii\_avg\_sim')
Bode_Darklines(3)

figure(3)
bode(Gid_avg_sim,Bode_O)
legend('Gid\_avg\_sim')
Bode_Darklines(3)

figure(4)
bode(Gvd_avg_sim,Bode_O)
legend('Gvd\_avg\_sim')
Bode_Darklines(3)

figure(5)
bode(Gig_avg_sim,Bode_O)
legend('Gig\_avg\_sim')
Bode_Darklines(3)

figure(6)
bode(Gvg_avg_sim,Bode_O)
legend('Gvg\_avg\_sim')
Bode_Darklines(3)