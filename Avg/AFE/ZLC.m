%% dq average model and analytical model of VSI Copyright ?2013 Boeing. All rights reserved.
clc
clear all;
close all

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
Bode_O.PhaseWrapping='on';

s = tf([1 0],[0 1]);
omega = 2*pi*400;%2*pi*60;%
%% impedance of the cable (1km):
% low voltage: R: 0.642Ohm/km X: 0.083Ohm/km
RL = 0.642;
L = 0.083/(omega);
% medium voltage: R: 0.161Ohm/km X: 0.190Ohm/km
RL = 0.161;
L = 0.19/omega;
% % high voltage: R: 0.06Ohm/km X: 0.191Ohm/km
% RL = 0.06;
% L = 0.191/omega
% C = 10e-6;
% RC = 2*0.000004/(omega*C);
s = tf([1 0],[0 1]);
ZLdq = [L*s + RL -L*omega;L*omega L*s + RL];

% YCdq = [C*s -C*60*2*pi;C*60*2*pi C*s];
% ZCdq = [RC 0;0 RC]+1/YCdq;
% ZLCdq = 1/(1/ZLdq+1/ZCdq);

n=1e2;
f=logspace(-1,4,n);
w=f*2*pi;
% 
% Zdd = freqresp(ZLCdq(1,1),w);
% Zdq = freqresp(ZLCdq(1,2),w);
% Zqd = freqresp(ZLCdq(2,1),w);
% Zqq = freqresp(ZLCdq(2,2),w);
% 
% Zdd_frd = frd(Zdd,w);
% Zdq_frd = frd(Zdq,w);
% Zqd_frd = frd(Zqd,w);
% Zqq_frd = frd(Zqq,w);
% Z_LC_dq_frd = [Zdd_frd Zdq_frd;Zqd_frd Zqq_frd];

Zdd = freqresp(ZLdq(1,1),w);
Zdq = freqresp(ZLdq(1,2),w);
Zqd = freqresp(ZLdq(2,1),w);
Zqq = freqresp(ZLdq(2,2),w);

Zdd_frd = frd(Zdd,w);
Zdq_frd = frd(Zdq,w);
Zqd_frd = frd(Zqd,w);
Zqq_frd = frd(Zqq,w);
Z_L_dq_frd = [Zdd_frd Zdq_frd;Zqd_frd Zqq_frd];

% figure(1)
% bode(Z_LC_dq_frd(1,1),Z_LC_dq_frd(2,1),Bode_O)
AZ = 1/sqrt(2)*[1 1i;1 -1i];
% ZLC_pn = AZ*Z_LC_dq_frd*inv(AZ);
% figure(2)
% bode(ZLC_pn(1,1),ZLC_pn(2,2),Bode_O)
    
fig=figure(3)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(ZLdq(1,1),ZLdq(2,1),Bode_O)
Bode_Darklines(3)

ZRL_pn = AZ*Z_L_dq_frd*inv(AZ);
fig=figure(4)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(ZRL_pn(1,1),ZRL_pn(2,2),Bode_O)
Bode_Darklines(3)