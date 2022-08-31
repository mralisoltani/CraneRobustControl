% #### Author : Ali Soltani #### %
% #### Data/Time : 5/8/2019  22:35 #### %
% #### Different Weights From Article

clear all ;
close all ;
clc ;
warning off;

%% System Equations :

ra1 = 10 ;
ra2 = 100 ;
kt1 = 10 ;
kt2 = 10 ;
r1 = 0.1 ;
r2 = 0.05 ;
mt = 50 ;
ml = 5 ;
j1 = 0.04 ;
j2 = 0.01 ;
ke1 = 3 ;
ke2 = 1 ;
g = 9.81 ;
c1 = kt1 / ra1 ;
c1 = 1 ;
c2 = ke1 ;
c2 = 3 ;
c3 = kt2/ra2 ;
c3 = 0.1 ;
c4 = ke2 ;
c4 = 1 ;
a1 = j1/r1 + (ml+mt)*r1 ;
a2 = ml*r1 ;
a3 = c1*c2/r1 ;
a4 = ml*r2 ;
a5 = j2/r2 + ml*r2 ;
a6 = ml*g*r2 ;
a7 = c3*c4/r2 ;
a8 = g ;

syms x1 x2 x3 x4 x5 x6 u1 u2 P L


x1_dot = x4 ;
x2_dot = x5 ;
x3_dot = x6 ;

b1 = a1 - a2*cos(x3)^2 ;
b2 = a2*sin(x3) ;
b3 = a4*sin(x3) ;
b4 = a5 ;
pho1 = c1*u1 - a3*x1_dot + a2*a8*sin(x3)*cos(x3) + a2*x2*(x3_dot)^2*sin(x3) ;
pho2 = c3*u2 - a7*x2_dot + a4*x2*(x3_dot)^2 - a6*(2-cos(x3)) ;


x4_dot = (b4*pho1 - b2*pho2) / (b1*b4 - b2*b3) ;
x5_dot = (-b3*pho1 + b1*pho2) / (b1*b4 - b2*b3) ;
x6_dot = (-x4_dot*cos(x3) - 2*x5*x6 -a8*sin(x3)) / x2 ;

% Linearization about x10 = P, x20 = L, x30=x40=x50=x60 = 0, u10 = 0, u20 =
% a6/a3


A = jacobian([x1_dot;x2_dot;x3_dot;x4_dot;x5_dot;x6_dot],[x1,x2,x3,x4,x5,x6]) ;
B = jacobian([x1_dot;x2_dot;x3_dot;x4_dot;x5_dot;x6_dot],[u1,u2]) ;

P = 3; L = 1;
x1 = P; x2 = L; x3 = 0 ; x4 = 0 ; x5 = 0 ; x6 = 0 ; u1 = 0 ; u2 = a6/a3 ;

A = vpa(subs(A)) ;
B = vpa(subs(B)) ;
A_new = zeros(6,6) ;
B_new = zeros(6,2) ;
for i = 1:36
    q = sym2poly(A(i)) ;
    A_new(i) = q ;
end
for i = 1:12
    q = sym2poly(B(i)) ;
    B_new(i) = q ;
end
C_new = [1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0] ;
D_new = zeros(3,2) ;

Crane_Dynamic = ss(A_new,B_new,C_new,D_new) ;

% K = [264.153 -0.753 -444.030 278.453 -0.233 685.454;...
%     -536.062 60.015 -272.155 -476.293 13.260 -1306.792] ;
% P = [-3.5+2.5*1i -3.5-2.5*1i -3+5.5*1i -3+-5.5*1i -4.5+3.2*1i -4.5-3.2*1i] ;
% K = place(A_new,B_new,P) ;
% A_cl = A_new - B_new * K ;
% Crane_Dynamic = ss(A_new,B_new,C_new,D_new) ;
% Crane_Dynamic_cl = ss(A_cl,B_new,C_new,D_new) ;
% 
% A_1 =[0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1;...
%     0 0 (a2*a8)/(a1-a2) -a3/(a1-a2) 0 0;...
%     0 0 0 0 -a7/a5 0;0 0 -a1*a8/(L*(a1-a2)) -a3/(L*(a1-a2)) 0 0] ;
% B_1 = [0 0;0 0;0 0;c1/(a1-a2) 0;0 c3/a5;c1/(L*(a1-a2)) 0] ;
% Crane_Dynamic_1 = ss(A_1,B_1,C_new,D_new) ;
% A_cl_1 = A_1 - B_1*K ;
% Crane_Dynamic_cl_1 = ss(A_cl_1,B_1,C_new,D_new) ;
% 
% % another model :
% l = 1 ;
% AA = [0 1 0 0 0 0;0 -5.73 1.014 0 0 0;...
%     0 0 0 1 0 0;0 5.73/l -10.81/l 0 0 0;...
%     0 0 0 0 0 1;0 0 0 0 0 -2.229];
%     %x1 = cart position x2 = cart velocity x3 = swing angle x4 = swing angle velocity 
%     % x5 = rope lenght x6 = rope velocity
% BB = [0 0;0.7376 0;0 0;-7.7376/l 0;0 0;0 0.11] ;
% CC = [0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 0 1 0] ;
% DD = [0 0;0 0;0 0] ;
% Crane_Dynamics = ss(AA,BB,CC,DD) ; 
% step(Crane_Dynamics) ;

%% Nominal Plant :

[n1,d1] = ss2tf(A_new,B_new,C_new,D_new,1) ;
sys11 = tf(n1(1,:),d1) ;
sys13 = tf(n1(3,:),d1) ;

[n2,d2] = ss2tf(A_new,B_new,C_new,D_new,2) ;
sys22 = tf(n2(2,:),d2) ;

bode(sys11,sys13,sys22,{1e-1,1e3}) ;
legend('1 to cart position','1 to swing angle','2 to rope position') ;



%% Wighting functions :

s = tf('s') ;
Wa1 = (10*s^2+5*s+1)/(s^2+5*s+10) ;
Wa2 = Wa1 ;
W_P_perf = 11/(100*s+1) ; % 11/(1000*s+1) everything is ok!
W_L_perf = 11/(150*s+1) ; % 11/(1500*s+1)
W_d_1 = 2/(1000*s+1) ;
W_d_2 = 4/(1000*s+1) ;
W_sa_perf = 20/(0.1*s+1) ; %20*s/s ;
W_noise1 = 0.01*s/s ;
W_noise2 = 0.03*s/s ;
figure
 bodemag(Wa1,Wa2,W_P_perf,W_L_perf,W_d_1,W_d_2,W_sa_perf,W_noise1,W_noise2) ;
W_s_mag.LineWidth = 3 ;
grid minor
legend('Wa1','Wa2','W-P-perf','W-L-perf','W-d-1','W-d-2','W-sa-perf','W-noise1','W-noise2') ;


%% Controller Design :
clc ;
[A_P,B_P,C_P,D_P] = linmod('Nominal_Plant') ;
P = ss(A_P,B_P,C_P,D_P) ;
Iz = [1:2]' ;
Iv = [1:2]' ;
Iw = [3:8]' ;
Ie = [3:5]' ;
Iy = [6:9]' ;
Iu = [9:10]' ;
nz = 2;
nv = 2;
nw = 6;
ne = 3;
ny = 4;
nu = 2;
nctrl = length(Iu) ;
nmeas = length(Iy) ;
Pnomdesign = P([Ie;Iy],[Iw;Iu]) ;
Pnomdesign = minreal(Pnomdesign) ;
[Knom,Gnom,gamma,info] = hinfsyn(Pnomdesign,nmeas,nctrl,'METHOD','ric','TOLGAM',0.1) ;
gamma

K = Knom ;
[Aclp,Bclp,Cclp,Dclp] = linmod('Unweighted_Plant_1') ; %not considered SA as measurement
Gclp = ss(Aclp,Bclp,Cclp,Dclp) ;
Gclp_nom = Gclp([3:5],[3:8]) ;
Gclp_nom = minreal(Gclp_nom) ;

[n1,d1] = ss2tf(Gclp_nom.A,Gclp_nom.B,Gclp_nom.C,Gclp_nom.D,5) ;
[n2,d2] = ss2tf(Gclp_nom.A,Gclp_nom.B,Gclp_nom.C,Gclp_nom.D,6) ;
G11 = tf(n1(1,:),d1) ; 
G13 = tf(n1(3,:),d1) ;
G22 = tf(n2(2,:),d2) ;
bodemag(G11,G13,G22)
legend('1 to cart position','1 to swing angle','2 to rope position') ;



Grob = lft(P,Knom) ;
omega = logspace(4,-2,300) ;
Grob_f = frd(Grob,omega) ; %frequency response
RS_blk = [1 0;1 0] ;
RP_blk = [1 0;1 0;6 3] ;
muRS = mussv(Grob_f(Iz,Iv),RS_blk) ; %robust stability
muNP = svd(Grob_f(Ie,Iw)) ; % nominal performance
[muRP,muinfo0] = mussv(Grob_f,RP_blk) ; % robust performance

figure
muRS_plt = semilogx(muRS(1),'r') ;
muRS_plt.LineWidth = 3 ;
hold on
muNP_plt = semilogx(muNP(1),'c') ;
muNP_plt.LineWidth = 1 ;
muRP_plt = semilogx(muRP(1)) ;
muRP_plt.LineWidth = 1.5 ;
legend('Robust Stability','Nominal Performance','Robust Performance');
xlabel('Fequency [rad/s]');
grid on


%% D-K iteration :

%% 1-st iteration :

[VDelta0,VSigma0] = mussvextract(muinfo0) ;
D10 = VSigma0.DLeft ;
D0_perf = D10(3,3) ;
D0_1 = D10(1,1)/D0_perf ;
D0_2 = D10(2,2)/D0_perf ;

D0_1a = fitfrd(genphase(D0_1),0) ;
D0_1b = fitfrd(genphase(D0_1),1) ;
D0_1c = fitfrd(genphase(D0_1),2) ;
D0_1d = fitfrd(genphase(D0_1),3) ;
D0_1e = fitfrd(genphase(D0_1),7) ;
figure
bodemag(D0_1,D0_1a,D0_1b,D0_1c,D0_1d,D0_1e) ;
legend('Original','0','first','second','third','6-th') ; 


D0_2a = fitfrd(genphase(D0_2),0) ;
D0_2b = fitfrd(genphase(D0_2),1) ;
D0_2c = fitfrd(genphase(D0_2),2) ;
D0_2d = fitfrd(genphase(D0_2),3) ;
D0_2e = fitfrd(genphase(D0_2),10) ;
figure
bodemag(D0_2,D0_2a,D0_2b,D0_2c,D0_2d) ;
legend('Original','0','first','second','third') ;

% Pmu1design = [D0_hat,zeros(nz,ne*nmeas);zeros(ne,nz),eye(ne,ne),...
%     zeros(ne,nmeas);zeros(nmeas,nz+ne),eye(nmeas,nmeas)]...
%     * P * ...
%     [inv(D0_hat),zeros(nv,nw+nctrl);zeros(nw,nv),eye(nw,nw),...
%     zeros(nw,nctrl);zeros(nctrl,nv+nw),eye(nctrl,nctrl)] ;

% selecting fitted curves :
clear('D0_hat') ;
D0_1a = frd(D0_1a,omega) ; %0-th order is ok
D0_2a = frd(D0_2,omega) ; %0-th order is ok
clear('D0_hat') ;
D0_hat(2,2) = D0_2a ; % 
D0_hat(1,1) = D0_1a ; % 10 & 10000 ;
Pmu1design_test = blkdiag(D0_hat,eye(3)) * Grob_f * blkdiag(inv(D0_hat),eye(6)) ;
[muRP_test,muinfo0_test] = mussv(Pmu1design_test,RP_blk) ; % robust performance
figure
semilogx(muRP_test(1)) ;
hold on
semilogx(muRP(1)) ;
legend('fitted','original') ;
grid on


% after evaluating sensivities, 0-th order is selected for both
clear('D0_hat')
D0_hat(1,1) = D0_1d ;
D0_hat(2,2) = D0_1d ;
Pmu1design = blkdiag(D0_hat,eye(7)) * P * blkdiag(inv(D0_hat),eye(8)) ;
[Kmu1,Gmu1,gamma1,info1] = hinfsyn(minreal(Pmu1design),nmeas,nctrl,'METHOD','ric'); %,'TOLGAM',0.1) ; 
gamma1

Gmu1 = lft(P,Kmu1) ;
Gmu1_f = frd(Gmu1,omega) ; %frequency response
muRS1 = mussv(Gmu1_f(Iz,Iv),RS_blk) ; %robust stability
muNP1 = svd(Gmu1_f(Ie,Iw)) ; % nominal performance
[muRP1,muinfo1] = mussv(Gmu1_f,RP_blk) ; % robust performance

figure
muRS1_plt = semilogx(muRS1(1),'r') ;
muRS1_plt.LineWidth = 3 ;
hold on
muNP1_plt = semilogx(muNP1(1),'c') ;
muNP1_plt.LineWidth = 1 ;
muRP1_plt = semilogx(muRP1(1)) ;
muRP1_plt.LineWidth = 1.5 ;
legend('Robust Stability','Nominal Performance','Robust Performance');
xlabel('Fequency [rad/s]');
grid on



%% 2-nd iteration :

[VDelta1,VSigma1] = mussvextract(muinfo1) ;
D11 = VSigma1.DLeft ;
D1_perf = D11(3,3) ;
D1_1 = D11(1,1)/D1_perf ;
D1_2 = D11(2,2)/D1_perf ;

D1_1a = fitfrd(genphase(D1_1),0) ;
D1_1b = fitfrd(genphase(D1_1),1) ;
D1_1c = fitfrd(genphase(D1_1),2) ;
D1_1d = fitfrd(genphase(D1_1),3) ;
D1_1e = fitfrd(genphase(D1_2),12) ;
figure
bodemag(D1_1,D1_1a,D1_1b,D1_1c,D1_1d) ;
legend('Original','0','first','second','third') ; 
% the frequencies between 1 - 1e1 are important due to robust stability plot
% so we select 2-nd order  for D1_1

D1_2a = fitfrd(genphase(D1_2),0) ;
D1_2b = fitfrd(genphase(D1_2),1) ;
D1_2c = fitfrd(genphase(D1_2),2) ;
D1_2d = fitfrd(genphase(D1_2),3) ;
D1_2e = fitfrd(genphase(D1_2),12) ;

figure
bodemag(D1_2,D1_2a,D1_2b,D1_2c,D1_2d,D1_2e) ;
legend('Original','0','first','second','third','fourth') ;
% 3-rd order selected for D1_2

% % selecting fitted curves :
% clear('D1_hat') ;
% D1_1a = frd(D1_1,omega) ; %0-th order is ok
% D1_2a = frd(D1_2,omega) ; %0-th order is ok
% D1_hat(2,2) = D1_2a ; % 
% D1_hat(1,1) = D1_1a ; % 10 & 10000 ;
% Pmu1design_test = blkdiag(D1_hat,eye(3)) * Gmu1_f * blkdiag(inv(D1_hat),eye(6)) ;
% [muRP_test,muinfo0_test] = mussv(Pmu1design_test,RP_blk) ; % robust performance
% 
% figure
% semilogx(muRP_test(1)) ;
% hold on
% semilogx(muRP1(1)) ;
% legend('fitted','original') ;
% grid on

clear('D1_hat')

D1_hat(1,1) = (5.894e-05* s^3 + 0.008507* s^2 + 1.2* s + 0.01)/(s^3 + 0.2377* s^2 + 10.23* s + 1) ; % 1.78 not sensitive to this one
D1_hat(2,2) = (5.894e-05* s^3 + 0.008507* s^2 + 1.2* s + 1)/(s^3 + 0.2377* s^2 + 10.23* s + 1); % 0.4284
Pmu2design = blkdiag(D1_hat,eye(7)) * P * blkdiag(inv(D1_hat),eye(8)) ;

[Kmu2,Gmu2,gamma2,info2] = hinfsyn(minreal(Pmu2design),nmeas,nctrl,'METHOD','ric','TOLGAM',0.1) ; 
gamma2


Gmu2 = lft(P,Kmu2) ;

Gmu2_f = frd(Gmu2,omega) ; %frequency response
muRS2 = mussv(Gmu2_f(Iz,Iv),RS_blk) ; %robust stability
muNP2 = svd(Gmu2_f(Ie,Iw)) ; % nominal performance
[muRP2,muinfo2] = mussv(Gmu2_f,RP_blk) ; % robust performance

figure
muRS2_plt = semilogx(muRS2(1),'r') ;
muRS2_plt.LineWidth = 3 ;
hold on
muNP2_plt = semilogx(muNP2(1),'c') ;
muNP2_plt.LineWidth = 1 ;
muRP2_plt = semilogx(muRP2(1)) ;
muRP2_plt.LineWidth = 1.5 ;
legend('Robust Stability','Nominal Performance','Robust Performance');
xlabel('Fequency [rad/s]');
grid on

% %% 3-rd iteration 
% 
% [VDelta2,VSigma2] = mussvextract(muinfo2) ;
% D12 = VSigma1.DLeft ;
% D2_perf = D12(3,3) ;
% D2_1 = D12(1,1)/D2_perf ;
% D2_2 = D12(2,2)/D2_perf ;
% 
% D2_1a = fitfrd(genphase(D2_1),0) ;
% D2_1b = fitfrd(genphase(D2_1),1) ;
% D2_1c = fitfrd(genphase(D2_1),2) ;
% D2_1d = fitfrd(genphase(D2_1),3) ;
% figure
% bodemag(D2_1,D2_1a,D2_1b,D2_1c,D2_1d) ;
% legend('Original','0','first','second','third') ; 
% % the frequencies between 1 - 1e1 are important due to robust stability plot
% % so we select 2-nd order  for D1_1
% 
% D2_2a = fitfrd(genphase(D2_2),0) ;
% D2_2b = fitfrd(genphase(D2_2),1) ;
% D2_2c = fitfrd(genphase(D2_2),2) ;
% D2_2d = fitfrd(genphase(D2_2),3) ;
% D2_2e = fitfrd(genphase(D2_2),4) ;
% 
% figure
% bodemag(D2_2,D2_2a,D2_2b,D2_2c,D2_2d,D2_2e) ;
% legend('Original','0','first','second','third','fourth') ;
% % 2-nd order selected for D2_2
% 
% % % selecting fitted curves :
% % D2_1a = frd(D2_1,omega) ; %0-th order is ok
% % D2_2a = frd(D2_2,omega) ; %0-th order is ok
% % D2_hat(2,2) = D2_2a ; % 
% % D2_hat(1,1) = D2_1a ; % 10 & 10000 ;
% % Pmu1design_test = blkdiag(D2_hat,eye(3)) * Gmu2_f * blkdiag(inv(D2_hat),eye(6)) ;
% % [muRP_test,muinfo0_test] = mussv(Pmu1design_test,RP_blk) ; % robust performance
% % figure
% % semilogx(muRP_test(1)) ;
% % hold on
% % semilogx(muRP2(1)) ;
% % legend('fitted','original') ;
% % grid on
% 
% clear('D2_hat') ;
% D2_hat(1,1) = (5.894e-05* s^3 + 0.008507* s^2 + 1.2* s + 0.01)/(s^3 + 0.2377* s^2 + 10.23* s + 1) ;
% D2_hat(2,2) =  (5.894e-05* s^3 + 0.008507* s^2 + 1.2* s + 0.01)/(s^3 + 0.2377* s^2 + 10.23* s + 1) ;
% Pmu3design = blkdiag(D2_hat,eye(7)) * P * blkdiag(inv(D2_hat),eye(8)) ;
% 
% [Kmu3,Gmu3,gamma3,info3] = hinfsyn(minreal(Pmu3design),nmeas,nctrl,'METHOD','ric','TOLGAM',0.1) ; %without gamma
% gamma3
% 
% Gmu3 = lft(P,Kmu3) ;
% Gmu3_f = frd(Gmu3,omega) ; %frequency response
% muRS3 = mussv(Gmu3_f(Iz,Iv),RS_blk) ; %robust stability
% muNP3 = svd(Gmu3_f(Ie,Iw)) ; % nominal performance
% [muRP3,muinfo3] = mussv(Gmu3_f,RP_blk) ; % robust performance
% 
% figure
% muRS3_plt = semilogx(muRS3(1),'r') ;
% muRS3_plt.LineWidth = 3 ;
% hold on
% muNP3_plt = semilogx(muNP3(1),'c') ;
% muNP3_plt.LineWidth = 1 ;
% muRP3_plt = semilogx(muRP3(1)) ;
% muRP3_plt.LineWidth = 1.5 ;
% legend('Robust Stability','Nominal Performance','Robust Performance');
% xlabel('Fequency [rad/s]');
% grid on


% %% 4-th iteration 
% 
% [VDelta3,VSigma3] = mussvextract(muinfo3) ;
% D13 = VSigma1.DLeft ;
% D3_perf = D13(3,3) ;
% D3_1 = D13(1,1)/D3_perf ;
% D3_2 = D13(2,2)/D3_perf ;
% 
% D3_1a = fitfrd(genphase(D3_1),0) ;
% D3_1b = fitfrd(genphase(D3_1),1) ;
% D3_1c = fitfrd(genphase(D3_1),2) ;
% D3_1d = fitfrd(genphase(D3_1),3) ;
% figure
% bodemag(D3_1,D3_1a,D3_1b,D3_1c,D3_1d) ;
% legend('Original','0','first','second','third') ; 
% % the frequencies between 1 - 1e1 are important due to robust stability plot
% % so we select 2-nd order  for D1_1
% 
% D3_2a = fitfrd(genphase(D3_2),0) ;
% D3_2b = fitfrd(genphase(D3_2),1) ;
% D3_2c = fitfrd(genphase(D3_2),2) ;
% D3_2d = fitfrd(genphase(D3_2),3) ;
% D3_2e = fitfrd(genphase(D3_2),4) ;
% 
% figure
% bodemag(D3_2,D3_2a,D3_2b,D3_2c,D3_2d,D3_2e) ;
% legend('Original','0','first','second','third','fourth') ;
% % 2-nd order selected for D2_2
% 
% % % selecting fitted curves :
% % clear('D3_hat')
% % D3_1a = frd(D3_1,omega) ; %0-th order is ok
% % D3_2a = frd(D3_2,omega) ; %0-th order is ok
% % D3_hat(2,2) = D3_2a ; % 
% % D3_hat(1,1) = D3_1a ; % 10 & 10000 ;
% % Pmu1design_test = blkdiag(D3_hat,eye(3)) * Gmu3_f * blkdiag(inv(D3_hat),eye(6)) ;
% % [muRP_test,muinfo0_test] = mussv(Pmu1design_test,RP_blk) ; % robust performance
% % 
% % semilogx(muRP_test(1)) ;
% % hold on
% % semilogx(muRP3(1)) ;
% % legend('fitted','original') ;
% % grid on
% 
% clear('D3_hat') ;
% D3_hat(1,1) = (s+1)/(s+10) ; %D3_1d ;
% D3_hat(2,2) = (s+1)/(s+10) ; %D3_2d;
% Pmu4design = blkdiag(D3_hat,eye(7)) * P * blkdiag(inv(D3_hat),eye(8)) ;
% 
% [Kmu4,Gmu4,gamma4,info4] = hinfsyn(minreal(Pmu4design),nmeas,nctrl,'METHOD','ric','TOLGAM',0.1) ; 
% gamma4
% 
% Gmu4 = lft(P,Kmu4) ;
% Gmu4_f = frd(Gmu4,omega) ; %frequency response
% muRS4 = mussv(Gmu4_f(Iz,Iv),RS_blk) ; %robust stability
% muNP4 = svd(Gmu4_f(Ie,Iw)) ; % nominal performance
% [muRP4,muinfo4] = mussv(Gmu4_f,RP_blk) ; % robust performance
% 
% figure
% muRS4_plt = semilogx(muRS4(1),'r') ;
% muRS4_plt.LineWidth = 3 ;
% hold on
% muNP4_plt = semilogx(muNP4(1),'c') ;
% muNP4_plt.LineWidth = 1 ;
% muRP4_plt = semilogx(muRP4(1)) ;
% muRP4_plt.LineWidth = 1.5 ;
% legend('Robust Stability','Nominal Performance','Robust Performance');
% xlabel('Fequency [rad/s]');
% grid on
% 
% %% 5-th iteration
% 
% [D14,Dr4] = mussvextract(muinfo3) ;
% D4_perf = D14(3,3) ;
% D4_1 = D14(1,1)/D4_perf ;
% D4_2 = D14(2,2)/D4_perf ;
% 
% D4_1a = fitfrd(genphase(D4_1),0) ;
% D4_1b = fitfrd(genphase(D4_1),1) ;
% D4_1c = fitfrd(genphase(D4_1),2) ;
% D4_1d = fitfrd(genphase(D4_1),3) ;
% figure
% bodemag(D4_1,D4_1a,D4_1b,D4_1c,D4_1d) ;
% legend('Original','0','first','second','third') ; 
% % the frequencies between 1 - 1e1 are important due to robust stability plot
% % so we select 2-nd order  for D1_1
% 
% D4_2a = fitfrd(genphase(D4_2),0) ;
% D4_2b = fitfrd(genphase(D4_2),1) ;
% D4_2c = fitfrd(genphase(D4_2),2) ;
% D4_2d = fitfrd(genphase(D4_2),3) ;
% D4_2e = fitfrd(genphase(D4_2),4) ;
% 
% figure
% bodemag(D4_2,D4_2a,D4_2b,D4_2c,D4_2d,D4_2e) ;
% legend('Original','0','first','second','third','fourth') ;
% % 2-nd order selected for D2_2
% 
% 
% D4_hat(1,1) = 10*D4_1a ;
% D4_hat(2,2) = D4_2a;
% Pmu5design = blkdiag(D4_hat,eye(6)) * Pnomdesign * blkdiag(inv(D4_hat),eye(6)) ;
% 
% [Kmu5,Gmu5,gamma5,info5] = hinfsyn(Pmu5design,nmeas,nctrl,'METHOD','lmi','TOLGAM',0.1) ; 
% gamma5
% 
% Gmu5 = lft(P,Kmu5) ;
% Gmu5_f = frd(Gmu5,omega) ; %frequency response
% muRS5 = mussv(Gmu5_f(Iz,Iv),RS_blk) ; %robust stability
% muNP5 = svd(Gmu5_f(Ie,Iw)) ; % nominal performance
% [muRP5,muinfo5] = mussv(Gmu5_f,RP_blk) ; % robust performance
% 
% figure
% muRS5_plt = semilogx(muRS5(1),'r') ;
% muRS5_plt.LineWidth = 3 ;
% hold on
% muNP5_plt = semilogx(muNP5(1),'c') ;
% muNP5_plt.LineWidth = 1 ;
% muRP5_plt = semilogx(muRP5(1)) ;
% muRP5_plt.LineWidth = 1.5 ;
% legend('Robust Stability','Nominal Performance','Robust Performance');
% xlabel('Fequency [rad/s]');
% grid on


%% Worst Case :

mudata = frdata(muRP); % extract data
maxmu = max(mudata); % find the peak
maxidx = find(maxmu == max(maxmu)); % find index for max over omega
maxidx = maxidx(1); % ensure only one frequency chosen.
Delta0 = mussvunwrap(muinfo0); % Delta from Khinf analysis
Delta0data = frdata(Delta0);
Delta0data_w = Delta0data(:,:,maxidx);
Delta0_wc = ss(zeros(nv,nz));
for i = 1:2
    delta_i = Delta0data_w(i,i);
    gamma = abs(delta_i);
    if imag(delta_i) > 0
        delta_i = -1*delta_i;
        gamma = -1*gamma;
    end
    x = real(delta_i)/abs(gamma); % fit a Pade with the
    tau = 2*omega(maxidx)*(sqrt((1+x)/(1-x))); % same phase
    Delta0_wc(i,i) = gamma*(-s + tau/2)/(s+tau/2);
end
nDelta = norm(Delta0data_w); % the size should be 1/mu.
Delta0_wc = Delta0_wc/nDelta; % scale the perturbation to be of size 1.

Ppert = lft(Delta0_wc,P([Iz;Iy],[Iv;Iu])) ;
Ppert = frd(Ppert,omega) ;
semilogx(Ppert)