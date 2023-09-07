addpath(genpath('C:\Users\exx241\Dropbox\Task\T4'));
clear all; close all; clc;

tic; sq0; disp('PAP2D model'); t1=toc
   
tic; sq;  disp('PAHM model'); t=toc

tic; sqp; disp('P2D model'); tp=toc

tic; sqh; disp('Homogenised model'); th=toc

save('myfile.mat','t1','t','tp','th');%

clearvars; 

load ss.mat; load sk.mat; % ss= pap2d, sk =pahm

g= readtable('ab1.csv'); g1= table2array(g(:,1)); g2=table2array(g(:,2));
g1=lpf(g1,7e-1);g2=lpf(g2,7e-1);g1=[1,g1];g2=[g2(1),g2];
load sh.mat; load sp.mat; load myfile; p=par; % sh= HM, sp= P2D
g3=1:g1(end); g4=interp1(g1,g2,g3);
tr=min([length(sk.t),length(ss.t),length(sp.t),length(sh.t),length(g3)]);
plott;
ess=sqrt( sum( (g4(1:tr)-ss.v(1:tr,3)').^2 )/tr )
eh=sqrt(  sum( (g4(1:tr)-sk.v(1:tr,3)').^2 )/tr )
esp=sqrt( sum( (g4(1:tr)-sp.v(1:tr,3)').^2 )/tr )
esh=sqrt( sum( (g4(1:tr)-sh.v(1:tr,3)').^2 )/tr )

% g3 = polyfit(g1,g2,6);g4 = 1:g1(end);
% g5 = g3(1)*g4.^6+g3(2)*g4.^5+g3(3)*g4.^4+g3(4)*g4.^3+g3(5)*g4.^2+g3(6)*g4+g3(7);