%%
g1=[0,16.31,44.852,73.394,110.092,142.712,171.254,207.951,236.493,260.958,293.578,330.275,362.895,403.67,444.444,472.987,505.607,534.149,558.614,591.233,623.853,656.473,693.17,721.713,750.255,791.03,823.649,856.269,892.966,929.664,974.516,1019.368,1072.375,1125.382,1174.312,1223.242,1264.016,1296.636,1341.488,1378.186,1427.115,1480.122,1545.362,1573.904,1614.679,1643.221,1692.151,1716.616,1753.313,1785.933,1818.552,1855.25,1896.024,1940.877,1993.884,2030.581,2071.356,2108.053,2140.673,2173.293,2214.067,2258.919,2307.849,2344.546,2389.399,2442.406,2503.568,2568.807,2617.737,2678.899,2727.829,2772.681,2817.533,2866.463,2943.935,2988.787,3033.639,3082.569,3139.653,3184.506,3229.358,3282.365,3343.527,3392.457,3441.386,3494.393,3539.246,3592.253,3628.95,3661.57,3690.112,3714.577,3730.887,3759.429,3779.817,3796.126,3820.591,3832.824,3849.134,3857.288,3869.521,3877.676,3883.792];
g2=[4.113,4.099,4.086,4.069,4.058,4.049,4.035,4.023,4.017,4.005,3.996,3.98,3.966,3.952,3.941,3.931,3.922,3.915,3.906,3.897,3.888,3.881,3.872,3.86,3.851,3.842,3.828,3.821,3.809,3.8,3.784,3.775,3.761,3.745,3.733,3.72,3.71,3.701,3.692,3.685,3.671,3.662,3.646,3.639,3.632,3.623,3.616,3.612,3.605,3.598,3.593,3.588,3.586,3.582,3.575,3.572,3.565,3.561,3.559,3.559,3.554,3.549,3.545,3.54,3.536,3.533,3.524,3.517,3.508,3.503,3.496,3.494,3.492,3.483,3.476,3.464,3.46,3.453,3.444,3.427,3.418,3.407,3.391,3.381,3.374,3.363,3.349,3.338,3.324,3.303,3.282,3.25,3.223,3.197,3.167,3.124,3.055,2.965,2.9,2.824,2.737,2.67,2.583];

% g= readtable('ab1.csv'); g1= table2array(g(:,1)); g2=table2array(g(:,2));
plot(g1,g2)
g1=lpf(g1,7e-1);g2=lpf(g2,7e-1);g1=[1,g1];g2=[g2(1),g2];

g3=1:g1(end); g4=interp1(g1,g2,g3);

tr=min([length(sk.t),length(ss.t),length(sp.t),length(sh.t),length(g3)]);
%tr=min([length(sk.t),length(sh.t),length(g3)]);
 e=(100/mean(g4)).*[sqrt( sum( (g4(1:tr)-ss.v(1:tr,3)').^2 )/tr ),sqrt( sum( (g4(1:tr)-sk.v(1:tr,3)').^2 )/tr ),...
     sqrt( sum( (g4(1:tr)-sp.v(1:tr,3)').^2 )/tr ),sqrt( sum( (g4(1:tr)-sh.v(1:tr,3)').^2 )/tr ),...
 sqrt( sum( (sh.v(1:tr,3)-sk.v(1:tr,3)).^2 )/tr ),sqrt( sum( (sh.v(1:tr,3)-sqq.v(1:tr,3)).^2 )/tr )];
% (100/mean(sh.v(1:tr,3)))*sqrt( sum( (sh.v(1:tr,3)-sk.v(1:tr,3)).^2 )/tr )
 fprintf('\rEsp2d\t=\t%5.5f\r',e(1));fprintf('\rEshm\t=\t%5.5f\r',e(2));
 fprintf('\rEp2d\t=\t%5.5f\r',e(3));   fprintf('\rEhm\t=\t%5.5f\r',e(4));
 fprintf('\rEsfhm\t=\t%5.5f\r',e(5));
% 
 fprintf('\rTsp2d\t=\t%5.5f\r',ss.Tsp2d);fprintf('\rTshm\t=\t%5.5f\r',sk.Tshm);
 fprintf('\rTp2d\t=\t%5.5f\r',sp.Tp);fprintf('\rThm\t=\t%5.5f\r',sh.Th);
 fprintf('\rTsfhm\t=\t%5.5f\r',sqq.T4p);
%%=
figure
plot(sh.t,sh.v(:,3),'LineWidth',2);hold on;
plot(sk.t,sk.v(:,3),'--','LineWidth',2);hold on;
plot(sp.t,sp.v(:,3),'LineWidth',2);hold on;
plot(ss.t,ss.v(:,3),'-.','LineWidth',2);hold on;
plot(sqq.t,sqq.v(:,3),'-.','LineWidth',2);hold on;
plot(g3,g4,'LineWidth',2);
%ax=gca; ax.FontSize = 14;
xlabel('Time[s]');ylabel('Voltage[V]','FontSize',14);
%title('Voltage Output ','FontSize',14);
legend('FHM Model' ,'Simplified FHM Model','P2D Model','Simplified P2D Model','Further Simplified FHM Model','Experiment','location','southwest','FontSize',12);
ylim([2.5 4.15])
xlim([0 4000])
grid on
%
figure
plot(sh.t,sh.v(:,3),'LineWidth',2);hold on;
plot(sk.t,sk.v(:,3),'--','LineWidth',2);hold on;
plot(sp.t,sp.v(:,3),'LineWidth',2);hold on;
plot(ss.t,ss.v(:,3),'-.','LineWidth',2);hold on;
plot(sqq.t,sqq.v(:,3),'-.','LineWidth',2);hold on;
%plot(g3,g4,'LineWidth',2);
%ax=gca; ax.FontSize = 14;
xlabel('Time[s]');ylabel('Voltage[V]','FontSize',14);
%title('Voltage Output ','FontSize',14);
legend('FHM Model' ,'Simplified FHM Model','P2D Model','Simplified P2D Model','Further Simplified FHM Model','location','southwest','FontSize',12);
ylim([2.5 4.15])
xlim([0 1000])
grid on
%%
figure
plot(sh.SOC(:,1),'LineWidth',2);hold on;
plot(sk.SOC(:,1),'--','LineWidth',2);hold on;
plot(sp.SOC(:,1),'LineWidth',2);hold on;
plot(ss.SOC(:,1),'-.','LineWidth',2);hold on;
plot(sqq.SOC(:,1),'-.','LineWidth',2);hold on;
xlabel('time[s]');ylabel('SOC [%]');
%title(' Anode SOC Output  ');
legend('FHM SOC','Simplified FHM SOC','P2D SOC','Simplified P2D SOC','Further Simplified FHM Model');
grid on
%%
% figure
% plot(sh.SOC(:,2));hold on;
% plot(sk.SOC(:,2));hold on;
% plot(sp.SOC(:,2));hold on;
% plot(ss.SOC(:,2));
% xlabel('time[s]');ylabel('SOC [%]');
% title('Cathode SOC Output ');
% legend('HM SOC','SHM SOC','P2D SOC','SP2D SOC');
%%
% figure
% plot(sh.t,sh.c(1:end,2));hold on;
% plot(sk.t,sk.c(1:end,2));hold on;
% plot(sp.t,sp.c(1:end,2));hold on;
% plot(ss.t,ss.c(2:end,2));hold on
% xlabel('time[s]');ylabel('charge[mol]');
% title(' charge Anode  ');
% legend('HM Voltage','SHM Voltage','P2D Voltage','SP2D Voltage','Experiment');
%%
% figure;pu=2;
% plot(sh.t,sh.v(:,pu));hold on;
% plot(sk.t,sk.v(:,pu));hold on;
% plot(sp.t,sp.v(:,pu));hold on;
% plot(ss.t,ss.v(:,pu));hold on;
% xlabel('time[s]');ylabel('Voltage[V]');
% title(' Voltage Anode  ');
% legend('HM Voltage','SHM Voltage','P2D Voltage','SP2D Voltage','Experiment');
%%
% figure;pu=1;
% plot(sh.t,sh.ocp(:,pu));hold on;
% plot(sk.t,sk.ocp(:,pu));hold on;
% plot(sp.t,sp.ocp(:,pu));hold on;
% plot(ss.t,ss.ocp(2:end,pu));hold on;
% xlabel('time[s]');ylabel('Voltage[V]');
% title(' Voltage Anode  ');
% legend('HM Voltage','SHM Voltage','P2D Voltage','SP2D Voltage','Experiment');
%%
p=par;
figure;
plot(p.zn', sh.c(4,1:100)','--','LineWidth',2);hold on;
plot(p.zn', sh.c(8,1:100)','--','LineWidth',2);hold on;
plot(p.zn', sh.c(12,1:100)','--','LineWidth',2);hold on;
plot(p.zn', sh.c(16,1:100)','--','LineWidth',2);hold on;
plot(p.zn', sk.c(4,1:100)','-.','LineWidth',2);hold on;
plot(p.zn', sk.c(8,1:100)','-.','LineWidth',2);hold on;
plot(p.zn', sk.c(12,1:100)','-.','LineWidth',2);hold on;
plot(p.zn', sk.c(16,1:100)','-.','LineWidth',2);hold on;
xlabel('Anode length[m]');ylabel('Li-ion Concentration in electrode [mol/m^3]');
legend('4s FHM model','8s FHM model','12s FHM model','16s FHM model','4s SFHM model','8s SFHM model','12s SFHM model','16s SFHM model','location','northeast');
grid on
p=par;
figure;
plot(p.zn', sk.c(4,1:100)','--','LineWidth',2);hold on;
plot(p.zn', sk.c(8,1:100)','-.','LineWidth',2);hold on;
plot(p.zn', sk.c(12,1:100)','-*','LineWidth',2);hold on;
plot(p.zn', sk.c(16,1:100)','-+','LineWidth',2);hold on;
xlabel('Anode length[m]');ylabel('Li-ion Concentration in electrode [mol/m^3]');
legend('4s','8s','12s','16s','location','northeast');
grid on
%%
% g3 = polyfit(g1,g2,6);g4 = 1:g1(end);
% g5 = g3(1)*g4.^6+g3(2)*g4.^5+g3(3)*g4.^4+g3(4)*g4.^3+g3(5)*g4.^2+g3(6)*g4+g3(7);
% figure
% plot(sh.ocp(:,end));hold on;plot(sk.ocp(:,end));hold on;
% plot(sp.ocp(:,end));hold on;plot(ss.ocp(:,end));
% xlabel('time[s]');ylabel('SOC [%]');
% title('Cathode OCP Output ');
% legend('HM SOC' ,'PAHM SOC','P2D SOC','PAP2D SOC');
% %%
% figure
% plot(sh.ocp(:,1));hold on;plot(sk.ocp(:,1));hold on;
% plot(sp.ocp(:,1));hold on;plot(ss.ocp(:,1));
% xlabel('time[s]');ylabel('Voltage [%]');
% title('Anode OCP Output ');
% legend('HM OCP' ,'PAHM OCP','P2D OCP','PAP2D OCP');

%%
% figure
% mesh(sk.phis(:,1:p.n)); hold on; mesh(ss.phis(:,1:p.n)); hold on;
% mesh(sh.phis(:,1:p.n)); hold on; mesh(sp.phis(:,1:p.n)); hold on;
% legend('PAHM phisn','PAP2D phisn','HM phisp','P2D phisp');