figure
plot(sh.t,sh.v(:,3));hold on;
plot(sk.t,sk.v(:,3));hold on;
plot(sp.t,sp.v(:,3));hold on;
plot(ss.t,ss.v(:,3));hold on
plot(g3,g4);
xlabel('time[s]');ylabel('Voltage[V]');
title(' Voltage Output  ');
% legend('P2D Voltage','HM Voltage','exp');
 legend('HM Voltage','SHM Voltage','P2D Voltage','SP2D Voltage','Experiment');
%%
figure
plot(sh.SOC(:,1));hold on;
plot(sk.SOC(:,1));hold on;
plot(sp.SOC(:,1));hold on;
plot(ss.SOC(:,1));hold on;
xlabel('time[s]');ylabel('SOC [%]');
title(' Anode SOC Output  ');
legend('HM SOC','SHM SOC','P2D SOC','SP2D SOC');
%%
figure
plot(sh.SOC(:,2));hold on;
plot(sk.SOC(:,2));hold on;
plot(sp.SOC(:,2));hold on;
plot(ss.SOC(:,2));
xlabel('time[s]');ylabel('SOC [%]');
title('Cathode SOC Output ');
legend('HM SOC','SHM SOC','P2D SOC','SP2D SOC');
% %%
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