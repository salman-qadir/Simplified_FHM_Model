function p = par
p.c=1.962;p.ln=53.2e-6; p.ls=24.7e-6; p.lp=39.9e-6;fg=99; 
p.zn=0:p.ln/fg:p.ln;p.zs=0:p.ls/fg:p.ls;p.zp=0:p.lp/fg:p.lp;
p.n=length(p.zn);p.s=length(p.zn);p.p=length(p.zp);p.x=p.n+p.s+p.p;
p.l=p.ln+p.ls+p.lp;p.zs1=p.ln:p.ls/fg:(p.ln+p.ls);p.zp1=(p.ln+p.ls):p.lp/fg:p.l;
p.z=0:p.l/(3*fg):p.l;p.nsn=0.626;p.nsp=0.574;p.nen=0.30;p.nes=0.35;p.nep=0.36;
p.xn0=1.2e-3;p.xn1=0.7916;p.xp1=.3494;p.xp0=.9886;p.csn=27088;p.csp=48700;
p.t=1;p.tp=0.4;p.f=96487;p.a=0.1037;p.r=8.314472;p.aa=0.5;p.rc=.027;p.ta=318;
p.kb=(p.aa*p.f)/(p.r*p.ta);p.ce=1e4;
% p.nsnn=p.nsn;% p.nsnn=.532347;%.55;% % p.nspp=p.nsp;% p.nspp=.5807;%.68;
%% FHM Parameters
p.ksn=113;p.ksp=113;p.de=1.18e-11;p.ke=0.060;
% p.dsn=8.8-12;%.5e-12;%%1e-10;%%.47e-12;%1.4e-10;.33e-10
% p.dsp=6e-11;%;.5e-12;%%2e-11;%%.43e-12;%4e-11;1e-11
 p.kn=155;
 p.kp=110;
p.dsn=1.4e-10;
p.dsp=4e-11;
%p.kn=243;
%p.kp=165;
%% P2D model parameters
p.rr=5;p.rs=5e-6;p.rx=(p.rs)/(p.rr):(p.rs)/(p.rr):p.rs;p.an=(3*p.nsn)/p.rs;
p.ap=(3*p.nsp)/p.rs;p.b=1.5;p.ksnp=113;p.kspp=113;p.dep=.99e-11;p.kep=.048;
% p.dsnp=1e-13;%2e-13;%%3.5e-14;% .99e-14;%.92e-10; solid diffusion
% p.dspp=.2e-13;%.7e-12;%;%1.42e-13;%.15e-12;%2.50e-14;1.42e-13; [m^2/s]
% p.knp=.2e-3;%2e-3;%;%2.98e-5;%%.2e-3;
% p.kpp=1e-3;%2e-3;%%1.72e-5;%.9e-3;%2.79e-11;%2e-7;
p.dsnp=2e-13;
p.dspp=.7e-12;
p.knp=2e-3;
p.kpp=2e-3;
%% homogenised model electrode matrices
crn=(p.dsn*p.t)/(p.zn(1,2)^2);p.crn=crn;
c1n=(-.5*crn*diag(ones(1,p.n-1),-1)+diag((1+crn)*ones(1,p.n),0)-.5*crn*diag(ones(1,p.n-1),1));
% c1n(1,2)=2*c1n(1,2);c1n(p.n,p.n-1)=2*c1n(p.n,p.n-1);
c1n(1,1)=1+.5*crn;   c1n(p.n,p.n)=1+.5*crn;
c2n=(.5*crn*diag(ones(1,p.n-1),-1)+diag((1-crn)*ones(1,p.n),0)+.5*crn*diag(ones(1,p.n-1),1));
% c2n(1,2)=2*c2n(1,2);c2n(p.n,p.n-1)=2*c2n(p.n,p.n-1);
c2n(1,1)=1-.5*crn;   c2n(p.n,p.n)=1-.5*crn;
%% Cathode solid
crp=(p.dsp*p.t)/(p.zp(1,2)^2);p.crp=crp;
c1p=(-.5*crp*diag(ones(1,p.p-1),-1)+diag((1+crp)*ones(1,p.p),0)-.5*crp*diag(ones(1,p.p-1),1));
% c1p(1,2)=2*c1p(1,2);c1p(p.p,p.p-1)=2*c1p(p.p,p.p-1);
c1p(1,1)=1+.5*crp;   c1p(p.p,p.p)=1+.5*crp;
c2p=(.5*crp*diag(ones(1,p.p-1),-1)+diag((1-crp)*ones(1,p.p),0)+.5*crp*diag(ones(1,p.p-1),1));
% c2p(1,2)=2*c2p(1,2);c2p(p.p,p.p-1)=2*c2p(p.p,p.p-1);
c2p(1,1)=1-.5*crp;   c2p(p.p,p.p)=1-.5*crp;
%% solid
p.c1n=c1n;p.c2n=c2n;p.c1p=c1p;p.c2p=c2p;
%% P2D model electrode Matrices
p.cnp=(p.dsnp*p.t)/(2*p.rx(1,1)^2);p.cn1p=(p.dsnp*p.t)./(p.rx(1,1).*p.rx);
c1np=(-(p.cnp+p.cn1p)'.*diag(ones(1,p.rr-1),1)+...
(1+2*p.cnp+p.cn1p)'.*diag(ones(1,p.rr),0)+...
-p.cnp*diag(ones(1,p.rr-1),-1));
c1np(1,1)=1+p.cnp+p.cn1p(1);c1np(p.rr,p.rr)=1+p.cnp;
% c1n(1,2)=-2*p.crn-p.crn1(1);c1n(p.rr,p.rr-1)=-2*p.crn-p.crn1(end);
c2np=((p.cnp+p.cn1p)'.*diag(ones(1,p.rr-1),1)+...
(1-2*p.cnp-p.cn1p)'.*diag(ones(1,p.rr),0)+...
p.cnp*diag(ones(1,p.rr-1),-1));
c2np(1,1)=1-p.cnp-p.cn1p(1);c2np(p.rr,p.rr)=1-p.cnp;
% c1n(1,2)=2*p.crn+p.crn1(1);c1n(p.rr,p.rr-1)=2*p.crn+p.crn1(end);
%% Cathode Solid Concentration
p.cpp=(p.dspp*p.t)/(2*p.rx(1,1)^2);p.cp1p=(p.dspp*p.t)./(p.rx(1,1).*p.rx);
c1pp=(-(p.cpp+p.cp1p)'.*diag(ones(1,p.rr-1),1)+...
(1+2*p.cpp+p.cp1p)'.*diag(ones(1,p.rr),0)+...
-p.cpp*diag(ones(1,p.rr-1),-1));
c1pp(1,1)=1+p.cpp+p.cp1p(1);c1pp(p.rr,p.rr)=1+p.cpp;
c2pp=((p.cpp+p.cp1p)'.*diag(ones(1,p.rr-1),1)+...
(1-2*p.cpp-p.cp1p)'.*diag(ones(1,p.rr),0)+...
p.cpp*diag(ones(1,p.rr-1),-1));
c2pp(1,1)=1-p.cpp-p.cp1p(1);c2pp(p.rr,p.rr)=1-p.cpp;
%% Solid Concentration 
p.c1np=kron(eye(p.n),c1np);p.c2np=kron(eye(p.n),c2np);
p.c1pp=kron(eye(p.p),c1pp);p.c2pp=kron(eye(p.p),c2pp);
p.b1p=repmat([zeros(p.rr-1,1);1],p.n,1);p.b2p=repmat([zeros(p.rr-1,1);1],p.p,1);
%  Temp dependent parameters i.e. diffusion and reaction rate
%%FHM ANODE DIFFUSION
x1=[270.147,271.144,272.199,273.49,274.78,275.601,277.126,277.83,278.886,279.589,280.762,281.818,283.226,284.751,286.393,287.683,288.974,290.499,292.141,293.666,294.487,295.66,296.598,298.006,298.827,299.531,300.587,301.994,303.167,303.988,304.809,306.334,307.155,308.094,308.915,309.619,310.44,311.378,312.199,313.021,314.076,314.897,315.601,316.07,316.774,317.243,317.83,318.651,319.12,319.824,320.293,320.762,321.466,321.818,322.522,323.109,323.695,324.516,324.985,325.572,326.041,326.628,327.214,327.683,328.152,328.739,329.208,329.443,329.971];
y1=1e-10.*[0.036,0.036,0.03,0.03,0.036,0.036,0.048,0.054,0.06,0.06,0.06,0.078,0.085,0.097,0.109,0.121,0.139,0.157,0.193,0.211,0.223,0.241,0.26,0.29,0.314,0.326,0.356,0.38,0.435,0.453,0.489,0.549,0.586,0.61,0.652,0.694,0.742,0.779,0.827,0.869,0.942,0.996,1.032,1.08,1.123,1.177,1.219,1.286,1.334,1.382,1.437,1.491,1.539,1.606,1.648,1.72,1.781,1.889,1.944,2.01,2.082,2.137,2.252,2.312,2.36,2.463,2.535,2.584,2.65];
% FHM ANODE  REACTION RATE
x2=[270.015,271.141,272.207,273.393,274.696,276.119,277.659,278.37,279.437,280.978,282.4,283.23,284.415,285.6,286.667,288.326,289.63,290.815,292.237,293.778,295.2,296.385,297.57,298.874,299.941,301.481,302.785,304.326,305.63,306.696,307.526,308.948,310.37,311.2,312.03,312.741,313.452,314.519,314.993,315.822,316.415,317.244,317.837,318.43,319.378,320.207,320.8,321.511,322.341,322.933,323.526,324.593,325.422,326.015,326.844,327.437,328.03,329.096,329.867];
y2=[84.337,86.747,89.759,91.566,93.976,97.59,101.205,103.012,105.422,108.434,111.446,115.06,117.47,121.084,124.096,128.313,132.53,136.145,139.759,143.976,148.795,151.807,154.217,159.036,163.855,168.675,172.892,177.108,183.133,185.542,189.157,194.578,198.193,201.205,204.819,209.036,210.843,214.458,216.867,219.277,222.289,224.699,228.313,231.325,233.735,236.747,240.361,243.373,246.386,249.398,253.012,256.024,260.241,263.253,265.663,268.675,272.289,275.904,278.313];
% FHM CATHODE DIFFUSION
x3=[271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330];
y3= 1e-11.*[0.099,0.128,0.128,0.138,0.155,0.173,0.19,0.203,0.227,0.257,0.273,0.292,0.31,0.33,0.369,0.379,0.426,0.464,0.5,0.535,0.563,0.608,0.669,0.715,0.767,0.824,0.887,0.959,1.037,1.106,1.183,1.276,1.356,1.444,1.541,1.653,1.758,1.879,2.013,2.139,2.273,2.427,2.594,2.719,2.916,3.089,3.294,3.459,3.711,3.948,4.153,4.445,4.68,5.007,5.253,5.611,5.927,6.268,6.623,6.959];
% FHM CATHODE REACTION RATE
x4=[270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329];
y4=[55.605,56.927,58.199,59.471,60.88,62.432,64.224,66.26,68.001,69.559,71.042,73.077,74.785,75.811,78.567,79.929,82.018,83.464,85.234,87.162,89.119,91.829,93.747,95.936,98.287,100.513,102.686,104.643,106.867,109.228,111.49,114.256,116.564,118.848,121.587,124.192,126.65,129.025,131.896,134.575,136.949,140.247,143.201,145.802,148.768,152.035,154.579,158.025,160.569,164.018,167.016,169.896,172.949,176.055,180.341,183.486,186.332,190.151,193.51,197.175];
% P2D ANODE DIFFUSION
x5=[271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329];
y5=1e-13.*[0.248,0.264,0.281,0.294,0.307,0.323,0.338,0.354,0.37,0.392,0.412,0.43,0.448,0.471,0.494,0.516,0.537,0.566,0.589,0.61,0.636,0.679,0.709,0.739,0.771,0.806,0.835,0.878,0.919,0.953,0.988,1.038,1.079,1.125,1.17,1.221,1.272,1.324,1.378,1.428,1.474,1.537,1.602,1.658,1.721,1.796,1.859,1.924,1.993,2.074,2.152,2.235,2.326,2.406,2.489,2.58,2.666,2.762,2.874];
% P2D ANODE  REACTION RATE
x6=[271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329];
y6=1e-3.*[0.015,0.018,0.014,0.01,0.01,0.01,0.016,0.027,0.037,0.035,0.03,0.034,0.039,0.039,0.039,0.039,0.048,0.06,0.064,0.069,0.086,0.102,0.116,0.132,0.152,0.169,0.179,0.191,0.215,0.239,0.268,0.297,0.327,0.372,0.42,0.467,0.514,0.574,0.641,0.704,0.761,0.839,0.948,1.039,1.136,1.25,1.394,1.514,1.665,1.857,2.039,2.211,2.447,2.685,2.962,3.214,3.557,3.872,4.211];
% P2D CATHODE DIFFUSION
x7=[271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329];
y7=1e-12.*[0.01,0.012,0.012,0.012,0.012,0.013,0.016,0.017,0.02,0.023,0.025,0.027,0.029,0.032,0.034,0.038,0.044,0.048,0.052,0.057,0.061,0.065,0.072,0.08,0.085,0.096,0.105,0.115,0.125,0.136,0.147,0.167,0.179,0.192,0.208,0.225,0.244,0.266,0.29,0.312,0.335,0.363,0.392,0.422,0.457,0.493,0.531,0.572,0.616,0.665,0.715,0.766,0.827,0.894,0.97,1.022,1.11,1.178,1.274];
% CATHODE ANODE  REACTION RATE
x8=[271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329];
y8=1e-3.*[0.249,0.259,0.274,0.29,0.301,0.311,0.33,0.348,0.365,0.381,0.395,0.409,0.428,0.453,0.479,0.5,0.521,0.544,0.57,0.595,0.623,0.652,0.677,0.702,0.727,0.761,0.795,0.833,0.873,0.907,0.937,0.97,1.004,1.048,1.094,1.139,1.188,1.234,1.277,1.32,1.372,1.433,1.484,1.533,1.588,1.657,1.724,1.787,1.849,1.919,1.994,2.058,2.117,2.22,2.293,2.363,2.442,2.53,2.621];
% figure;plot(x1,y1,'LineWidth',2);title(' FHM Anode Diffusion Parameter','FontSize',14);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Diffusion constant D_{s,n}[m^2s^{-1}]', 'Interpreter', 'tex','FontSize',14);grid on
% figure;plot(x2,y2,'LineWidth',2);title(' FHM Anode Reaction Rate Parameter','FontSize',14);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Reatction Rate k_{s,n}[Amol^{-1}]', 'Interpreter', 'tex','FontSize',14);grid on
% figure;plot(x3,y3,'LineWidth',2);title(' FHM Cathode Diffusion Parameter','FontSize',14);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Diffusion constant D_{s,p}[m^2s^{-1}]', 'Interpreter', 'tex','FontSize',14);grid on
% figure;plot(x4,y4,'LineWidth',2);title(' FHM Cathode Reaction Rate Parameter','FontSize',14);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Reatction Rate k_{s,p}[Amol^{-1}]', 'Interpreter', 'tex','FontSize',14);grid on
% figure;plot(x5,y5,'LineWidth',2);title(' P2D Anode Diffusion Parameter','FontSize',14);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Diffusion constant D_{s,n}[m^2s^{-1}]', 'Interpreter', 'tex','FontSize',14);grid on
% figure;plot(x6,y6,'LineWidth',2);title(' P2D Anode Reaction Rate Parameter','FontSize',14);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Reatction Rate k_{s,n}[Am^{2.5}mol^{-1.5}]', 'Interpreter', 'tex','FontSize',14);grid on
% figure;plot(x7,y7,'LineWidth',2);title(' P2D Cathode Diffusion Parameter','FontSize',14);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Diffusion constant D_{s,p}[m^2s^{-1}]', 'Interpreter', 'tex','FontSize',14);grid on
% figure;plot(x8,y8,'LineWidth',2);title(' P2D Cathode Reaction Rate Parameter','FontSize',14);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Reatction Rate k_{s,p}[Am^{2.5}mol^{-1.5}]', 'Interpreter', 'tex','FontSize',14);grid on
%% 
% figure; yyaxis left;plot(x1,y1,'--','LineWidth',2);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Diffusion constant D_{s,n}[m^2s^{-1}]', 'Interpreter', 'tex','FontSize',14);grid on; %legend('Anode');
% yyaxis right;plot(x3,y3,'-.','LineWidth',2);title(' FHM Diffusion Parameter','FontSize',14);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Diffusion constant D_{s,p}[m^2s^{-1}]', 'Interpreter', 'tex','FontSize',14);legend('Anode','Cathode','location','southeast');grid on;
% figure;yyaxis left;plot(x2,y2,'--','LineWidth',2);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Reatction Rate k_{s,n}[Amol^{-1}]', 'Interpreter', 'tex','FontSize',14);grid on; %legend('Anode');
% yyaxis right;plot(x4,y4,'-.','LineWidth',2);title(' FHM Reaction Rate Parameter','FontSize',14);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Reatction Rate k_{s,p}[Amol^{-1}]', 'Interpreter', 'tex','FontSize',14);legend('Anode','Cathode','location','southeast');grid on;
% figure;yyaxis left;plot(x5,y5,'--','LineWidth',2);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Diffusion constant D_{s,n}[m^2s^{-1}]', 'Interpreter', 'tex','FontSize',14);grid on; %legend('Anode');
% yyaxis right;plot(x7,y7,'-.','LineWidth',2);title(' P2D Diffusion Parameter','FontSize',14);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Diffusion constant D_{s,p}[m^2s^{-1}]', 'Interpreter', 'tex','FontSize',14);legend('Anode','Cathode','location','southeast');grid on;
% figure;yyaxis left;plot(x6,y6,'--','LineWidth',2);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Reatction Rate k_{s,n}[Am^{2.5}mol^{-1.5}]', 'Interpreter', 'tex','FontSize',14);grid on; %legend('Anode'); 
% yyaxis right;plot(x8,y8,'-.','LineWidth',2);title(' P2D Reaction Rate Parameter','FontSize',14);ax = gca;ax.FontSize = 14;xlabel('Temprature [K]','FontSize',14); ylabel('Reatction Rate k_{s,p}[Am^{2.5}mol^{-1.5}]', 'Interpreter', 'tex','FontSize',14);legend('Anode','Cathode','location','southeast');grid on;
%  
%%  Thermal dynamics
% p.ro = 1500;   % density, from  ref[14]
% p.cp = 700;   % Heat capacity[J/kg-K], from ref[14] 
% p.h1 = 0.36;  % Heat transfer coefficient, [W/K-m^2] 0
% %p.cp = 750; % Heat capacity [J/kg/K] p.ro  = 1626;  % Density [kg/m3]
% p.ht=65e-3;   % 18650 height [m]
% p.dm=18e-3;   % 18650 diameter [m]
% p.sv=4*(1 + p.dm/p.ht/2)/p.dm; % Surface area to volume ratio for an 18650 cell [m-1]
% p.vc=pi*(p.dm/2)^2*p.ht; % Volume of the cell
% p.h2=30; % Convection heat transfer coefficient [W/m2/K] 
%% Cutoff voltages
%p.rn = 1e-3;       % Resistivity of SEI layer, [Ohms*m^2]
%p.rp = 0;       % Resistivity of SEI layer, [Ohms*m^2]
% p.v_max = 4.2; %
% p.v_min = 2.5; 
%% Thermodynamic Params
% Thermal dynamics
% p.ro =1500;  % density, from  ref[14]
%  p.C_p = 700;   % Heat capacity[J/kg-K], from ref[14] 
%  p.h = 0.36;   % Heat transfer coefficient, [W/K-m^2] 0
% Activation Energies
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
% All units are [J/mol]
% p.E.kn = 37.48e3;
% p.E.kp = 39.57e3;
% p.E.Dsn = 42.77e3;
% p.E.Dsp = 18.55e3;
% p.E.De = 37.04e3;
% p.E.kappa_e = 34.70e3;
% Reference temperature
% % Specific interfacial surface area
% p.a_s_n = 3*p.e_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
% p.a_s_p = 3*p.e_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]
% p.L_ccn = 25e-6;    % Thickness of negative current collector [m]
% p.L_ccp = 25e-6;    % Thickness of negative current collector [m]
% % make element to caclulate phi_{s} by Saehong Park 
% p.epsilon_f_n = 0.1;  % Volume fraction of filler in neg. electrode
% p.epsilon_f_p = 0.2;  % Volume fraction of filler in pos. electrode
% % Mass densities
% p.rho_sn = 1800;    % Solid phase in negative electrode [kg/m^3]
% p.rho_sp = 5010;    % Solid phase in positive electrode [kg/m^3]
% p.rho_e =  1324;    % Electrolyte [kg/m^3]
% p.rho_f = 1800;     % Filler [kg/m^3]
% p.rho_ccn = 8954;   % Current collector in negative electrode
% p.rho_ccp = 2707;   % Current collector in positive electrode
% % Compute cell mass [kg/m^2]
% p.m_n = p.L_n * (p.rho_e*p.epsilon_e_n + p.rho_sn*p.epsilon_s_n + p.rho_f*p.epsilon_f_n);
% p.m_s = p.L_s * (p.rho_e*p.epsilon_e_n);
% p.m_p = p.L_p * (p.rho_e*p.epsilon_e_p + p.rho_sp*p.epsilon_s_p + p.rho_f*p.epsilon_f_p);
% p.m_cc = p.rho_ccn*p.L_ccn + p.rho_ccp*p.L_ccp;
% % Lumped density [kg/m^2]
% p.rho_avg = p.m_n + p.m_s + p.m_p + p.m_cc;
%%
