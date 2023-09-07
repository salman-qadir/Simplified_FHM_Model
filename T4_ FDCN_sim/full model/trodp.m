 function [cn,cp,csn,csp]=trodp(cnn,cpp,j1,j3,jn,jp,j11,j33,p)
%% Anode
jna=jn/(p.f*p.an);
bn=(p.b1p.*(j1+j11))*((p.cnp+p.cn1p(end))*p.rx(1))/(p.dsnp*p.f*p.an);
csn= ( (p.c1np\p.c2np)*cnn' - p.c1np\bn )';
cn=cnn(p.rr:p.rr:p.n*p.rr)-(p.rx(1)*jna)/p.dsnp;

%% CATHODE
jpa=jp/(p.f*p.ap);
bp=(p.b2p.*(p.cpp+p.cp1p(end))*p.rx(1).*(j3+j33))/(p.dspp*p.f*p.ap);
csp= ( (p.c1pp\p.c2pp)*cpp' - p.c1pp\bp )';
cp=cpp(p.rr:p.rr:p.p*p.rr)-(p.rx(1)*jpa)/p.dspp;

end
% plot((p.b1.*(j1+j11)))
%  b=9;

 

% cnn=(csn'+p.t*(p.A1*csn'+p.B1.*ja))';
% cpp=(csp'+p.t*(p.A3*csp'+p.B3.*jb))';
% cn=((p.C1*csn')+(p.D1'.*jna'))';cp=((p.C3*csp')+(p.D3'.*jpa'))';

% cnn = cnn1- (p1.t*3*jn)/(p1.rs*p1.an*p1.f);
% qsn = qsp1-(p1.rs^(-2))*p1.t*( 30*p1.dsn*qsn1+22.5*(jn/(p1.an*p1.f)));
% csn = cnn +(8*p1.rs*qsn)/35 - (p1.rs*jn)/(p1.dsn*35*p1.an*p1.f);
% 
% cpp = cpp1-(p1.t*3*jp)/(p1.rs*p1.ap*p1.f);
% qsp = qsp1-(p1.rs^(-2))*p1.t*( 30*p1.dsp*qsp1+22.5*(jp/(p1.ap*p1.f)));
% csp = cpp +(8*p1.rs*qsp)/35 - (p1.rs*jp)/(p1.dsp*35*p1.ap*p1.f);
