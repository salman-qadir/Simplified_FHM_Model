  function [csn,csp]=trodh(cnn,cpp,jn,jp,jn1,jp1,p)
%% ANODE
bn=-((jn+jn1)*p.t*.5)/(p.f);
bn(end)=bn(end)-(p.crn*p.zn(1,2)*.5*p.ln*(jn(end)+jn1(end)))/(3*p.nsn*p.f*p.dsn);
csn= ( (p.c1n\p.c2n)*cnn' + p.c1n\bn' )';
%% CATHODE
bp=-((jp+jp1)*p.t*.5)/(p.f);
bp(end)=bp(end)-(p.crp*p.zp(1,2)*.5*p.lp*(jp(end)+jp1(end)))/(3*p.nsp*p.f*p.dsp);
csp= ( (p.c1p\p.c2p)*cpp' + p.c1p\bp' )';
end