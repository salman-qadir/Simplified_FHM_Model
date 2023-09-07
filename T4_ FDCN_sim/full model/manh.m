function sh = manh(i,p,sh)
%% differential
[jn,jp,cnn,cpp,cen1,ces1,cep1,jn1,jp1,T1]=fomh(sh.y0,p);
[cen,ces,cep] = lytf(cen1,ces1,cep1,jn,jp,jn1,jp1,p,...
    de(cen1,p,p.ta),de(ces1,p,p.ta),de(cep1,p,p.ta));cew=[cen,ces,cep];
[cn,cp]=trodh(cnn,cpp,jn,jp,jn1,jp1,p);
%% algebraic
[~,~,~,~,phis1,phis3,~,~,~,phie1,phie2,phie3]=...
    phi(jn,jp,sh.u(i),cen,ces,cep,p,ke(cen,p,p.ta),ke(ces,p,p.ta),ke(cep,p,p.ta),p.ksn,p.ksp);
[ecd1,ecd3,~,~]=ecd(cn,cp,cen,cep,p);jn1=jn;jp1=jp;
[v,ficn,ficp,op1,op3,jn,jp,ocp1,ocp3] = volt(sh.u(i),ecd1,ecd3,phie1,(phie3),phis1,(phis3),cn',cp',p);
%T=temh(T1,op1,op3,dn,dp,jn,jp,isn,isp,ien,iep,phie1,phie3,s.u(i),p);
sh.y0=[jn',jp',cn,cp,cew,jn1,jp1,T1];
%% store
sh.v(i,:)=[ficn(1,1),ficp(1,1),v];
sh.ce(i+1,:)=[cen,ces,cep];sh.x(i,:)=[cn/p.csn,sh.zz,cp/p.csp];
sh.phis(i,:)=[phis1',sh.zz,phis3'];sh.op(i,:)=[op1',sh.zz,op3'];
sh.ocp(i,:)=[ocp1',sh.zz,ocp3'];sh.j(i+1,:)=[jn',sh.zz,jp'];
sh.ecd(i,:)=[ecd1',sh.zz,ecd3'];sh.phie(i,:)=[phie1',phie2',phie3'];
sh.SOC(i+1,:)=[100*(mean(cn/p.csn)-p.xn0)/(p.xn1-p.xn0),...
    100*(mean(cp/p.csp)-p.xp0)/(p.xp1-p.xp0)];sh.T(i)=(T1);
sh.c(i,:)=[cn,sh.zz,cp];sh.f(i,:)=[ficn',ficp'];
end
function keff = ke(~,p,~)
% y= .0911+(1.9101*cee)/1e3-1.052*(cee/1e3).^2+.1554*(cee/1e3).^3;
 keff = p.ke;%(1e-4*ce.*((-10.5+0.668*1e-3*ce+0.494*1e-6*ce.^2) +...
%                 (0.074  -1.78*1e-5*ce -8.86*1e-10*ce.^2).*T + (-6.96*1e-5+2.8*1e-8*ce).*T.^2).^2);
end
function Deff = de(ce,p,~)
%   Deff = .0911+(1.9101*ce)/1e3-1.052*(ce/1e3).^2+.1554*(ce/1e3).^3;
%  Deff1 = 1e-4*10.^((-4.43-54./(T-229-5e-3*ce)-0.22e-3*ce));
   Deff = p.de.*(ce./ce);%2.5e-11.*(ce./ce);%1e-4*10.^(-(4.43+54./(T-229-1e-3*ce)+0.22*1e-3*ce));
 end
% cn=c1nn1+(p.t/(p.f*p.ln))*(-jn/p.ln -s.u(i));
% c1n=c1nn1+(8*q1sn1*p.ln)/15 -jn/(p.f*p.dsn*15);
% q1sn = q1sn1+p.t*( -((jn(end)-jn(1))/(p.f*p.ln))+...
%     +(p.dsn*12*(7.5*c1nn1-7.5*c1n-5*q1sn1*p.ln))/p.ln^3);
% cp=c1pp1+(p.t/(p.f*p.lp))*(-jp/p.lp -s.u(i));
% c1p=c1pp1+(8*q1sp1*p.lp)/15 -jp/(p.f*p.dsp*15);
% q1sp = q1sp1+p.t*( -((jp(end)-jp(1))/(p.f*p.lp))+...
%     +(p.dsp*12*(7.5*c1pp1-7.5*c1p-5*q1sp1*p.lp))/p.lp^3);
% function [ye,yde] =df(s,cee)
% ye=diag(5.34e-8*exp(-0.65*cee/1e3));
% yde=(-0.65*ye)/1e3;
% % ye=(1e-4)*10^(-4.43+(54/(p.T_a-229+cee))+.22*cee);
% end
%%
% cnn=csn+s.ts*(s.c.A1*csn+s.c.B1.*jj1);
% cpp=csp+s.ts*(s.c.A3*csp+s.c.B3.*jj3);
% cn=(s.c.C1*csn)+(s.c.D1'.*j11);cp=(s.c.C3*csp)+(s.c.D3'.*j33);
% x1=(cn/p.csnm)';x3=(cp/p.cspm)';[ocp1,ocp3]=ocp(x1',x3');
% [de,dde]=df(s,c_e);
%  ce=c_e+s.ts*(de*s.c.ce1*c_e+dde*s.c.d*(s.c.ce2*c_e).^2+s.c.ce3.*j4);
