function sp = manp(i,p,sp)
%% differential
[jn,jp,cen1,ces1,cep1,csn,csp,j1,j3,j11,j33,jn1,jp1,T1]=fomp(sp.y0,p);
[cen,ces,cep] = lytf(cen1,ces1,cep1,jn,jp,jn1,jp1,p,...
    de(cen1,p,p.ta),de(ces1,p,p.ta),de(cep1,p,p.ta));cew=[cen,ces,cep];
[cn,cp,cnn,cpp]=trodp(csn,csp,j1,j3,jn,jp,j11,j33,p);
%% algebraic
[~,~,~,~,phis1,phis3,phie1,phie2,phie3,~,~,~]=...
    phi( jn,jp,sp.u(i),cen,ces,cep,...
    p,ke(cen,p,p.ta),ke(ces,p,p.ta),ke(cep,p,p.ta),p.ksnp,p.kspp);
[~,~,ecd1,ecd3] = ecd(cn,cp,cen,cep,p);jn1=jn;jp1=jp;
[v,ficn,ficp,op1,op3,jn,jp,ocp1,ocp3] = volt(sp.u(i),ecd1,ecd3,phie1,phie3,phis1,phis3,cn',cp',p);
% T=temp(T1,op1,op3,dn,dp,jn,jp,isn,isp,ien,iep,phie1,phie3,s.u(i),p);
sp.y0=[jn',jp',cew,cnn,cpp,jn1,jp1,T1];
%% store
sp.v(i,:)=[ficn(1,1),ficp(1,1),v];sp.ce(i,:)=cew;sp.x(i,:)=[cn,sp.zz,cp];
sp.phis(i,:)=[phis1',sp.zz,phis3'];sp.op(i,:)=[op1',sp.zz,op3'];
sp.ocp(i,:)=[ocp1',sp.zz,ocp3'];sp.j(i+1,:)=[jn',sp.zz,jp'];
sp.ecd(i,:)=[ecd1',sp.zz,ecd3'];sp.phie(i,:)=[phie1',phie2',phie3'];
sp.SOC(i+1,:)=[100*((mean(cn/p.csn))-p.xn0)/(p.xn1-p.xn0),...
    100*((mean(cp/p.csp))-p.xp0)/(p.xp1-p.xp0)];
sp.c(i,:)=[cnn,sp.zz,cpp];sp.T(i)=(T1);
end
%%
function keff = ke(~,p,~,~,~)
% y= .0911+(1.9101*cee)/1e3-1.052*(cee/1e3).^2+.1554*(cee/1e3).^3;
 keff = p.kep;%e^b*(1e-4*ce.*((-10.5+0.668*1e-3*ce+0.494*1e-6*ce.^2) +...
%                 (0.074  -1.78*1e-5*ce -8.86*1e-10*ce.^2).*T + (-6.96*1e-5+2.8*1e-8*ce).*T.^2).^2);
end
 function Deff = de(ce,p,~)
%   Deff = .0911+(1.9101*ce)/1e3-1.052*(ce/1e3).^2+.1554*(ce/1e3).^3;
%  Deff1 = 1e-4*10.^((-4.43-54./(T-229-5e-3*ce)-0.22e-3*ce));
   Deff = p.dep.*(ce./ce);%2.5e-11.*(ce./ce);%1e-4*10.^(-(4.43+54./(T-229-1e-3*ce)+0.22*1e-3*ce));
 end
%% Temperature
% Equilibrium Potential and Gradient wrt bulk concentration
% theta_avg_n = c_avg_n / p.c_s_n_max;
% theta_avg_p = c_avg_p / p.c_s_p_max;
% [Unb,~,dUnbdT] = refPotentialAnode(p, theta_avg_n);
% [Upb,~,dUpbdT] = refPotentialCathode(p, theta_avg_p); 
% % Heat generated from intercalation (w/o boundaries for NOW)
% Q_nx = p.a_s_n*p.Faraday * jn .* (Unb - T*dUnbdT);
% Q_n = sum(Q_nx) * p.delta_x_n * p.L_n; 
% Q_px = p.a_s_p*p.Faraday * jp .* (Upb - T*dUpbdT);
% Q_p = sum(Q_px) * p.delta_x_p * p.L_p;
% Q_inter = Q_n + Q_p; 
% % Temperature ODE
% T_dot = (p.h*(p.T_amb - T) - Cur*Volt - Q_inter) / (p.rho_avg*p.C_p);
%%
% k1=real((s.u(i))./(s.a1*ecd1));k3=real((-s.u(i))./(s.a3*ecd3));
% ficn=real(sp.phie(i,1:p.n)'+sp.ocp(i,1:p.n)'-sp.phis(i,1:p.n)'+asinh(k1)/s.kb);
% ficp=real(sp.phie(i,p.n+p.p+1:p.x)'+sp.ocp(i,p.n+p.p+1:p.x)'-sp.phis(i,p.n+p.p+1:p.x)'+asinh(k3)/s.kb);

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
