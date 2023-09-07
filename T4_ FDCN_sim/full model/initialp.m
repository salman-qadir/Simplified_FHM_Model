function sp=initialp(p,cr)
sp.cr=cr;sp.tsp=0:p.t:5e3/abs(sp.cr);sp.zz=ones(1,p.s)*nan;
sp.u=((sp.cr*p.c)/p.a).*ones(1,length(sp.tsp));

jn=(sp.u(1)/(p.ln)).*ones(1,length(p.zn));jn1=jn;
jp=-(sp.u(1)/(p.lp)).*ones(1,length(p.zp));jp1=jp;
cnn=p.xn1*p.csn*ones(1,p.n*p.rr);
cpp=p.xp1*p.csp*ones(1,p.p*p.rr);

sp.ce(1,:)=p.ce*ones(1,p.x);sp.j(1,:)=[jn,sp.zz,jp];
cen=p.ce*ones(1,p.n);ces=p.ce*ones(1,p.s);cep=p.ce*ones(1,p.p);
ceq=[cen,ces,cep];

sp.SOC(1,:)=[100*(mean(cnn/p.csn)-p.xn0)/(p.xn1-p.xn0),...
    100*(mean(cpp/p.csp)-p.xp0)/(p.xp1-p.xp0)];

sp.y0=[jn,jp,ceq,cnn,cpp,jn1,jp1,p.ta];
end
%%
% x1=-1:.01:2;x=x1.^1;fs=x.^2/2;w=1;
% xa=w*fs(1);xb=w*fs(end);
% ff=5.392*sin(2*pi*x);
% % [qw,qe,xj]=intg(length(x));
% % xf=x'-xa*(-.5)-xb*(0.5);
% % fd=qw.*[xa;x';xb];
% fr=cumsimpsum(x)+fs(1);
% fr1=cumtrapz(x)*.01+fs(1);
% % fr1=-.5*x(1)*g;
%  f=(deff(length(x))*fs')/.01;
% % f1=f*ff/.1;
% % f2=g*ff*.1;
%  plot(x)
% %  figure
% % sum(ans(:,1))
% hold on
% plot(fr)
% %   figure
% hold on
%  plot(fr1)
% % figure
% % plot(x)
%  hold on
%  plot(fs)
%%
% csn1=cnn1-(p.rs*jn)/(p.dsn*35*p.an*p.f);
% csp1=cpp1-(p.rs*jp)/(p.dsp*35*p.ap*p.f);csn=csn1;csp=csp1;
% cn  = cnn1 +(8*p.rs*qsn1)/35 - (p.rs*jn)/(p.dsn*35*p.an*p.f);
% cp  = cpp1 +(8*p.rs*qsp1)/35 - (p.rs*jp)/(p.dsp*35*p.ap*p.f);
% ecd1=p.kn*p.f*real(((p.csn*ones(p.n,1)-cn').*cn'.*cen').^p.aa);
% ecd3=p.kp*p.f*real(((p.csp*ones(p.p,1)-cp').*cp'.*cep').^p.aa);
% [phis1,phis3,phie1,phie2,phie3]=phie(jna,jpa,sp.u(1),cen,ces,cep,p);
% k1=real(sp.u(1)./(p.ln*p.f*p.an*sp.w1'.*ecd1));
% k3=real(-sp.u(1)./(p.lp*p.f*p.ap*sp.w3'.*ecd3));
% ficn=real(phie1'+un1-phis1'+(2*p.r*p.ta*asinh(p.f*k1'))/p.f );
% ficp=real(phie3'+up1-phis3'+(2*p.r*p.ta*asinh(p.f*k3'))/p.f );
% [un,up]=ocp(cn/p.csnm,cp/p.cspm);
% u=up-un;
%  r=lyte(ce0,p,j1,j3)
%%
