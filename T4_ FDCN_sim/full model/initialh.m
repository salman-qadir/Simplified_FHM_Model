function sh=initialh(p,cr)
sh.cr=cr;sh.tsp=0:p.t:5e3/abs(sh.cr);
sh.zz=ones(1,p.s)*nan;
sh.u=((sh.cr*p.c)/p.a).*ones(1,length(sh.tsp));

jn=(sh.u(1)/p.ln).*ones(1,length(p.zn));
jp=-(sh.u(1)/p.lp).*ones(1,length(p.zp));
cn=p.xn1*p.csn*ones(1,length(p.zn));
cp=p.xp1*p.csp*ones(1,length(p.zp));

sh.ce(1,:)=p.ce*ones(1,p.x);sh.j(1,:)=[jn,sh.zz,jp];
sh.x(1,:)=[cn/p.csn,sh.zz,cp/p.csp];jn1=jn;jp1=jp;

sh.SOC(1,:)=[100*(mean(cn/p.csn)-p.xn0)/(p.xn1-p.xn0),...
    100*(mean(cp/p.csp)-p.xp0)/(p.xp1-p.xp0)];

sh.y0=[jn,jp,cn,cp,sh.ce(1,:),jn1,jp1,p.ta];%ceq=sh.ce(1,:)
end
%%
% csn1=cnn1-(p.rs*jn)/(p.dsn*35*p.an*p.f);
% csp1=cpp1-(p.rs*jp)/(p.dsp*35*p.ap*p.f);csn=csn1;csp=csp1;
% cn  = cnn1 +(8*p.rs*qsn1)/35 - (p.rs*jn)/(p.dsn*35*p.an*p.f);
% cp  = cpp1 +(8*p.rs*qsp1)/35 - (p.rs*jp)/(p.dsp*35*p.ap*p.f);
% ecd1=p.kn*p.f*real(((p.csn*ones(p.n,1)-cn').*cn'.*cen').^p.aa);
% ecd3=p.kp*p.f*real(((p.csp*ones(p.p,1)-cp').*cp'.*cep').^p.aa);
% [phis1,phis3,phie1,phie2,phie3]=phie(jna,jpa,sh.u(1),cen,ces,cep,p);
% k1=real(sh.u(1)./(p.ln*p.f*p.an*sh.w1'.*ecd1));
% k3=real(-sh.u(1)./(p.lp*p.f*p.ap*sh.w3'.*ecd3));
% ficn=real(phie1'+un1-phis1'+(2*p.r*p.ta*asinh(p.f*k1'))/p.f );
% ficp=real(phie3'+up1-phis3'+(2*p.r*p.ta*asinh(p.f*k3'))/p.f );
% [un,up]=ocp(cn/p.csnm,cp/p.cspm);
% u=up-un;
%  r=lyte(ce0,p,j1,j3)
%%
