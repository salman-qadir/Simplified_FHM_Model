function [isn,isp,ien,iep,phis1,phis3,phi1,phi2,phi3,phin,phis,phip]=phi(j11,j33,iap,cen,ces,cep,p,ken,kes,kep,ksn,ksp)
ien=cumsp(p.zn',j11');
iep=cumsp(p.zp',j33')+iap;
ie2=repmat(iap,p.s,1);
isn=iap - ien; isp=iap - iep;
phis1=-(cumsp(p.zn',isn))/ksn;
y4=cumsp(p.zp',isp);
phis3=-( (y4-y4(end))*p.zp(1,2) )/ksp; %% check for chnaging sign if required

ki=(2*(1-p.tp)*p.r*p.ta)/p.f;
y3=cumsp(p.zp',wrev((iep)./ken));
phi3=wrev(real( y3 + wrev(ki*log(cep)')- ki*log(cep(1,end)) ));

y2=cumsp(p.zs',wrev((ie2)./kes));
phi2=wrev(real( y2+wrev(ki*log(ces)')-ki*log(ces(1,end))+phi3(1,1))) ;

y1=cumsp(p.zn',wrev((ien)./kep));
phi1=wrev(real( y1+ wrev(ki*log(cen)')-ki*log(cen(1,end))+phi2(1,1) ) );
%%
k1=(p.tp*p.r*p.ta)/p.f;
y3=cumsp(p.zp',wrev((iep)./ken));
phip=wrev(real( y3 + wrev(k1*log(cep)')- k1*log(cep(1,end)) ));

y2=cumsp(p.zs',wrev((ie2)./kes));
phis=wrev(real( y2+wrev(k1*log(ces)')-k1*log(ces(1,end))+phip(1,1))) ;

y1=cumsp(p.zn',wrev((ien)./kep));
phin=wrev(real( y1+ wrev(k1*log(cen)')-k1*log(cen(1,end))+phis(1,1) ) );
end

%%
% ien=cumsim(j11)*p.zn(1,2);
% iep=cumsim(j33)*p.zp(1,2)+iap;
% ie2=repmat(iap,p.s,1);
% isn=iap - ien; isp=iap - iep;
% phis1=-(cumsim(isn)*p.zn(1,2))/p.ksn;
% y4=cumsim(isp);
% phis3=-( (y4-y4(end))*p.zp(1,2) )/p.ksp; %% check for chnaging sign if required
% 
% ki=(2*(1-p.tp)*p.r*p.t)/p.f;
% y3=cumsim(wrev((iep)./ke(cep,p.ta)'))*p.zp(1,2);
% phi3=wrev(real( y3 + wrev(ki*log(cep)')- ki*log(cep(1,end)) ));
% 
% y2=cumsim(wrev(ie2./ke(ces,p.ta)'))*p.zs(1,2);
% phi2=wrev(real( y2+wrev(ki*log(ces)')-ki*log(ces(1,end))+phi3(1,1))) ;
% 
% y1=cumsim(wrev(ien./ke(cen,p.ta)'))*p.zn(1,2);
% phi1=wrev(real( y1+ wrev(ki*log(cen)')-ki*log(cen(1,end))+phi2(1,1) ) );
%%
% ie1=.5*p.ln*cumsummat(p.n)*j11';
% ie3=.5*p.lp*cumsummat(p.p)*j33'+iap;
% ie2=repmat(iap,p.s,1);
% is1=iap - ie1;is3=iap - ie3;
% phis1=(-0.5*p.ln*cumsummat(p.n)*is1)/p.ksn;
% phis3=(-0.5*p.lp*flip(cumsummat(p.p))*is3)/p.ksp;
% %[phis1,phis3] = fis(oc,j11,j33,p,iap);
% kfi=(2*(p.tp)*p.r*p.t)/p.f;c3=-real(kfi*log(cep(1,end)));
% phie3=real((p.lp*.5*flip(cumsummat(p.p)))*(ie3./ke(cep)')+c3+kfi*log(cep)');
% 
% c2=-real(kfi*log(ces(1,end)))+phie3(1,1);
% phie2=real((p.ls*.5*flip(cumsummat(p.s)))*(ie2./ke(ces)')+c2+kfi*log(ces)');
% 
% c1=-real(kfi*log(cen(1,end)))+phie2(1,1);
% phie1=real((p.ln*.5*flip(cumsummat(p.n)))*(ie1./ke(cen)')+c1+kfi*log(cen)');
% end
% 
% function y = ke(cee)
% y= .0911+(1.9101*cee)/1e3-1.052*(cee/1e3).^2+.1554*(cee/1e3).^3;
% end
% function[phis1,phis3] = fis(oc,jna,jpa,p,u)
% a=(diag(ones(p.n-3,1),-1)-2*diag(ones(p.n-2,1),0)+diag(ones(p.n-3,1),1));
% hn=p.ln/(p.n-1);
% a(end,end)=-1;
% fn=(a*p.ksn)\((hn.^2)*jna(2:p.n-1)');
% phis1=[0;fn;fn(end)];
% 
% hp=p.lp/(p.p-1);
% gp=(p.ksp)\(jpa(2:p.p-1)'*hp.^2);gp(end)=gp(end)+(u*hp)/p.ksp;
% gp(1,1)=gp(1,1)-oc;
% a(end,end)=-1;
% a(1,1)=-2.000001;fp=a\gp;
% phis3=[fp(1);fp;fp(end)-(u*hp)/p.ksp];
% end

