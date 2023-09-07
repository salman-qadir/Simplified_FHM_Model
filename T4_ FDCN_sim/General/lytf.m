 function [cen,ces,cep] = lytf(cen1,ces1,cep1,jn,jp,jn1,jp1,p,den,des,dep)
 bn=jn'*((1-p.tp)/(p.nen*p.f));bp=jp'*((1-p.tp)/(p.nep*p.f));
 bn1=jn1'*((1-p.tp)/(p.nen*p.f));bp1=jp1'*((1-p.tp)/(p.nep*p.f));
 z =met(den,des,dep,p);ceq=[cen1,ces1,cep1]; 
b=[bn;0*ces1';bp];b1=[bn1;0*ces1';bp1];
ce1=(z.a1\z.a2 )*ceq';
ce2=(z.a1)\(p.t*.5*b);
ce3=(z.a1)\(p.t*.5*b1);
ce=ce1+ce2+ce3;cer = ce';
cen=cer(1,1:p.n);ces=cer(1,p.n+1:p.n+p.s);
cep=cer(1,p.n+p.s+1:p.x);
%%
% cen = ( ( z.a1n\z.a2n )*cen1' + (z.a1n)\(p.t*bn) )';
% ces = ( ( z.a1s\z.a2s )*ces1' )';
% cep = ( ( z.a1p\z.a2p )*cep1' + (z.a1p)\(p.t*bp) )';
%%

% cen= (cen1'+  m.a1n*[cen1';ces1(1,1)] +m.a2n*[cen2';ces2(1,1)] + p.t*bn )';
% ces= (ces1'+ m.a1s*[cen1(1,p.n);ces1';cep1(1,1)]+m.a2s*[cen2(1,p.n);ces2';cep2(1,1)] )';
% cep= (cep1'+  m.a1p*[ces1(1,1);cep1'] + m.a2p*[ces2(1,1);cep2'] + p.t*bp  )';


% cen= (cen1'+  electrolyteDe(cen1)'.*(m.a1n*cen1') + p.t*bn )';
% ces= (ces1'+  electrolyteDe(ces1)'.*(m.a1s*ces1') )';
% cep= (cep1'+  electrolyteDe(cep1)'.*(m.a1p*cep1' ) + p.t*bp  )';


% cen=ce(1,1:p.n);ces=ce(1,p.n+1:p.n+p.s);
% 
% cep=ce(1,p.n+p.s+1:p.x);
%plot([cen,ces,cep]')

 end
 

 function z =met(den,des,dep,p) 
 %% Anode Electrolyte Concentration

rn=(p.t*den')./(p.nen*p.zn(1,2)^2);kn=(p.zs(1,2))/(p.zn(1,2));kn1=1/(1+kn); kn2=kn/(1+kn);
a1n=(-.5*rn.*diag(ones(1,p.n-1),-1)+(1+rn).*diag(ones(1,p.n),0)+-.5*rn.*diag(ones(1,p.n-1),1));
a2n=(.5*rn.*diag(ones(1,p.n-1),-1)+(1-rn).*diag(ones(1,p.n),0)+.5*rn.*diag(ones(1,p.n-1),1));

%  a1n(1,1)=.5*a1n(1,1);a2n(1,1)=.5*a2n(1,1);  % Collector
   a1n(1,2)=2*a1n(1,2);a2n(1,2)=2*a2n(1,2);  %anod
%   a1n(p.n,p.n)=.5*a1n(p.n,p.n);a2n(p.n,p.n)=.5*a2n(p.n,p.n);  % sep
%%
%% Separator  Electrolyte Concentration
rs=(p.t*des')./(p.nes*p.zs(1,2)^2);kp=(p.zp(1,2))/(p.zs(1,2));kp1=1/(1+kp); kp2=kp/(1+kp);
a1s=(-.5*rs.*diag(ones(1,p.s-1),-1)+(1+rs).*diag(ones(1,p.s),0)+-.5*rs.*diag(ones(1,p.s-1),1));
a2s=(.5*rs.*diag(ones(1,p.s-1),-1)+(1-rs).*diag(ones(1,p.s),0)+.5*rs.*diag(ones(1,p.s-1),1));
%   a1s(1,1)=.5*a1s(1,1);a2s(1,1)=.5*a2s(1,1);  %anod
%     a1s(1,2)=2*a1s(1,2);a2s(1,2)=2*a2s(1,2);  %anod
%    a1s(p.s,p.s)=.5*a1s(p.s,p.s);a2s(p.s,p.s)=.5*a2s(p.s,p.s);  % Cathod
%% Cathode Electrolyte Concentration
rp=(p.t*dep')./(p.nep*p.zp(1,2)^2);
a1p=(-.5*rp.*diag(ones(1,p.p-1),-1)+(1+rp).*diag(ones(1,p.p),0)+-.5*rp.*diag(ones(1,p.p-1),1));
a2p=(.5*rp.*diag(ones(1,p.p-1),-1)+(1-rp).*diag(ones(1,p.p),0)+.5*rp.*diag(ones(1,p.p-1),1));

% a1p(1,1)=.5*a1p(1,1);a2p(1,1)=.5*a2p(1,1);  % 
%   a1p(1,2)=2*a1p(1,2);a2p(1,2)=2*a2p(1,2);  %anod
%  a1p(p.p,p.p)=.5*a1p(p.p,p.p);a2p(p.p,p.p)=.5*a2p(p.p,p.p);% Col
  a1p(p.p,p.p-1)=2*a1p(p.p,p.p-1);a2p(p.p,p.p-1)=2*a2p(p.p,p.p-1);% Col
 %%
 a1=blkdiag(a1n,a1s,a1p);a2=blkdiag(a2n,a2s,a2p);
 %% anode separator side
 a1(p.n,p.n)=a1(p.n,p.n)-.5*rn(p.n)*kn2; % cn N-1
 a2(p.n,p.n)=a2(p.n,p.n)+.5*rn(p.n)*kn2; 
 a1(p.n,p.n+1)=-.5*rn(p.n)*kn1; % cs 0
 a2(p.n,p.n+1)=.5*rn(p.n)*kn1;
 %% separator anode side
 a1(p.n+1,p.n+1)=a1(p.n+1,p.n+1)-.5*rs(1)*kn1; % cs 0
 a2(p.n+1,p.n+1)=a2(p.n+1,p.n+1)+.5*rs(1)*kn1;
 a1(p.n+1,p.n)=-.5*rs(1)*kn2; % cn  N-1
 a2(p.n+1,p.n)=.5*rs(1)*kn2;
 %% separator cathode side
 a1(p.n+p.s,p.n+p.s)=a1(p.n+p.s,p.n+p.s)-.5*rs(p.s)*kp2; % cs  N-1
 a2(p.n+p.s,p.n+p.s)=a2(p.n+p.s,p.n+p.s)+.5*rs(p.s)*kp2;
 a1(p.n+p.s,p.n+p.s+1)=-.5*rs(p.s)*kp1; % cp 0
 a2(p.n+p.s,p.n+p.s+1)=.5*rs(p.s)*kp1;
 %% cathode separator side
 a1(p.n+p.s+1,p.n+p.s+1)=a1(p.n+p.s+1,p.n+p.s+1)-.5*rp(1)*kp1; % cp 0
 a2(p.n+p.s+1,p.n+p.s+1)=a2(p.n+p.s+1,p.n+p.s+1)+.5*rp(1)*kp1; 
 a1(p.n+p.s+1,p.n+p.s)=-.5*rp(1)*kp2; % cs  N-1
 a2(p.n+p.s+1,p.n+p.s)=.5*rp(1)*kp2; 
 %%
  z.a1n=a1n;z.a2n=a2n;z.a1s=a1s;z.a2s=a2s;z.a1p=a1p;z.a2p=a2p;
 z.a1=a1;z.a2=a2;
 jk=9;
 end
 
 
