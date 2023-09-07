 function [cen,ces,cep,qf,g1] = lyte(u,p,qb,g,de)
qf= qb+p.t*([ (u-u*p.tp)/p.f + de*2*g(1)*p.ln,...
    de*2*g(3)*(p.ls),-(u-u*p.tp)/p.f + de*2*g(6)*(p.lp)]);  
a=[2*de*p.ln,0,-2*de*p.ln,-de,0,0,0;...
     0,0,de*2*(p.ln+p.ls),de,0,2*de*p.lp,0;...
   p.ln^2,1,-p.ln^2,-p.ln,-1,0,0;...
   0,0,(p.ln+p.ls)^2,(p.ln+p.ls),1,-p.lp^2,-1;...
   (1/3)*p.nen*p.ln^3,p.nen*p.ln,0,0,0,0,0;...
0,0,(1/3)*p.nes*((p.ln+p.ls)^3-p.ln^3),.5*p.nes*((p.ln+p.ls)^2-p.ln^2)...
,p.nes*p.ls,0,0;0,0,0,0,0,(1/3)*p.nep*p.lp^3,p.nep*p.lp];
b=[0;0;0;0;qb'];
%   g11=(a\b)'; 
a1=a*2e9; b1=b*2e9;
g1=(a1\b1)';
%  lk=g1-g11; kj=[g11;g1];
cen=g(1)*((p.zn).^2)+g(2)*ones(1,length(p.zn));
ces=g(3)*((p.zs1).^2)+g(4)*p.zs1+g(5)*ones(1,length(p.zs1));
cep=g(6)*((p.l-p.zp1).^2)+g(7)*ones(1,length(p.zp1));
 end
 %   qf= qb+p.t*([ (u-u*p.tp)/p.f+de*2*g(2)*p.ln,...
%  de*2*g(5)*p.ls,-(u-u*p.tp)/p.f+de*2*g(7)*p.lp]);
 
% a=[1,p.ln^2,-1,0,0,0,0;0,0,-1,-p.ls,-p.ls^2,1,p.lp^2;...
%    0,2*p.ln*de,0,-de,0,0,0;0,0,0,de,2*de*p.ls,0,2*p.lp*de;...
%    p.nen*p.ln,(1/3)*p.nen*p.ln^3,0,0,0,0,0;...
%    0,0,p.nes*p.ls,(1/2)*p.nes*(p.ls^2),(1/3)*p.nes*(p.ls^3),0,0;...
%    0,0,0,0,0,p.nep*p.lp,(1/3)*p.nep*p.lp^3];


 % cep=g(6)*((p.zp1).^2)+g(7)*ones(1,length(p.zp1));
% cen=g1(2)*((p.zn).^2)+g1(1);
% ces=g1(5)*((p.zs1).^2)+g1(4)*p.zs1+g1(3);
% cep=flip(g1(7)*((p.zp1).^2)+g1(6));
% ce=[cen,ces,cep];

% cen1=g(2)*((p.zn).^2)+g(1);
% ces1=g(5)*((p.zs).^2)+g(4)*p.zs+g(3);
% cep1=g(7)*((p.zp).^2)+g(6);
% ce=[cen1,ces1,cep1];



% a=[2*de*p.ln,0,-2*de*p.ln,-de,0,0,0;...
%      0,0,de*2*(p.ln+p.ls),de,0,-2*de*(p.ln+p.ls),0;...
%    p.ln^2,1,-p.ln^2,-p.ln,-1,0,0;...
%    0,0,(p.ln+p.ls)^2,(p.ln+p.ls),1,-(p.ln+p.ls)^2,-1;...
%    (1/3)*p.nen*p.ln^3,p.nen*p.ln,0,0,0,0,0;...
% 0,0,(1/3)*p.nes*((p.ln+p.ls)^3-p.ln^3),.5*p.nes*((p.ln+p.ls)^2-p.ln^2)...
% ,p.nes*p.ls,0,0;0,0,0,0,0,(1/3)*p.nep*(p.l^3-(p.ln+p.ls)^3),p.nep*p.lp];

% de=-84.3-2.2*(1e-3).*ce -540./(p.ta-229-5*(1e-3)*ce);
% a=[1,p.ln^2,-1,0,0,0,0;0,0,-1,-p.ls,-p.ls^2,1,p.lp^2;
%    0,2*p.ln*den,0,-desn,0,0,0;0,0,0,desp,2*desp*p.ls,0,2*p.lp*dep;
%    p.nen*p.ln,1/3*p.nen*p.ln^3,0,0,0,0,0;
%    0,0,p.nes*p.ls,1/2*p.nes*p.ls^2,1/3*p.nes*p.ls^3,0,0;
%    0,0,0,0,0,p.nep*p.lp,1/3*p.nep*p.lp^3];
% den=(p.nen^p.b)*de(1,length(p.zn));desn=(p.nes^p.b)*de(1,length(p.zn)+1);
% desp=(p.nes^p.b)*de(1,length(p.zn)+length(p.zp));
% dep=(p.nep^p.b)*de(1,length(p.zn)+length(p.zp)+1);
% 
% dn=(p.nen^p.b)*mean(de(1,1:length(p.zn)));
% ds=(p.nes^p.b)*mean(de(1,length(p.zn)+1:length(p.zn)+length(p.zs)));
% dp=(p.nep^p.b)*mean(de(1,length(p.zn)+length(p.zs)+1:p.x));
% qf= qb+p.t*[ (u-u*p.tp)/p.f+den*2*g(2)*p.ln,...
%     desn*2*g(5)*p.ls,-(u-u*p.tp)/p.f+dep*2*g(7)*p.lp];