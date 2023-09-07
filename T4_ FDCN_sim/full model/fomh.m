function [jn,jp,cn,cp,cen,ces,cep,jn1,jp1,T]=fomh(y,p)

jn=y(1:p.n);
jp=y(p.n+1:p.n+p.p);
cn=y(p.n+p.p+1:2*p.n+p.p);
cp=y(2*p.n+p.p+1:2*(p.n+p.p));
cen=y(2*(p.n+p.p)+1:2*(p.n+p.p)+p.n);
ces=y(2*(p.n+p.p)+p.n+1:2*(p.n+p.p)+p.n+p.s);
cep=y(2*(p.n+p.p)+p.n+p.s+1:2*(p.n+p.p)+p.x);
jn1=y(2*(p.n+p.p)+p.x+1:2*(p.n+p.p)+p.x+p.n);
jp1=y(2*(p.n+p.p)+p.x+p.n+1:3*(p.n+p.p)+p.x);
T=y(p.x+3*(p.n+p.p)+1);

% cnn1=y(1:p.n);
% cpp1=y(p.n+1:p.n+p.p);
% qsn1=y(p.n+p.p+1:2*p.n+p.p);
% qsp1=y(2*p.n+p.p+1:2*p.n+2*p.p);
% 
% jn=y(2*(p.n+p.p)+1:2*(p.n+p.p)+p.n);
% jp=y(2*(p.n+p.p)+p.n+1:3*(p.n+p.p));
% % fic1=y(3*(p.n+p.p)+1:3*(p.n+p.p)+p.n);
% % fic3=y(3*(p.n+p.p)+p.n+1:4*(p.n+p.p));
% qb=y(3*(p.n+p.p)+1:3*(p.n+p.p)+3);
% g2=y(3*(p.n+p.p)+3+1:3*(p.n+p.p)+3+7);
% % c_e=y(p.r*(p.n+p.p)+1:p.r*(p.n+p.p)+p.x-3);
% % cen=c_e(1:p.n);ces=c_e(p.n+1:p.n+p.s-3);
% % cep=c_e(p.n+p.s-3+1:p.x-3);
% % [csn,csp,c_e,cen,ces,cep,j11,j33,j4,jj1,jj3] = fom(y,p)
end

