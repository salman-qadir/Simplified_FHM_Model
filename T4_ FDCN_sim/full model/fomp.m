function [jn,jp,cen,ces,cep,csn,csp,j1,j3,j11,j33,jn1,jp1,T]=fomp(y,p)
% out=[jn',jp',ceq,cnn,cpp];

jn=y(1:p.n);
jp=y(p.n+1:p.n+p.p);
cen=y(p.n+p.p+1:2*p.n+p.p);
ces=y(2*p.n+p.p+1:2*p.n+p.p+p.s);
cep=y(2*p.n+p.p+p.s+1:p.n+p.p+p.x);
csn=y(p.n+p.p+p.x+1:p.n+p.p+p.x+p.rr*p.n);
csp=y(p.n+p.p+p.x+p.rr*p.n+1:p.x+(p.rr+1)*(p.n+p.p));

jn1=y(p.x+(p.rr+1)*(p.n+p.p)+1:p.x+(p.rr+1)*(p.n+p.p)+p.n );
jp1=y(p.x+(p.rr+1)*(p.n+p.p)+p.n+1:p.x+(p.rr+2)*(p.n+p.p));
T=y(p.x+(p.rr+2)*(p.n+p.p)+1);

j1=kron(jn,ones(1,p.rr))';j3=kron(jp,ones(1,p.rr))';
j11=kron(jn1,ones(1,p.rr))';j33=kron(jp1,ones(1,p.rr))'; 

% csn=y(p.n+p.p+p.x+1:p.n+p.p+p.x+p.rr);
% csp=y(p.n+p.p+p.x+p.rr+1:p.x+2*p.rr+p.n+p.p);
% 
% j11=y(p.x+2*p.rr+p.n+p.p+1:p.x+p.rr*2+p.n+p.p+p.n );
% j33=y(p.x+p.rr*2+p.n+p.p+p.n+1:p.x+p.rr*2+2*(p.n+p.p) );
%%
% csn=y(1:p.r*p.n);
% csp=y(p.r*p.n+1:p.r*(p.n+p.p));
% 
% c_e=y(p.r*(p.n+p.p)+1:p.r*(p.n+p.p)+p.x-3);
% cen=c_e(1:p.n);ces=c_e(p.n+1:p.n+p.s-3);
% cep=c_e(p.n+p.s-3+1:p.x-3);

% j11=y(p.r*(p.n+p.p)+p.x-3+1:p.r*(p.n+p.p)+p.x-3+p.n);
% j33=y(p.r*(p.n+p.p)+p.x-3+1+p.n:(p.r+1)*(p.n+p.p)+p.x-3);
% j4=[j11;zeros(p.s-3,1);j33];

% we=repmat(jn,1,p.rr);[q1,q2]=size(we);j1=reshape(we',q1*q2,1); 
% we=repmat(jp,1,p.rr);[q1,q2]=size(we);j3=reshape(we',q1*q2,1);
% we=repmat(jnn,1,p.rr);[q1,q2]=size(we);j11=reshape(we',q1*q2,1); 
% we=repmat(jpp,1,p.rr);[q1,q2]=size(we);j33=reshape(we',q1*q2,1);
end

