function [ecd1,ecd3,ecdn,ecdp] = ecd(cn,cp,cen,cep,p)

ecd1=p.kn*real(((1-cn'/p.csn).*cn'.*cen').^p.aa);
ecd3=p.kp*real(((1-cp'/p.csp).*cp'.*cep').^p.aa);

ecdn=p.an*p.kn*real(((p.csn*ones(p.n,1)-cn').*cn'.*cen').^p.aa);
ecdp=p.ap*p.kp*real(((p.csp*ones(p.p,1)-cp').*cp'.*cep').^p.aa);

end

