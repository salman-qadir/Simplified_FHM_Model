function [c,f,s] = hyl(x,t,u,dudx,jf,p)
c = 1;
f = p*dudx;
s=-jf;
% s =-jf*ones(1,length(x));
end

