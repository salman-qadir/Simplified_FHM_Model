 function [csn,csp]=pdepem(cn,cp,jn,jp,p)

x = linspace(0,p.ln,50);
t = linspace(0,3.6e2,3.6e2);
jf=jn/p.f;
bc = @(xl,ul,xr,ur,t) hbc(xl,ul,xr,ur,t,jf,cn,p.ln);
ic = @(x) hic(x,cn);
yl=@(x,t,u,dudx) hyl(x,t,u,dudx,jf,p.dsn);
m = 0;
csn = pdepe(m,yl,ic,bc,x,t);

jf=jp/p.f;
x = linspace(0,p.lp,50);
bc = @(xl,ul,xr,ur,t) hbc(xl,ul,xr,ur,t,jf,cp,p.lp);
ic = @(x) hic(x,cp);
yl=@(x,t,u,dudx) hyl(x,t,u,dudx,jf,p.dsp);
m = 0;
csp = pdepe(m,yl,ic,bc,x,t);
end

% % u = sol(:,:,1);
% surf(x,t,sol)
% xlabel('x')
% ylabel('t')
% zlabel('u(x,t)')