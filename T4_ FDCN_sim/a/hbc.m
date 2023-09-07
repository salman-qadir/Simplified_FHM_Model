function [pl,ql,pr,qr] = hbc(xl,ul,xr,ur,t,jf,cnn,ln)
pl = 0;%ul-cnn;
ql = 1;
pr = jf*ln;
qr = 1;
end